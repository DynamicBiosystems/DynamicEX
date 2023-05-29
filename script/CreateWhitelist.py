import time
import sys
import logging
import pandas as pd
from pathlib import Path
from pandarallel import pandarallel
###barcode
##import umi_tools.network as network
from umi_tools._dedup_umi import edit_distance
###fastq
import fastqandfurious.fastqandfurious as fqf
from fastqandfurious.fastqandfurious import entryfunc

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--expect-cells',required=True)
parser.add_argument('--nthreads',required=True,help="maximum number of threads to use")
parser.add_argument('--inputdir',required=True,help="inputdir")
parser.add_argument('--bc-pattern',required=True)
args = parser.parse_args()

cutoff = int(int(args.expect_cells)*1.5)
threads = int(args.nthreads)
inputdir = args.inputdir
bcpattern = args.bc_pattern

logger = logging.getLogger()
logger.setLevel("INFO")

###并行操作
pandarallel.initialize(nb_workers=threads)

###
dd={}
d={}

bufsize = 30000
start = time.time()
barcode_len = bcpattern.count("C")

def barcode_counts(file):
    counts = {}
    with open(file, "rb") as fh:
        it = fqf.readfastq_iter(fh,bufsize,entryfunc)
        for entry in it:
            barcode = entry[1][:barcode_len]
            counts.setdefault(barcode,0)
            counts[barcode]+=1
    return counts
    
logging.info(f"Loading barcode")
split_file = list(Path(inputdir).glob("xa*_whitelist"))
df_bc = pd.DataFrame(split_file,columns=['files'])
#print(df_bc)

df_bc['barcode_counts'] = df_bc["files"].parallel_apply(barcode_counts)

logging.info(f"Statistics of barcode completed,time consuming:{(time.time()-start)/60:.2f} min")
logging.info(f"Merge barcode")
start2 = time.time()
#for info in df_bc['barcode_counts'].values:
#    for key in info:
#        dd.setdefault(key,set())
#        dd[key].update(info.get(key))
#d = {key:len(dd.get(key)) for key in dd}

for info in df_bc['barcode_counts'].values:
    for key in info:
        d.setdefault(key,0)
        d[key]+=info.get(key)


logging.info(f"Merge barcode complete,time consuming:{(time.time()-start2):.2f} sec")

logging.info(f"Total barcode:{len(d):,}")
logging.info(f"Started barcode cut off")
df = pd.DataFrame(d.items(),columns=['barcode','counts'])
df = df.sort_values(by=['counts'],ascending=False)
df['part'] = ["part1" if i < cutoff else "part2" for i in range(len(df))]

part1 = df.query("part=='part1'")
part2 = df.query("part=='part2'")

def edit_distance1(i):
    for barcode in part1['barcode'].values:
        if edit_distance(i,barcode) <= 1:
            return True
    return False

start3 = time.time()
part2["match"]=part2["barcode"].parallel_apply(edit_distance1)
part2 = part2.query("match == False")

logging.info(f"Finished barcode cut off,time consuming:{(time.time()-start3)/60:.02f} min")

logging.info(f"Remaining barcode:{len(part2):,}")
logging.info(f"Saving whitelist")

with open(str(Path(inputdir).parent/"barcode_whitelist.txt"),'w') as out:
    for i in part1['barcode'].values:
        out.write(i.decode("utf-8")+"\n")
    for i in part2['barcode'].values:
        out.write(i.decode("utf-8")+"\n")

logging.info(f"Finished successfully,time consuming:{(time.time()-start)/60:.02f} min")
