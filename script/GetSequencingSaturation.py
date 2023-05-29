import sys
import pandas as pd
from pathlib import Path
import numpy as np
import pysam
from pandarallel import pandarallel
import os
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="ticks")

parser = argparse.ArgumentParser()
parser.add_argument('--bam',required=True,help="bam file")
parser.add_argument('--barcodes',required=True,help="barcodes.tsv")
parser.add_argument('--summary',required=True,help="star summary")
parser.add_argument('--nthreads',required=True,help="maximum number of threads to use")
parser.add_argument('--ouputdir',required=True,help="outputdir")
parser.add_argument('--reads',action='store_true',help="Calculate reads saturation")

args = parser.parse_args()
###并行操作
pandarallel.initialize(nb_workers=int(args.nthreads))
bamfile = args.bam
barcodes = args.barcodes
summary = args.summary

df = pd.read_csv(summary,index_col=0,header=None)

Number_of_Reads = df.loc['Number of Reads'].values[0]
mgene = pd.read_csv(summary,index_col=0,header=None)
if 'Median Gene per Cell' in df.index:
    mgene = df.loc['Median Gene per Cell'].values[0]
else:
    mgene = df.loc['Median GeneFull per Cell'].values[0]

output = Path(args.ouputdir)

temp_sat = output/".saturation"
temp_sat.mkdir(exist_ok=True)


barcode = {line.strip() for line in open(barcodes)}
bam = pysam.AlignmentFile(bamfile, "rb")

def counts(ref):
    bam = pysam.AlignmentFile(bamfile, "rb")
    out = open(temp_sat/f"{ref}.txt","w",encoding="utf-8")
    l=[]
    for read in bam.fetch(reference=ref):
        CB = read.get_tag("CB") ###矫正后的Barcode
        UB = read.get_tag("UB") ###矫正后的UMI
        GX = read.get_tag("GX") ###Geneid
        if all([i!="-" for i in (CB,UB,GX)]):
            out.write(f"{CB}\t{UB}\t{GX}\n")


chrom = [bam.get_reference_name(i) for i in range(bam.nreferences)]
df = pd.DataFrame(chrom,columns=['chrom'])


df_chr = df[df.chrom.str.contains("chr")]
df_other = df[~df.chrom.str.contains("chr")]

if len(df_chr) > 0:
    df_chr['counts'] = df_chr['chrom'].parallel_apply(counts)
if len(df_other) > 0:
    df_other['counts'] = df_other['chrom'].parallel_apply(counts)

os.system(f"cat {temp_sat}/* > {output/'bam2df.txt'}")

df = pd.read_csv(str(output/'bam2df.txt'),sep='\t',header=None,dtype="category",engine="c")

def downsample(frac):
    if frac == 0:
        return 0,0,0
    d = {}
    n_deduped_reads = 0
    N_umis = 0
    N_reads = 0
    dt = df.sample(frac=frac,random_state=0)
    medain_gene = {}
    cells = set(df[0].values)&barcode
    total_reads = frac*Number_of_Reads
    for value in dt.values:
        bc,umi,gene = value
        if bc in cells:
            medain_gene.setdefault(value[0],{})
            medain_gene[value[0]][value[-1]] = 1
            key = "".join(value)
            d.setdefault(key,0)
            d[key]+=1
    for key in d:
        if d[key] == 1:
            n_deduped_reads+=1
        N_umis+=1
        N_reads+=d[key]
    umis_saturation = (1-n_deduped_reads/N_umis)
    reads_saturation = (1-n_deduped_reads/N_reads)
    medain_gene2 = {bc:len(medain_gene[bc]) for bc in medain_gene}
    Median_Gene_Per_Cell = pd.DataFrame(medain_gene2.items())[1].median()
    return (umis_saturation,total_reads/len(cells),int(Median_Gene_Per_Cell))



plot = pd.DataFrame([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],columns=["frac"])
plot['param'] = plot['frac'].parallel_apply(downsample)
ll = []
for i in plot['param'].values:
    ll.append(list(i))
    
dft = pd.DataFrame(ll,columns=['saturation','Mean_Reads_per_cell',"Median_Gene_Per_Cell"])
diffgene = mgene - list(dft['Median_Gene_Per_Cell'])[-1]
diffgene = round(diffgene)
def correct_mgene(elem):
    if elem == 0:
        return elem
    return elem+diffgene
    
dft['Median_Gene_Per_Cell'] = dft['Median_Gene_Per_Cell'].apply(correct_mgene)
dft.to_csv("saturation.csv",index=False)
#pd.DataFrame(ll,columns=['saturation','Mean_Reads_per_cell',"Median_Gene_Per_Cell"]).to_csv("saturation.csv",index=False)
#plot.to_csv("saturation.csv")


yl,xl = [],[]
median_gene = []
for y,x,y2 in dft.values:
    yl.append(y)
    xl.append(x)
    median_gene.append(y2)

plt.figure(figsize=(10,8))
plt.plot(xl,median_gene,lw=2,color='blue')
plt.ylabel("Median Gene Per Cell",fontsize=12,fontweight="bold")
plt.xlabel("Mean Reads per Cell",fontsize=12,fontweight="bold")
plt.xlim(0,)
plt.ylim(0,)
#plt.yticks([0,0.2,0.4,0.6,0.8,1.0])

#####坐标轴
plt.gca().spines["right"].set_alpha(.0)
plt.gca().spines["top"].set_alpha(.0)

plt.savefig(str(output/"Median_Gene_Per_Cell.pdf"))

plt.figure(figsize=(10,8))
plt.plot(xl,yl,lw=2,color='blue')
plt.ylabel("Sequencing Saturation",fontsize=12,fontweight="bold")
plt.xlabel("Mean Reads per Cell",fontsize=12,fontweight="bold")
plt.xlim(0,)
plt.ylim(0,)
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])

#####坐标轴
plt.gca().spines["right"].set_alpha(.0)
plt.gca().spines["top"].set_alpha(.0)
plt.savefig(str(output/"SequencingSaturation.pdf"))


df = pd.read_csv(summary,index_col=0,header=None)
df.loc['Sequencing Saturation'][1] = yl[-1]
df.to_csv(summary,header=False)


