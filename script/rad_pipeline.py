import re
import os
import json
import logging
import argparse
from glob import glob
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('-d',required=True,help="Data directory")
parser.add_argument('-s',required=True,help="samples")

args = parser.parse_args()

dataDir = str(args.d).split(",")
samples = str(args.s).split(",")

gzfile = set()
for i in dataDir:
    for sample in Path(i).rglob("**/*gz"):
        gzfile.add(sample)

d={}
for ID in samples:
    d.setdefault(ID,{})
    d[ID].setdefault("R1",[])
    d[ID].setdefault("R2",[])
    for sample in gzfile:
        sample = str(sample)
        if re.match(f"*{ID}*_R1*f*q.gz",sample):
            d[ID]['R1'].append(sample)
            R2 = sample.replace("_R1","_R2")
            if Path(R2).exists():
                d[ID]['R2'].append(R2)

for i in d:
    d[i]['R1'].sort()
    d[i]['R2'].sort()
df = pd.DataFrame(d)
print(df)    
