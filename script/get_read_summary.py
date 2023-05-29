import sys
import pandas as pd
from pathlib import Path

outname = sys.argv[1]
l=[]
files = list(Path().cwd().rglob("**/*_reads_summary.csv"))

for file in files:
    name = file.name.split("_")[0]
    df = pd.read_csv(file,index_col=0)

    df = df.sum().to_frame(name)
    total = df.loc['total'].values[0]
    df[name] = df[name].apply(lambda e:f"{e:,}({e/total:.2%})")
    df = df.T
    df.columns=['总reads数','有效reads数','不含ployT','不符合固定序列','插入片段过短']
    l.append(df)
pd.concat(l).to_excel(outname)

