import pandas as pd
import sys
import os
from pathlib import Path

starsummary,wpssummary,threads,pigz = sys.argv[1:5]

df = pd.read_csv(starsummary,index_col=0,header=None)
df = df.dropna()
df.to_csv(wpssummary,header=False)

filtered = Path(starsummary).parent/"filtered"
raw = Path(starsummary).parent/"raw"

if not list(filtered.glob("*gz")):
    os.system(f"{pigz} -p {threads} {filtered}/*")
    
if not list(raw.glob("*gz")):
    os.system(f"{pigz} -p {threads} {raw}/*")