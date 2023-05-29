#0d08ef48983e

import pandas as pd
import scanpy as sc
import numpy as np
from pathlib import Path
import argparse
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser()
parser.add_argument('--inputdir',required=True,help="10X martix dir")

args = parser.parse_args()

inputdir = args.inputdir
basedir = Path(inputdir).parent/"analysis/clustering"

adata = sc.read_10x_mtx(inputdir)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)
sc.tl.umap(adata)
sc.tl.leiden(adata,resolution=0.9)

leiden_df = adata.obs['leiden'].to_frame(name="Cluster")
leiden_df.index.name = "Barcode"
leiden_df.to_csv(basedir/"graphclust/clusters.csv")

umap_df = pd.DataFrame(adata.obsm['X_umap'],index=adata.obs_names,columns=['umap1','umap2'])
umap_df.index.name = "Barcode"
###Barcode,Cluster
for n_clusters in range(2,11):
    kmeans = KMeans(n_clusters=n_clusters)
    (basedir/f'kmeans_{n_clusters}_clusters').mkdir(parents=True, exist_ok=True)
    umap_df['Cluster'] = kmeans.fit_predict(umap_df.values)
    temp = umap_df['Cluster'].to_frame(name="Cluster")
    temp.index.name = "Barcode"
    temp.to_csv(basedir/f'kmeans_{n_clusters}_clusters/clusters.csv')






