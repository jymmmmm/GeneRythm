# -*- coding: utf-8 -*-
"""
@File   ： clustering_gene_embedding.py
@Time   ： 2023/10/17 18:14
@Author ： Jia Yiming
"""
from __future__ import print_function

import mygene
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline
from sklearn.cluster import KMeans
from scipy.fftpack import fft

"""
计算p-value判断结果表现
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import stats
import seaborn as sns



gene_info = pd.read_csv("human_fetal_immune/gene_meta_3kgenes_final.csv")
non_zero_index = np.load('non_zero_index1.npy')
data = np.load('mu1.npy')
adata = ad.AnnData(data)
adata.obs_names = gene_info.loc[non_zero_index, "gene_id"]
adata.var_names = ["C" + str(i) for i in range(20)]
adata.obs['indexa'] = range(len(non_zero_index))

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)

sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden'], legend_loc='on data', legend_fontsize='x-large', save='test.pdf')

for i in range(len(gene_info)):
    gene_info.loc[i, "gene_id"] = gene_info.loc[i, "gene_id"][:15]

adata.obs_names = gene_info.loc[non_zero_index, "gene_id"]