# -*- coding: utf-8 -*-
"""
@File   ： ablation study.py
@Time   ： 2022/7/31 09:41
@Author ： Jia Yiming
"""
from __future__ import print_function

import math

from matplotlib import pyplot as plt

"""
计算p-value判断结果表现
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import mygene
import json
import requests

gene_info = pd.read_csv("human_fetal_immune/gene_meta_3kgenes_final.csv")
non_zero_index = np.load('non_zero_index1.npy')
data0 = np.load('mu1.npy')
data1 = np.load('mu_t1.npy')
data2 = np.load('mu_f.npy')

adata0 = ad.AnnData(data0)
adata0.obs_names = gene_info.loc[non_zero_index, "gene_id"]
adata0.var_names = ["C" + str(i) for i in range(20)]
adata0.obs['indexa'] = range(len(non_zero_index))
sc.pp.neighbors(adata0, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata0)
sc.tl.leiden(adata0)
sc.pl.umap(adata0, color=['leiden'])
for i in range(len(gene_info)):
    gene_info.loc[i, "gene_id"] = gene_info.loc[i, "gene_id"][:15]
adata0.obs_names = gene_info.loc[non_zero_index, "gene_id"]

adata1 = ad.AnnData(data1)
adata1.obs_names = gene_info.loc[non_zero_index, "gene_id"]
adata1.var_names = ["C" + str(i) for i in range(20)]
adata1.obs['indexa'] = range(len(non_zero_index))
sc.pp.neighbors(adata1, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata1)
sc.tl.leiden(adata1)
sc.pl.umap(adata1, color=['leiden'])
for i in range(len(gene_info)):
    gene_info.loc[i, "gene_id"] = gene_info.loc[i, "gene_id"][:15]
adata1.obs_names = gene_info.loc[non_zero_index, "gene_id"]

adata2 = ad.AnnData(data2)
adata2.obs_names = gene_info.loc[non_zero_index, "gene_id"]
adata2.var_names = ["C" + str(i) for i in range(10)]
adata2.obs['indexa'] = range(len(non_zero_index))
sc.pp.neighbors(adata2, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata2)
sc.tl.leiden(adata2)
sc.pl.umap(adata2, color=['leiden'])
for i in range(len(gene_info)):
    gene_info.loc[i, "gene_id"] = gene_info.loc[i, "gene_id"][:15]
adata1.obs_names = gene_info.loc[non_zero_index, "gene_id"]
mg = mygene.MyGeneInfo()

p_values0 = []
p_values_log0 = []
p_values1 = []
p_values_log1 = []
p_values2 = []
p_values_log2 = []

for i in range(len(adata0.obs['leiden'].unique())):
    ids = []
    gene_ids = mg.getgenes(adata0.obs[adata0.obs['leiden'] == str(i)].index, 'name, symbol, entrezgene', as_dataframe=True)
    gene_ids.index.name = "UNIPROT"
    gene_ids.reset_index(inplace=True)
    for i in list(gene_ids.loc[:, "_id"]):
        if isinstance(i, str) and not i.startswith("E"):
            ids.append(int(i))
    data = {"Genes": ids, "Categories": [
        {
            "Type": "GeneOntologyMolecularFunction",
            "PValue": 0.05,
            "MinGenes": 1,
            "MaxGenes": 1500,
            "MaxResults": 50,
            "Correction": "FDR"
        },
        {
            "Type": "GeneOntologyBiologicalProcess",
            "PValue": 0.05,
            "MinGenes": 1,
            "MaxGenes": 1500,
            "MaxResults": 50,
            "Correction": "FDR"
        },
        {
            "Type": "GeneOntologyCellularComponent",
            "PValue": 0.05,
            "MinGenes": 1,
            "MaxGenes": 1500,
            "MaxResults": 50,
            "Correction": "FDR"
        }
    ]}
    j_str = json.dumps(data)

    res = requests.post(url='https://toppgene.cchmc.org/API/enrich',
                        headers={"Content-Type": "text/json"},
                        data=j_str)

    data = json.loads(res.text)
    if data["Annotations"] is None:
        continue
    for i in data["Annotations"]:
        p_values_log0.append(math.log10(i["PValue"]))
        p_values0.append(i["PValue"])


for i in range(len(adata1.obs['leiden'].unique())):
    ids = []
    if len(adata1.obs[adata1.obs['leiden'] == str(i)]) == 0:
        continue
    gene_ids = mg.getgenes(adata1.obs[adata1.obs['leiden'] == str(i)].index, 'name, symbol, entrezgene', as_dataframe=True)
    gene_ids.index.name = "UNIPROT"
    gene_ids.reset_index(inplace=True)
    for i in list(gene_ids.loc[:, "_id"]):
        if isinstance(i, str) and not i.startswith("E"):
            ids.append(int(i))
    data = {"Genes": ids, "Categories": [
        {
            "Type": "GeneOntologyMolecularFunction",
            "PValue": 0.05,
            "MinGenes": 1,
            "MaxGenes": 1500,
            "MaxResults": 50,
            "Correction": "FDR"
        },
        {
            "Type": "GeneOntologyBiologicalProcess",
            "PValue": 0.05,
            "MinGenes": 1,
            "MaxGenes": 1500,
            "MaxResults": 50,
            "Correction": "FDR"
        },
        {
            "Type": "GeneOntologyCellularComponent",
            "PValue": 0.05,
            "MinGenes": 1,
            "MaxGenes": 1500,
            "MaxResults": 50,
            "Correction": "FDR"
        }
    ]}
    j_str = json.dumps(data)

    res = requests.post(url='https://toppgene.cchmc.org/API/enrich',
                        headers={"Content-Type": "text/json"},
                        data=j_str)
    data = json.loads(res.text)
    if data["Annotations"] is None:
        continue
    for i in data["Annotations"]:
        p_values_log1.append(math.log10(i["PValue"]))
        p_values1.append(i["PValue"])

for i in range(len(adata2.obs['leiden'].unique())):
    ids = []
    if len(adata2.obs[adata1.obs['leiden'] == str(i)]) == 0:
        continue
    gene_ids = mg.getgenes(adata2.obs[adata2.obs['leiden'] == str(i)].index, 'name, symbol, entrezgene', as_dataframe=True)
    gene_ids.index.name = "UNIPROT"
    gene_ids.reset_index(inplace=True)
    for i in list(gene_ids.loc[:, "_id"]):
        if isinstance(i, str) and not i.startswith("E"):
            ids.append(int(i))
    data = {"Genes": ids, "Categories": [
        {
            "Type": "GeneOntologyMolecularFunction",
            "PValue": 0.05,
            "MinGenes": 1,
            "MaxGenes": 1500,
            "MaxResults": 50,
            "Correction": "FDR"
        },
        {
            "Type": "GeneOntologyBiologicalProcess",
            "PValue": 0.05,
            "MinGenes": 1,
            "MaxGenes": 1500,
            "MaxResults": 50,
            "Correction": "FDR"
        },
        {
            "Type": "GeneOntologyCellularComponent",
            "PValue": 0.05,
            "MinGenes": 1,
            "MaxGenes": 1500,
            "MaxResults": 50,
            "Correction": "FDR"
        }
    ]}
    j_str = json.dumps(data)

    res = requests.post(url='https://toppgene.cchmc.org/API/enrich',
                        headers={"Content-Type": "text/json"},
                        data=j_str)
    data = json.loads(res.text)
    if data["Annotations"] is None:
        continue
    for i in data["Annotations"]:
        p_values_log2.append(math.log10(i["PValue"]))
        p_values2.append(i["PValue"])


x_norm0=np.array(p_values_log0)
x_norm1=np.array(p_values_log1)
x_norm2=np.array(p_values_log2)
plt.rcParams['axes.unicode_minus']=False#显示负号\n",
plt.figure(figsize=(6,4))## 设置画布\n",
plt.hist(x_norm1,bins=50,color='b')
plt.hist(x_norm0,bins=50,color='r')
plt.hist(x_norm2,bins=50,color='g')


plt.show()

print(len(x_norm0))
print(len(adata0.obs['leiden'].unique()))
print(len(x_norm0)/len(adata0.obs['leiden'].unique()))
print(len(x_norm1))
print(len(adata1.obs['leiden'].unique()))
print(len(x_norm1)/len(adata1.obs['leiden'].unique()))
print(len(x_norm2))
print(len(adata2.obs['leiden'].unique()))
print(len(x_norm2)/len(adata2.obs['leiden'].unique()))

np.save("p_value", p_values0)
np.save("p_value_t", p_values1)
np.save("p_value_f", p_values2)