# -*- coding: utf-8 -*-
"""
@File   ： build_graph_mouse.py
@Time   ： 2022/10/25 09:44
@Author ： Jia Yiming
"""
"""
生成knn，用于之后的gcn
"""

import numpy as np
import pandas as pd
import mygene

graph_df = pd.read_csv('mippie_subnet3.txt', sep='\t', header=None)

graph_df = graph_df.iloc[:,:3]

print(graph_df)
a = pd.DataFrame(graph_df.iloc[:,0])
b = pd.DataFrame(graph_df.iloc[:,1])
a.columns = ['names']
b.columns = ['names']

a =  pd.Series(a.iloc[1:,0])
b =  pd.Series(b.iloc[1:,0])
graph_genes = pd.unique(pd.concat([a,b]))

print(graph_genes)

gene_info = pd.read_csv("/Users/jiayiming/Desktop/麦吉尔暑研/scSTEM_sample_data/mouse_embryo_neural/gene_meta_5kgenes_final.csv")
non_zero_index = np.load('non_zero_index3.npy')

for i in range(len(gene_info)):
    gene_info.loc[i,"gene_id"] = gene_info.loc[i,"gene_id"][:18]
print(gene_info.head())

mg = mygene.MyGeneInfo()
gene_id = mg.getgenes(gene_info.loc[:,"gene_id"], fields='_id')
ids = [gene_id[i].get('_id') for i in range(len(gene_id))]

id_idx = {}
for i in range(len(gene_id)):
    if ids[i] is None or ids[i].startswith('ENS'):
        continue
    id_idx[ids[i]] = i
print(id_idx)

edge = []
for i in range(len(graph_df)):
    if id_idx.__contains__(graph_df.iloc[i][0]) and id_idx.__contains__(graph_df.iloc[i][1]):
        node1 = id_idx[graph_df.iloc[i][0]]
        node2 = id_idx[graph_df.iloc[i][1]]
        score = graph_df.iloc[i][2]
        if node1 < node2:
            edge.append([node1, node2, score])
        elif node1 > node2:
            edge.append([node2, node1, score])


node_id = []
for i in graph_genes:
    if id_idx.__contains__(i):
        node_id.append(id_idx[i])

graph_table = {}
for i in node_id:
    graph_table[i] = []

for i in edge:
    if i[0] != i[1]:
        graph_table[i[0]].append(i[2])
        graph_table[i[1]].append(i[2])
for i in node_id:
    graph_table[i].sort(reverse=True)

k = 10
edge_threshold = {}
for i in node_id:
    if len(graph_table[i]) == 0:
        edge_threshold[i] = 1
    elif len(graph_table[i]) < k:
        edge_threshold[i] = graph_table[i][-1]
    else:
        edge_threshold[i] = graph_table[i][k-1]

index = []
for i in range(len(edge)):

    if edge[i][2] >= min(edge_threshold[edge[i][0]],edge_threshold[edge[i][1]]):
        index.append([edge[i][0], edge[i][1]])
        index.append([edge[i][1], edge[i][0]])

graph_index = np.array(index)    #存储knn
print(graph_index.shape)

print(graph_index)
np.save("graph_index3", graph_index)