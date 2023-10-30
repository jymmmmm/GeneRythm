# -*- coding: utf-8 -*-
"""
@File   ： build_graph_human.py
@Time   ： 2022/8/7 17:22
@Author ： Jia Yiming
"""

"""
生成knn，用于之后的gcn
"""

import numpy as np
import pandas as pd
import mygene

graph_df = pd.read_csv('out.txt', sep='\t', header=None)

graph_df = graph_df.iloc[:,:5]

a = pd.DataFrame(graph_df.iloc[:,0])
b = pd.DataFrame(graph_df.iloc[:,2])
a.columns = ['names']
b.columns = ['names']

a =  pd.Series(a.iloc[:,0])
b =  pd.Series(b.iloc[:,0])

graph_genes = pd.unique(pd.concat([a,b]))

gene_info = pd.read_csv("/Users/jiayiming/Desktop/麦吉尔暑研/scSTEM_sample_data/human_fetal_immune/gene_meta_3kgenes_final.csv")
non_zero_index = np.load('non_zero_index1.npy')
# non_zero_index = np.array(range(3000))

for i in range(len(gene_info)):
    gene_info.loc[i,"gene_id"] = gene_info.loc[i,"gene_id"][:15]

mg = mygene.MyGeneInfo()
gene_name =  mg.getgenes(gene_info.loc[:,"gene_id"], fields='symbol')
names = [gene_name[i].get('symbol') for i in range(len(gene_name))]

names_dic = {}
for i in range(len(non_zero_index)):
    if names[non_zero_index[i]] is not None:
        names_dic[names[non_zero_index[i]]] = i

edge = []
for i in range(len(graph_df)):
    if names_dic.__contains__(graph_df.iloc[i][0][:-6]) and names_dic.__contains__(graph_df.iloc[i][2][:-6]):
        node1 = names_dic[graph_df.iloc[i][0][:-6]]
        node2 = names_dic[graph_df.iloc[i][2][:-6]]
        score = graph_df.iloc[i][4]
        if node1 < node2:
            edge.append([node1, node2, score])
        elif node1 > node2:
            edge.append([node2, node1, score])

node_id = []
for i in graph_genes:
    if names_dic.__contains__(i[:-6]):
        node_id.append(names_dic[i[:-6]])

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
np.save("graph_index", graph_index)