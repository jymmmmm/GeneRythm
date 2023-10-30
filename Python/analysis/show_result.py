# -*- coding: utf-8 -*-
"""
@File   ： check_result.py
@Time   ： 2022/7/31 09:41
@Author ： Jia Yiming
"""
from __future__ import print_function

import mygene
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline
from sklearn.cluster import KMeans
from scipy.fftpack import fft

"""
calculate and check the result
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


trajectory_info = pd.read_csv("/Users/jiayiming/Desktop/麦吉尔暑研/summer_intern/human_fetal_immune.csv")
a = pd.Index(trajectory_info.loc["path"] == 5)  # dataset1
y = np.array(trajectory_info.iloc[:3000, a])
x = np.array(trajectory_info.loc["time", a])

x.sort()
x_smooth = np.linspace(start=x.min(), stop=x.max(), num=100)
y_smooth = make_interp_spline(x, y.transpose())(x_smooth)
yf_smooth = fft(y_smooth.transpose())

yf_info = np.array([abs(yf_smooth[i, :50]) for i in range(yf_smooth.shape[0])])

for i in range(len(yf_info)):
    if yf_info[i].max() == 0:
        continue
    yf_info[i] = yf_info[i] / yf_info[i].max()

for i in range(len(y)):
    y[i, :] = y[i, :] - y[i, 0]


cluster = adata.obs[adata.obs['leiden'] == '8']['indexa']


# np.save('cluster1_8.npy', cluster)
# cluster = np.array(cluster)
# np.save('cluster_index.npy',cluster)


def movingaverage(data, window_size):
    window = np.ones(int(window_size)) / float(window_size)
    a = np.convolve(data[window_size:], window, 'valid')
    return np.concatenate([np.zeros(window_size), a])


for j in [8]:
    cluster_num = str(j)
    cluster = adata.obs[adata.obs['leiden'] == cluster_num]['indexa']

    # for s in adata.obs_names[adata.obs['leiden'] == cluster_num]:
    #     print(s)
    # print("********************")
    error_t = []
    for i in range(len(cluster)):
        # print(np.max(y[cluster[i],:] - np.min(y[cluster[i],:])))
        error_t.append(np.max(y[cluster[i], :] - np.min(y[cluster[i], :])))

    error_index_t = np.argsort(error_t)
    error_t = np.array(error_t)
    num = 0
    cluster_filter = []
    yf_info_filter = []
    cluater_index = []
    name = []
    for i in error_index_t[:]:
        # a = movingaverage(y[cluster[i],:],5)
        a = movingaverage(y[cluster[i], :], 5)
        if abs(np.max(abs(a))) < 0.1 or abs(np.max(abs(a))) > 5:
            continue
        num = num + 1
        cluster_filter.append(a)
        cluater_index.append(cluster[i])
        yf_info_filter.append(yf_info[cluster[i], :])
        name.append(adata.obs_names[cluster[i]])
    print(len(name))
    cluster_filter = np.array(cluster_filter)
    cluater_index = np.array(cluater_index)
    yf_info_filter = np.array(yf_info_filter)


    color_map = ['b', 'r', 'y', 'c', 'g']
    n_cluster = 5
    index = [[] for i in range(n_cluster)]

    if len(cluster_filter) >= 10:
        y_pred = KMeans(n_clusters=n_cluster, random_state=9).fit_predict(cluster_filter)
        color = []
        number_of_genes = np.zeros([n_cluster])
        for i in range(len(y_pred)):
            index[y_pred[i]].append(i)
        for i in y_pred:
            number_of_genes[i] = number_of_genes[i] + 1
            color.append(color_map[i])
        for i in range(len(cluster_filter)):
            if number_of_genes[y_pred[i]] < 5:
                continue
            # plt.plot(x[:-4], cluster_filter[i], c=color[i])
        plt.xlabel("Pseudotime", fontdict={'family': 'Arial', 'size': 18})
        plt.ylabel("Expression log2 fold change", fontdict={'family': 'Arial', 'size': 18})
        plt.title('cluster' + str(j), fontdict={'family': 'Arial', 'size': 18})
        filename = 'figure of ds1/time data' + cluster_num + '.pdf'
        # plt.savefig(filename, bbox_inches="tight")
        # plt.show()
        plt.close()
        print(y_pred)

        avg = np.zeros([n_cluster, cluster_filter.shape[1]])
        std = np.zeros([n_cluster, cluster_filter.shape[1]])
        yf_avg = np.zeros([n_cluster, yf_info_filter.shape[1]])
        for i in range(n_cluster):
            avg[i] = np.mean(cluster_filter[index[i]], axis=0)
            std[i] = np.std(cluster_filter[index[i]], axis=0)
            yf_avg[i] = np.mean(yf_info_filter[index[i]], axis=0)
        m = 0
        for i in range(n_cluster):
            if number_of_genes[i] < 5:
                continue
            # plt.plot(range(cluster_filter.shape[1]), avg[i], c = color_map[i])
            plt.errorbar(x[:-4], avg[i], yerr=std[i], c=color_map[i], label = 'path'+ str(m))
            m = m+1
        plt.xlabel("Pseudotime", fontdict={'family': 'Arial', 'size': 18})
        plt.ylabel("Expression log2 fold change", fontdict={'family': 'Arial', 'size': 18})
        plt.title('cluster' + str(j), fontdict={'family': 'Arial', 'size': 18})
        plt.legend()
        filename = 'figure of ds1/avg time data' + cluster_num + '.pdf'
        # plt.savefig(filename, bbox_inches="tight")
        plt.show()
        plt.close()
        avgf_smooth = fft(avg.transpose())
        print(avgf_smooth.shape)

        print('cluster' + str(j))
        m = 0
        for i in range(n_cluster):
            if number_of_genes[i] < 5:
                continue
            # plt.plot(range(cluster_filter.shape[1]), avg[i], c = color_map[i])
            plt.plot(range(len(yf_avg[i])), yf_avg[i], c=color_map[i], label = 'path'+str(m))
            m = m+1
        plt.xlabel("Frequency", fontdict={'family': 'Arial', 'size': 18})
        plt.ylabel("Amplitude", fontdict={'family': 'Arial', 'size': 18})
        plt.title('cluster' + str(j), fontdict={'family': 'Arial', 'size': 18})
        plt.legend()
        filename = 'figure of ds1/frequency' + cluster_num + '.pdf'
        # plt.savefig(filename, bbox_inches="tight")
        plt.show()
        plt.close()
        plt.gcf().clear()
        print(yf_info_filter.shape)


    # mg = mygene.MyGeneInfo()
    # gene_ids = mg.getgenes(name, 'name, symbol, entrezgene', as_dataframe=True)
    # gene_ids.index.name = "UNIPROT"
    # gene_ids.reset_index(inplace=True)
    # for index, row in gene_ids.iterrows():
    #     # print(row['symbol'])
    #     if pd.isna(row['symbol']):
    #         gene_ids.at[index, 'symbol'] = name[index]
    #         print(row['entrezgene'])
    #
    #
    # yff = pd.DataFrame(data = yf_info_filter, index=gene_ids['symbol']).iloc[:, :15]
    # print(yff)
    # sns.heatmap(data = yff, cmap = 'Reds')
    # filename = 'figure of ds1/heatmap' + cluster_num + '.pdf'
    # # plt.savefig(filename, bbox_inches="tight")
    # plt.show()






