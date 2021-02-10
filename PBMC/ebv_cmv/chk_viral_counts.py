
import os
import time
import datetime
import glob
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scprep
import graphtools as gt
import phate
from scipy import sparse
from scipy.stats import zscore 
from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from scipy.stats import mannwhitneyu, tiecorrect, rankdata
from statsmodels.stats.multitest import multipletests
import warnings
import sys
sys.path.append('/home/ngr4/project/scripts/')
import utils

plt.rc('font', size = 8)
plt.rc('font', family='sans serif')
plt.rcParams['pdf.fonttype']=42
plt.rcParams['ps.fonttype']=42
plt.rcParams['legend.frameon']=False
sns.set_style("ticks")

dfp = '/ycga-gpfs/sequencers/pacbio/gw92/10x/data_sent/cl934/073120/*/filtered_feature_bc_matrix/'
data_files = glob.glob(dfp)

pfp = '/home/ngr4/project/collabs/carrie_lucas/results'
pdfp = '/home/ngr4/project/collabs/carrie_lucas/data/processed/'

data_files = glob.glob(dfp)

if True:
    # first load
    adatas = {}
    for i, file in enumerate(data_files):
        if i==0:
            adata = sc.read_10x_mtx(file)
            batch_key = file.split('/filtered_')[0].split('/')[-1].split('_cite_EBV')[0]
            adata.var_names_make_unique()
        else:
            adatas[file.split('/filtered_')[0].split('/')[-1].split('_cite_EBV')[0]] = sc.read_10x_mtx(file)
            adatas[file.split('/filtered_')[0].split('/')[-1].split('_cite_EBV')[0]].var_names_make_unique()

    adata = adata.concatenate(*adatas.values(), batch_categories=[batch_key]+list(adatas.keys()))
    del adatas
    
ebv_genes = adata.var_names.to_list()[-85:]

adata.obs['comb_ebv_raw_counts'] = adata[:, ebv_genes].X.todense().sum(axis=1)

print(adata.obs.groupby('batch').sum())


sc.pp.calculate_qc_metrics(adata,inplace=True)
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['pmito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
print('Ncells=%d have >10percent mt expression' % np.sum(adata.obs['pmito']>0.1))
print('Ncells=%d have <200 genes expressed' % np.sum(adata.obs['n_genes_by_counts']<200))

print(adata.obs.groupby('batch').mean())