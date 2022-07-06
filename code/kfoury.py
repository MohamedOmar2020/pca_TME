import os
from pathlib import Path
import scanpy as sc
from matplotlib.pyplot import ion
from scutils.figures.base import basics
from scutils.figures.prostate import annotate_cell_types_prostate
from scutils.qc import PreprocessRNA
import pandas as pd
import scipy.sparse as sp
import loompy as lp
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import anndata as ad
import sc_toolbox as sct

# solid metastatic tissue (Tumor)
# liquid BM at the vertebral level of spinal cord compression (Involved)
# liquid BM from a different vertebral body distant from the tumor site but within the surgical field (Distal)

###
metadata = pd.read_csv('data/kfoury/GSE143791_cell.annotation.human.csv', sep=',', header=0)
##
kfoury_dir = 'data/kfoury/GSE143791_RAW'
filenames = [sample for sample in os.listdir(kfoury_dir) if sample.endswith('.csv.gz')]
paths = [os.path.join(kfoury_dir, file) for file in filenames]
adatas = [sc.read_csv(filename).transpose() for filename in paths]

for i in adatas:
    i.var['gene'] = i.var_names

# reset index
for i in adatas:
    i.var.reset_index(inplace=True, drop=True)


adata = ad.concat(adatas=adatas, join = 'inner', label = 'sample', merge='first')
adata
adata.var_names = adata.var['gene']
adata.var_names

# Write
adata.write('data/kfoury/kfoury_raw.h5ad')

## read the adata
adata = sc.read('data/kfoury/kfoury_raw.h5ad')


#######################################################
# keep only bone mets samples
adata.obs['cellID'] = adata.obs_names

adata.obs['organ'] = adata.obs['cellID'].str.split('_', expand = True)[0]
#data.obs['organ'] = adata.obs['organ'].str.replace('.', '-')
#adata.obs['condition'] = adata.obs['organ'].str.split('-', expand = True)[0]
adata.obs['organ'].value_counts()

adata_BoneMet = adata[adata.obs['organ'].str.startswith('BMET'),:]
adata_BoneMet = adata_BoneMet[adata_BoneMet.obs['organ'].str.endswith('Tumor'),:]

adata_BoneMet.obs['organ'].value_counts()

# Write
adata_BoneMet.write('data/kfoury/adata_BoneMet.h5ad')

## read the adata
adata_BoneMet = sc.read('data/kfoury/adata_BoneMet.h5ad')

#####################################################
## Preprocess

adata_BoneMet.var['mt'] = adata_BoneMet.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata_BoneMet, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


sc.pp.filter_cells(adata_BoneMet, min_genes=200)
sc.pp.filter_genes(adata_BoneMet, min_cells=3)

sc.pl.scatter(adata_BoneMet, x="total_counts", y="pct_counts_mt", save="qc2.png")
sc.pl.scatter(adata_BoneMet, x="total_counts", y="n_genes_by_counts", save="qc3.png")

adata_BoneMet = adata_BoneMet[adata_BoneMet.obs.n_genes_by_counts < 4000, :]

sc.pp.normalize_total(adata_BoneMet, target_sum=1e4)
sc.pp.log1p(adata_BoneMet)
sc.pp.highly_variable_genes(adata_BoneMet, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_BoneMet.raw = adata_BoneMet

# Regress out effects of total counts per cell
sc.pp.regress_out(adata_BoneMet, ['total_counts'])

# scale the data to unit variance.
sc.pp.scale(adata_BoneMet, max_value=10)

# pca
sc.tl.pca(adata_BoneMet, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_BoneMet, log=True)

# Computing the neighborhood graph
sc.pp.neighbors(adata_BoneMet)

# umap
sc.tl.umap(adata_BoneMet)
sc.pl.umap(adata_BoneMet, color = ['sample'])

# clustering
sc.tl.leiden(adata_BoneMet, resolution=0.3)
adata_BoneMet.obs['leiden'].value_counts()
sc.pl.umap(adata_BoneMet, color = ['leiden'])

# DE genes
sc.tl.rank_genes_groups(adata_BoneMet, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata_BoneMet, n_genes=25, sharey=False)

# Write
adata_BoneMet.write('data/kfoury/adata_BoneMet_proc.h5ad')

########################
adata_BoneMet.obs_names
adata_BoneMet.obs['sample'].value_counts()

#adata.obs['cellID'] = adata.obs_names

#adata.obs['cellID'] = adata.obs['cellID'].str.replace('.', '-')
#adata.obs['organ'] = adata.obs['cellID'].str.split('-', expand = True)[0]
#adata.obs['condition'] = adata.obs['organ'].str.split('-', expand = True)[1]
#adata.obs['organ'].value_counts()
#adata.obs['condition'].value_counts()

########################################
# subset to cells present in metadata
adata = adata[metadata['barcode'],:]
adata

if (adata.obs_names == metadata['barcode']).all():
    print("cells match")

adata.obs['barcode'] = adata.obs_names

# add cell type column to adata
adata.obs = adata.obs.merge(metadata, left_on='barcode', right_on='barcode')

adata.obs['cells'].value_counts()

# Write
adata.write('data/kfoury/kfoury_fil_annot.h5ad')

######################################################################################
## read the adata
adata = sc.read('data/kfoury/kfoury_fil_annot.h5ad')







