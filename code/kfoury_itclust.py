
import tensorflow as tf
import ItClust as ic
import scanpy.api as sc
import os
from numpy.random import seed
from tensorflow import set_random_seed
import pandas as pd
import numpy as np
import warnings
import anndata as ad
#os.environ["CUDA_VISIBLE_DEVICES"]="1"
warnings.filterwarnings("ignore")
#import sys
#!{sys.executable} -m pip install 'scanpy==1.4.4.post1'
#Set seeds
seed(20180806)
np.random.seed(10)
set_random_seed(20180806) # on GPU may be some other default



###
metadata = pd.read_csv('data/kfoury/GSE143791_cell.annotation.human.csv', sep=',', header=0)
metadata['cells'].value_counts()
##

##
kfoury_dir = 'data/kfoury/GSE143791_RAW'
filenames = [sample for sample in os.listdir(kfoury_dir) if sample.endswith('.csv.gz')]
filenames_bmet = [sample for sample in filenames if 'BMET' in sample]
filenames_bmet = [sample for sample in filenames if 'Tumor' in sample]

paths = [os.path.join(kfoury_dir, file) for file in filenames_bmet]
adatas = [sc.read_csv(filename).transpose() for filename in paths]

for i in adatas:
    i.var_names_make_unique()
    i.var['gene'] = i.var_names

sc.pl.violin(adatas[8], ['POSTN'])

# reset index
for i in adatas:
    i.var.reset_index(inplace=True, drop=True)


adata = ad.AnnData.concatenate(adatas, join = 'inner', batch_key = 'sample')
adata

adata.var_names = adata.var['gene']
adata.var_names

sc.pl.violin(adata, ['APC'])

matching = [gene for gene in adata.var['gene'] if "MKI67" in gene]
matching


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


#####################################################
## Preprocess

adata_BoneMet.var['mt'] = adata_BoneMet.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
adata_BoneMet.var['mt'].value_counts()
sc.pp.calculate_qc_metrics(adata_BoneMet, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_cells(adata_BoneMet, min_genes=200)
#sc.pp.filter_genes(adata_BoneMet, min_cells=3)

matching = [gene for gene in adata_BoneMet.var['gene'] if "MKI67" in gene]
matching


sc.pl.scatter(adata_BoneMet, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata_BoneMet, x="total_counts", y="n_genes_by_counts")

adata_BoneMet = adata_BoneMet[adata_BoneMet.obs.n_genes_by_counts < 4000, :]
adata_BoneMet = adata_BoneMet[adata_BoneMet.obs.pct_counts_mt < 15, :]

sc.pp.normalize_total(adata_BoneMet, target_sum=1e4)
adata_BoneMet.layers["counts"] = adata_BoneMet.X.copy()

sc.pp.log1p(adata_BoneMet)
adata_BoneMet.raw = adata_BoneMet

# ComBat batch correction
sc.pp.combat(adata_BoneMet, key='sample')

sc.pp.highly_variable_genes(
    adata_BoneMet,
    flavor="seurat_v3",
    n_top_genes=4000,
    layer="counts",
    subset=True,
)

# Regress out effects of total counts per cell
sc.pp.regress_out(adata_BoneMet, ['total_counts', 'pct_counts_mt'])

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


########################
# check if bone mets cells are present in the metadata file
adata_BoneMet.obs['cellID'].isin(metadata["barcode"]).value_counts()

# subset the metadata cells to those present in adata
metadata_BoneMets = pd.merge(adata_BoneMet.obs, metadata, left_on='cellID', right_on='barcode')

# reinsert into adata_BoneMets
adata_BoneMet.obs = adata_BoneMet.obs.merge(metadata, left_on='cellID', right_on='barcode', how='inner')
adata_BoneMet.obs['cells'].value_counts()

# DE genes
sc.tl.rank_genes_groups(adata_BoneMet, 'cells', method='t-test')
sc.pl.rank_genes_groups(adata_BoneMet, n_genes=25, sharey=False)

# plot all celltypes
sc.pl.umap(adata_BoneMet, color = ['cells'], save='_kfoury_boneMets_celltypes.png')

# plot stromal cells
sc.pl.umap(adata_BoneMet, color = ['cells'], groups = ['Progenitors', 'Osteoblasts', 'Osteoclasts', 'Endothelial', 'Pericytes'], save='_kfoury_boneMets_stroma.png')


##################################################################################
## subset to the stromal cells ???
adata_BoneMet.obs_names = adata_BoneMet.obs['cellID']
adata_BoneMet_stroma = adata_BoneMet[adata_BoneMet.obs['cells'].isin(['Osteoblasts', 'Osteoclasts', 'Endothelial', 'Pericytes'],)]


##################################################################################
## subset to the stromal cells ???
adata_BoneMet.obs_names = adata_BoneMet.obs['cellID']
adata_BoneMet_stroma = adata_BoneMet[adata_BoneMet.obs['cells'].isin(['Osteoblasts', 'Osteoclasts', 'Endothelial', 'Pericytes'],)]







