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
metadata['cells'].value_counts()
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

matching = [gene for gene in adata.var['gene'] if "MKI67" in gene]
matching

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
adata_BoneMet.var['mt'].value_counts()
sc.pp.calculate_qc_metrics(adata_BoneMet, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


sc.pp.filter_cells(adata_BoneMet, min_genes=200)
#sc.pp.filter_genes(adata_BoneMet, min_cells=3)

matching = [gene for gene in adata_BoneMet.var['gene'] if "MKI67" in gene]
matching


sc.pl.scatter(adata_BoneMet, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata_BoneMet, x="total_counts", y="n_genes_by_counts")

adata_BoneMet = adata_BoneMet[adata_BoneMet.obs.n_genes_by_counts < 4000, :]

sc.pp.normalize_total(adata_BoneMet, target_sum=1e4)
sc.pp.log1p(adata_BoneMet)
sc.pp.highly_variable_genes(adata_BoneMet, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_BoneMet.raw = adata_BoneMet
adata_BoneMet.layers['counts'] = adata_BoneMet.X

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

# read
adata_BoneMet = sc.read_h5ad('data/kfoury/adata_BoneMet_proc.h5ad')

########################
# check if bone mets cells are present in the metadata file
adata_BoneMet.obs['cellID'].isin(metadata["barcode"]).value_counts()

# subset the metadata cells to those present in adata
metadata_BoneMets = pd.merge(adata_BoneMet.obs, metadata, left_on='cellID', right_on='barcode')

# reinsert into adata_BoneMets
adata_BoneMet.obs = adata_BoneMet.obs.merge(metadata, left_on='cellID', right_on='barcode', how='inner')
adata_BoneMet.obs['cells'].value_counts()

# plot all celltypes
sc.pl.umap(adata_BoneMet, color = ['cells'], save='_kfoury_boneMets_celltypes.png')

# plot stromal cells
sc.pl.umap(adata_BoneMet, color = ['cells'], groups = ['Progenitors', 'Osteoblasts', 'Endothelial', 'Pericytes'], save='_kfoury_boneMets_stroma.png')

# Write
adata_BoneMet.write('data/kfoury/adata_BoneMet_proc.h5ad')

# read
adata_BoneMet = sc.read('data/kfoury/adata_BoneMet_proc.h5ad')

##################################################################################
## subset to the stromal cells ???
adata_BoneMet.obs_names = adata_BoneMet.obs['cellID']
adata_BoneMet_stroma = adata_BoneMet[adata_BoneMet.obs['cells'].isin(['Progenitors', 'Osteoblasts', 'Endothelial', 'Pericytes'],)]

# Write
adata_BoneMet_stroma.write('data/kfoury/adata_BoneMet_stroma.h5ad')

# read
adata_BoneMet_stroma = sc.read('data/kfoury/adata_BoneMet_stroma.h5ad')

####
adata_BoneMet_stroma.obs['cells'].value_counts()
sc.pl.umap(adata_BoneMet_stroma, color = ['cells'], save='_kfoury_boneMets_stroma.png')


##################################################################################
##################################################################################
## Mapping

# read the human bone mets stroma data
adata_BoneMet_stroma = sc.read('data/kfoury/adata_BoneMet_stroma.h5ad')

matching = [gene for gene in adata_BoneMet_stroma.var['gene'] if "APC" in gene]
matching

# read the mouse data
adata_mouse = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)
adata_mouse.obs['cluster'] = adata_mouse.obs['cluster'].astype('str')
adata_mouse.obs['cluster'].value_counts()

# filter the na
adata_mouse = adata_mouse[adata_mouse.obs['cluster'] != 'nan', :]
adata_mouse.obs['cluster'].value_counts()

# change the human var_names to match mouse gene symbols
adata_BoneMet_stroma.var_names = [gene.title() for gene in adata_BoneMet_stroma.var_names]

# subset also the adata_BoneMet_stroma.raw
tempAdata = adata_BoneMet_stroma.raw.to_adata()
tempAdata.var_names = [gene.title() for gene in tempAdata.var_names]
adata_BoneMet_stroma.raw = tempAdata

adata_BoneMet_stroma.var_names
adata_BoneMet_stroma.raw.var_names

# find common genes
var_names = adata_mouse.var_names.intersection(adata_BoneMet_stroma.var_names)
len(var_names)

# subset
adata_mouse = adata_mouse[:, var_names]
adata_BoneMet_stroma = adata_BoneMet_stroma[:, var_names]

# recompute the neighbors for the adata_mouse (ref dataset)
sc.pp.pca(adata_mouse)
sc.pp.neighbors(adata_mouse)
sc.tl.umap(adata_mouse)
sc.pl.umap(adata_mouse, color='cluster')

# mapping the clusters from mouse to human using ingest
sc.tl.ingest(adata_BoneMet_stroma, adata_mouse, obs='cluster')

adata_BoneMet_stroma.obs['cluster'].value_counts()
adata_BoneMet_stroma.uns['cluster_colors'] = adata_mouse.uns['cluster_colors']  # fix colors


# remove 0,5,2,4 since they have 1/0 samples
adata_BoneMet_stroma = adata_BoneMet_stroma[~adata_BoneMet_stroma.obs['cluster'].isin(['0', '2','4','5'])]
adata_BoneMet_stroma.obs['cluster'].value_counts()

sc.pl.umap(adata_BoneMet_stroma, color='cluster', save = '_kfoury_projectedClusters.png')


# plot the inferred clusters in human
#sc.pp.pca(adata_BoneMet_stroma)
#sc.pp.neighbors(adata_BoneMet_stroma)
#sc.tl.umap(adata_BoneMet_stroma)
#sc.pl.umap(adata_BoneMet_stroma, color='cluster')

# human
sc.tl.rank_genes_groups(adata_BoneMet_stroma, 'cluster', pts=True, use_raw = False, reference='1', method = 't-test_overestim_var')
sc.pl.rank_genes_groups(adata_BoneMet_stroma, n_genes=25, sharey=True)

####################################################################
## get cluster markers

# mouse
sc.tl.rank_genes_groups(adata_mouse, 'cluster', pts=True, use_raw = False, reference='1',  method = 't-test_overestim_var')

markers_mouse_c0 = sc.get.rank_genes_groups_df(adata_mouse, group = '0')
markers_mouse_c1 = sc.get.rank_genes_groups_df(adata_mouse, group = '1')
markers_mouse_c2 = sc.get.rank_genes_groups_df(adata_mouse, group = '2')
markers_mouse_c3 = sc.get.rank_genes_groups_df(adata_mouse, group = '3')
markers_mouse_c4 = sc.get.rank_genes_groups_df(adata_mouse, group = '4')
markers_mouse_c5 = sc.get.rank_genes_groups_df(adata_mouse, group = '5')
markers_mouse_c6 = sc.get.rank_genes_groups_df(adata_mouse, group = '6')
markers_mouse_c7 = sc.get.rank_genes_groups_df(adata_mouse, group = '7')

########
# human
#sc.tl.rank_genes_groups(adata_BoneMet_stroma, 'cluster', pts=True, use_raw = False)

markers_human_c1 = sc.get.rank_genes_groups_df(adata_BoneMet_stroma, group = '1')
markers_human_c3 = sc.get.rank_genes_groups_df(adata_BoneMet_stroma, group = '3')
markers_human_c6 = sc.get.rank_genes_groups_df(adata_BoneMet_stroma, group = '6')
markers_human_c7 = sc.get.rank_genes_groups_df(adata_BoneMet_stroma, group = '7')

# fix gene symbols
#markers_human_c1['names'] = markers_human_c1['names'].str.title()
#markers_human_c3['names'] = markers_human_c3['names'].str.title()
#markers_human_c6['names'] = markers_human_c6['names'].str.title()
#markers_human_c7['names'] = markers_human_c7['names'].str.title()

####################
## get the top markers

# mouse
top100_mouse_c0 = markers_mouse_c0[0:100]
top100_mouse_c1 = markers_mouse_c1[0:100]
top100_mouse_c2 = markers_mouse_c2[0:100]
top100_mouse_c3 = markers_mouse_c3[0:100]
top100_mouse_c4 = markers_mouse_c4[0:100]
top100_mouse_c5 = markers_mouse_c5[0:100]
top100_mouse_c6 = markers_mouse_c6[0:100]
top100_mouse_c7 = markers_mouse_c7[0:100]

# human
top100_human_c1 = markers_human_c1[0:100]
top100_human_c3 = markers_human_c3[0:100]
top100_human_c6 = markers_human_c6[0:100]
top100_human_c7 = markers_human_c7[0:100]

#######################
## get the genes in common
intersection_c1 = pd.merge(top100_mouse_c1, top100_human_c1, how='inner', on=['names'])
intersection_c3 = pd.merge(top100_mouse_c3, top100_human_c3, how='inner', on=['names'])
intersection_c6 = pd.merge(top100_mouse_c6, top100_human_c6, how='inner', on=['names'])
intersection_c7 = pd.merge(top100_mouse_c7, top100_human_c7, how='inner', on=['names'])





###############################
# re-cap the gene symbols
adata_BoneMet_stroma.var_names = [gene.upper() for gene in adata_BoneMet_stroma.var_names]
# subset also the adata_human.raw
tempAdata = adata_BoneMet_stroma.raw.to_adata()
tempAdata.var_names = [gene.upper() for gene in tempAdata.var_names]
adata_BoneMet_stroma.raw = tempAdata

dp2 = sc.pl.DotPlot(adata_BoneMet_stroma,
                                var_names = ['SFRP2', 'WNT5A', 'LGR5', 'APC', 'WNT4',
                                             'WNT6', 'NOTUM', 'WIF1', 'NKD1', 'FZD1',
                                             'WNT2', 'WNT10A', 'DKK2', 'RORB', 'CXXC4',
                                             'NFAT5', #'APOE',
                                             'DACT1', 'CTNNB1', 'LEF1',
                                             'TCF4', 'MYC', 'MKI67', 'H2AFX', 'TOP2A',
                                             'CCNB1', 'CCNB2', 'STMN1', 'PTN', 'MDK',
                                             'TUBB3', 'MRC2', 'FN1', 'TNC', 'COL12A1',
                                             'COL14A1', 'COL16A1', 'MMP19', 'CTHRC1',
                                             #'WISP1',
                                             'FZD1', 'FZD2', 'SFRP4', 'BMP1',
                                             'TLE3', #'TGFB1',
                                             'POSTN'],
                                categories_order = ['1','3','6','7'],
                                groupby='cluster', cmap = 'Reds', use_raw=True,
                                )
dp2.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/kfoury_dotplot.png')


########################
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







