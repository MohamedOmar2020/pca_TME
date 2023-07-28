
import sys

sys.path.insert(0, "/Users/mohamedomar/Documents/Research/Projects/PCa_TME")

from pathlib import Path

import scanpy as sc
import scvi
from matplotlib.pyplot import ion
import matplotlib as plt
from scutils.figures.base import basics
from scutils.figures.prostate import annotate_cell_types_prostate
from scutils.qc import PreprocessRNA
import anndata as ad
import squidpy as sq
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import SpatialDE
import seaborn as sns

# set global paths
basepath = "/Users/mohamedomar/Documents/Research/Projects/PCa_TME"
outspath = basepath + "coculture/outs"


#############################################
# set figure parameters
sc.settings.figdir = 'figures/coculture2'
sc.set_figure_params(dpi_save = 300, transparent = False, fontsize =7, format='tiff')
plt.rcParams["font.family"] = "arial"
plt.rcParams['font.size'] = 9
plt.rcParams['font.style'] = 'italic'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['legend.fontsize'] = 9
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 9
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9

#############################################################################################
############
# function to load and merge the adata objects

def load():
    ion()
    ## make sure to re-check batches

    # fibro_only
    fibro_only = sc.read_h5ad(basepath + "/coculture/kb_counts/07_07/fibro_only/adata.h5ad")
    #fibro_only.var_names = fibro_only.var['gene_name']
    fibro_only.var_names_make_unique()
    fibro_only.obs["batch"] = 0
    fibro_only.obs["cells"] = 'fibro_only'

    ##########################################
    ## PRN_fibro

    PRN_fibro = sc.read_h5ad(basepath + "/coculture/kb_counts/07_07/PRN_fibro/adata.h5ad")
    #PRN_fibro.var_names = PRN_fibro.var['gene_name']
    PRN_fibro.var_names_make_unique()
    PRN_fibro.obs["batch"] = 0
    PRN_fibro.obs["cells"] = 'PRN_fibro'

    ##########################################
    ## TRG_fibro

    TRG_fibro = sc.read_h5ad(basepath + "/coculture/kb_counts/07_07/TRG_fibro/adata.h5ad")
    #TRG_fibro.var_names = TRG_fibro.var['gene_name']
    TRG_fibro.var_names_make_unique()
    TRG_fibro.obs["batch"] = 0
    TRG_fibro.obs["cells"] = 'TRG_fibro'

    ############################################
    # Merge all adatas
    adata = fibro_only.concatenate(
             TRG_fibro, PRN_fibro, join="outer",
    )
    adata.layers["counts"] = adata.X.copy()
    return adata




def preprocess(adata):
    ion()
    sc.set_figure_params(dpi=400, dpi_save=400, figsize=(4, 4))


    adata.var_names = adata.var['gene_name'].astype(str)
    adata.var_names_make_unique()

    # get mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save="qc1.png",
    )
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save="_qc2.png")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save="_qc3.png")
    """
    min_counts=400
    max_counts=500000
    min_genes=300
    max_genes=8000
    max_percent_mito=0.15
    """
    pre = PreprocessRNA(subset_hvg=False, )
    adata = pre.preprocess(adata)
    print('finished processing adata: ', adata)
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=4000,
        layer="counts",
        #batch_key="batch",
        subset=True,
    )
    print('finished detecting HVGs: ', adata)

    #scvi.data.setup_anndata(adata, layer="counts", batch_key="batch")
    #scvi.model.SCVI.setup_anndata(adata, layer='counts', batch_key='batch')
    #print('finished setting up anndata with scvi: ', adata)
    #vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
    #print('training: ', vae)
    #vae.train()
    #print('finished training ', adata)
    #adata.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.pca(adata, n_comps=10)
    print('computing neighbors')
    sc.pp.neighbors(adata)
    print('computing UMAP')
    sc.tl.umap(adata)
    print('computing leiden clusters')
    sc.tl.leiden(adata, resolution=1.0)
    basics(adata)
    print('saving adata: ', adata)
    adata.write_h5ad("coculture/outs/07_07/adata_fibroblasts_07_07_processed.h5ad")
    return adata


##########################################################################
## Load and merge all the adata objects
##########################################################################
adata_coculture_07_07 = load()
adata_coculture_07_07.var_names = adata_coculture_07_07.var['gene_name']


preprocess(adata_coculture_07_07)

##########################################################################
# load the processed adata
##########################################################################
adata_coculture_07_07 = sc.read_h5ad("coculture/outs/07_07/adata_fibroblasts_07_07_processed.h5ad")

adata_coculture_07_07.var_names_make_unique()
adata_coculture_07_07.uns['log1p']["base"] = None

# make var names unique in .raw
adata_coculture_07_07_raw = adata_coculture_07_07.raw.to_adata()
adata_coculture_07_07_raw.var_names = adata_coculture_07_07_raw.var_names.astype('str')
adata_coculture_07_07_raw.var_names_make_unique()
adata_coculture_07_07.raw = adata_coculture_07_07_raw


##########################
# re-cluster
#########################
sc.tl.leiden(adata_coculture_07_07, resolution=1.0)
sc.pl.umap(adata_coculture_07_07, color = ['leiden'], save = '_leiden')

adata_coculture_07_07.obs['leiden'].value_counts()

##################################################
## marker genes
##################################################
# fibro clusters per genotype
adata_coculture_07_07.obs['cells'].value_counts()

# pd.crosstab(adata_coculture_07_07.obs['leiden'], adata_coculture_07_07.obs['cells'])
#
# # marker genes per genotype
# sc.tl.rank_genes_groups(adata_coculture_07_07, "cells", method="t-test", pts=True, use_raw = True)
# sc.pl.rank_genes_groups(adata_coculture_07_07, groups=None, n_genes=25, groupby="cells", sharey=False, save='_topGenes_fibro_07_07')
#
# markers_fibro_genotype = sc.get.rank_genes_groups_df(adata_coculture_07_07, group=None)
# markers_fibro_genotype.to_csv('coculture/outs/07_07/markers_per_genotype.csv')

# marker genes per cluster
sc.tl.rank_genes_groups(adata_coculture_07_07, "leiden", method="t-test", pts=True, use_raw = True)
sc.pl.rank_genes_groups(adata_coculture_07_07, groups=None, n_genes=25, groupby="cells", sharey=False, save='_topGenes_fibro_07_07_leiden')

sc.tl.rank_genes_groups(adata_coculture_07_07, "cells", method="t-test", pts=True, use_raw = True, reference='fibro_only')
sc.pl.rank_genes_groups(adata_coculture_07_07, groups=None, n_genes=25, groupby="cells", sharey=False, save='_topGenes_fibro_07_07_genotype')


markers_fibro_cluster = sc.get.rank_genes_groups_df(adata_coculture_07_07, group=None)
markers_fibro_cluster.to_csv('coculture/outs/07_07/markers_per_cluster.csv')

##########################################################################



####################################################################################################################
adata_mouse = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)
adata_mouse.obs['cluster'] = adata_mouse.obs['cluster'].astype('str')
adata_mouse.obs['cluster'].value_counts()

adata_mouse.obs['condition'].value_counts()

#############################
# filter the na
adata_mouse = adata_mouse[adata_mouse.obs['cluster'] != 'nan', :]
adata_mouse.obs['cluster'].value_counts()

#############################
# add c
adata_mouse.obs['cluster'] = 'c' + adata_mouse.obs['cluster']
adata_mouse.obs['cluster'].value_counts()

# get raw data
adata_mouse_raw = adata_mouse.raw.to_adata()
adata_mouse_raw.var_names_make_unique()

##########################################################
# fix GEMMs names
##########################################################
adata_mouse.obs['key_new'] = adata_mouse.obs['key']
adata_mouse.obs['key_new'].value_counts()
adata_mouse.obs['key_new'].replace('terg', 'T-ERG', inplace=True)
adata_mouse.obs['key_new'].replace('himyc', 'Hi-MYC', inplace=True)
adata_mouse.obs['key_new'].replace('fvbn', 'FVBN', inplace=True)
adata_mouse.obs['key_new'].replace('pten', 'NP', inplace=True)
adata_mouse.obs['key_new'].replace('mycn', 'PRN', inplace=True)
adata_mouse.obs['key_new'].replace('129b6', 'B6.129', inplace=True)
adata_mouse.obs['key_new'].replace('b6', 'B6', inplace=True)
adata_mouse.obs['key_new'].replace('129b6_pten', 'WT for PN', inplace=True)
adata_mouse.obs['key_new'].replace('129b6_mycn', 'WT for PRN', inplace=True)




##########################################################################
# load the processed co-culture adata
##########################################################################
adata_coculture_07_07 = sc.read_h5ad("coculture/outs/07_07//adata_fibroblasts_07_07_processed.h5ad")

adata_coculture_07_07.var_names_make_unique()
adata_coculture_07_07.uns['log1p']["base"] = None

# make var names unique in .raw
adata_coculture_07_07_raw = adata_coculture_07_07.raw.to_adata()
adata_coculture_07_07_raw.var_names = adata_coculture_07_07_raw.var_names.astype('str')
adata_coculture_07_07_raw.var_names_make_unique()
adata_coculture_07_07.raw = adata_coculture_07_07_raw

##########################################################################
# separate into TRG and PRN
##########################################################################
adata_coculture_TRG = adata_coculture_07_07[adata_coculture_07_07.obs['cells'].isin(['TRG_fibro'])]
adata_coculture_PRN = adata_coculture_07_07[adata_coculture_07_07.obs['cells'].isin(['PRN_fibro'])]
adata_coculture_normal = adata_coculture_07_07[adata_coculture_07_07.obs['cells'].isin(['fibro_only'])]

##########################################################################
# find common genes
##########################################################################
var_names = adata_mouse.var_names.intersection(adata_coculture_07_07.var_names)
len(var_names)


# subset
adata_mouse = adata_mouse[:, var_names]
adata_coculture_07_07 = adata_coculture_07_07[:, var_names]
adata_coculture_TRG = adata_coculture_TRG[:, var_names]
adata_coculture_PRN = adata_coculture_PRN[:, var_names]
adata_coculture_normal = adata_coculture_normal[:, var_names]

##########################################################################
# recompute the neighbors for the adata_mouse (ref dataset)
##########################################################################
sc.pp.pca(adata_mouse)
sc.pp.neighbors(adata_mouse)
sc.tl.umap(adata_mouse)
sc.pl.umap(adata_mouse, color = 'cluster')
adata_coculture_07_07.uns['cluster_colors'] = adata_mouse.uns['cluster_colors']  # fix colors


##########################################################################
# Ingest all
##########################################################################
sc.tl.ingest(adata_coculture_07_07, adata_mouse, obs='cluster')
adata_coculture_07_07.obs.cluster.value_counts()

# remove rare cells
adata_coculture_07_07 = adata_coculture_07_07[~adata_coculture_07_07.obs['cluster'].isin(['c0', 'c2', 'c4', 'c5'])]

##########################################################################
# Ingest using adata_mouse_terg as reference and fibro_TRG as the target
##########################################################################

# get T-erg model
adata_mouse_terg = adata_mouse[adata_mouse.obs.key_new=='T-ERG']

#recompute the neighbors for the adata_mouse_terg (ref dataset)
sc.pp.pca(adata_mouse_terg)
sc.pp.neighbors(adata_mouse_terg)
sc.tl.umap(adata_mouse_terg)

#Ingest using adata_mouse_terg as reference and adata_human_ERGpos as the target
sc.tl.ingest(adata_coculture_TRG, adata_mouse_terg, obs='cluster')

#cell frequency in each cluster
adata_coculture_TRG.obs.cluster.value_counts()

# remove empty clusters
adata_coculture_TRG.obs['cluster'].value_counts()
adata_coculture_TRG = adata_coculture_TRG[~adata_coculture_TRG.obs['cluster'].isin(['c0', 'c2', 'c4', 'c6', 'c7'])]


# ##########################################################################
# # Ingest using adata_mouse_PRN as reference and fibro_PRN as the target
# ##########################################################################
#
# Get the PRN model
adata_mouse_PRN = adata_mouse[adata_mouse.obs.key_new=='PRN']

# recompute the neighbors for the adata_mouse_PRN (ref dataset)
sc.pp.pca(adata_mouse_PRN)
sc.pp.neighbors(adata_mouse_PRN)
sc.tl.umap(adata_mouse_PRN)

# Ingest using adata_mouse_terg as reference and adata_human_ERGpos as the target
sc.tl.ingest(adata_coculture_PRN, adata_mouse_PRN, obs='cluster')

# cell frequency in each cluster
adata_coculture_PRN.obs.cluster.value_counts()

# # remove empty clusters (cells < 5)
adata_coculture_PRN = adata_coculture_PRN[~adata_coculture_PRN.obs['cluster'].isin(['c0', 'c1', 'c2', 'c3', 'c4', 'c5'])]

##########################################################################
# Ingest using adata_mouse_normal as reference and fibro_normal as the target
##########################################################################

# get all wt models
adata_mouse_normal = adata_mouse[adata_mouse.obs.condition=='wildtype']

#recompute the neighbors for the adata_mouse_terg (ref dataset)
sc.pp.pca(adata_mouse_normal)
sc.pp.neighbors(adata_mouse_normal)
sc.tl.umap(adata_mouse_normal)

#Ingest using adata_mouse_terg as reference and adata_human_ERGpos as the target
sc.tl.ingest(adata_coculture_normal, adata_mouse_normal, obs='cluster')

#cell frequency in each cluster
adata_coculture_normal.obs.cluster.value_counts()

# remove empty clusters
adata_coculture_normal.obs['cluster'].value_counts()
adata_coculture_normal = adata_coculture_normal[adata_coculture_normal.obs['cluster'].isin(['c1', 'c3', 'c4'])]


# ##########################################################################
# # merge the ingested co-culture data
# ##########################################################################
adata_coculture_list = [adata_coculture_normal, adata_coculture_TRG, adata_coculture_PRN]
adata_coculture_sep = ad.concat(adata_coculture_list, label = 'model', join = 'outer')

adata_coculture_sep.obs['model'].replace('0', 'control', inplace=True)
adata_coculture_sep.obs['model'].replace('1', 'TRG', inplace=True)
adata_coculture_sep.obs['model'].replace('2', 'PRN', inplace=True)

adata_coculture_sep.obs['model'].value_counts()

# ##########################################################################
# # save
# ##########################################################################
adata_coculture_ingested_together = adata_coculture_07_07.copy()
adata_coculture_ingested_sep = adata_coculture_sep.copy()


adata_coculture_ingested_together.write("coculture/outs/07_07/adata_coculture_ingested_together.h5ad")
adata_coculture_ingested_sep.write("coculture/outs/07_07/adata_coculture_ingested_sep.h5ad")


# read
adata_coculture_ingested_together = sc.read_h5ad("coculture/outs/07_07/adata_coculture_ingested_together.h5ad")
adata_coculture_ingested_sep = sc.read_h5ad("coculture/outs/07_07/adata_coculture_ingested_sep.h5ad")

# ####################################################################################################################################################
#
# ##########################################################################
# # Downstream
# ##########################################################################
#
# # remove empty clusters
# adata_coculture.obs['cluster'].value_counts()
# adata_coculture = adata_coculture[~adata_coculture.obs['cluster'].isin(['c0', 'c5', 'c2', 'c4'])]

# ######################################
# ## confusion matrix
#
# # first combine adata mouse with adata human
# #adata_CAFs.obs['cluster_predicted'] = adata_CAFs.obs['cluster']
# #del adata_CAFs.obs['cluster']
# adata_coculture.obs['cluster'] = 'predicted_' + adata_coculture.obs['cluster'].astype(str)
# adata_coculture.obs['cluster'].value_counts()
#
# # same for mouse
# adata_mouse.obs['cluster'] = 'original_' + adata_mouse.obs['cluster'].astype(str)
# adata_mouse.obs['cluster'].value_counts()
#
#
# # make var and obs names unique
# adata_coculture.obs_names_make_unique()
# adata_mouse.obs_names_make_unique()
#
# adata_coculture.var_names_make_unique()
# adata_mouse.var_names_make_unique()
#
# tempAdata = adata_coculture.raw.to_adata()
# tempAdata.var_names_make_unique()
# adata_coculture.raw = tempAdata
#
# tempAdata_mouse = adata_mouse.raw.to_adata()
# tempAdata_mouse.var_names_make_unique()
# adata_mouse.raw = tempAdata_mouse
#
# #######
# # join
# adata_list = [adata_mouse, adata_coculture]
# adata_all = ad.concat(adata_list,  join="outer")
#
# adata_all.obs['cluster'] = adata_all.obs['cluster'].astype('category')
# adata_all.obs['cluster'].value_counts()
#
# sc.pl.correlation_matrix(adata_all, 'cluster', save='_corr.png')



##############################################################################
# figures
##############################################################################
# umap
sc.pl.umap(adata_coculture_ingested_together, color = ['cells', 'cluster'], save = '_coculture_model_cluster_together')
sc.pl.umap(adata_coculture_ingested_sep, color = ['cells', 'cluster'], save = '_coculture_model_cluster_sep')


sc.pl.umap(adata_coculture_ingested_together, color=['cluster', 'Ar', 'Postn'], save = '_coculture_clusters_Ar_Postn_together')
sc.pl.umap(adata_coculture_ingested_sep, color=['cluster', 'Ar', 'Postn'], save = '_coculture_clusters_Ar_Postn_sep')

###################
# marker genes
sc.tl.rank_genes_groups(adata_coculture_ingested_together, 'cluster', pts=True, use_raw = False)
sc.pl.rank_genes_groups(adata_coculture_ingested_together, n_genes=25, sharey=False, save = '_coculture_markers_together.png')

sc.tl.rank_genes_groups(adata_coculture_ingested_sep, 'cluster', pts=True, use_raw = False)
sc.pl.rank_genes_groups(adata_coculture_ingested_sep, n_genes=25, sharey=False, save = '_coculture_markers_sep.png')
markers_fibro_clusters = sc.get.rank_genes_groups_df(adata_coculture_ingested_sep, group=None)
markers_fibro_clusters.to_csv('coculture/outs/07_07/markers_fibro_clusters.csv')

sc.tl.rank_genes_groups(adata_coculture_ingested_sep, 'cells', pts=True, use_raw = False, reference='fibro_only')
sc.pl.rank_genes_groups(adata_coculture_ingested_sep, n_genes=25, sharey=False, save = '_coculture_markers_genotype.png')
markers_fibro_genotypes = sc.get.rank_genes_groups_df(adata_coculture_ingested_sep, group=None)
markers_fibro_genotypes.to_csv('coculture/outs/07_07/markers_fibro_genotypes.csv')

pd.crosstab(adata_coculture_ingested_sep.obs['cluster'], adata_coculture_ingested_sep.obs['cells'])






# dotplot of GEMM-specific clusters
dp2 = sc.pl.DotPlot(adata_coculture_ingested_sep, var_names = ['Ar', 'Sfrp2', 'Wnt5a', 'Lgr5', 'Apc',
                                                          'Wnt4', 'Wnt6', 'Notum', 'Wif1',
                                                          'Nkd1', 'Wnt2', 'Wnt10a',
                                                          'Dkk2', 'Rorb', 'Cxxc4', 'Nfat5',
                                                          'Apoe', 'Dact1', 'Ctnnb1', 'Lef1',
                                                          'Tcf4', 'Myc', 'Mki67', 'H2afx',
                                                          'Top2a', 'Ccnb1', 'Ccnb2', 'Stmn1',
                                                          'Ptn', 'Mdk', 'Tubb3', 'Mrc2', 'Fn1',
                                                          'Tnc', 'Col12a1', 'Col14a1', 'Col16a1',
                                                          'Mmp19', 'Cthrc1', 'Wisp1', 'Fzd1', 'Fzd2',
                                                          'Sfrp4', 'Bmp1', 'Tle3', 'Tgfb1', 'Postn'],
                                #categories_order = ['c0','c1','c2','c3','c4','c5','c6','c7'],
                                groupby='cluster', cmap = 'Reds', figsize=[15, 3],
                                )

dp2.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/coculture2/dotplot_coculture_clusters.tiff')







################################################################
## RF
################################################################

# Get common genes
common_genes = adata_mouse.var_names.intersection(adata_coculture_07_12.var_names)

# Subset the AnnData objects to include only the common genes
adata_mouse_common = adata_mouse[:, common_genes]
adata_coculture_07_12_common = adata_coculture_07_12[:, common_genes]

print(adata_mouse_common.X.shape)
print(adata_coculture_07_12_common.X.shape)


# Then proceed with the RandomForestClassifier
rf = RandomForestClassifier()
rf.fit(adata_mouse_common.X, adata_mouse_common.obs['cluster'])
adata_coculture_07_12_common.obs['cluster'] = rf.predict(adata_coculture_07_12_common.X)


adata_coculture_07_12_common.obs['cluster'].value_counts()

adata_coculture_07_12_common = adata_coculture_07_12_common[~adata_coculture_07_12_common.obs['cluster'].isin(['c4'])]

sc.tl.rank_genes_groups(adata_coculture_07_12_common, 'cluster', pts=True, use_raw = False)
sc.pl.rank_genes_groups(adata_coculture_07_12_common, n_genes=25, sharey=False)

sc.pl.violin(adata_coculture_07_12_common, keys=['Ar', 'Postn'], groupby='cluster')