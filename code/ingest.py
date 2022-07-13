
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

sc.settings.figdir = 'figures'
sc.set_figure_params(dpi_save = 300)
#plt.rcParams.update({'xtick.labelsize' : '50'})
#plt.rcParams.update({'ytick.labelsize' : '50'})

adata_mouse = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)
adata_mouse.obs['cluster'] = adata_mouse.obs['cluster'].astype('str')
adata_mouse.obs['cluster'].value_counts()

# filter the na
adata_mouse = adata_mouse[adata_mouse.obs['cluster'] != 'nan', :]
adata_mouse.obs['cluster'].value_counts()

# Process the human data using ingest
# load the human data
adata_human = sc.read_h5ad('human/h5ads/erg_fibroblasts_scvi_v6_regulons.h5ad', chunk_size=100000)

adata_human.obs.erg.value_counts()

# change the var_names to match mouse gene symbols
adata_human.var_names = [gene.title() for gene in adata_human.var_names]

# subset also the adata_human.raw
tempAdata = adata_human.raw.to_adata()
tempAdata.var_names = [gene.title() for gene in tempAdata.var_names]
adata_human.raw = tempAdata

adata_human.var_names
adata_human.raw.var_names

# find common genes
var_names = adata_mouse.var_names.intersection(adata_human.var_names)
len(var_names)

# subset
adata_mouse = adata_mouse[:, var_names]
adata_human = adata_human[:, var_names]

# recompute the neighbors for the adata_mouse (ref dataset)
sc.pp.pca(adata_mouse)
sc.pp.neighbors(adata_mouse)
sc.tl.umap(adata_mouse)
sc.pl.umap(adata_mouse, color='cluster', save = 'mouse_umap.png')

# mapping the clusters from mouse to human using ingest
sc.tl.ingest(adata_human, adata_mouse, obs='cluster')

adata_human

# divide the human data into ERG pos and neg
adata_human_ERGpos = adata_human[adata_human.obs.erg=='positive']
adata_human_ERGneg = adata_human[adata_human.obs.erg=='negative']

# Number of cells in each cluster: mouse vs human
print(adata_mouse.obs['cluster'].value_counts())
print(adata_human.obs['cluster'].value_counts())
print(adata_human_ERGpos.obs['cluster'].value_counts())
print(adata_human_ERGneg.obs['cluster'].value_counts())

adata_human.uns['cluster_colors'] = adata_mouse.uns['cluster_colors']  # fix colors

# Ingest using adata_mouse_terg as reference and adata_human_ERGpos as the target
# mouse T-erg model
print(adata_mouse.obs.key.value_counts())
adata_mouse_terg = adata_mouse[adata_mouse.obs.key=='terg']
adata_mouse_nonterg = adata_mouse[adata_mouse.obs.key!='terg']
print(adata_mouse_terg.obs.cluster.value_counts())
print(adata_mouse_nonterg.obs.cluster.value_counts())

# recompute the neighbors for the adata_mouse_terg (ref dataset)
sc.pp.pca(adata_mouse_terg)
sc.pp.neighbors(adata_mouse_terg)
sc.tl.umap(adata_mouse_terg)
sc.pl.umap(adata_mouse_terg, color='cluster', save = 'mouse_terg_umap.png')

# recompute the neighbors for the adata_mouse_nonterg (ref dataset)
sc.pp.pca(adata_mouse_nonterg)
sc.pp.neighbors(adata_mouse_nonterg)
sc.tl.umap(adata_mouse_nonterg)
sc.pl.umap(adata_mouse_nonterg, color='cluster', save = 'mouse_nonterg_umap.png')

# Ingest using adata_mouse_terg as reference and adata_human_ERGpos as the target
sc.tl.ingest(adata_human_ERGpos, adata_mouse_terg, obs='cluster')

# Ingest using adata_mouse_nonterg as reference and adata_human_ERGneg as the target
sc.tl.ingest(adata_human_ERGneg, adata_mouse_nonterg, obs='cluster')

print(len(adata_mouse_terg.obs_names))
print(len(adata_human_ERGpos.obs_names))

# cell frequency in each cluster
print(adata_mouse_terg.obs.cluster.value_counts())
print(adata_human_ERGpos.obs.cluster.value_counts())

# cell frequency in each cluster
print(adata_mouse_nonterg.obs.cluster.value_counts())
print(adata_human_ERGneg.obs.cluster.value_counts())

# concat adata_human_ERGpos and adata_human_ERGneg into a new adata_human
adata_human_list = [adata_human_ERGpos, adata_human_ERGneg]
adata_human_new = ad.concat(adata_human_list, label = 'erg', join = 'outer')

# save the human count matrix with the cluster info
adata_human.write('human/h5ads/erg_fibroblasts_scvi_v6_regulons_annot.h5ad')
adata_human_new.write('human/h5ads/erg_fibroblasts_scvi_v6_regulons_annot_new.h5ad')
adata_human_ERGpos.write('human/h5ads/erg_fibroblasts_scvi_v6_regulons_ERGpos_annot.h5ad')

###########################################################################
# load the annot human data
#%%
adata_human_new = sc.read_h5ad('human/erg_fibroblasts_scvi_v6_regulons_annot_new.h5ad', chunk_size=100000)

# re-cap the gene symbols
adata_human_new.var_names = [gene.upper() for gene in adata_human_new.var_names]
# subset also the adata_human.raw
tempAdata = adata_human_new.raw.to_adata()
tempAdata.var_names = [gene.upper() for gene in tempAdata.var_names]
adata_human_new.raw = tempAdata

#######################
## umap clusters
# mouse
sc.pl.umap(adata_mouse, color='cluster', save= '_mouse.png')

# human
sc.pl.umap(adata_human_new, color='cluster', save = '_human.png')

##########################
# Acta2 umap
sc.pl.umap(adata_mouse, color=['cluster', 'Acta2'], save = '_mouse_ACTA2.png')

sc.pl.umap(adata_human_new, color=['cluster', 'ACTA2'], save = '_human_ACTA2.png')

##########################
## cluster markers

# mouse
sc.tl.rank_genes_groups(adata_mouse, 'cluster', pts=True, use_raw = False)
sc.pl.rank_genes_groups(adata_mouse, n_genes = 25, sharey=False, save = '_mouseMarkers_notRaw.png')
markers_mouse_c3 = sc.get.rank_genes_groups_df(adata_mouse, group = '3')

# human
sc.tl.rank_genes_groups(adata_human_new, 'cluster', pts=True, use_raw = False)
sc.pl.rank_genes_groups(adata_human_new, n_genes=25, sharey=False, save = '_humanMarkers_notRaw.png')

###########################
# Plot Ar mouse
sc.pl.violin(adata_mouse, ['Ar'], groupby = 'cluster', save = '_AR_mouse.png')
sc.pl.umap(adata_mouse, color='Ar', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Ar_mouse_all.png')
############################
# plot Soat1 mouse
sc.pl.violin(adata_mouse, ['Soat1'], groupby = 'cluster', save = '_Soat1_mouse.png')
sc.pl.umap(adata_mouse, color='Soat1', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Soat1_mouse_all.png')
############################
# plot Acat1 mouse
sc.pl.violin(adata_mouse, ['Acat1'], groupby = 'cluster', save = '_Acat1_mouse.png')
sc.pl.umap(adata_mouse, color='Acat1', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Acat1_mouse_all.png')
############################
# plot AR human
sc.pl.violin(adata_human_new, ['AR'], groupby = 'cluster', save = '_AR_human.png')
sc.pl.umap(adata_human_new, color=['AR'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Ar_human_all.png')
#############################
# plot SOAT1 human
sc.pl.violin(adata_human_new, ['SOAT1'], groupby = 'cluster', save = '_SOAT1_human.png')
sc.pl.umap(adata_human_new, color='SOAT1', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_SOAT1_human_all.png')
##############################
# plot ACAT1 human
sc.pl.violin(adata_human_new, ['ACAT1'], groupby = 'cluster', save = '_ACAT1_human.png')
sc.pl.umap(adata_human_new, color='ACAT1', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_ACAT1_human_all.png')
###############################
# violin for Periostin mouse
sc.pl.violin(adata_mouse, ['Postn'], groupby = 'cluster', save = '_Postn_mouse.png')
# umap for Periostin mouse
sc.pl.umap(adata_mouse, color= 'Postn', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_mouse_all.png')
###############################

# dotplot for periostin per different mouse models
dp = sc.pl.DotPlot(adata_mouse, var_names = 'Postn', groupby = 'key', cmap = 'Reds')
dp.legend(width=2.5).savefig('figures/dotplot_Postn_mouse.png')

################
# dotplot for Ar per different mouse models
dp = sc.pl.DotPlot(adata_mouse, var_names = 'Ar', groupby = 'key', cmap = 'Reds')
dp.legend(width=2.5).savefig('figures/dotplot_Ar_mouse.png')

###############################
# violin for Periostin human
sc.pl.violin(adata_human_new, ['POSTN'], groupby = 'cluster', save = '_Postn_human.png')

# umap for Periostin human
sc.pl.umap(adata_human_new, color='POSTN', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_human_all.png')


####################################################################
#################################################################
## plot PdgfrA and PdgfrB

## plot Pdgfra
# all clusters
sc.pl.umap(adata_mouse, color=['cluster', 'Pdgfra', 'Pdgfrb'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Pdgfra_Pdgfrb_mouse.png')

# dotplot for Pdgfra per different mouse models
dp = sc.pl.DotPlot(adata_mouse, var_names = ['Pdgfra', 'Pdgfrb'], groupby = 'cluster', cmap = 'Reds')
dp.legend(width=2.5).savefig('figures/dotplot_Pdgfra_Pdgfrb_mouse.png')


####################################################################
####################################################################

# c0 markers
sc.pl.rank_genes_groups_dotplot(adata_mouse, var_names = ['Acta2', 'Myl9', 'Myh11', 'Tagln', 'Pdgfra', 'Mustn1', 'Angpt2', 'Notch3'], save = '_c0_mouse.png')

sc.pl.rank_genes_groups_dotplot(adata_human_new, var_names = ['ACTA2', 'MYL9', 'MYH11', 'TAGLN', 'PDGFRA', 'MUSTN1', 'ANGPT2', 'NOTCH3'], save = '_c0_human.png')

###############################
# mouse
dp = sc.pl.DotPlot(adata_mouse, var_names = ['Acta2', 'Myl9', 'Myh11', 'Tagln', 'Pdgfra',
                                             'Mustn1', 'Angpt2', 'Notch3', 'Sfrp1', 'Gpx3',
                                             'C3', 'C7', 'Cfh', 'Ccl11', 'Cd55', 'Ptx3',
                                             'Thbd', 'Ifi204', 'Ifi205', 'Ifi207', 'Jun',
                                             'Junb', 'Jund', 'Fos', 'Fosb', 'Fosl2', 'Atf3',
                                             'Mafb', 'Maff', 'Nek2', 'Id1', 'Id3', 'Btg2',
                                             'Gadd45a', 'Hes1', 'Bcl3', 'Socs1', 'Socs3',
                                             'Il6', 'Irf1', 'Map3k8', 'Gadd45b', 'Gadd45g',
                                             'Dusp1', 'Dusp6', 'Klf4'],
                   groupby='cluster', categories_order = ['0','1','2','3','4','5','6','7'], cmap = 'Reds'
                   )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot1.png')
################################
# human
dp2 = sc.pl.DotPlot(adata_mouse, var_names = ['Sfrp2', 'Wnt5a', 'Lgr5', 'Apc',
                                                          'Wnt4', 'Wnt6', 'Notum', 'Wif1',
                                                          'Nkd1', 'Fzd1', 'Wnt2', 'Wnt10a',
                                                          'Dkk2', 'Rorb', 'Cxxc4', 'Nfat5',
                                                          'Apoe', 'Dact1', 'Ctnnb1', 'Lef1',
                                                          'Tcf4', 'Myc', 'Mki67', 'H2afx',
                                                          'Top2a', 'Ccnb1', 'Ccnb2', 'Stmn1',
                                                          'Ptn', 'Mdk', 'Tubb3', 'Mrc2', 'Fn1',
                                                          'Tnc', 'Col12a1', 'Col14a1', 'Col16a1',
                                                          'Mmp19', 'Cthrc1', 'Wisp1', 'Fzd1', 'Fzd2',
                                                          'Sfrp4', 'Bmp1', 'Tle3', 'Tgfb1', 'Tgfb1', 'Postn'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )

dp2.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot2.png')

###############################
dp = sc.pl.DotPlot(adata_human_new, var_names = ['ACTA2', 'MYL9', 'MYH11', 'TAGLN',
                                                              'PDGFRA', 'MUSTN1', 'ANGPT2', 'NOTCH3',
                                                              'SFRP1', 'GPX3', 'C3', 'C7', 'CFH', 'CCL11',
                                                              'CD55', 'PTX3', 'THBD', 'IFI16', 'JUN', 'JUNB',
                                                              'JUND', 'FOS', 'FOSB', 'FOSL2', 'ATF3',
                                                              'MAFB', 'MAFF', 'NEK2', 'ID1', 'ID3', 'BTG2',
                                                              'GADD45A', 'HES1', 'BCL3', 'SOCS1', 'SOCS3',
                                                              'IL6', 'IRF1', 'MAP3K8', 'GADD45B', 'GADD45G',
                                                              'DUSP1', 'DUSP6', 'KLF4'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot1.png')

###############################
dp2 = sc.pl.DotPlot(adata_human_new,
                                var_names = ['SFRP2', 'WNT5A', 'LGR5', 'APC', 'WNT4',
                                             'WNT6', 'NOTUM', 'WIF1', 'NKD1', 'FZD1',
                                             'WNT2', 'WNT10A', 'DKK2', 'RORB', 'CXXC4',
                                             'NFAT5', 'APOE', 'DACT1', 'CTNNB1', 'LEF1',
                                             'TCF4', 'MYC', 'MKI67', 'H2AFX', 'TOP2A',
                                             'CCNB1', 'CCNB2', 'STMN1', 'PTN', 'MDK',
                                             'TUBB3', 'MRC2', 'FN1', 'TNC', 'COL12A1',
                                             'COL14A1', 'COL16A1', 'MMP19', 'CTHRC1',
                                             'WISP1', 'FZD1', 'FZD2', 'SFRP4', 'BMP1',
                                             'TLE3', 'TGFB1', 'TGFB1', 'POSTN'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp2.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot2.png')

#################################################################
## which cells in mouse c5,6,7 co-express Postn plus Mki67, and Ar?

# subset to c5, c6, and c7
adata_mouse_c5c6c7 = adata_mouse[adata_mouse.obs['cluster'].isin(['5','6','7'])]
adata_mouse_c5c6c7.obs['cluster'].value_counts()

# plot Postn
sc.pl.violin(adata_mouse_c5c6c7, ['Postn'], groupby = 'cluster')

## plot Mki67
# all clusters
sc.pl.violin(adata_mouse, ['Mki67'], groupby = 'cluster', save='_Mki67_mouse.png')
sc.pl.umap(adata_mouse, color=['Postn', 'Mki67'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Mki67_mouse.png')

# just c5, c6, and c7
sc.pl.violin(adata_mouse_c5c6c7, ['Mki67'], groupby = 'cluster', save='_Postn+Mki67+.png')
sc.pl.umap(adata_mouse_c5c6c7, color=['Postn', 'Mki67'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Mki67_mouse_c5c6c7.png')

##########
## get cells that are Mki67+
countMat_c5c6c7 = pd.concat([adata_mouse_c5c6c7.to_df(), adata_mouse_c5c6c7.obs], axis=1)

mki67_cond = (countMat_c5c6c7.Mki67 >= countMat_c5c6c7['Mki67'].quantile(0.9))

adata_mouse_c5c6c7.obs['Mki67_status'] = 'negative'
adata_mouse_c5c6c7.obs.loc[mki67_cond, 'Mki67_status'] = 'positive'
adata_mouse_c5c6c7.obs['Mki67_status'].value_counts()

## N of Mki67+ cells in Postn+ clusters
pd.crosstab(adata_mouse_c5c6c7.obs['cluster'], adata_mouse_c5c6c7.obs['Mki67_status'])

##############
## do the same for all adata
countMat_all = pd.concat([adata_mouse.to_df(), adata_mouse.obs], axis=1)

mki67_cond_all = (countMat_all.Mki67 >= countMat_all['Mki67'].quantile(0.9))

adata_mouse.obs['Mki67_status'] = 'negative'
adata_mouse.obs.loc[mki67_cond_all, 'Mki67_status'] = 'positive'
adata_mouse.obs['Mki67_status'].value_counts()

## N of Mki67+ cells in all clusters
pd.crosstab(adata_mouse.obs['cluster'], adata_mouse.obs['Mki67_status'])

###########
## plot cells by Mki67 status

# only in c5 c6 and c7
sc.pl.umap(adata_mouse_c5c6c7, color=['Mki67', 'Mki67_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_c5c6c7_Mki67Expression_Mki67status.png')

# all clusters
sc.pl.umap(adata_mouse, color=['Mki67', 'Mki67_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_mouseAll_Mki67Expression_Mki67status.png')

Mki67_relativeFrequency_all = sct.tools.relative_frequency_per_cluster(adata_mouse, group_by='cluster', xlabel='Mki67_status', condition=None)
Mki67_relativeFrequency_all['cluster'] = 'c'+Mki67_relativeFrequency_all['cluster']
sct.plot.cluster_composition_stacked_barplot(Mki67_relativeFrequency_all, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_mouse.uns['Mki67_status_colors'], save = 'figures/Mki67_status.png')

#################################################################
## which cells in mouse co-express Postn plus Lgr5

# subset to c5 c6 and c7
adata_mouse_c5c6c7 = adata_mouse[adata_mouse.obs['cluster'].isin(['5','6','7'])]
adata_mouse_c5c6c7.obs['cluster'].value_counts()

# plot Postn
sc.pl.violin(adata_mouse_c5c6c7, ['Postn'], groupby = 'cluster')

## plot Mki67
# all clusters
sc.pl.violin(adata_mouse, ['Lgr5'], groupby = 'cluster', save='_Lgr5_mouse.png')
sc.pl.umap(adata_mouse, color=['Postn', 'Lgr5'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Lgr5_mouse.png')

# just c5, c6, and c7
sc.pl.violin(adata_mouse_c5c6c7, ['Lgr5'], groupby = 'cluster', save='_Postn+Lgr5+.png')
sc.pl.umap(adata_mouse_c5c6c7, color=['Postn', 'Lgr5'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Lgr5_mouse_c5c6c7.png')

##########
## get cells that are Lgr5+
countMat_c5c6c7 = pd.concat([adata_mouse_c5c6c7.to_df(), adata_mouse_c5c6c7.obs], axis=1)

Lgr5_cond = (countMat_c5c6c7.Lgr5 >= countMat_c5c6c7['Lgr5'].quantile(0.9))

adata_mouse_c5c6c7.obs['Lgr5_status'] = 'negative'
adata_mouse_c5c6c7.obs.loc[Lgr5_cond, 'Lgr5_status'] = 'positive'
adata_mouse_c5c6c7.obs['Lgr5_status'].value_counts()

## N of Mki67+ cells in Postn+ clusters
pd.crosstab(adata_mouse_c5c6c7.obs['cluster'], adata_mouse_c5c6c7.obs['Lgr5_status'])

############
## do the same for all adata
countMat_all = pd.concat([adata_mouse.to_df(), adata_mouse.obs], axis=1)

Lgr5_cond_all = (countMat_all.Lgr5 >= countMat_all['Lgr5'].quantile(0.9))

adata_mouse.obs['Lgr5_status'] = 'negative'
adata_mouse.obs.loc[Lgr5_cond_all, 'Lgr5_status'] = 'positive'
adata_mouse.obs['Lgr5_status'].value_counts()

## N of Mki67+ cells in all clusters
pd.crosstab(adata_mouse.obs['cluster'], adata_mouse.obs['Lgr5_status'])

############################
## plot cells by Lgr5 status

# only in c5 c6 and c7
sc.pl.umap(adata_mouse_c5c6c7, color=['Lgr5', 'Lgr5_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_c5c6c7_Lgr5Expression_Lgr5status.png')

# all clusters
sc.pl.umap(adata_mouse, color=['Lgr5', 'Lgr5_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_mouseAll_Lgr5Expression_Lgr5status.png')

Lgr5_relativeFrequency_all = sct.tools.relative_frequency_per_cluster(adata_mouse, group_by='cluster', xlabel='Lgr5_status', condition=None)
Lgr5_relativeFrequency_all['cluster'] = 'c'+Lgr5_relativeFrequency_all['cluster']
sct.plot.cluster_composition_stacked_barplot(Lgr5_relativeFrequency_all, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_mouse.uns['Lgr5_status_colors'], save = 'figures/Lgr5_status.png')

##################################################################
## Postn+ cells that express Lef1

## plot Lef1
# all clusters
sc.pl.violin(adata_mouse, ['Lef1'], groupby = 'cluster', save='_Lef1_mouse.png')
sc.pl.umap(adata_mouse, color=['Postn', 'Lef1'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Lef1_mouse.png')

# just c5, c6, and c7
sc.pl.violin(adata_mouse_c5c6c7, ['Lef1'], groupby = 'cluster', save='_Postn+Lef1+.png')
sc.pl.umap(adata_mouse_c5c6c7, color=['Postn', 'Lef1'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Lef1_mouse_c5c6c7.png')

##########
## get cells that are Lef1+
# Lef1 is not present in the normalized count matrix, present only in the raw one
# We get the raw then normalize it to match that done in Mki67
adata_mouse_c5c6c7_raw = adata_mouse_c5c6c7.raw.to_adata()
#sc.pp.log1p(adata_mouse_c5c6c7_raw)
#sc.pp.scale(adata_mouse_c5c6c7_raw, max_value=10)
countMat_c5c6c7['Lef1'].mean()
sc.pl.violin(adata_mouse_c5c6c7_raw, ['Lef1'], groupby = 'cluster')

countMat_c5c6c7 = pd.concat([adata_mouse_c5c6c7_raw.to_df(), adata_mouse_c5c6c7_raw.obs], axis=1)

Lef1_cond = (countMat_c5c6c7.Lef1 >= countMat_c5c6c7['Lef1'].quantile(0.962))

adata_mouse_c5c6c7.obs['Lef1_status'] = 'negative'
adata_mouse_c5c6c7.obs.loc[Lef1_cond, 'Lef1_status'] = 'positive'
adata_mouse_c5c6c7.obs['Lef1_status'].value_counts()

## N of Lef1+ cells in Postn+ clusters
pd.crosstab(adata_mouse_c5c6c7.obs['cluster'], adata_mouse.obs['Lef1_status'])

############
## do the same for all adata
adata_mouse_raw = adata_mouse.raw.to_adata()
countMat_all = pd.concat([adata_mouse_raw.to_df(), adata_mouse_raw.obs], axis=1)

Lef1_cond_all = (countMat_all.Lef1 >= countMat_all['Lef1'].quantile(0.962))

adata_mouse.obs['Lef1_status'] = 'negative'
adata_mouse.obs.loc[Lef1_cond_all, 'Lef1_status'] = 'positive'
adata_mouse.obs['Lef1_status'].value_counts()

###########
## plot cells by Lef1 status

# only in c5 c6 and c7
sc.pl.umap(adata_mouse_c5c6c7, color=['Lef1', 'Lef1_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_c5c6c7_Lef1Expression_Lef1status.png')

# all clusters
sc.pl.umap(adata_mouse, color=['Lef1', 'Lef1_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_mouseAll_Lef1Expression_Lef1status.png')

## N of Lef1+ cells in all clusters
pd.crosstab(adata_mouse.obs['cluster'], adata_mouse.obs['Lef1_status'])


####################################################################
#################################################################
## which cells in mouse co-express Postn plus Nrg1

# subset to c5 c6 and c7
adata_mouse_c5c6c7 = adata_mouse[adata_mouse.obs['cluster'].isin(['5','6','7'])]
adata_mouse_c5c6c7.obs['cluster'].value_counts()

## plot Nrg1
# all clusters
sc.pl.violin(adata_mouse, ['Nrg1'], groupby = 'cluster', save='_Nrg1_mouse.png')
sc.pl.umap(adata_mouse, color=['Postn', 'Nrg1'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Nrg1_mouse.png')

# just c5, c6, and c7
sc.pl.violin(adata_mouse_c5c6c7, ['Nrg1'], groupby = 'cluster', save='_Postn+Nrg1+.png')
sc.pl.umap(adata_mouse_c5c6c7, color=['Postn', 'Nrg1'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Nrg1_mouse_c5c6c7.png')

##########
## get cells that are Nrg1+
countMat_c5c6c7 = pd.concat([adata_mouse_c5c6c7.to_df(), adata_mouse_c5c6c7.obs], axis=1)

Nrg1_cond = (countMat_c5c6c7.Nrg1 >= countMat_c5c6c7['Nrg1'].quantile(0.9))

adata_mouse_c5c6c7.obs['Nrg1_status'] = 'negative'
adata_mouse_c5c6c7.obs.loc[Nrg1_cond, 'Nrg1_status'] = 'positive'
adata_mouse_c5c6c7.obs['Nrg1_status'].value_counts()

## N of Mki67+ cells in Postn+ clusters
pd.crosstab(adata_mouse_c5c6c7.obs['cluster'], adata_mouse_c5c6c7.obs['Nrg1_status'])

############
## do the same for all adata
countMat_all = pd.concat([adata_mouse.to_df(), adata_mouse.obs], axis=1)

Nrg1_cond_all = (countMat_all.Nrg1 >= countMat_all['Nrg1'].quantile(0.9))

adata_mouse.obs['Nrg1_status'] = 'negative'
adata_mouse.obs.loc[Nrg1_cond_all, 'Nrg1_status'] = 'positive'
adata_mouse.obs['Nrg1_status'].value_counts()

## N of Mki67+ cells in all clusters
pd.crosstab(adata_mouse.obs['cluster'], adata_mouse.obs['Nrg1_status'])

############################
## plot cells by Lgr5 status

# only in c5 c6 and c7
sc.pl.umap(adata_mouse_c5c6c7, color=['Nrg1', 'Nrg1_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_c5c6c7_Nrg1Expression_Nrg1status.png')

# all clusters
sc.pl.umap(adata_mouse, color=['Nrg1', 'Nrg1_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_mouseAll_Nrg1Expression_Nrg1status.png')

Nrg1_relativeFrequency_all = sct.tools.relative_frequency_per_cluster(adata_mouse, group_by='cluster', xlabel='Nrg1_status', condition=None)
Nrg1_relativeFrequency_all['cluster'] = 'c'+Nrg1_relativeFrequency_all['cluster']
sct.plot.cluster_composition_stacked_barplot(Nrg1_relativeFrequency_all, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_mouse.uns['Nrg1_status_colors'], save = 'figures/Nrg1_status.png')

####################################################################
#################################################################
## which cells in mouse co-express Postn plus Nrg2

# subset to c5 c6 and c7
adata_mouse_c5c6c7 = adata_mouse[adata_mouse.obs['cluster'].isin(['5','6','7'])]
adata_mouse_c5c6c7.obs['cluster'].value_counts()

## plot Nrg1
# all clusters
sc.pl.violin(adata_mouse, ['Nrg2'], groupby = 'cluster', save='_Nrg2_mouse.png')
sc.pl.umap(adata_mouse, color=['Postn', 'Nrg2'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Nrg2_mouse.png')

# just c5, c6, and c7
sc.pl.violin(adata_mouse_c5c6c7, ['Nrg2'], groupby = 'cluster', save='_Postn+Nrg2+.png')
sc.pl.umap(adata_mouse_c5c6c7, color=['Postn', 'Nrg2'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Nrg2_mouse_c5c6c7.png')

##########
## get cells that are Nrg2+

############
## do the same for all adata
adata_mouse_raw = adata_mouse.raw.to_adata()
countMat_all = pd.concat([adata_mouse_raw.to_df(), adata_mouse_raw.obs], axis=1)

Nrg2_cond_all = (countMat_all.Nrg2 >= countMat_all['Nrg2'].quantile(0.94))

adata_mouse_raw.obs['Nrg2_status'] = 'negative'
adata_mouse_raw.obs.loc[Nrg2_cond_all, 'Nrg2_status'] = 'positive'
adata_mouse_raw.obs['Nrg2_status'].value_counts()

## N of Mki67+ cells in all clusters
pd.crosstab(adata_mouse_raw.obs['cluster'], adata_mouse_raw.obs['Nrg2_status'])

############################
## plot cells by Lgr5 status

# all clusters
sc.pl.umap(adata_mouse_raw, color=['Nrg2', 'Nrg2_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_mouseAll_Nrg2Expression_Nrg2status.png')

Nrg2_relativeFrequency_all = sct.tools.relative_frequency_per_cluster(adata_mouse_raw, group_by='cluster', xlabel='Nrg2_status', condition=None)
Nrg2_relativeFrequency_all['cluster'] = 'c'+Nrg2_relativeFrequency_all['cluster']
sct.plot.cluster_composition_stacked_barplot(Nrg2_relativeFrequency_all, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_mouse_raw.uns['Nrg2_status_colors'], save = 'figures/Nrg2_status.png')

####################################################################
#################################################################
## which cells in mouse express Oncecut2

# subset to c5 c6 and c7
adata_mouse_c5c6c7 = adata_mouse[adata_mouse.obs['cluster'].isin(['5','6','7'])]
adata_mouse_c5c6c7.obs['cluster'].value_counts()

## plot Nrg1
# all clusters
sc.pl.violin(adata_mouse, ['Onecut2'], groupby = 'cluster', save='_Onecut2_mouse.png')
sc.pl.umap(adata_mouse, color=['Postn', 'Onecut2'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Onecut2_mouse.png')

# just c5, c6, and c7
sc.pl.violin(adata_mouse_c5c6c7, ['Onecut2'], groupby = 'cluster', save='_Postn+Onecut2+.png')
sc.pl.umap(adata_mouse_c5c6c7, color=['Postn', 'Onecut2'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Onecut2_mouse_c5c6c7.png')

##########
## get cells that are Nrg1+
countMat_c5c6c7 = pd.concat([adata_mouse_c5c6c7.to_df(), adata_mouse_c5c6c7.obs], axis=1)

Onecut2_cond = (countMat_c5c6c7.Onecut2 >= countMat_c5c6c7['Onecut2'].quantile(0.9))

adata_mouse_c5c6c7.obs['Onecut2_status'] = 'negative'
adata_mouse_c5c6c7.obs.loc[Onecut2_cond, 'Onecut2_status'] = 'positive'
adata_mouse_c5c6c7.obs['Onecut2_status'].value_counts()

## N of Mki67+ cells in Postn+ clusters
pd.crosstab(adata_mouse_c5c6c7.obs['cluster'], adata_mouse_c5c6c7.obs['Onecut2_status'])

############
## do the same for all adata
countMat_all = pd.concat([adata_mouse.to_df(), adata_mouse.obs], axis=1)

Onecut2_cond_all = (countMat_all.Onecut2 >= countMat_all['Onecut2'].quantile(0.98))

adata_mouse.obs['Onecut2_status'] = 'negative'
adata_mouse.obs.loc[Onecut2_cond_all, 'Onecut2_status'] = 'positive'
adata_mouse.obs['Onecut2_status'].value_counts()

## N of Mki67+ cells in all clusters
pd.crosstab(adata_mouse.obs['cluster'], adata_mouse.obs['Onecut2_status'])

############################
## plot cells by OnceCut2 status

# all clusters
sc.pl.umap(adata_mouse, color=['Onecut2', 'Onecut2_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_mouseAll_Onecut2Expression_Onecut2status.png')

Onecut2_relativeFrequency_all = sct.tools.relative_frequency_per_cluster(adata_mouse, group_by='cluster', xlabel='Nrg1_status', condition=None)
Onecut2_relativeFrequency_all['cluster'] = 'c'+Onecut2_relativeFrequency_all['cluster']
sct.plot.cluster_composition_stacked_barplot(Onecut2_relativeFrequency_all, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_mouse.uns['Onecut2_status_colors'], save = 'figures/Onecut2_status.png')

####################################################################
#################################################################
## which cells in mouse co-express Ccl11 plus Ccl20

# subset to c5 c6 and c7
adata_mouse_c0c1 = adata_mouse[adata_mouse.obs['cluster'].isin(['0','1'])]
adata_mouse_c0c1.obs['cluster'].value_counts()

## plot Nrg1
# all clusters
sc.pl.violin(adata_mouse, ['Ccl11', 'Ccl20'], groupby = 'cluster', save='_Ccl11_Ccl20_mouse.png')
sc.pl.umap(adata_mouse, color=['Ccl11', 'Ccl20'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Ccl11_Ccl20_mouse.png')

# just c0, c6, and c7
sc.pl.umap(adata_mouse_c0c1, color=['Ccl11', 'Ccl20'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Ccl11_Ccl20_mouse_c0c1.png')

##########
## get cells that are CCl11+

countMat_all = pd.concat([adata_mouse.to_df(), adata_mouse.obs], axis=1)

Ccl11_cond_all = (countMat_all.Ccl11 >= countMat_all['Ccl11'].quantile(0.9))

adata_mouse.obs['Ccl11_status'] = 'negative'
adata_mouse.obs.loc[Ccl11_cond_all, 'Ccl11_status'] = 'positive'
adata_mouse.obs['Ccl11_status'].value_counts()

## N of Mki67+ cells in all clusters
pd.crosstab(adata_mouse.obs['cluster'], adata_mouse.obs['Ccl11_status'])

########
## get cells that are CCl20+

countMat_all = pd.concat([adata_mouse.to_df(), adata_mouse.obs], axis=1)

Ccl20_cond_all = (countMat_all.Ccl20 >= countMat_all['Ccl20'].quantile(0.99))

adata_mouse.obs['Ccl20_status'] = 'negative'
adata_mouse.obs.loc[Ccl20_cond_all, 'Ccl20_status'] = 'positive'
adata_mouse.obs['Ccl20_status'].value_counts()

## N of Mki67+ cells in all clusters
pd.crosstab(adata_mouse.obs['cluster'], adata_mouse.obs['Ccl20_status'])

############################
## plot cells by Ccl11 and Ccl20 status

# Ccl11
sc.pl.umap(adata_mouse, color=['Ccl11', 'Ccl11_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_mouseAll_Ccl11Expression_Ccl11status.png')

Ccl11_relativeFrequency_all = sct.tools.relative_frequency_per_cluster(adata_mouse, group_by='cluster', xlabel='Ccl11_status', condition=None)
Ccl11_relativeFrequency_all['cluster'] = 'c'+Ccl11_relativeFrequency_all['cluster']
sct.plot.cluster_composition_stacked_barplot(Ccl11_relativeFrequency_all, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_mouse.uns['Ccl11_status_colors'], save = 'figures/Ccl11_status.png')

#######
# Ccl20
sc.pl.umap(adata_mouse, color=['Ccl20', 'Ccl20_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_mouseAll_Ccl20Expression_Ccl20status.png')

Ccl20_relativeFrequency_all = sct.tools.relative_frequency_per_cluster(adata_mouse, group_by='cluster', xlabel='Ccl20_status', condition=None)
Ccl20_relativeFrequency_all['cluster'] = 'c'+Ccl20_relativeFrequency_all['cluster']
sct.plot.cluster_composition_stacked_barplot(Ccl20_relativeFrequency_all, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_mouse.uns['Ccl20_status_colors'], save = 'figures/Ccl20_status.png')


########################################################################################################################################
########################################################################################################################################
# umaps for Notch1, Notch2, Notch3, Hes1Tubb3, Soat1, Acat1 and Lef1
sc.pl.umap(adata_mouse, color=['Notch1', 'Notch2', 'Notch3', 'Hes1', 'Tubb3', 'Soat1', 'Acat1', 'Lef1'], color_map = 'RdBu_r', vmin='p1', vmax='p99', ncols=3, save = '_Notch_Hes1_Tubb3_Soat1_Acat1_Lef1_mouse.png')

# dotplot
dp = sc.pl.DotPlot(adata_mouse,
                             var_names = ['Notch1', 'Notch2', 'Notch3', 'Hes1', 'Tubb3', 'Soat1', 'Acat1', 'Lef1'],
                             categories_order = ['0','1','2','3','4','5','6','7'],
                             groupby='cluster', cmap = 'Reds'
                            )
dp.savefig('figures/DotPlot_Notch_Hes1_Tubb3_Soat1_Acat1_Lef1_mouse.png')

####################################################################
## get cluster markers

# mouse
sc.tl.rank_genes_groups(adata_mouse, 'cluster', pts=True, use_raw = False)

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
sc.tl.rank_genes_groups(adata_human_new, 'cluster', pts=True, use_raw = False)

markers_human_c0 = sc.get.rank_genes_groups_df(adata_human_new, group = '0')
markers_human_c1 = sc.get.rank_genes_groups_df(adata_human_new, group = '1')
markers_human_c2 = sc.get.rank_genes_groups_df(adata_human_new, group = '2')
markers_human_c3 = sc.get.rank_genes_groups_df(adata_human_new, group = '3')
markers_human_c4 = sc.get.rank_genes_groups_df(adata_human_new, group = '4')
markers_human_c5 = sc.get.rank_genes_groups_df(adata_human_new, group = '5')
markers_human_c6 = sc.get.rank_genes_groups_df(adata_human_new, group = '6')
markers_human_c7 = sc.get.rank_genes_groups_df(adata_human_new, group = '7')

# fix gene symbols
markers_human_c0['names'] = markers_human_c0['names'].str.title()
markers_human_c1['names'] = markers_human_c1['names'].str.title()
markers_human_c2['names'] = markers_human_c2['names'].str.title()
markers_human_c3['names'] = markers_human_c3['names'].str.title()
markers_human_c4['names'] = markers_human_c4['names'].str.title()
markers_human_c5['names'] = markers_human_c5['names'].str.title()
markers_human_c6['names'] = markers_human_c6['names'].str.title()
markers_human_c7['names'] = markers_human_c7['names'].str.title()

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
top100_human_c0 = markers_human_c0[0:100]
top100_human_c1 = markers_human_c1[0:100]
top100_human_c2 = markers_human_c2[0:100]
top100_human_c3 = markers_human_c3[0:100]
top100_human_c4 = markers_human_c4[0:100]
top100_human_c5 = markers_human_c5[0:100]
top100_human_c6 = markers_human_c6[0:100]
top100_human_c7 = markers_human_c7[0:100]

#######################
## get the genes in common
intersection_c0 = pd.merge(top100_mouse_c0, top100_human_c0, how='inner', on=['names'])
intersection_c1 = pd.merge(top100_mouse_c1, top100_human_c1, how='inner', on=['names'])
intersection_c2 = pd.merge(top100_mouse_c2, top100_human_c2, how='inner', on=['names'])
intersection_c3 = pd.merge(top100_mouse_c3, top100_human_c3, how='inner', on=['names'])
intersection_c4 = pd.merge(top100_mouse_c4, top100_human_c4, how='inner', on=['names'])
intersection_c5 = pd.merge(top100_mouse_c5, top100_human_c5, how='inner', on=['names'])
intersection_c6 = pd.merge(top100_mouse_c6, top100_human_c6, how='inner', on=['names'])
intersection_c7 = pd.merge(top100_mouse_c7, top100_human_c7, how='inner', on=['names'])

##########################
## Plot the genes in common

# dotplot 1 (c0,c1,c2): common markers only
c0c1c2_common = [intersection_c0['names'][1:20], intersection_c1['names'][1:20], intersection_c2['names'][1:20]]
c0c1c2_common = [x for xs in c0c1c2_common for x in xs]

# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = c0c1c2_common,
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c0c1c2_CommonOnly.png')

#####
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in c0c1c2_common],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c0c1c2_CommonOnly.png')

####################
# dotplot 1 (c3,c4, c5, c6, c7): common markers only
c3c4c5c6c7_common = [intersection_c3['names'][1:20], intersection_c4['names'][1:20], intersection_c5['names'][1:20], intersection_c6['names'][1:20], intersection_c7['names'][1:20]]
c3c4c5c6c7_common = [x for xs in c3c4c5c6c7_common for x in xs]

# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = c3c4c5c6c7_common,
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c3c4c5c6c7_CommonOnly.png')

######
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in c3c4c5c6c7_common],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c3c4c5c6c7_CommonOnly.png')

##############################################################################
## Plot the genes in common: individual clusters

# dotplot c0: common markers only
# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = intersection_c0['names'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c0_CommonOnly.png')

#######
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in intersection_c0['names']],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c0_CommonOnly.png')

###########################

# dotplot c1: common markers only
# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = intersection_c1['names'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c1_CommonOnly.png')

#####
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in intersection_c1['names']],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c1_CommonOnly.png')

####################

# dotplot c2: common markers only
# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = intersection_c2['names'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c2_CommonOnly.png')

#####
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in intersection_c2['names']],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c2_CommonOnly.png')

########################
# dotplot c3: common markers only
# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = intersection_c3['names'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c3_CommonOnly.png')

########
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in intersection_c3['names']],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c3_CommonOnly.png')

##########################
# dotplot c4: common markers only
# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = intersection_c4['names'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c4_CommonOnly.png')

#######
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in intersection_c4['names']],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c4_CommonOnly.png')

########################
# dotplot c5: common markers only
# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = intersection_c5['names'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c5_CommonOnly.png')

#####
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in intersection_c5['names']],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c5_CommonOnly.png')

########################
# dotplot c6: common markers only
# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = intersection_c6['names'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c6_CommonOnly.png')

#####
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in intersection_c6['names']],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c6_CommonOnly.png')

############################
# dotplot c7: common markers only
# mouse
dp = sc.pl.DotPlot(adata_mouse,
                                var_names = intersection_c7['names'],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/mouse_dotplot_c7_CommonOnly.png')

#######
# human
dp = sc.pl.DotPlot(adata_human_new,
                                var_names = [gene.upper() for gene in intersection_c7['names']],
                                categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_dotplot_c7_CommonOnly.png')

###################################

#################################################################
## which cells in human c5,6,7 co-express Postn plus Mki67, and Ar?

# subset to c5, c6, and c7
adata_human_new_c5c6c7 = adata_human_new[adata_human_new.obs['cluster'].isin(['5','6','7'])]
adata_human_new_c5c6c7.obs['cluster'].value_counts()

# plot Postn
sc.pl.violin(adata_human_new_c5c6c7, ['POSTN'], groupby = 'cluster')

## plot Mki67
# all clusters
sc.pl.violin(adata_human_new, ['MKI67'], groupby = 'cluster', save='_Mki67_human.png')
sc.pl.umap(adata_human_new, color=['POSTN', 'MKI67'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Mki67_human.png')

# just c5, c6, and c7
sc.pl.violin(adata_human_new_c5c6c7, ['MKI67'], groupby = 'cluster', save='_human_Postn+Mki67+.png')
sc.pl.umap(adata_human_new_c5c6c7, color=['POSTN', 'MKI67'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Mki67_human_c5c6c7.png')

##########
## get cells that are Mki67+
countMat_c5c6c7_human = pd.concat([adata_human_new_c5c6c7.to_df(), adata_human_new_c5c6c7.obs], axis=1)

mki67_cond_human = (countMat_c5c6c7_human.MKI67 >= countMat_c5c6c7_human['MKI67'].quantile(0.9))

adata_human_new_c5c6c7.obs['MKI67_status'] = 'negative'
adata_human_new_c5c6c7.obs.loc[mki67_cond_human, 'MKI67_status'] = 'positive'
adata_human_new_c5c6c7.obs['MKI67_status'].value_counts()

## N of Mki67+ cells in Postn+ clusters
pd.crosstab(adata_human_new_c5c6c7.obs['cluster'], adata_human_new_c5c6c7.obs['MKI67_status'])

############
## do the same for all adata
countMat_all_human = pd.concat([adata_human_new.to_df(), adata_human_new.obs], axis=1)

mki67_cond_all_human = (countMat_all_human.MKI67 >= countMat_all_human['MKI67'].quantile(0.9))

adata_human_new.obs['MKI67_status'] = 'negative'
adata_human_new.obs.loc[mki67_cond_all_human, 'MKI67_status'] = 'positive'
adata_human_new.obs['MKI67_status'].value_counts()

## N of Mki67+ cells in all clusters
pd.crosstab(adata_human_new.obs['cluster'], adata_human_new.obs['MKI67_status'])

###########
## plot cells by Mki67 status

# only in c5 c6 and c7
sc.pl.umap(adata_human_new_c5c6c7, color=['MKI67', 'MKI67_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_c5c6c7_MKI67Expression_MKI67status.png')

# all clusters
sc.pl.umap(adata_human_new, color=['MKI67', 'MKI67_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_humanAll_MKI67Expression_MKI67status.png')

Mki67_relativeFrequency_all_human = sct.tools.relative_frequency_per_cluster(adata_human_new, group_by='cluster', xlabel='MKI67_status', condition=None)
Mki67_relativeFrequency_all_human['cluster'] = 'c'+Mki67_relativeFrequency_all_human['cluster']
sct.plot.cluster_composition_stacked_barplot(Mki67_relativeFrequency_all_human, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_human_new.uns['MKI67_status_colors'], save = 'figures/MKI67_status_human.png')

#################################################################
## which cells in human co-express Postn plus Lgr5

adata_human_new_raw = adata_human_new.raw.to_adata()

# subset to c5 c6 and c7
adata_human_new_raw_c5c6c7 = adata_human_new_raw[adata_human_new_raw.obs['cluster'].isin(['5','6','7'])]
adata_human_new_raw_c5c6c7.obs['cluster'].value_counts()

# plot Postn
sc.pl.violin(adata_human_new_raw_c5c6c7, ['POSTN'], groupby = 'cluster')

## plot Mki67
# all clusters
sc.pl.violin(adata_human_new_raw, ['LGR5'], groupby = 'cluster', save='_Lgr5_human.png')
sc.pl.umap(adata_human_new_raw, color=['POSTN', 'LGR5'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Lgr5_human.png')

# just c5, c6, and c7
sc.pl.violin(adata_human_new_raw_c5c6c7, ['LGR5'], groupby = 'cluster', save='_Postn+Lgr5+_human.png')
sc.pl.umap(adata_human_new_raw_c5c6c7, color=['POSTN', 'LGR5'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_Lgr5_human_c5c6c7.png')


##########
## get cells that are Lgr5+

countMat_c5c6c7_human = pd.concat([adata_human_new_raw_c5c6c7.to_df(), adata_human_new_raw_c5c6c7.obs], axis=1)

Lgr5_cond_human = (countMat_c5c6c7_human.LGR5 >= countMat_c5c6c7_human['LGR5'].quantile(0.999))

adata_human_new_raw_c5c6c7.obs['LGR5_status'] = 'negative'
adata_human_new_raw_c5c6c7.obs.loc[Lgr5_cond_human, 'LGR5_status'] = 'positive'
adata_human_new_raw_c5c6c7.obs['LGR5_status'].value_counts()

## N of Mki67+ cells in Postn+ clusters
pd.crosstab(adata_human_new_raw_c5c6c7.obs['cluster'], adata_human_new_raw_c5c6c7.obs['Lgr5_status'])

############
## do the same for all adata
countMat_all_human = pd.concat([adata_human_new_raw.to_df(), adata_human_new_raw.obs], axis=1)

Lgr5_cond_all_human = (countMat_all_human.LGR5 >= countMat_all_human['LGR5'].quantile(0.9985))

adata_human_new_raw.obs['LGR5_status'] = 'negative'
adata_human_new_raw.obs.loc[Lgr5_cond_all_human, 'LGR5_status'] = 'positive'
adata_human_new_raw.obs['LGR5_status'].value_counts()

## N of Mki67+ cells in all clusters
pd.crosstab(adata_human_new_raw.obs['cluster'], adata_human_new_raw.obs['LGR5_status'])

############
## plot cells by Lgr5 status

# only in c5, c6, and c7
sc.pl.umap(adata_mouse_c5c6c7, color=['LGR5', 'LGR5_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_c5c6c7_Lgr5Expression_Lgr5status.png')

# all clusters
sc.pl.umap(adata_human_new_raw, color=['LGR5', 'LGR5_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_humanAll_Lgr5Expression_Lgr5status.png')

Lgr5_relativeFrequency_all_human = sct.tools.relative_frequency_per_cluster(adata_human_new_raw, group_by='cluster', xlabel='LGR5_status', condition=None)
Lgr5_relativeFrequency_all_human['cluster'] = 'c'+Lgr5_relativeFrequency_all_human['cluster']
sct.plot.cluster_composition_stacked_barplot(Lgr5_relativeFrequency_all_human, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_human_new_raw.uns['LGR5_status_colors'], save = 'figures/Lgr5_status_human.png')

##############################################################################################

markers_human_c0.to_csv('markers_human_c0.csv')
markers_human_c1.to_csv('markers_human_c1.csv')
markers_human_c2.to_csv('markers_human_c2.csv')
markers_human_c3.to_csv('markers_human_c3.csv')
markers_human_c4.to_csv('markers_human_c4.csv')
markers_human_c5.to_csv('markers_human_c5.csv')
markers_human_c6.to_csv('markers_human_c6.csv')
markers_human_c7.to_csv('markers_human_c7.csv')


#####################################
# umap by genotype
adata_mouse.obs.condition.value_counts()
# mutant
sc.pl.umap(
    adata_mouse[adata_mouse.obs["condition"] == "mutant"],
    color="key", size = 10,
    save = '_mouse_models_mutant.png'
)
# Wildtype
sc.pl.umap(
    adata_mouse[adata_mouse.obs["condition"] == "wildtype"],
    color="key", size = 10,
    save = '_mouse_models_wildtype.png'
)

#####################################
# umapo by n_genes
sc.pl.umap(
    adata_mouse,
    color="n_genes", size = 10,
    title = 'Genes',
    save = '_mouse_genes.png'
)

# umap by n_counts
sc.pl.umap(
    adata_mouse,
    color="n_counts", size = 10,
    title = 'Counts',
    save = '_mouse_counts.png'
)

#####################################
# mouse smooth muscles

# umap just c0
sc.pl.umap(
    adata_mouse,
    color="cluster", size = 10,
    groups = ['0'],
    title = 'c0',
    palette = 'turbo',
    na_in_legend = False,
    #add_outline = True, outline_width = [0.05,0.005],
    save = '_mouse_smooth_muscle.png'
)


adata_mouse_myo = sc.read_h5ad('outs/h5ads/fapcm_myo.h5ad', chunk_size=100000)

mapdict = {
    "0": "c0.1",
    "1": "c0.2",
    "2": "c0.1",
    "3": "c0.2",
}
adata_mouse_myo.obs["cluster"] = adata_mouse_myo.obs["leiden"].map(mapdict)
adata_mouse_myo.obs["cluster"] = adata_mouse_myo.obs["cluster"].astype("category")

adata_mouse_myo.obs.cluster.value_counts()

sc.tl.dendrogram(adata_mouse_myo, groupby="cluster")
sc.pl.dotplot(
    adata_mouse_myo,
    var_names=[
        "Rgs5",
        "Mef2c",
        "Vtn",
        "Cygb",
        "Pdgfrb",
        "Rras",
        "Rasl12",
        "Rspo3",
        "Nrg2",
        "Hopx",
        "Actg2",
    ],
    groupby="cluster",
    cmap="Reds",
    vmax=6,
    save="_myo_markers.png",
)

# umap myo and pericytes
sc.pl.umap(
    adata_mouse_myo,
    color="cluster", size = 50,
    color_map = 'viridis',
    title = 'c0 subtypes',
    save = '_smoothMuscle_subtypes.png'
)

# myo mouse models
sc.pl.umap(
    adata_mouse_myo,
    color="key", size = 50,
    title = 'mouse models',
    save = '_smoothMuscle_models.png'
)

######
# myo markers
sc.pl.umap(
    adata_mouse_myo,
    color=[
        "Acta2",
        "Myl9",
        "Myh11",
        "Tagln",
        "Rgs5",
        "Mef2c",
        "Pdgfrb"
    ],
    cmap="RdBu_r",
    vmax=5,
    save="_myo_markers.png"
)

## divide by c0.1 and c0.2
# c0.1
sc.pl.umap(
    adata_mouse_myo[adata_mouse_myo.obs['cluster']=='c0.1'],
    color=[
        "Acta2",
        "Myl9",
        "Myh11",
        "Tagln",
        "Rgs5",
        "Mef2c",
        "Pdgfrb"
    ],
    cmap="RdBu_r",
    vmax=5,
    save="_myo_c0.1_markers.png"
)

# c0.2
sc.pl.umap(
    adata_mouse_myo[adata_mouse_myo.obs['cluster']=='c0.2'],
    color=[
        "Acta2",
        "Myl9",
        "Myh11",
        "Tagln",
        "Rgs5",
        "Mef2c",
        "Pdgfrb"
    ],
    cmap="RdBu_r",
    vmax=5,
    save="_myo_c0.2_markers.png"
)
##############################






top100_mouse_c0 = markers_mouse_c0[0:100]
top100_mouse_c1 = markers_mouse_c1[0:100]
top100_mouse_c2 = markers_mouse_c2[0:100]
top100_mouse_c3 = markers_mouse_c3[0:100]
top100_mouse_c4 = markers_mouse_c4[0:100]
top100_mouse_c5 = markers_mouse_c5[0:100]
top100_mouse_c6 = markers_mouse_c6[0:100]
top100_mouse_c7 = markers_mouse_c7[0:100]

top100_human_c0 = markers_human_c0[0:100]
top100_human_c1 = markers_human_c1[0:100]
top100_human_c2 = markers_human_c2[0:100]
top100_human_c3 = markers_human_c3[0:100]
top100_human_c4 = markers_human_c4[0:100]
top100_human_c5 = markers_human_c5[0:100]
top100_human_c6 = markers_human_c6[0:100]
top100_human_c7 = markers_human_c7[0:100]

intersection_c0 = pd.merge(top100_mouse_c0, top100_human_c0, how='inner', on=['names'])
intersection_c1 = pd.merge(top100_mouse_c1, top100_human_c1, how='inner', on=['names'])
intersection_c2 = pd.merge(top100_mouse_c2, top100_human_c2, how='inner', on=['names'])
intersection_c3 = pd.merge(top100_mouse_c3, top100_human_c3, how='inner', on=['names'])
intersection_c4 = pd.merge(top100_mouse_c4, top100_human_c4, how='inner', on=['names'])
intersection_c5 = pd.merge(top100_mouse_c5, top100_human_c5, how='inner', on=['names'])
intersection_c6 = pd.merge(top100_mouse_c6, top100_human_c6, how='inner', on=['names'])
intersection_c7 = pd.merge(top100_mouse_c7, top100_human_c7, how='inner', on=['names'])

markers_mouse_c5 = sc.get.rank_genes_groups_df(adata_mouse, group = '5')



