
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


sc.pl.umap(adata_mouse, color='cluster', save= '_mouse.png')


sc.pl.umap(adata_human_new, color='cluster', save = '_human.png')

sc.pl.umap(adata_mouse, color=['cluster', 'Acta2'], save = '_mouse_ACTA2.png')

sc.pl.umap(adata_human_new, color=['cluster', 'ACTA2'], save = '_human_ACTA2.png')

sc.tl.rank_genes_groups(adata_mouse, 'cluster', pts=True, use_raw = False)
sc.pl.rank_genes_groups(adata_mouse, n_genes = 25, sharey=False, save = '_mouseMarkers.png')
markers_mouse_c3 = sc.get.rank_genes_groups_df(adata_mouse, group = '3')


sc.tl.rank_genes_groups(adata_human_new, 'cluster', pts=True, use_raw = False)
sc.pl.rank_genes_groups(adata_human_new, n_genes=25, sharey=False, save = '_humanMarkers.png')

sc.pl.violin(adata_mouse, ['Ar'], groupby = 'cluster', save = '_AR_mouse.png')

sc.pl.umap(adata_mouse, color='Ar', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Ar_mouse_all.png')


sc.pl.violin(adata_human_new, ['AR'], groupby = 'cluster', save = '_AR_human.png')

sc.pl.umap(adata_human_new, color=['AR'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Ar_human_all.png')

# violin for Periostin mouse
sc.pl.violin(adata_mouse, ['Postn'], groupby = 'cluster', save = '_Postn_mouse.png')

# umap for Periostin mouse
sc.pl.umap(adata_mouse, color= 'Postn', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_mouse_all.png')

# dotplot for periostin per different mouse models
dp = sc.pl.DotPlot(adata_mouse, var_names = 'Postn', groupby = 'key', cmap = 'Reds')
dp.legend(width=2.5).savefig('figures/dotplot_Postn_mouse.png')

# dotplot for Ar per different mouse models
dp = sc.pl.DotPlot(adata_mouse, var_names = 'Ar', groupby = 'key', cmap = 'Reds')
dp.legend(width=2.5).savefig('figures/dotplot_Ar_mouse.png')

# violin for Periostin human
sc.pl.violin(adata_human_new, ['POSTN'], groupby = 'cluster', save = '_Postn_human.png')

# umap for Periostin human
sc.pl.umap(adata_human_new, color='POSTN', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_human_all.png')

sc.pl.rank_genes_groups_dotplot(adata_mouse, var_names = ['Acta2', 'Myl9', 'Myh11', 'Tagln', 'Pdgfra', 'Mustn1', 'Angpt2', 'Notch3'], save = '_c0_mouse.png')

sc.pl.rank_genes_groups_dotplot(adata_human_new, var_names = ['ACTA2', 'MYL9', 'MYH11', 'TAGLN', 'PDGFRA', 'MUSTN1', 'ANGPT2', 'NOTCH3'], save = '_c0_human.png')

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


adata_human_raw = adata_human.raw.to_adata()
adata_human_raw.X = sp.csr_matrix.todense(adata_human_raw.X)
adata_human_raw.X = adata_human_raw.to_df()
adata_human_raw.write('outs/forCellChat/adata_human_norm.h5ad')

sc.tl.rank_genes_groups(adata_human_new, 'cluster', method='t-test')
markers_human_c0 = sc.get.rank_genes_groups_df(adata_human_new, group = '0')
markers_human_c1 = sc.get.rank_genes_groups_df(adata_human_new, group = '1')
markers_human_c2 = sc.get.rank_genes_groups_df(adata_human_new, group = '2')
markers_human_c3 = sc.get.rank_genes_groups_df(adata_human_new, group = '3')
markers_human_c4 = sc.get.rank_genes_groups_df(adata_human_new, group = '4')
markers_human_c5 = sc.get.rank_genes_groups_df(adata_human_new, group = '5')
markers_human_c6 = sc.get.rank_genes_groups_df(adata_human_new, group = '6')
markers_human_c7 = sc.get.rank_genes_groups_df(adata_human_new, group = '7')

#markers_human_c0.sort_values(by = ['logfoldchanges'], inplace=True, ascending = False)
#markers_human_c1.sort_values(by = ['logfoldchanges'], inplace=True, ascending = False)
#markers_human_c2.sort_values(by = ['logfoldchanges'], inplace=True, ascending = False)
#markers_human_c3.sort_values(by = ['logfoldchanges'], inplace=True, ascending = False)
#markers_human_c4.sort_values(by = ['logfoldchanges'], inplace=True, ascending = False)
#markers_human_c5.sort_values(by = ['logfoldchanges'], inplace=True, ascending = False)
#markers_human_c6.sort_values(by = ['logfoldchanges'], inplace=True, ascending = False)
#markers_human_c7.sort_values(by = ['logfoldchanges'], inplace=True, ascending = False)

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
        #"TangI",
        "Rgs5",
        "Mef2c",
        "Pdgfrb"
    ],
    cmap="RdBu_r",
    vmax=5,
    save="_myo_markers.png",
)

markers_human_c3['names'] = markers_human_c3['names'].str.title()

markers_mouse_c0 = sc.get.rank_genes_groups_df(adata_mouse, group = '0')
markers_mouse_c1 = sc.get.rank_genes_groups_df(adata_mouse, group = '1')
markers_mouse_c2 = sc.get.rank_genes_groups_df(adata_mouse, group = '2')
markers_mouse_c3 = sc.get.rank_genes_groups_df(adata_mouse, group = '3')
markers_mouse_c4 = sc.get.rank_genes_groups_df(adata_mouse, group = '4')
markers_mouse_c5 = sc.get.rank_genes_groups_df(adata_mouse, group = '5')
markers_mouse_c6 = sc.get.rank_genes_groups_df(adata_mouse, group = '6')
markers_mouse_c7 = sc.get.rank_genes_groups_df(adata_mouse, group = '7')

markers_human_c3['names'] = markers_human_c3['names'].str.title()

top50_mouse = markers_mouse_c3[0:100]
top50_human = markers_human_c3[0:100]

s1 = pd.merge(top50_mouse, top50_human, how='inner', on=['names'])
s1

markers_mouse_c5 = sc.get.rank_genes_groups_df(adata_mouse, group = '5')



