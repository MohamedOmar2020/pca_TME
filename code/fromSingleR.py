
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

sc.settings.figdir = 'singleR'
sc.set_figure_params(dpi_save = 300)
#plt.rcParams.update({'xtick.labelsize' : '50'})
#plt.rcParams.update({'ytick.labelsize' : '50'})

adata_mouse = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)
adata_mouse.obs['cluster'] = adata_mouse.obs['cluster'].astype('str')
adata_mouse.obs['cluster'].value_counts()

# filter the na
adata_mouse = adata_mouse[adata_mouse.obs['cluster'] != 'nan', :]
adata_mouse.obs['cluster'].value_counts()


###########################################################################
# load the annot human data
#%%
adata_human_singleR = sc.read_h5ad('human/human_pred.h5ad')

# fix the var_names
adata_human_singleR.var_names = adata_human_singleR.var['features']

adata_human_singleR_raw = adata_human_singleR.raw.to_adata()
adata_human_singleR_raw.X = sp.csr_matrix.todense(adata_human_singleR_raw.X)
adata_human_singleR_raw.X = adata_human_singleR_raw.to_df()
del adata_human_singleR_raw.obs['leiden']
#adata_human_raw.write('outs/forCellChat/adata_human_norm.h5ad')

#sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata_human_singleR_raw, min_cells=3)
sc.pp.normalize_total(adata_human_singleR_raw, target_sum=1e4)
sc.pp.log1p(adata_human_singleR_raw)
sc.pp.highly_variable_genes(adata_human_singleR_raw, min_mean=0.0125, max_mean=3, min_disp=0.5)
#adata_human_singleR_raw = adata_human_singleR_raw[:, adata_human_singleR_raw.var.highly_variable]
sc.pp.scale(adata_human_singleR_raw, max_value=10)


#sc.pp.combat(adata_human_singleR, key='case', covariates=['condition', 'erg'])

# PCA and umap
sc.tl.pca(adata_human_singleR_raw, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_human_singleR_raw, log=True)
sc.pp.neighbors(adata_human_singleR_raw, n_neighbors=10, n_pcs=40)
sc.tl.paga(adata_human_singleR_raw, groups='labels')
sc.pl.paga(adata_human_singleR_raw, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata_human_singleR_raw, init_pos='paga')

sc.pl.umap(adata_human_singleR_raw, color='labels', save = '_human.png')


sc.pl.umap(adata_mouse, color=['cluster', 'Acta2'], save = '_mouse_ACTA2.png')

sc.pl.umap(adata_human_singleR, color=['cluster', 'ACTA2'], save = '_human_ACTA2.png')

sc.tl.rank_genes_groups(adata_mouse, 'cluster', pts=True, use_raw = False)
sc.pl.rank_genes_groups(adata_mouse, n_genes = 25, sharey=False, save = '_mouseMarkers_notRaw.png')
markers_mouse_c3 = sc.get.rank_genes_groups_df(adata_mouse, group = '3')


sc.tl.rank_genes_groups(adata_human_singleR, 'cluster', pts=True, use_raw = False)
sc.pl.rank_genes_groups(adata_human_singleR, n_genes=25, sharey=False, save = '_humanMarkers_notRaw.png')

sc.pl.violin(adata_mouse, ['Ar'], groupby = 'cluster', save = '_AR_mouse.png')

sc.pl.umap(adata_mouse, color='Ar', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Ar_mouse_all.png')


sc.pl.violin(adata_human_singleR, ['AR'], groupby = 'cluster', save = '_AR_human.png')

sc.pl.umap(adata_human_singleR, color=['AR'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Ar_human_all.png')

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
sc.pl.violin(adata_human_singleR, ['POSTN'], groupby = 'cluster', save = '_Postn_human.png')

# umap for Periostin human
sc.pl.umap(adata_human_singleR, color='POSTN', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_Postn_human_all.png')

sc.pl.rank_genes_groups_dotplot(adata_mouse, var_names = ['Acta2', 'Myl9', 'Myh11', 'Tagln', 'Pdgfra', 'Mustn1', 'Angpt2', 'Notch3'], save = '_c0_mouse.png')

sc.pl.rank_genes_groups_dotplot(adata_human_singleR, var_names = ['ACTA2', 'MYL9', 'MYH11', 'TAGLN', 'PDGFRA', 'MUSTN1', 'ANGPT2', 'NOTCH3'], save = '_c0_human.png')

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


dp = sc.pl.DotPlot(adata_human_singleR, var_names = ['ACTA2', 'MYL9', 'MYH11', 'TAGLN',
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


dp2 = sc.pl.DotPlot(adata_human_singleR,
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

#sc.tl.rank_genes_groups(adata_human_singleR, 'cluster', method='t-test')
markers_human_c0 = sc.get.rank_genes_groups_df(adata_human_singleR, group = '0')
markers_human_c1 = sc.get.rank_genes_groups_df(adata_human_singleR, group = '1')
markers_human_c2 = sc.get.rank_genes_groups_df(adata_human_singleR, group = '2')
markers_human_c3 = sc.get.rank_genes_groups_df(adata_human_singleR, group = '3')
markers_human_c4 = sc.get.rank_genes_groups_df(adata_human_singleR, group = '4')
markers_human_c5 = sc.get.rank_genes_groups_df(adata_human_singleR, group = '5')
markers_human_c6 = sc.get.rank_genes_groups_df(adata_human_singleR, group = '6')
markers_human_c7 = sc.get.rank_genes_groups_df(adata_human_singleR, group = '7')

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

markers_human_c0['names'] = markers_human_c0['names'].str.title()
markers_human_c1['names'] = markers_human_c1['names'].str.title()
markers_human_c2['names'] = markers_human_c2['names'].str.title()
markers_human_c3['names'] = markers_human_c3['names'].str.title()
markers_human_c4['names'] = markers_human_c4['names'].str.title()
markers_human_c5['names'] = markers_human_c5['names'].str.title()
markers_human_c6['names'] = markers_human_c6['names'].str.title()
markers_human_c7['names'] = markers_human_c7['names'].str.title()

markers_mouse_c0 = sc.get.rank_genes_groups_df(adata_mouse, group = '0')
markers_mouse_c1 = sc.get.rank_genes_groups_df(adata_mouse, group = '1')
markers_mouse_c2 = sc.get.rank_genes_groups_df(adata_mouse, group = '2')
markers_mouse_c3 = sc.get.rank_genes_groups_df(adata_mouse, group = '3')
markers_mouse_c4 = sc.get.rank_genes_groups_df(adata_mouse, group = '4')
markers_mouse_c5 = sc.get.rank_genes_groups_df(adata_mouse, group = '5')
markers_mouse_c6 = sc.get.rank_genes_groups_df(adata_mouse, group = '6')
markers_mouse_c7 = sc.get.rank_genes_groups_df(adata_mouse, group = '7')


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



x = pd.merge(top100_mouse_c0, top100_human_c3, how='inner', on=['names'])
