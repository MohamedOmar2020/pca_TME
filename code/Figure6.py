# load libraries
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
import plotly.express as px



#############################################
# set figure parameters
sc.settings.figdir = 'figures/figures_cell'
sc.set_figure_params(dpi_save = 300, transparent = False, fontsize =9, format='tiff')
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

###########################
# load the human primary data
###########################
adata_human_mesenchyme = sc.read_h5ad('data/for_human/adata_human.h5ad', chunk_size=100000)
adata_human_mesenchyme.obs['cluster'] = adata_human_mesenchyme.obs['cluster'].astype('str')

###########################
# load the human bone metastatic data
adata_human_bone = sc.read_h5ad('data/for_human/adata_human_BoneMet_stroma.h5ad', chunk_size=100000)
adata_human_bone.obs['cluster'] = adata_human_bone.obs['cluster'].astype('str')

###########################
# 6A
###########################
# parallel categories plot for clusters per ERG status
adata_human_mesenchyme.obs['cluster'] = adata_human_mesenchyme.obs['cluster'].astype('category')
adata_human_mesenchyme.obs.sort_values(by = ['cluster'], inplace=True)

color_map = dict(
    zip(adata_human_mesenchyme.obs['cluster'].cat.categories.tolist(), sc.pl.palettes.vega_20_scanpy)
)
adata_human_mesenchyme.obs["color"] = adata_human_mesenchyme.obs['cluster'].map(color_map)

fig = px.parallel_categories(
    adata_human_mesenchyme.obs,
    dimensions=['cluster', 'erg'],
    color="color",
    labels={'erg': "ERG status", 'cluster': "cluster"},
)
fig.update_layout(autosize=True, width =400, height=600, font_size = 9, font_family='Arial Black')
fig.update_yaxes(tickfont_family="Arial Black")
fig.write_image("figures/figures_cell/human_parallel_categories_cluster.svg", scale = 2)

# umap of clusters
sc.pl.umap(adata_human_mesenchyme, color=['cluster'], save = '_human_clusters')

# umap of AR expression
sc.pl.umap(adata_human_mesenchyme, color='AR', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_AR')

###########################
# 6B
###########################
# dotplot common clusters in human
dp = sc.pl.DotPlot(adata_human_mesenchyme, var_names = ['ACTA2', 'MYL9', 'MYH11', 'TAGLN',
                                                              'PDGFRA', 'MUSTN1', 'ANGPT2', 'NOTCH3',
                                                              'SFRP1', 'GPX3', 'C3', 'C7', 'CFH', 'CCL11',
                                                              'CD55', 'PTX3', 'THBD', 'IFI16', 'JUN', 'JUNB',
                                                              'JUND', 'FOS', 'FOSB', 'FOSL2', 'ATF3',
                                                              'MAFB', 'MAFF', 'NEK2', 'ID1', 'ID3', 'BTG2',
                                                              'GADD45A', 'HES1', 'BCL3', 'SOCS1', 'SOCS3',
                                                              'IL6', 'IRF1', 'MAP3K8', 'GADD45B', 'GADD45G',
                                                              'DUSP1', 'DUSP6', 'KLF4'],
                                categories_order = ['c0','c1','c2','c3','c4','c5','c6','c7'],
                                groupby='cluster', cmap = 'Reds', figsize = [15, 3]
                                )

dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).legend(width = 1.2).savefig('figures/figures_cell/human_dotplot_common.tiff')

###########################
# 6C
###########################
## umap of transferred clusters and original cell types
sc.pl.umap(adata_human_bone, color='cluster', save = '_kfoury_projectedClusters')
sc.pl.umap(adata_human_bone, color='cells',  save = '_kfoury_originalCells')

## umap of BGN expression
sc.pl.umap(adata_human_bone, color='BGN', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_kfoury_BGN')

###########################
# 6D
###########################
# violin plot of MKI67
sc.pl.violin(adata_human_bone, ['MKI67'], groupby = 'cluster', use_raw=False, save='_kfoury_MKI67')

# violin plot of POSTN
sc.pl.violin(adata_human_bone, ['POSTN'], groupby = 'cluster', use_raw=False, save='_kfoury_POSTN')
