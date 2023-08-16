
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
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cmasher as cmr



#############################################
# set figure parameters
sc.settings.figdir = 'figures/figures_cell'
sc.set_figure_params(dpi_save = 300, transparent = False, fontsize =9, format='png')
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
# load the data
###########################
adata_mouse_mesenchyme = sc.read_h5ad('data/for_mouse/adata_mouse.h5ad', chunk_size=100000)
adata_mouse_mesenchyme.obs['cluster'] = adata_mouse_mesenchyme.obs['cluster'].astype('str')

##########################################################
# fix GEMMs names
adata_mouse_mesenchyme.obs['key_new'] = adata_mouse_mesenchyme.obs['key']
adata_mouse_mesenchyme.obs['key_new'].value_counts()
adata_mouse_mesenchyme.obs['key_new'].replace('terg', 'T-ERG', inplace=True)
adata_mouse_mesenchyme.obs['key_new'].replace('himyc', 'Hi-MYC', inplace=True)
adata_mouse_mesenchyme.obs['key_new'].replace('fvbn', 'FVBN', inplace=True)
adata_mouse_mesenchyme.obs['key_new'].replace('pten', 'NP', inplace=True)
adata_mouse_mesenchyme.obs['key_new'].replace('mycn', 'PRN', inplace=True)
adata_mouse_mesenchyme.obs['key_new'].replace('129b6', 'B6.129', inplace=True)
adata_mouse_mesenchyme.obs['key_new'].replace('b6', 'B6', inplace=True)
adata_mouse_mesenchyme.obs['key_new'].replace('129b6_pten', 'WT for NP', inplace=True)
adata_mouse_mesenchyme.obs['key_new'].replace('129b6_mycn', 'WT for PRN', inplace=True)

###########################
# 1A
###########################
# umap clusters
sc.pl.umap(adata_mouse_mesenchyme, color = 'cluster', save = '_mouse_clusters', title = '')


#############
## umap by genotype
adata_mouse_mesenchyme.obs.condition.value_counts()

# plot mutants
sc.pl.umap(
    adata_mouse_mesenchyme[adata_mouse_mesenchyme.obs["condition"] == "mutant"],
    color="key_new", size = 10, title = '', palette = 'viridis',
    save = '_mouse_models_mutants'
)

# plot Wildtype

sc.pl.umap(
    adata_mouse_mesenchyme[adata_mouse_mesenchyme.obs["condition"] == "wildtype"],
    color="key_new", size = 10, title = '',
    save = '_mouse_models_wildtype'
)


###########################
# 1B
###########################
# parallel categories mouse models

# Arrange the cell annotation dataframe by the clusters
adata_mouse_mesenchyme.obs.sort_values(by = ['cluster'], inplace=True)

adata_mouse_mesenchyme.obs['cluster'] = adata_mouse_mesenchyme.obs['cluster'].astype('category')
color_map = dict(
    zip(adata_mouse_mesenchyme.obs['cluster'].cat.categories.tolist(), sc.pl.palettes.vega_20_scanpy)
)
adata_mouse_mesenchyme.obs["color"] = adata_mouse_mesenchyme.obs['cluster'].map(color_map)
fig = px.parallel_categories(
    adata_mouse_mesenchyme.obs,
    dimensions=['cluster', 'key_new'],
    color="color",
    labels={'key_new': "model", 'cluster': "cluster"},
)
fig.update_layout(autosize=False, width=400, height=600, font_size = 9, font_family="Arial Black")
fig.update_yaxes(tickfont_family="Arial Black")
fig.write_image("figures/figures_cell/parallel_categories_cluster.png", scale = 2)

###############################
# 1C-D
###############################
# refer to the R script


#############################################