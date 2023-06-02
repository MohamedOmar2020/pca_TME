
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

import warnings
warnings.filterwarnings('ignore')

#############################################
# set figure parameters
sc.settings.figdir = 'figures/cellchat_allCompartments/new'
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
adata_all = sc.read_h5ad('data/mouse_all.h5ad')
adata_all.obs['cluster'].value_counts()

##########################################################
# dotplot common clusters in human
EpithC0C1C2_path = [
    # EGF
    'Areg', 'Hbegf', 'Ereg', 'Egfr',
    # EPHA
    'Efna1', 'Efna5', 'Epha2', 'Epha3', 'Epha4', 'Epha7', 'Ephb2',
    # EPHB
    'Efnb1', 'Efnb2', 'Epha4', 'Ephb2', 'Ephb4',
    # MIF
    'Mif', 'Ackr3', 'Cd74', 'Cxcr4', 'Cd44',
    # PDGF
    'Pdgfa', 'Pdgfc', 'Pdgfd', 'Pdgfra', 'Pdgfrb',
    # THBS
    'Thbs1', 'Thbs2', 'Thbs3', 'Sdc1', 'Sdc4', 'Cd47', 'Itga3', 'Itgav', 'Itgb1', 'Itgb3',
    # WNT
    'Wnt10a', 'Wnt2', 'Wnt4', 'Wnt6', 'Fzd1', 'Fzd2', 'Lrp6',
    # GAS
    'Gas6', 'Axl'
]

# dp = sc.pl.DotPlot(adata_all[adata_all.obs['cluster'].isin(['epithelium', 'c0', 'c1', 'c2'])], var_names = EpithC0C1C2_path,
#                                 categories_order = ['epithelium', 'c0', 'c1', 'c2'],
#                                 groupby='cluster', cmap = 'Reds', figsize = [16, 5]
#                                 )
#
# dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).legend(width = 1.2).savefig('figures/cellchat_allCompartments/new/EpithC0C1C2_path.png')

sc.pl.dotplot(adata_all[adata_all.obs['cluster'].isin(['epithelium', 'c0', 'c1', 'c2'])], EpithC0C1C2_path, 'cluster', swap_axes = False, dendrogram=False, save='EpithC0C1C2_path.png')

marker_genes_dict = {
    'B-cell': ['CD79A', 'MS4A1'],
    'Dendritic': ['FCER1A', 'CST3'],
    'Monocytes': ['FCGR3A'],
    'NK': ['GNLY', 'NKG7'],
    'Other': ['IGLL1'],
    'Plasma': ['IGJ'],
    'T-cell': ['CD3D'],
}

ax = sc.pl.heatmap(pbmc, marker_genes_dict, groupby='clusters', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r', dendrogram=True, swap_axes=True, figsize=(11,4))
