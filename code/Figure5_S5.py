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
# load the mouse data
###########################
adata_mouse_mesenchyme = sc.read_h5ad('data/for_mouse/adata_mouse.h5ad', chunk_size=100000)
adata_mouse_mesenchyme.obs['cluster'] = adata_mouse_mesenchyme.obs['cluster'].astype('str')

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
# 5B
###########################

# Postn expression / PRN clusters
ax = sc.pl.umap(adata_mouse_mesenchyme, color=['cluster', 'Postn', 'Ar'], title = ['PRN clusters', 'Postn expression', 'Ar epression'], groups = ['c5', 'c6', 'c7'], na_in_legend = False, color_map = 'RdBu_r', vmin='p1', vmax='p99', show = False)
ax[0].title.set_size(12)
ax[1].title.set_size(12)
ax[2].title.set_size(12)
ax[0].title.set_fontweight('bold')
ax[1].title.set_fontweight('bold')
ax[2].title.set_fontweight('bold')
plt.savefig('figures/figures_cell/umap_Postn_expression_PRN_clusters.tiff', dpi=200, bbox_inches='tight')

# dotplot for periostin per different mouse models
dp = sc.pl.DotPlot(adata_mouse_mesenchyme, var_names = 'Postn', groupby = 'key_new', cmap = 'Reds')
dp.legend(width=2.5).savefig('figures/figures_cell/dotplot_Postn_mouseModels.tiff')

#######
# dotplot for Ar per different mouse models
dp = sc.pl.DotPlot(adata_mouse_mesenchyme, var_names = 'Ar', groupby = 'key_new', cmap = 'Reds')
dp.legend(width=2.5).savefig('figures/figures_cell/dotplot_Ar_mouseModels.tiff')

############################################################################################################
# Figure S5

###########################
# S5A
###########################

# violin plots
ax = sc.pl.violin(adata_mouse_mesenchyme, ['C1qa', 'C1qb', 'C1qc'], groupby = 'cluster', use_raw=True, show=False)
ax[0].title.set_size(12)
ax[1].title.set_size(12)
ax[2].title.set_size(12)
ax[0].title.set_fontweight('bold')
ax[1].title.set_fontweight('bold')
ax[2].title.set_fontweight('bold')

plt.savefig('figures/figures_cell/violin_mouse_complement.tiff', dpi=400, bbox_inches='tight')

# dotplot
dp = sc.pl.DotPlot(adata_mouse_mesenchyme, var_names = ['C1qa', 'C1qb', 'C1qc'], groupby = 'key_new', cmap = 'Reds')
dp.legend(width=2.5).savefig('figures/figures_cell/C1Q_GEMMs.png')

###########################
# S5B
###########################

ax = sc.pl.violin(adata_human_mesenchyme, ['C1QA', 'C1QB', 'C1QC'], groupby = 'cluster', use_raw=True, show = False)

ax[0].title.set_size(12)
ax[1].title.set_size(12)
ax[2].title.set_size(12)

ax[0].title.set_fontweight('bold')
ax[1].title.set_fontweight('bold')
ax[2].title.set_fontweight('bold')

plt.savefig('figures/figures_cell/violin_human_complement.tiff', dpi=400, bbox_inches='tight')

###########################
# S5C
###########################
ax = sc.pl.violin(adata_human_bone, ['C1QA', 'C1QB', 'C1QC'], groupby = 'cluster', use_raw=True, show=False)

ax[0].title.set_size(12)
ax[1].title.set_size(12)
ax[2].title.set_size(12)

ax[0].title.set_fontweight('bold')
ax[1].title.set_fontweight('bold')
ax[2].title.set_fontweight('bold')

plt.savefig('figures/figures_cell/violin_BoneMets_complement.tiff', dpi=400, bbox_inches='tight')
