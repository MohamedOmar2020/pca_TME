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
from scipy import stats


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
adata_human_mesenchyme.obs['cluster'].value_counts()

######################################################
# Identify cells that are annotated as c5, c6, or c7
######################################################
selected_clusters = ['c5', 'c6', 'c7']
selected_cells_mask = adata_human_mesenchyme.obs['cluster'].isin(selected_clusters)

# Calculate the median expression for JUN and FOS across all cells
jun_median = np.median(adata_human_mesenchyme[:, 'JUN'].X)
fos_median = np.median(adata_human_mesenchyme[:, 'FOS'].X)

# Get the JUN and FOS expression for the selected cells
jun_values_selected = adata_human_mesenchyme[selected_cells_mask, 'JUN'].X
fos_values_selected = adata_human_mesenchyme[selected_cells_mask, 'FOS'].X

# Identify cells within the selected clusters that have JUN and FOS expression above median values
above_median_jun_fos_mask = (jun_values_selected > jun_median) & (fos_values_selected > fos_median)

# Combine the masks to identify the final selection of cells to re-assign
final_selection_mask = selected_cells_mask.copy()
final_selection_mask[selected_cells_mask] = above_median_jun_fos_mask

# Re-assign these cells to cluster c2
adata_human_mesenchyme.obs.loc[final_selection_mask, 'cluster'] = 'c2'
