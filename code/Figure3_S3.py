
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
from openpyxl import load_workbook



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
# load the mouse stroma data
###########################
adata_mouse_mesenchyme = sc.read_h5ad('data/for_mouse/adata_mouse.h5ad', chunk_size=100000)
adata_mouse_mesenchyme.obs['cluster'] = adata_mouse_mesenchyme.obs['cluster'].astype('str')

############################
# load the mouse immune data
############################
mouse_immune = sc.read_h5ad('outs/h5ads/mouse_immune.h5ad', chunk_size=100000)


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
# 3A
###########################
# dotplot of common clusters
VarsOfInterest = ['Acta2', 'Myl9', 'Myh11', 'Tagln', 'Pdgfra',
                                             'Mustn1', 'Angpt2', 'Notch3', 'Sfrp1', 'Gpx3',
                                             'C3', 'C7', 'Cfh', 'Ccl11', 'Cd55', 'Ptx3',
                                             'Thbd', 'Ifi204', 'Ifi205', 'Ifi207', 'Jun',
                                             'Junb', 'Jund', 'Fos', 'Fosb', 'Fosl2', 'Atf3',
                                             'Mafb', 'Maff', 'Nek2', 'Id1', 'Id3', 'Btg2',
                                             'Gadd45a', 'Hes1', 'Bcl3', 'Socs1', 'Socs3',
                                             'Il6', 'Irf1', 'Map3k8', 'Gadd45b', 'Gadd45g',
                                             'Dusp1', 'Dusp6', 'Klf4'
                  ]

dp = sc.pl.DotPlot(adata_mouse_mesenchyme, var_names = VarsOfInterest,
                   figsize = [15, 3],
                   groupby='cluster', categories_order = ['c0','c1','c2','c3','c4','c5','c6','c7'], cmap = 'Reds'
                   )
dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/figures_cell/mouse_dotplot_common.tiff')


###########################
# 3b histology
###########################

###########################
# 3c and e: Refere to the R script
###########################

##########################################################
# Figure 3d stacked violin plot for immune cell types per mouse model

mouse_immune.obs['cluster'] = mouse_immune.obs['immune']
mouse_immune.obs['cluster'] = mouse_immune.obs['cluster'].astype('category')
mouse_immune.obs['cluster'].value_counts()

colMap = ['powderblue', 'plum', 'maroon', 'grey', 'darkseagreen', 'salmon']
relativeFrequency_immune = sct.tools.relative_frequency_per_cluster(mouse_immune, group_by='key_new', xlabel='cluster')
sct.plot.cluster_composition_stacked_barplot(relativeFrequency_immune, xlabel='key_new', figsize=(50, 40), width=0.7, margins=(0.02, 0.02), label_size=0, tick_size=50, colors=colMap, save = 'figures/immune_frequency_Per_GEMM.png')

# save to source data
with pd.ExcelWriter('tables/Source_data.xlsx', engine='openpyxl', mode='a') as writer:
    relativeFrequency_immune.to_excel(writer, sheet_name='3d')

