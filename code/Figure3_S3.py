
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

################################
# violin plots
########################################
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import numpy as np

################
# for c0

c0 = ['Acta2', 'Myl9', 'Myh11', 'Tagln', 'Pdgfra', 'Mustn1', 'Angpt2', 'Notch3']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c0), figsize=(35, 5))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c0)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_mouse_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_mouse_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_mouse_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c0_expression = df[df['cluster'] == 'c0'][gene]
    other_expression = df[df['cluster'] != 'c0'][gene]
    t_stat, p_value = stats.ttest_ind(c0_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    x_c0 = df['cluster'].unique().tolist().index('c0')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    else:
        ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=14)
    ax.set_ylabel(gene, fontsize=14)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/mouse_violin_c0_small.png')

############
## For c1

c1 = ['Sfrp1', 'Gpx3', 'C3', 'C7', 'Cfh', 'Ccl11', 'Cd55', 'Ptx3', 'Thbd', 'Ifi204', 'Ifi205', 'Ifi207']
c1_small = ['Sfrp1', 'Gpx3', 'C3', 'C7', 'Ccl11', 'Cd55', 'Ptx3', 'Ifi204']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c1_small), figsize=(35, 5))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c1_small)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_mouse_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_mouse_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_mouse_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c1_expression = df[df['cluster'] == 'c1'][gene]
    other_expression = df[df['cluster'] != 'c1'][gene]
    t_stat, p_value = stats.ttest_ind(c1_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    x_c1 = df['cluster'].unique().tolist().index('c1')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    else:
        ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=14)
    ax.set_ylabel(gene, fontsize=14)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)


# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/mouse_violin_c1_small.png')

############
## For c2

c2 = [
    'Jun', 'Junb', 'Jund', 'Fos', 'Fosb', 'Fosl2', 'Atf3',
        'Mafb', 'Maff', 'Nek2', 'Id1', 'Id3', 'Btg2',
        'Gadd45a', 'Hes1', 'Bcl3', 'Socs1', 'Socs3',
        'Il6', 'Irf1', 'Map3k8', 'Gadd45b', 'Gadd45g',
        'Dusp1', 'Dusp6', 'Klf4'
      ]

c2_small = ['Jun', 'Fos', 'Atf3', 'Mafb', 'Id1', 'Gadd45g', 'Dusp1', 'Klf4']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c2_small), figsize=(35, 5))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c2_small)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_mouse_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_mouse_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_mouse_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c2_expression = df[df['cluster'] == 'c2'][gene]
    other_expression = df[df['cluster'] != 'c2'][gene]
    t_stat, p_value = stats.ttest_ind(c2_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    x_c2 = df['cluster'].unique().tolist().index('c2')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    else:
        ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=14)
    ax.set_ylabel(gene, fontsize=14)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/mouse_violin_c2_small.png')

















###########################
# 3B-C: histology
###########################

###########################
# 3D-E: Refere to the R script
###########################
