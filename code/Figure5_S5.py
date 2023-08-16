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
################################
# violin plots
########################################
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import numpy as np

################################################################
# for c0
################################################################

# mouse
common = ['Acta2', 'Myl9', 'Myh11', 'Tagln', 'Jun', 'Fos']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(common), figsize=(25, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, common)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_mouse_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_mouse_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_mouse_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    #c0_expression = df[df['cluster'] == 'c0'][gene]
    #other_expression = df[df['cluster'] != 'c0'][gene]
    #t_stat, p_value = stats.ttest_ind(c0_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    #x_c0 = df['cluster'].unique().tolist().index('c0')
    #y_max = df[gene].max()
    # Add the p-value to the plot
    #if p_value < 0.0001:
    #    ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    #else:
    #    ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    #ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=16)
    ax.set_ylabel(gene, fontsize=16)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/mouse_violin_commonClusters.png')

######################
# human
common_human = ['ACTA2', 'MYL9', 'MYH11', 'TAGLN', 'JUN', 'FOS']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(common_human), figsize=(25, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, common_human)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_human_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    #c0_expression = df[df['cluster'] == 'c0'][gene]
    #other_expression = df[df['cluster'] != 'c0'][gene]
    #t_stat, p_value = stats.ttest_ind(c0_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    #x_c0 = df['cluster'].unique().tolist().index('c0')
    #y_max = df[gene].max()
    # Add the p-value to the plot
    #if p_value < 0.0001:
    #    ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    #else:
    #    ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    #ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=16)
    ax.set_ylabel(gene, fontsize=16)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_violin_commonClusters.png')

#################
## human bone metastasis
# Create a figure and axes objects
fig, axs = plt.subplots(1, len(common_human), figsize=(25, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, common_human)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_bone.raw[:, gene].X,
                                   index=adata_human_bone.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_bone.obs)

    # Perform a statistical test (e.g., t-test)
    #c0_expression = df[df['cluster'] == 'c0'][gene]
    #other_expression = df[df['cluster'] != 'c0'][gene]
    #t_stat, p_value = stats.ttest_ind(c0_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    #x_c0 = df['cluster'].unique().tolist().index('c0')
    #y_max = df[gene].max()
    # Add the p-value to the plot
    #if p_value < 0.0001:
    #    ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    #else:
    #    ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    #ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=16)
    ax.set_ylabel(gene, fontsize=16)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_boneMets_violin_common_clusters.png')


######################################################
## For c3-c4
######################################################

# mouse
c3c4 = ['Sfrp2', 'Notum', 'Wnt4', 'Wnt5a', 'Rorb', 'Tcf4']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c3c4), figsize=(25, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c3c4)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_mouse_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_mouse_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_mouse_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    #c2_expression = df[df['cluster'] == 'c2'][gene]
    #other_expression = df[df['cluster'] != 'c2'][gene]
    #t_stat, p_value = stats.ttest_ind(c2_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    #x_c2 = df['cluster'].unique().tolist().index('c2')
    #y_max = df[gene].max()
    # Add the p-value to the plot
    #if p_value < 0.0001:
    #    ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    #else:
    #    ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=16)
    ax.set_ylabel(gene, fontsize=16)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/mouse_violin_c3c4_clusters.png')

#############
# human: primary tumor:
c3c4_human = ['SFRP2', 'NOTUM', 'WNT4', 'WNT5A', 'RORB', 'TCF4']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c3c4), figsize=(25, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c3c4_human)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_human_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    #c2_expression = df[df['cluster'] == 'c2'][gene]
    #other_expression = df[df['cluster'] != 'c2'][gene]
    #t_stat, p_value = stats.ttest_ind(c2_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    #x_c2 = df['cluster'].unique().tolist().index('c2')
    #y_max = df[gene].max()
    # Add the p-value to the plot
    #if p_value < 0.0001:
    #    ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    #else:
    #    ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=16)
    ax.set_ylabel(gene, fontsize=16)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_violin_c3c4_clusters.png')

#############
# human: bone mets

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c3c4_human), figsize=(25, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c3c4_human)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_bone.raw[:, gene].X,
                                   index=adata_human_bone.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_bone.obs)

    # Perform a statistical test (e.g., t-test)
    #c2_expression = df[df['cluster'] == 'c2'][gene]
    #other_expression = df[df['cluster'] != 'c2'][gene]
    #t_stat, p_value = stats.ttest_ind(c2_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    #x_c2 = df['cluster'].unique().tolist().index('c2')
    #y_max = df[gene].max()
    # Add the p-value to the plot
    #if p_value < 0.0001:
    #    ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    #else:
    #    ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=16)
    ax.set_ylabel(gene, fontsize=16)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/humanBoneMets_violin_c3c4_clusters.png')



######################################################
## For c5-c7
######################################################

# mouse
c5c6c7 = ['Postn', 'Fn1', 'Sfrp4', 'Tnc', 'Tgfb1', 'Cthrc1']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c5c6c7), figsize=(25, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c5c6c7)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_mouse_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_mouse_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_mouse_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    #c2_expression = df[df['cluster'] == 'c2'][gene]
    #other_expression = df[df['cluster'] != 'c2'][gene]
    #t_stat, p_value = stats.ttest_ind(c2_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, cut = 0.5, ax=ax, scale='width')

    # Add the p-value to the plot
    #x_c2 = df['cluster'].unique().tolist().index('c2')
    #y_max = df[gene].max()
    # Add the p-value to the plot
    #if p_value < 0.0001:
    #    ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    #else:
    #    ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=16)
    ax.set_ylabel(gene, fontsize=16)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/mouse_violin_c5c6c7_clusters.png')

#############
# human: primary tumor:
c5c6c7_human = ['POSTN', 'FN1', 'SFRP4', 'TNC', 'TGFB1', 'CTHRC1']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c3c4), figsize=(25, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c5c6c7_human)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_human_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    #c2_expression = df[df['cluster'] == 'c2'][gene]
    #other_expression = df[df['cluster'] != 'c2'][gene]
    #t_stat, p_value = stats.ttest_ind(c2_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], cut = 0.5, data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    #x_c2 = df['cluster'].unique().tolist().index('c2')
    #y_max = df[gene].max()
    # Add the p-value to the plot
    #if p_value < 0.0001:
    #    ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    #else:
    #    ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=16)
    ax.set_ylabel(gene, fontsize=16)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_violin_c5c6c7_clusters.png')

#############
# human: bone mets

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c5c6c7_human), figsize=(25, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c3c4_human)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_bone.raw[:, gene].X,
                                   index=adata_human_bone.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_bone.obs)

    # Perform a statistical test (e.g., t-test)
    #c2_expression = df[df['cluster'] == 'c2'][gene]
    #other_expression = df[df['cluster'] != 'c2'][gene]
    #t_stat, p_value = stats.ttest_ind(c2_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    #x_c2 = df['cluster'].unique().tolist().index('c2')
    #y_max = df[gene].max()
    # Add the p-value to the plot
    #if p_value < 0.0001:
    #    ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 14)
    #else:
    #    ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=16)
    ax.set_ylabel(gene, fontsize=16)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/humanBoneMets_violin_c5c6c7_clusters.png')
