
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
# load the data
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
# 4A
###########################
# dotplot of GEMM-specific clusters

GEMMs_markers = ['Sfrp2', 'Wnt5a', 'Lgr5', 'Apc',
                'Wnt4', 'Wnt6', 'Notum', 'Wif1',
                'Nkd1', 'Wnt2', 'Wnt10a',
                'Dkk2', 'Rorb', 'Cxxc4', 'Nfat5',
                'Apoe', 'Dact1', 'Ctnnb1', 'Lef1',
                'Tcf4', 'Myc', 'Mki67', 'H2afx',
                'Top2a', 'Ccnb1', 'Ccnb2', 'Stmn1',
                'Ptn', 'Mdk', 'Tubb3', 'Mrc2', 'Fn1',
                'Tnc', 'Col12a1', 'Col14a1', 'Col16a1',
                'Mmp19', 'Cthrc1', 'Wisp1', 'Fzd1', 'Fzd2',
                'Sfrp4', 'Bmp1', 'Tle3', 'Tgfb1', 'Postn'
                 ]

dp2 = sc.pl.DotPlot(adata_mouse_mesenchyme, var_names = GEMMs_markers,
                                categories_order = ['c0','c1','c2','c3','c4','c5','c6','c7'],
                                groupby='cluster', cmap = 'Reds', figsize=[15, 3],
                                )

dp2.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/figures_cell/mouse_dotplot_specific.tiff')



################################
# violin plots
########################################
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import numpy as np

################
# for c3 and c4

c3c4 = ['Sfrp2', 'Lgr5', 'Wnt2', 'Wnt4', 'Wnt5a', 'Wnt6', 'Notum', 'Wif1']
c3c4_small = ['Sfrp2', 'Wif1']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c3c4_small), figsize=(12, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c3c4_small)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_mouse_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_mouse_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_mouse_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c3c4_expression = df[df['cluster'].isin(['c3','c4'])][gene]
    other_expression = df[~df['cluster'].isin(['c3','c4'])][gene]
    t_stat, p_value = stats.ttest_ind(c3c4_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    x_c3c4 = df['cluster'].unique().tolist().index('c3')
    y_max = df[gene].max()
    #Add the p-value to the plot
    if p_value < 0.0001:
       ax.text(3.5, y_max+0.2, "p < 0.0001", ha='center', fontsize = 14)
    else:
       ax.text(3.5, y_max+0.2, f"p = {p_value:.4f}", ha='center', fontsize = 14)

    #Adjust the y-axis limit to accommodate the p-value text
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
plt.savefig('figures/figures_cell/mouse_violin_c3c4_small.png')


################
# for c5-c7

c5c6c7 = ['Tgfb1', 'Postn']

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c5c6c7), figsize=(12, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c5c6c7)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_mouse_mesenchyme.raw[:, gene].X.todense(),
                                   index=adata_mouse_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_mouse_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c5c6c7_expression = df[df['cluster'].isin(['c5','c6', 'c7'])][gene]
    other_expression = df[~df['cluster'].isin(['c5','c6', 'c7'])][gene]
    t_stat, p_value = stats.ttest_ind(c5c6c7_expression, other_expression)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, scale='width')

    # Add the p-value to the plot
    x_c3c4 = df['cluster'].unique().tolist().index('c3')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(3.5, y_max+0.2, "p < 0.0001", ha='center', fontsize = 14)
    else:
        ax.text(3.5, y_max+0.2, f"p = {p_value:.4f}", ha='center', fontsize = 14)

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
plt.savefig('figures/figures_cell/mouse_violin_c5c6c7_small.png')



###########################
# 4B-C: histology
###########################

###########################
# 4D
###########################

# subset to c5, c6, and c7
# adata_mouse_c5c6c7 = adata_mouse_mesenchyme[adata_mouse_mesenchyme.obs['cluster'].isin(['c5','c6','c7'])]
# adata_mouse_c5c6c7.obs['cluster'].value_counts()
#
# ## plot Mki67 expression in the PRN clusters
# sc.pl.violin(adata_mouse_c5c6c7, ['Mki67'], groupby = 'cluster', save='_Postn+Mki67+.png')
#
#
# ## get cells that are Mki67+
# countMat_all = pd.concat([adata_mouse_mesenchyme.to_df(), adata_mouse_mesenchyme.obs], axis=1)
#
# mki67_cond_all = (countMat_all.Mki67 >= countMat_all['Mki67'].quantile(0.9))
#
# adata_mouse_mesenchyme.obs['Mki67_status'] = 'negative'
# adata_mouse_mesenchyme.obs.loc[mki67_cond_all, 'Mki67_status'] = 'positive'
# adata_mouse_mesenchyme.obs['Mki67_status'].value_counts()
#
# ## plot cells by Mki67 status
# sc.pl.umap(adata_mouse_mesenchyme, color=['Mki67', 'Mki67_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_all_Mki67Expression_Mki67status.png')
# sc.pl.umap(adata_mouse_c5c6c7, color=['Mki67', 'Mki67_status'], color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_c5c6c7_Mki67Expression_Mki67status.png')
#
# # Relative frequency barplot
# adata_mouse_mesenchyme.obs['cluster'] = adata_mouse_mesenchyme.obs['cluster'].astype('category')
# adata_mouse_mesenchyme.obs['Mki67_status'] = adata_mouse_mesenchyme.obs['Mki67_status'].astype('category')
#
# Mki67_relativeFrequency_all = sct.tools.relative_frequency_per_cluster(adata_mouse_mesenchyme, group_by='cluster', xlabel='Mki67_status', condition=None)
# sct.plot.cluster_composition_stacked_barplot(Mki67_relativeFrequency_all, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_mouse_mesenchyme.uns['Mki67_status_colors'], save = 'figures/figures_cell/Mki67_status.png')

############################