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
sc.set_figure_params(dpi_save = 300, transparent = False, fontsize =7, format='tiff')
plt.rcParams["font.family"] = "arial"
plt.rcParams['font.size'] = 7
plt.rcParams['font.style'] = 'italic'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 7
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

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
# 7A
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
# 7B
###########################
# dotplot common clusters in human
dp = sc.pl.DotPlot(adata_human_mesenchyme, var_names = ['ACTA2', 'MYL9', 'JUN', 'JUNB', 'FOS', 'DUSP1',
                                                        'SFRP2', 'RORB', 'TCF4', 'WNT6',
                                                        'POSTN', 'FN1', 'SFRP4', 'TNC', 'CTHRC1'
                                                        ],
                                categories_order = ['c0','c1','c2','c3','c4','c5','c6','c7'],
                                groupby='cluster', cmap = 'Reds', figsize = [15, 3],
                                #log=True,
                                standard_scale = 'var',
                                #vmin = 0.15, vmax = 0.9,
                                expression_cutoff = 0.5,
                                layer = 'denoised',
                                use_raw=True,
                                )

dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5, dot_max = 1).legend(width = 1.2).savefig('figures/figures_cell/human_dotplot_common2.tiff')


markers = ['ACTA2', 'MYL9', 'JUN', 'JUNB', 'FOS', 'DUSP1',
            'SFRP2', 'RORB', 'TCF4', 'WNT6', 'POSTN', 'FN1', 'SFRP4', 'TNC', 'CTHRC1']

sc.pl.matrixplot(adata_human_mesenchyme,
                 markers,
                 groupby='cluster',
                 dendrogram=False,
                 categories_order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'],
                 standard_scale = 'var',
                 layer = 'denoised',
                 save = 'humanMarkers')


##################
# violin plots
human_markers = ['ACTA2', 'MYL9', 'GPX3', 'C7', 'JUN', 'FOS',
           'WNT5A', 'RORB', 'POSTN', 'SFRP4']

for i in human_markers:
    sc.pl.violin(adata_human_mesenchyme, i, groupby='cluster', use_raw=True, save=f'_human_common_{i}')

###########################
# DE results
sc.tl.rank_genes_groups(adata_human_mesenchyme, 'cluster', method='t-test')
results = adata_human_mesenchyme.uns['rank_genes_groups']

# violin plots with p-value
c0_markers = ['ACTA2', 'MYL9']
c1_markers = ['GPX3', 'C7']
c2_markers = ['JUN', 'FOS']
c3c4_markers = ['WNT4', 'RORB']
c5c6c7_markers = ['POSTN', 'SFRP4']

#######
## c0
#####
# check if the variance is normal
from scipy.stats import levene
gene_expression_ACTA2 = pd.DataFrame(adata_human_mesenchyme.raw[:, 'ACTA2'].X.toarray(),
                               index=adata_human_mesenchyme.obs_names,
                               columns=['ACTA2'])

df_ACTA2 = gene_expression_ACTA2.join(adata_human_mesenchyme.obs)

c0_expression_ACTA2 = df_ACTA2[df_ACTA2['cluster'].isin(['c0'])]['ACTA2']
other_expression_ACTA2 = df_ACTA2[~df_ACTA2['cluster'].isin(['c0'])]['ACTA2']

stat_ACTA2, p_value_ACTA2 = levene(c0_expression_ACTA2, other_expression_ACTA2)

# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c0_markers), figsize=(7, 4))
# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c0_markers)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_mesenchyme.raw[:, gene].X.toarray(),
                                   index=adata_human_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c0_expression = df[df['cluster'].isin(['c0'])][gene]
    other_expression = df[~df['cluster'].isin(['c0'])][gene]
    t_stat, p_value = stats.ttest_ind(c0_expression, other_expression, equal_var=False)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, cut=0.1, scale='width')

    # Add the p-value to the plot
    x_c0 = df['cluster'].unique().tolist().index('c0')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(0.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 7)
    else:
        ax.text(0.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 7)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=7)
    ax.set_ylabel(gene, fontsize=7)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_violin_c0Markers.png')

############
## c1
# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c1_markers), figsize=(7, 4))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c1_markers)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_mesenchyme.raw[:, gene].X.toarray(),
                                   index=adata_human_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c1_expression = df[df['cluster'].isin(['c1'])][gene]
    other_expression = df[~df['cluster'].isin(['c1'])][gene]
    t_stat, p_value = stats.ttest_ind(c1_expression, other_expression, equal_var=False)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, cut=0.1, scale='width')

    # Add the p-value to the plot
    x_c1 = df['cluster'].unique().tolist().index('c1')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(1, y_max+0.1, "p < 0.0001", ha='center', fontsize = 7)
    else:
        ax.text(1, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 7)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=7)
    ax.set_ylabel(gene, fontsize=7)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_violin_c1Markers.png')

############
## c2
# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c2_markers), figsize=(7, 4))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c2_markers)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_mesenchyme.raw[:, gene].X.toarray(),
                                   index=adata_human_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c2_expression = df[df['cluster'].isin(['c2'])][gene]
    other_expression = df[~df['cluster'].isin(['c2'])][gene]
    t_stat, p_value = stats.ttest_ind(c2_expression, other_expression, equal_var=False)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, cut=0.1, scale='width')

    # Add the p-value to the plot
    x_c2 = df['cluster'].unique().tolist().index('c2')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(2, y_max+0.1, "p < 0.0001", ha='center', fontsize = 7)
    else:
        ax.text(2, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 7)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=7)
    ax.set_ylabel(gene, fontsize=7)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_violin_c2Markers.png')

############
## c3-c4
# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c3c4_markers), figsize=(7, 4))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c3c4_markers)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_mesenchyme.raw[:, gene].X.toarray(),
                                   index=adata_human_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c3c4_expression = df[df['cluster'].isin(['c3', 'c4'])][gene]
    other_expression = df[~df['cluster'].isin(['c3', 'c4'])][gene]
    t_stat, p_value = stats.ttest_ind(c3c4_expression, other_expression, equal_var=False)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, cut=0.1, scale='width')

    # Add the p-value to the plot
    x_c3 = df['cluster'].unique().tolist().index('c3')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(4, y_max+0.1, "p < 0.0001", ha='center', fontsize = 7)
    else:
        ax.text(4, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 7)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=7)
    ax.set_ylabel(gene, fontsize=7)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_violin_c3c4Markers.png')

############
## c5-c7
# Create a figure and axes objects
fig, axs = plt.subplots(1, len(c5c6c7_markers), figsize=(7, 4))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, c5c6c7_markers)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_mesenchyme.raw[:, gene].X.toarray(),
                                   index=adata_human_mesenchyme.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_mesenchyme.obs)

    # Perform a statistical test (e.g., t-test)
    c5c6c7_expression = df[df['cluster'].isin(['c5', 'c6', 'c7'])][gene]
    other_expression = df[~df['cluster'].isin(['c5', 'c6', 'c7'])][gene]
    t_stat, p_value = stats.ttest_ind(c5c6c7_expression, other_expression, equal_var=False)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, cut=0.1, scale='width')

    # Add the p-value to the plot
    x_c5 = df['cluster'].unique().tolist().index('c5')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(6, y_max+0.1, "p < 0.0001", ha='center', fontsize = 7)
    else:
        ax.text(6, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 7)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=7)
    ax.set_ylabel(gene, fontsize=7)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_violin_c5c6c7Markers.png')

###########################
# 7C
###########################
## umap of transferred clusters and original cell types
sc.pl.umap(adata_human_bone, color='cluster', save = '_kfoury_projectedClusters')
sc.pl.umap(adata_human_bone, color='cells',  save = '_kfoury_originalCells')

## umap of BGN expression
#sc.pl.umap(adata_human_bone, color='BGN', color_map = 'RdBu_r', vmin='p1', vmax='p99', save = '_kfoury_BGN')

###########################
# 7D
###########################
# Violin plots of marker genes in the kfoury dataset
Bone_markers = ['MKI67', 'POSTN', 'BGN', 'SFRP4', 'TNFSF11',
                'BMP6', 'BMP7', 'RUNX2', 'PTHLH', 'EDN1',
                'CXCR4', 'CXCL12', 'TGFB1', 'SPP1',
                'WNT1', 'WNT3A', 'WNT5A', 'WNT7B',
                'ITGAV', 'ITGB3', 'VEGFA', 'CDH1', 'CDH2']

Bone_markers_small = ['POSTN', 'BGN', 'RUNX2', 'SPP1']

for i in Bone_markers:
    sc.pl.violin(adata_human_bone, i, groupby='cluster', use_raw=True, save=f'_kfoury_{i}')




# Create a figure and axes objects
fig, axs = plt.subplots(1, len(Bone_markers_small), figsize=(20, 6))

# Loop over the genes
for i, (ax, gene) in enumerate(zip(axs, Bone_markers_small)):
    # Get the expression data for the gene and add it to the data frame
    gene_expression = pd.DataFrame(adata_human_bone.raw[:, gene].X,
                                   index=adata_human_bone.obs_names,
                                   columns=[gene])

    df = gene_expression.join(adata_human_bone.obs)

    # Perform a statistical test (e.g., t-test)
    c5c6c7_expression = df[df['cluster'].isin(['c5','c6', 'c7'])][gene]
    other_expression = df[~df['cluster'].isin(['c5','c6', 'c7'])][gene]
    t_stat, p_value = stats.ttest_ind(c5c6c7_expression, other_expression, equal_var=False)

    # Create the violin plot for this gene
    sns.violinplot(x='cluster', y=gene, order=['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7'], data=df, ax=ax, cut=0.1, scale='width')

    # Add the p-value to the plot
    x_c5 = df['cluster'].unique().tolist().index('c5')
    y_max = df[gene].max()
    # Add the p-value to the plot
    if p_value < 0.0001:
        ax.text(3.5, y_max+0.1, "p < 0.0001", ha='center', fontsize = 7)
    else:
        ax.text(3.5, y_max+0.1, f"p = {p_value:.4f}", ha='center', fontsize = 7)

    # Adjust the y-axis limit to accommodate the p-value text
    ax.set_ylim([df[gene].min(), y_max + 0.5])

    # Set the title for this subplot
    #ax.set_title(gene, fontsize=16)

    # Set the label for the x and y axes
    ax.set_xlabel('cluster', fontsize=7)
    ax.set_ylabel(gene, fontsize=7)

    # Increase the size of the tick labels
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)

# Save the figure
plt.tight_layout()
plt.savefig('figures/figures_cell/human_violin_bone_markers.png')
