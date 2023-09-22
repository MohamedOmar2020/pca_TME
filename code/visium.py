
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import squidpy as sq
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
import SpatialDE
import seaborn as sns

#############################################
# set figure parameters
sc.settings.figdir = 'figures/visium'
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

################################################
# read the data
################################################
# read the spaceranger-processed visium data
adata_PRN1 = sq.read.visium(path = 'visium/outs/PRN1/outs', counts_file = 'filtered_feature_bc_matrix.h5')
adata_PRN1.obs_names
adata_PRN1.var_names
adata_PRN1.var_names_make_unique()

adata_PRN1_wt = sq.read.visium(path = 'visium/outs/PRN1_wt/outs', counts_file = 'filtered_feature_bc_matrix.h5')
adata_PRN1_wt.obs_names
adata_PRN1_wt.var_names
adata_PRN1_wt.var_names_make_unique()

###########################
# read and merge the initial clusters
###########################
PRN1_clusters = pd.read_csv('visium/outs/PRN1/outs/analysis/clustering/gene_expression_graphclust/clusters.csv')
PRN1_wt_clusters = pd.read_csv('visium/outs/PRN1_wt/outs/analysis/clustering/gene_expression_graphclust/clusters.csv')

# check if rownames/cells/barcodes are the same
PRN1_clusters.index = PRN1_clusters['Barcode']
print(adata_PRN1.obs.index.equals(PRN1_clusters.index))

PRN1_wt_clusters.index = PRN1_wt_clusters['Barcode']
print(adata_PRN1_wt.obs.index.equals(PRN1_wt_clusters.index))


# add the clusters to the adata object
adata_PRN1.obs = pd.merge(adata_PRN1.obs, PRN1_clusters, left_index=True, right_index=True)
adata_PRN1.obs['Cluster'] = adata_PRN1.obs['Cluster'].astype('str')
adata_PRN1.obs['Cluster'].value_counts()

adata_PRN1_wt.obs = pd.merge(adata_PRN1_wt.obs, PRN1_wt_clusters, left_index=True, right_index=True)
adata_PRN1_wt.obs['Cluster'] = adata_PRN1_wt.obs['Cluster'].astype('str')
adata_PRN1_wt.obs['Cluster'].value_counts()

####
# add c
adata_PRN1.obs['Cluster'] = 'PRN_c' + adata_PRN1.obs['Cluster']
adata_PRN1.obs['Cluster'].value_counts()

adata_PRN1_wt.obs['Cluster'] = 'wt_c' + adata_PRN1_wt.obs['Cluster']
adata_PRN1_wt.obs['Cluster'].value_counts()

############################################
## combine
############################################

adata_list = [adata_PRN1, adata_PRN1_wt]
adata_all = ad.concat(adata_list, join="outer", label = 'model', uns_merge='unique')
adata_all.obs['model'].value_counts()
adata_all.obs['model'].replace('0', 'PRN', inplace=True)
adata_all.obs['model'].replace('1', 'WT', inplace=True)

################################################
# QC and preprocessing
################################################
# combined adata
sc.pp.normalize_total(adata_all, inplace=True)
adata_all.layers["counts"] = adata_all.X.copy()
sc.pp.log1p(adata_all)
sc.pp.highly_variable_genes(adata_all, flavor="seurat", n_top_genes=4000)
adata_all.raw = adata_all

# PRN1
sc.pp.normalize_total(adata_PRN1, inplace=True)
adata_PRN1.layers["counts"] = adata_PRN1.X.copy()
sc.pp.log1p(adata_PRN1)
sc.pp.highly_variable_genes(adata_PRN1, flavor="seurat", n_top_genes=4000)
adata_PRN1.raw = adata_PRN1


# PRN1 wt
sc.pp.normalize_total(adata_PRN1_wt, inplace=True)
adata_PRN1_wt.layers["counts"] = adata_PRN1_wt.X.copy()
sc.pp.log1p(adata_PRN1_wt)
sc.pp.highly_variable_genes(adata_PRN1_wt, flavor="seurat", n_top_genes=4000)
adata_PRN1_wt.raw = adata_PRN1_wt

################################################
# UMAP and Clustering
################################################
# On PRN1
sc.pp.pca(adata_PRN1, n_comps=10)
sc.pp.neighbors(adata_PRN1)
sc.tl.umap(adata_PRN1)
sc.tl.leiden(adata_PRN1, key_added="leiden", resolution=0.5)

# On PRN1 wt
sc.pp.pca(adata_PRN1_wt, n_comps=10)
sc.pp.neighbors(adata_PRN1_wt)
sc.tl.umap(adata_PRN1_wt)
sc.tl.leiden(adata_PRN1_wt, key_added="leiden", resolution=0.5)

# on combined adata
sc.pp.pca(adata_all, n_comps=10)
sc.pp.neighbors(adata_all)
sc.tl.umap(adata_all)
sc.tl.leiden(adata_all, key_added="leiden", resolution=0.5)
adata_all.obs['leiden'].value_counts()

################################################
# save
################################################
adata_all.write('objs/visium_adata_all.h5ad')

# load
adata_all = sc.read_h5ad('objs/visium_adata_all.h5ad', chunk_size=100000)
adata_all.obs_names_make_unique()
adata_all.uns['log1p']["base"] = None

################################################
# plots
################################################

## spatial
# PRN
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['PRN'])], img_key="hires", color=["leiden"], library_id = 'PRN1', size=1, save = '_PRN1_leiden')
# WT
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['WT'])], img_key="hires", color=["leiden"], library_id = 'PRN1_wt', size=1, save = '_WT_leiden')

##############
## umap
# PRN
sc.pl.umap(adata_all[adata_all.obs['model'].isin(['PRN'])], color=["leiden"], size=10, legend_fontsize=6, save = '_PRN1_leiden')
# WT
sc.pl.umap(adata_all[adata_all.obs['model'].isin(['WT'])], color=["leiden"], size=10, legend_fontsize=6, save = '_WT_leiden')


#################
# expression of interesting genes

# PRN
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['PRN'])], img_key="hires", color=["Postn", 'Ar', 'Bgn'], library_id = 'PRN1', alpha=0.7, cmap = 'viridis', save = 'PRN1_Postn_Ar_Bgn')
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['PRN'])], img_key="hires", color=["Postn"], library_id = 'PRN1', alpha=0.8, cmap = 'viridis', save = 'PRN1_Postn')
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['PRN'])], img_key="hires", color=["Ar"], library_id = 'PRN1', alpha=0.8, cmap = 'viridis', save = 'PRN1_Ar')
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['PRN'])], img_key="hires", color=["Bgn"], library_id = 'PRN1', alpha=0.8, cmap = 'viridis', save = 'PRN1_Bgn')

sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['WT'])], img_key="hires", color=["Postn", 'Ar', 'Bgn'], library_id = 'PRN1_wt', alpha=0.7, cmap = 'viridis', save = 'PRN1_wt_Postn_Ar_Bgn')
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['WT'])], img_key="hires", color=["Postn"], library_id = 'PRN1_wt', alpha=0.8, cmap = 'viridis', save = 'PRN1_wt_Postn')
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['WT'])], img_key="hires", color=["Ar"], library_id = 'PRN1_wt', alpha=0.8, cmap = 'viridis', save = 'PRN1_wt_Ar')
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['WT'])], img_key="hires", color=["Bgn"], library_id = 'PRN1_wt', alpha=0.8, cmap = 'viridis', save = 'PRN1_wt_Bgn')

################
# Violin plots for c5:c7 markers
sc.pl.violin(adata_all, ['Postn', 'Ar'], groupby = 'model', use_raw=True)
sc.pl.violin(adata_all, ['Bgn'], groupby = 'model', use_raw=True, save='_Bgn')

sc.pl.violin(adata_all, ['Sfrp4', 'Mki67', 'Fn1', 'Tnc', 'Col12a1', 'Fzd1', 'Tgfb1'], groupby = 'model', use_raw=True, save='_PRNclusterMarkers_comparison')

# same but spatial
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['WT'])], img_key="hires", color=['Sfrp4', 'Mki67', 'Fn1', 'Tnc', 'Col12a1', 'Fzd1', 'Tgfb1'], library_id = 'PRN1_wt', alpha=0.8, cmap = 'viridis', save = 'PRN1_wt_PRNclustersMarkers')
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['PRN'])], img_key="hires", color=['Sfrp4', 'Mki67', 'Fn1', 'Tnc', 'Col12a1', 'Fzd1', 'Tgfb1'], library_id = 'PRN1', alpha=0.8, cmap = 'viridis', save = 'PRN1_PRNclustersMarkers')


##################################################
## marker genes
##################################################
# clusters by geneotype
pd.crosstab(adata_all.obs['leiden'], adata_all.obs['model'])

# by clusters
sc.tl.rank_genes_groups(adata_all, "leiden", method="t-test", pts=True, use_raw = True)
sc.pl.rank_genes_groups(adata_all, groups=None, n_genes=25, groupby="leiden", sharey=False, save='_topGenes_leiden')
sc.pl.rank_genes_groups_heatmap(adata_all, n_genes=3, use_raw=False, swap_axes=True, vmin=-5, vmax=5, cmap='bwr', figsize=(10,7), save='_topGenes_leiden')

# by model
sc.tl.rank_genes_groups(adata_all, "model", method="wilcoxon", pts=True, use_raw = True)
sc.pl.rank_genes_groups_dotplot(adata_all, groupby='model', n_genes=15, save='topGenesByModel')

##################################################
## Cell types annotation
##################################################
# For myofibroblasts
# ACTA2
# DES
# MYH11
# TAGLN

# For CAFs
# FAP
# PDGFRB
# S100A4
# ACTA2
# COL1A1


adata_all.obs["myofibroblasts"] = (
    adata_all.obs["leiden"].isin(['0', '6']).astype("category")
)

adata_all.obs["fibroblasts"] = (
    adata_all.obs["leiden"].isin(['11']).astype("category")
)

adata_all.obs["tumor epithelium"] = (
    adata_all.obs["leiden"].isin(['1', '2', '3']).astype("category")
)

adata_all.obs["normal epithelium"] = (
    adata_all.obs["leiden"].isin(['4', '9', '13', '15']).astype("category")
)

adata_all.obs["CAFs"] = (
    adata_all.obs["leiden"].isin(['5', '7', '14']).astype("category")
)

# adata_all.obs["myofibroblasts"] = (
#     adata_all.obs["leiden"].isin(['6']).astype("category")
# )

adata_all.obs["myeloid cells"] = (
    adata_all.obs["leiden"].isin(['8']).astype("category")
)

adata_all.obs["SV epithelial cells"] = (
    adata_all.obs["leiden"].isin(['10', '12']).astype("category")
)



adata_all.obs["cell types"] = "Other"
#adata_all.obs.loc[adata_all.obs["normal smooth muscles"] == True, "cell types"] = "normal smooth muscles"
adata_all.obs.loc[adata_all.obs['fibroblasts'] == True, "cell types"] = "fibroblasts"
adata_all.obs.loc[adata_all.obs['tumor epithelium'] == True, "cell types"] = "tumor epithelium"
adata_all.obs.loc[adata_all.obs['normal epithelium'] == True, "cell types"] = "normal epithelium"
adata_all.obs.loc[adata_all.obs['CAFs'] == True, "cell types"] = "CAFs"
adata_all.obs.loc[adata_all.obs['myofibroblasts'] == True, "cell types"] = "myofibroblasts"
adata_all.obs.loc[adata_all.obs['myeloid cells'] == True, "cell types"] = "myeloid cells"
adata_all.obs.loc[adata_all.obs['SV epithelial cells'] == True, "cell types"] = "SV epithelial cells"


adata_all.obs['cell types'] = adata_all.obs['cell types'].astype('category')
adata_all.obs['cell types'].value_counts()

# add a key for the compartment
adata_all.obs["stroma"] = (
    adata_all.obs["cell types"].isin(['CAFs', 'myofibroblasts', 'fibroblasts']).astype("category")
)

adata_all.obs["epithelium"] = (
    adata_all.obs["cell types"].isin(['normal epithelium', 'tumor epithelium']).astype("category")
)

adata_all.obs["compartment"] = "Other"
adata_all.obs.loc[adata_all.obs["stroma"] == True, "compartment"] = "stroma"
adata_all.obs.loc[adata_all.obs['epithelium'] == True, "compartment"] = "epithelium"

adata_all.obs['compartment'] = adata_all.obs['compartment'].astype('category')
adata_all.obs['compartment'].value_counts()

# plot top genes by cell type
sc.tl.rank_genes_groups(adata_all, "cell types", method="t-test", pts=True, use_raw = True)
sc.pl.rank_genes_groups(adata_all, groups=None, n_genes=25, groupby="cell types", sharey=False, save='_topGenes_celltypes')

##################################################
# some plots for the annotated cell types
##################################################
## spatial plots for annotated cell types

def generate_distinct_colors(n):
    return plt.cm.get_cmap('tab20b', n)

# Get the unique cell types across all data
all_cell_types = np.unique(adata_all.obs['cell types'])

# Generate a color palette with the same number of colors as cell types
cmap = generate_distinct_colors(len(all_cell_types))
color_palette = [cmap(i) for i in range(len(all_cell_types))]

# Create a dictionary mapping cell types to colors
cell_type_colors = dict(zip(all_cell_types, color_palette))

# Now, when we plot, the colors should be consistent
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['PRN'])], img_key="hires", color=["cell types"], library_id = 'PRN1', palette = cell_type_colors, size=1, save = '_PRN1_celltypes')
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['WT'])], img_key="hires", color=["cell types"], library_id = 'PRN1_wt', palette = cell_type_colors, size=1, save = '_WT_celltypes')

## spatial plots for compartment
# PRN
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['PRN'])], img_key="hires", color=["compartment"], library_id = 'PRN1', size=1, save = '_PRN1_compartment')
# WT
sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['WT'])], img_key="hires", color=["compartment"], library_id = 'PRN1_wt', size=1, save = '_WT_compartment')


##################################################
# plots for c5:c7 markers in the stroma
##################################################
sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])], ['Postn', 'Ar'], groupby = 'model', use_raw=True, save='_Postn_Ar_stroma')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])], ['Bgn'], groupby = 'model', use_raw=True, save='_Bgn_stroma')

# PRN vs wildtype
sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Sfrp4', 'Mki67', 'Fn1', 'Tnc', 'Col12a1', 'Fzd1', 'Tgfb1', 'Top2a', 'Col12a1', 'Col14a1', 'Col16a1'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma')

# individual plots
sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Mki67'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Mki67')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Postn'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Postn')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Fn1'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Fn1')


sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Bgn'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Bgn')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Tnc'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Tnc')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Col12a1'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Col12a1')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Fzd1'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Fzd1')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Tgfb1'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Tgfb1')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Top2a'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Top2a')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Col16a1'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Col16a1')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Vim'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Vim')

sc.pl.violin(adata_all[adata_all.obs['compartment'].isin(['stroma'])],
             ['Ar'],
             groupby = 'model',
             use_raw=True,
             save='_PRNclusterMarkers_comparison_stroma_Ar')


sc.pl.violin(adata_all[adata_all.obs['model'].isin(['PRN']) & adata_all.obs['compartment'].isin(['stroma'])],
             ['Ar', 'Postn'],
             groupby = 'model',
             use_raw=True,
             save='_Ar_Postn_PRN_stroma')

sc.pl.violin(adata_all[adata_all.obs['model'].isin(['PRN']) & adata_all.obs['compartment'].isin(['stroma'])],
             ['Ar'],
             groupby = 'model',
             use_raw=True,
             save='_Ar_PRN_stroma')

sc.pl.violin(adata_all[adata_all.obs['model'].isin(['PRN']) & adata_all.obs['compartment'].isin(['stroma'])],
             ['Postn'],
             groupby = 'model',
             use_raw=True,
             save='_Postn_PRN_stroma')

########################################
from scipy import stats

## compare Ar and Postn in PRN
PRN_stroma = adata_all[adata_all.obs['model'].isin(['PRN']) & adata_all.obs['compartment'].isin(['stroma'])]

##########
# p-value

ar_expression_PRN = PRN_stroma.raw[:, 'Ar'].X
postn_expression_PRN = PRN_stroma.raw[:, 'Postn'].X
ar_expression_PRN_dense = np.array(ar_expression_PRN.todense())
postn_expression_PRN_dense = np.array(postn_expression_PRN.todense())

# Then, perform a statistical test. Here's an example using a t-test from scipy:
t_stat, p_value = stats.ttest_ind(ar_expression_PRN_dense, postn_expression_PRN_dense)

#########
# Create a DataFrame from the expression data
df = pd.DataFrame({
    'Ar': np.array(PRN_stroma.raw[:, 'Ar'].X.todense()).ravel(),
    'Postn': np.array(PRN_stroma.raw[:, 'Postn'].X.todense()).ravel(),
})

# "Melt" the dataset to have genes and their values in two separate columns
df_melted = df.melt(var_name='Gene', value_name='Expression')

# Create the violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x='Gene', y='Expression', data=df_melted, scale='width', inner=None)
# Add significance line
#y_max = df_melted['Expression'].max()  # find maximum y value
#plt.plot([0, 1], [y_max + 0.45, y_max + 0.45], lw=1.5, color='black')  # draw horizontal line
#p_value_scalar = p_value[0] if isinstance(p_value, np.ndarray) else p_value
#plt.text(0.5, y_max + 0.5, f"p-value = {p_value_scalar:.2e}", ha='center')  # add p-value text

# Add title and labels with larger font size
plt.xlabel('Gene', size=16)
plt.ylabel('Expression', size=16)

# Increase the size of the tick labels
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

# save the plot
plt.savefig('figures/Visium/violin_Ar_Postn_PRN_stroma.png')

########################################
## compare Ar and Postn in WT
WT_stroma = adata_all[adata_all.obs['model'].isin(['WT']) & adata_all.obs['compartment'].isin(['stroma'])]

##########
# p-value

ar_expression_WT = WT_stroma.raw[:, 'Ar'].X
postn_expression_WT = WT_stroma.raw[:, 'Postn'].X
ar_expression_WT_dense = np.array(ar_expression_WT.todense())
postn_expression_WT_dense = np.array(postn_expression_WT.todense())

# t-test
t_stat_WT, p_value_WT = stats.ttest_ind(ar_expression_WT_dense, postn_expression_WT_dense)

######
# Create a DataFrame from the expression data
df_WT = pd.DataFrame({
    'Ar': np.array(WT_stroma.raw[:, 'Ar'].X.todense()).ravel(),
    'Postn': np.array(WT_stroma.raw[:, 'Postn'].X.todense()).ravel(),
})

# "Melt" the dataset to have genes and their values in two separate columns
df_WT_melted = df_WT.melt(var_name='Gene', value_name='Expression')

# Create the violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x='Gene', y='Expression', data=df_WT_melted, scale='width', inner=None)
# Add significance line
#y_max = df_WT_melted['Expression'].max()  # find maximum y value
#plt.plot([0, 1], [y_max + 0.45, y_max + 0.45], lw=1.5, color='black')  # draw horizontal line
#p_value_WT_scalar = p_value_WT[0] if isinstance(p_value_WT, np.ndarray) else p_value_WT
#plt.text(0.5, y_max + 0.5, f"p-value = {p_value_WT_scalar:.2e}", ha='center')  # add p-value text

# Add title and labels with larger font size
plt.xlabel('Gene', size=16)
plt.ylabel('Expression', size=16)

# Increase the size of the tick labels
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

# save the plot
plt.savefig('figures/Visium/violin_Ar_Postn_WT_stroma.png')



###################
# Compare AR in PRN and WT stroma
from scipy import stats
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Subset data for PRN stroma
PRN_stroma = adata_all[adata_all.obs['model'].isin(['PRN']) & adata_all.obs['compartment'].isin(['stroma'])]
ar_expression_PRN_dense = np.array(PRN_stroma.raw[:, 'Ar'].X.todense()).ravel()

# Subset data for WT stroma
WT_stroma = adata_all[adata_all.obs['model'].isin(['WT']) & adata_all.obs['compartment'].isin(['stroma'])]
ar_expression_WT_dense = np.array(WT_stroma.raw[:, 'Ar'].X.todense()).ravel()

# Perform t-test
t_stat, p_value = stats.ttest_ind(ar_expression_PRN_dense, ar_expression_WT_dense)

# Create a DataFrame for visualization
labels_PRN = ['PRN'] * len(ar_expression_PRN_dense)
labels_WT = ['WT'] * len(ar_expression_WT_dense)

df = pd.DataFrame({
    'Model': labels_PRN + labels_WT,
    'Ar Expression': np.concatenate([ar_expression_PRN_dense, ar_expression_WT_dense])
})

# Create the violin plot
plt.figure(figsize=(10, 6))
sns.violinplot(x='Model', y='Ar Expression', data=df, scale='width')

# Add significance line and p-value text (optional)
y_max = df['Ar Expression'].max()
plt.plot([0, 1], [y_max + 0.45, y_max + 0.45], lw=1.5, color='black')
plt.text(0.5, y_max + 0.5, f"p-value = {p_value:.2e}", ha='center')

# Add title and labels
plt.xlabel('Model', size=16)
plt.ylabel('Ar Expression', size=16)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

# Save the plot
plt.savefig('figures/Visium/violin_Ar_PRN_WT_stroma.png')


### Look deeper
print("PRN Ar Expression:")
print(f"Mean: {ar_expression_PRN_dense.mean()}")
print(f"Standard Deviation: {ar_expression_PRN_dense.std()}")
print(f"Sample Size: {len(ar_expression_PRN_dense)}\n")

print("WT Ar Expression:")
print(f"Mean: {ar_expression_WT_dense.mean()}")
print(f"Standard Deviation: {ar_expression_WT_dense.std()}")
print(f"Sample Size: {len(ar_expression_WT_dense)}")

print("Median of PRN Ar Expression:", np.median(ar_expression_PRN_dense))
print("Median of WT Ar Expression:", np.median(ar_expression_WT_dense))


from scipy.stats import mannwhitneyu

stat, p_mannwhitney = mannwhitneyu(ar_expression_PRN_dense, ar_expression_WT_dense)
print(f"P-value for Mann-Whitney U test: {p_mannwhitney}")

def cohen_d(group1, group2):
    diff = group1.mean() - group2.mean()
    pooled_var = ((len(group1) - 1) * group1.var() + (len(group2) - 1) * group2.var()) / (len(group1) + len(group2) - 2)
    d = diff / np.sqrt(pooled_var)
    return d

d = cohen_d(ar_expression_PRN_dense, ar_expression_WT_dense)
print(f"Cohen's d: {d}")


plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.hist(ar_expression_PRN_dense, bins=30, alpha=0.7, label='PRN', color='blue')
plt.title('PRN Ar Expression')
plt.xlabel('Expression level')
plt.ylabel('Frequency')

plt.subplot(1, 2, 2)
plt.hist(ar_expression_WT_dense, bins=30, alpha=0.7, label='WT', color='green')
plt.title('WT Ar Expression')
plt.xlabel('Expression level')
plt.ylabel('Frequency')

plt.tight_layout()
plt.savefig('figures/Visium/Ar_PRN_WT_stroma_histogram.png')

#######################################################################
# Calculate the Proportion of Zeros:

prn_zeros = (ar_expression_PRN_dense == 0).sum() / len(ar_expression_PRN_dense)
wt_zeros = (ar_expression_WT_dense == 0).sum() / len(ar_expression_WT_dense)

print(f"Proportion of zeros in PRN: {prn_zeros:.2f}")
print(f"Proportion of zeros in WT: {wt_zeros:.2f}")

from scipy.stats import chi2_contingency

# Create a contingency table
contingency_table = [
    [(ar_expression_PRN_dense == 0).sum(), (ar_expression_PRN_dense != 0).sum()],  # PRN: zeros and non-zeros
    [(ar_expression_WT_dense == 0).sum(), (ar_expression_WT_dense != 0).sum()]   # WT: zeros and non-zeros
]

chi2, p_val, _, _ = chi2_contingency(contingency_table)
print(f"Chi-squared p-value: {p_val:.2e}")

import matplotlib.pyplot as plt

labels = ['PRN', 'WT']
zeros = [prn_zeros, wt_zeros]
non_zeros = [1 - prn_zeros, 1 - wt_zeros]

x = range(len(labels))
plt.bar(x, zeros, width=0.4, label='Zeros', color='blue')
plt.bar(x, non_zeros, width=0.4, bottom=zeros, label='Non-Zeros', color='red')

plt.xlabel('Group')
plt.ylabel('Proportion')
plt.title('Proportion of Zeros vs Non-Zeros in PRN and WT')
plt.xticks(x, labels)
plt.legend(loc="upper left")

plt.show()


adata_stroma = adata_all[adata_all.obs['compartment'] == 'stroma'].copy()
# Extract Ar expression values for PRN and WT
ar_expression_PRN = adata_stroma[adata_stroma.obs['model'] == 'PRN', 'Ar'].X.toarray().flatten()
ar_expression_WT = adata_stroma[adata_stroma.obs['model'] == 'WT', 'Ar'].X.toarray().flatten()

# Perform t-test
t_stat, p_value_ttest = ttest_ind(ar_expression_PRN, ar_expression_WT)

# Perform Mann-Whitney U test
stat, p_value_mannwhitney = mannwhitneyu(ar_expression_PRN, ar_expression_WT)

print(f"P-value for t-test: {p_value_ttest:.2e}")
print(f"P-value for Mann-Whitney U test: {p_value_mannwhitney:.2e}")


#######################################################################
### remove the zeros
prn_non_zero = ar_expression_PRN_dense[ar_expression_PRN_dense > 0]
wt_non_zero = ar_expression_WT_dense[ar_expression_WT_dense > 0]

print("Non-zero PRN Ar Expression:")
print(f"Mean: {prn_non_zero.mean()}")
print(f"Median: {np.median(prn_non_zero)}")
print(f"Sample Size: {len(prn_non_zero)}\n")

print("Non-zero WT Ar Expression:")
print(f"Mean: {wt_non_zero.mean()}")
print(f"Median: {np.median(wt_non_zero)}")
print(f"Sample Size: {len(wt_non_zero)}")

plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.hist(prn_non_zero, bins=30, alpha=0.7, color='blue')
plt.title('Non-zero PRN Ar Expression')
plt.xlabel('Expression level')
plt.ylabel('Frequency')

plt.subplot(1, 2, 2)
plt.hist(wt_non_zero, bins=30, alpha=0.7, color='green')
plt.title('Non-zero WT Ar Expression')
plt.xlabel('Expression level')
plt.ylabel('Frequency')

plt.tight_layout()
plt.savefig('figures/Visium/Ar_PRN_WT_stroma_histogram.png')







sc.pl.spatial(adata_all[adata_all.obs['model'].isin(['PRN'])], img_key="hires", color=["Ar", 'Postn'], library_id = 'PRN1', size=1, alpha=0.7, cmap = 'viridis',  save = '_PRN1_Ar_Postn')



# PRN stroma vs PRN epithelium
sc.pl.violin(adata_all[adata_all.obs['model'].isin(['PRN'])],
             ['Sfrp4', 'Mki67', 'Fn1', 'Tnc', 'Col12a1', 'Fzd1', 'Tgfb1', 'Top2a', 'Col12a1', 'Col14a1', 'Col16a1'],
             groupby = 'compartment',
             use_raw=True,
             save='_PRNclusterMarkers_CompartmentComparison_inPRN')



##################################################
## Spatially variable genes
##################################################
#counts = pd.DataFrame(adata_PRN1.X.todense(), columns=adata_PRN1.var_names, index=adata_PRN1.obs_names)
#coord = pd.DataFrame(adata_PRN1.obsm['spatial'], columns=['x_coord', 'y_coord'], index=adata_PRN1.obs_names).to_numpy()
#results = SpatialDE.run(coord, counts)

#results.index = results["g"]
#adata_PRN1.var = pd.concat([adata_PRN1.var, results.loc[adata_PRN1.var.index.values, :]], axis=1)

# inspect significant genes with varying expression in space
#results.sort_values("qval").head(10)
#sc.pl.spatial(adata_PRN1, img_key="hires", color=["Olfm1", "Wdr5"], alpha=0.7)


#####################################################
# Image features
#####################################################
# load the img
# adata_PRN1.obsm['spatial'].shape
# adata_PRN1.obsm['spatial'].min()
# adata_PRN1.obsm['spatial'].max()
# adata_PRN1.uns['spatial']
#
# scale = adata_PRN1.uns['spatial']['PRN1']['scalefactors']['tissue_hires_scalef']
#
# #img = sq.im.ImageContainer(adata_PRN1.uns['spatial']['PRN1']['images']['hires'], scale=scale)
#
# img_path = 'visium/data/EC-HP-7056_CytAssit&Microscope_Images/Cytassist_MT.TIF'
# img = sq.im.ImageContainer(img_path, layer="image", scale=scale)
# img
# img.show(layer="image")
# plt.show()
#
#
# sq.im.calculate_image_features(
#     adata_PRN1,
#     img.compute(),
#     features="summary",
#     key_added=feature_name,
#     n_jobs=6,
#     #scale=scale
# )
#
#
# # combine features in one dataframe
# adata.obsm["features"] = pd.concat(
#     [adata.obsm[f] for f in adata.obsm.keys() if "features_summary" in f], axis="columns"
# )
# # make sure that we have no duplicated feature names in the combined table
# adata.obsm["features"].columns = ad.utils.make_index_unique(adata.obsm["features"].columns)

#####################################################
# Neighborhood enrichment
#####################################################
adata_all[adata_all.obs['model'].isin(['PRN'])] = adata_all[~adata_all.obs['cell types'].isin(['normal epithelium'])]



sq.gr.spatial_neighbors(adata_all[adata_all.obs['model'].isin(['PRN'])])
sq.gr.nhood_enrichment(adata_all[adata_all.obs['model'].isin(['PRN'])], cluster_key="cell types")
sq.pl.nhood_enrichment(adata_all[adata_all.obs['model'].isin(['PRN'])], cluster_key="cell types", figsize = [10,6])
plt.savefig('figures/visium/PRN_neighborhoods.png', dpi=200)

sq.gr.spatial_neighbors(adata_all)
sq.gr.nhood_enrichment(adata_all, cluster_key="cell types")
sq.pl.nhood_enrichment(adata_all[adata_all.obs['model'].isin(['PRN'])], cluster_key="cell types", figsize = [10,6])
plt.savefig('figures/visium/PRN_neighborhoods.png', dpi=200)


#####################################################
# co-occurence probability
#####################################################

sq.gr.co_occurrence(adata_PRN1, cluster_key='cell types')
sq.pl.co_occurrence(adata_PRN1, cluster_key="cell types", clusters="NE epithelium", figsize = [10,6])
plt.savefig('figures/visium/PRN1_NE_cooccurence.png', dpi=200)

#####################################################
# LR
#####################################################
sq.gr.ligrec(
        adata_PRN1,
        n_perms=1000,
        cluster_key="cell types",
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        seed = 123, show_progress_bar = True, corr_method='fdr_bh', alpha = 0.05,
)

sq.pl.ligrec(
        adata_PRN1,
        cluster_key="cell types",
        source_groups=["NE epithelium", 'Fn1+Tnc+Postn+ cells', "Wnt+ stroma", 'immune-reactive stroma'],
        target_groups=["NE epithelium", 'Fn1+Tnc+Postn+ cells', "Wnt+ stroma", 'immune-reactive stroma'],
        means_range=(0.3, np.inf),
        alpha=0.05,
        pvalue_threshold=0.05,
        swap_axes=True,
        kwargs={"width": 20},
        #figsize= (150,3),
        dpi = 200,
        save="lr.png"
)













