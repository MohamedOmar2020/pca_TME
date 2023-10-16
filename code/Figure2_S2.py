
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
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.plotting import plot_binarization, plot_rss
from pyscenic.rss import regulon_specificity_scores
from pyscenic.transform import df2regulons
from pyscenic.utils import load_motifs
from matplotlib.pyplot import gcf


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

# load the mouse data
adata_mouse_mesenchyme = sc.read_h5ad('data/for_mouse/adata_mouse.h5ad', chunk_size=100000)
adata_mouse_mesenchyme.obs['cluster'] = adata_mouse_mesenchyme.obs['cluster'].astype('str')

# load the mouse mensenchyme myofibroblast data
adata_mouse_mesenchyme_myo = sc.read_h5ad('data/for_mouse/adata_mouse_myo.h5ad', chunk_size=100000)
adata_mouse_mesenchyme_myo.obs['cluster'] = adata_mouse_mesenchyme_myo.obs['cluster'].astype('str')


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

adata_mouse_mesenchyme_myo.obs['key_new'] = adata_mouse_mesenchyme.obs['key']
adata_mouse_mesenchyme_myo.obs['key_new'].value_counts()
adata_mouse_mesenchyme_myo.obs['key_new'].replace('terg', 'T-ERG', inplace=True)
adata_mouse_mesenchyme_myo.obs['key_new'].replace('himyc', 'Hi-MYC', inplace=True)
adata_mouse_mesenchyme_myo.obs['key_new'].replace('fvbn', 'FVBN', inplace=True)
adata_mouse_mesenchyme_myo.obs['key_new'].replace('pten', 'NP', inplace=True)
adata_mouse_mesenchyme_myo.obs['key_new'].replace('mycn', 'PRN', inplace=True)
adata_mouse_mesenchyme_myo.obs['key_new'].replace('129b6', 'B6.129', inplace=True)
adata_mouse_mesenchyme_myo.obs['key_new'].replace('b6', 'B6', inplace=True)
adata_mouse_mesenchyme_myo.obs['key_new'].replace('129b6_pten', 'WT for NP', inplace=True)
adata_mouse_mesenchyme_myo.obs['key_new'].replace('129b6_mycn', 'WT for PRN', inplace=True)


###########################
# 2A
###########################
# umap showing just c0
sc.pl.umap(
    adata_mouse_mesenchyme,
    color="cluster", size = 10,
    groups = ['c0'],
    title = '',
    #palette = 'turbo',
    na_in_legend = False,
    #add_outline = True, outline_width = [0.05,0.005],
    save = '_mouse_smooth_muscle'
)

# umap myo and pericytes
sc.pl.umap(
    adata_mouse_mesenchyme_myo,
    color="cluster", size = 50,
    #palette = 'Set2',
    title = '',
    save = '_smoothMuscle_subtypes'
)

# myo mouse models
sc.pl.umap(
    adata_mouse_mesenchyme_myo,
    color="key_new", size = 50,
    palette = 'viridis',
    title = '',
    save = '_smoothMuscle_models'
)

###############################
# 2B
###############################
# myo markers
ax = sc.pl.umap(
    adata_mouse_mesenchyme_myo,
    color=[
        "Acta2",
        "Myl9",
        "Myh11",
        "Tagln",
        "Rgs5",
        "Mef2c",
        "Pdgfrb"
    ],
    cmap="RdBu_r",
    vmax=5,
    legend_fontsize=9, show=False,
    #save="_myo_markers"
)
ax[0].title.set_size(12)
ax[1].title.set_size(12)
ax[2].title.set_size(12)
ax[3].title.set_size(12)
ax[4].title.set_size(12)
ax[5].title.set_size(12)
ax[6].title.set_size(12)

ax[0].title.set_fontweight('bold')
ax[1].title.set_fontweight('bold')
ax[2].title.set_fontweight('bold')
ax[3].title.set_fontweight('bold')
ax[4].title.set_fontweight('bold')
ax[5].title.set_fontweight('bold')
ax[6].title.set_fontweight('bold')

plt.savefig('figures/figures_cell/umap_myo_markers.tiff', dpi=300, bbox_inches='tight')

###############################
# 2C
###############################
# dotplot
sc.pl.dotplot(
    adata_mouse_mesenchyme_myo,
    var_names=[
        "Rgs5",
        "Mef2c",
        "Vtn",
        "Cygb",
        "Pdgfrb",
        "Rras",
        "Rasl12",
        "Rspo3",
        "Nrg2",
        "Hopx",
        "Actg2",
    ],
    groupby="cluster",
    cmap="Reds",
    #vmax=3,
    save="myo_markers",
)

###############################
# Figure S2
###############################

###########################
# S2A
###########################
# umap by n_genes
sc.pl.umap(
    adata_mouse_mesenchyme,
    color="n_genes", size = 10,
    title = '', cmap = 'viridis',
    save = '_mouse_genes'
)

# umap by n_counts
sc.pl.umap(
    adata_mouse_mesenchyme,
    color="n_counts", size = 10,
    title = '', cmap = 'viridis',
    save = '_mouse_counts'
)


###########################
# S2B
###########################
## heatmap of fraction of cells of each stromal cluster by mouse models
adata_mouse_mesenchyme.obs['cluster'] = adata_mouse_mesenchyme.obs['cluster'].astype('category')

def distribution(
    adata,
    partition="cluster",
    labels="key_new",
    modelorder=['B6', 'B6.129', 'FVBN', 'WT for NP', 'WT for PRN', 'T-ERG', 'NP', 'Hi-MYC', 'PRN'],
    partitionorder=None,
    figsize=(10, 7),
):
    """
    plots for distribution of cells by key labels
    """

    models = []
    for name in adata.obs[labels].cat.categories.tolist():
        model = adata[adata.obs[labels] == name]
        models.append((name, model))

    # generate percentage dict by model
    bymodeldict = {}
    for name, model in models:
        total = len(model.obs)
        modeldict = {}
        for j in adata.obs[partition].values.categories:
            modeldict[j] = np.sum(model.obs[partition] == j) / total
        bymodeldict[name] = modeldict

    # subset wt and mut
    mutant = adata[adata.obs["condition"] == "mutant"]
    wildtype = adata[adata.obs["condition"] == "wildtype"]
    mutwt = [("mutant", mutant), ("wildtype", wildtype)]

    # compute composition of partition by keys
    byclusterdict = {}
    for cluster in adata.obs[partition].cat.categories:
        clusterdata = adata[adata.obs[partition] == cluster]
        total = np.sum(adata.obs[partition] == cluster)
        clusterdict = {}
        for key in adata.obs[labels].cat.categories.tolist():
            clusterdict[key] = np.sum(clusterdata.obs[labels] == key) / total
        byclusterdict[cluster] = clusterdict
    byclusterdf = pd.DataFrame.from_dict(byclusterdict, orient="index")
    byclusterdf.to_csv("modelcellcomposition2.csv")

    # generate percentage dict by mutwt
    bymutwtdict = {}
    for name, model in mutwt:
        total = len(model.obs)
        modeldict = {}
        for j in adata.obs[partition].values.categories:
            modeldict[j] = np.sum(model.obs[partition] == j) / total
        bymutwtdict[name] = modeldict

    # generate df wildtype/mutant
    lollipopdict = {}
    for i in adata.obs[partition].values.categories:
        lollipopdict[i] = bymutwtdict["wildtype"][i] / bymutwtdict["mutant"][i]
    lollipopdf = pd.DataFrame.from_dict(lollipopdict, orient="index")
    lollipopdf = lollipopdf.rename(columns={0: "ratio of percentage (wt/mut)"})
    if partitionorder:
        lollipopdf = lollipopdf.reindex(partitionorder)

    # plot heatmapall
    import seaborn as sns

    sns.set(rc={"figure.figsize": figsize})
    bymodeldf = pd.DataFrame.from_dict(bymodeldict, orient="index")
    if modelorder:
        bymodeldf = bymodeldf.reindex(modelorder)
    if partitionorder:
        bymodeldf = bymodeldf[partitionorder]
    bymodeldf.to_csv("modelcellcomposition.csv")
    ax = sns.heatmap(
        bymodeldf, cmap="RdBu_r", center=0.05, annot=True, fmt=".1%", cbar=False
    )
    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0)
    ax.get_figure().savefig(f"percentages.png", dpi=500, bbox_inches="tight")


distribution(adata_mouse_mesenchyme)

###########################
# S2C: scenic heatmap
###########################
# load the scenic output
auc_mtx = pd.read_csv('data/for_mouse/scenic_mesenchyme_auc.csv', index_col=0)
bin_mtx = pd.read_csv('data/for_mouse/scenic_mesenchyme_binary.csv', index_col=0)
thresholds = pd.read_csv('data/for_mouse/scenic_mesenchyme_thresholds.csv', index_col=0).threshold

# Arrange the cell annotation dataframe by the clusters
adata_mouse_mesenchyme.obs.sort_values(by = ['cluster'], inplace=True)

# Re-arrange the cells in auc_mtx to match those in cellAnnotation
auc_mtx['CellID'] = auc_mtx.index
cell_index = pd.DataFrame(adata_mouse_mesenchyme.obs_names)
cell_index = cell_index.rename(columns = {0:'CellID'})
auc_mtx_ord = cell_index.merge(auc_mtx, on = "CellID")
auc_mtx_ord.index = auc_mtx_ord['CellID']
del auc_mtx_ord['CellID']
#%%
# calc regulon specificity score
rss_cellType = regulon_specificity_scores( auc_mtx_ord, adata_mouse_mesenchyme.obs['cluster'] )


# Select the top 5 regulons from each cell type
cats = sorted(list(set(adata_mouse_mesenchyme.obs['cluster'])))
topreg = []
for i,c in enumerate(cats):
    topreg.extend(
        list(rss_cellType.T[c].sort_values(ascending=False)[:5].index)
    )
topreg = list(set(topreg))

#%%
# Re-arrange the cells in bin_mtx to match those in cellAnnotation
bin_mtx['CellID'] = bin_mtx.index
cell_index = pd.DataFrame(adata_mouse_mesenchyme.obs_names)
cell_index = cell_index.rename(columns = {0:'CellID'})
bin_mtx_ord = cell_index.merge(bin_mtx, on = "CellID")
bin_mtx_ord.index = bin_mtx_ord['CellID']
del bin_mtx_ord['CellID']

#%%
# Filter to the top 10 regulons per cluster
bin_mtx_top10 = bin_mtx_ord[topreg]

# Transpose
bin_mtx_top10_transposed = bin_mtx_top10.transpose()

#%%
adata_mouse_mesenchyme.obs.cluster = adata_mouse_mesenchyme.obs.cluster.astype('category')
adata_mouse_mesenchyme.obs.cluster.value_counts()

def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f
#%%
colors = sns.color_palette('bright',n_colors=len(cats) )
colorsd = dict( zip( cats, colors ))
colormap = [ colorsd[x] for x in adata_mouse_mesenchyme.obs['cluster'] ]

##########################
# binary heatmap for all regulons in all clusters
N_COLORS = len(adata_mouse_mesenchyme.obs.cluster.dtype.categories)
COLORS = [color['color'] for color in mpl.rcParams["axes.prop_cycle"]]

cell_type_color_lut = dict(zip(adata_mouse_mesenchyme.obs.cluster.dtype.categories, COLORS))
cell_type_color_lut = dict(zip(adata_mouse_mesenchyme.obs.cluster.dtype.categories, adata_mouse_mesenchyme.uns['cluster_colors']))
#cell_id2cell_type_lut = anndata.var.set_index('cell_id').cell_type.to_dict()
bw_palette = sns.xkcd_palette(["white", "black"])

sns.set()
sns.set(font_scale=7.0)
sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 1})
g = sns.clustermap(bin_mtx_ord.T,
                col_colors=bin_mtx_ord.index.map(adata_mouse_mesenchyme.obs['cluster'].to_dict()).map(cell_type_color_lut),
                yticklabels=True, col_cluster = False, row_cluster = True, dendrogram_ratio=(.05, .05),
                cmap=bw_palette, figsize=(100,100), linewidths=0, xticklabels=False)

# for the legend
for cluster in cats:
    g.ax_col_dendrogram.bar(0, 0, color=colorsd[cluster],
                            label=cluster, linewidth=0.1)
g.ax_col_dendrogram.legend(loc="center", ncol=8, bbox_transform=gcf().transFigure)
col = g.ax_col_dendrogram.get_position()
g.ax_col_dendrogram.set_position([col.x0, col.y0, col.width, col.height])

g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontweight = 'normal')
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_xlabel('Cells')
g.ax_heatmap.set_ylabel('Regulons')
g.ax_col_colors.set_yticks([0.5])
g.ax_col_colors.set_yticklabels(['Cluster'])
g.cax.set_visible(False)
g.ax_col_dendrogram.set_visible(True)
g.ax_row_dendrogram.set_visible(False)
plt.tight_layout()
g.fig.savefig('figures/figures_cell/scenic_mouse_binary.tif', dpi=300)




