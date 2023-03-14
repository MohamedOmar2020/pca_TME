
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



#############################################
# set figure parameters
sc.settings.figdir = 'figures/figures_cell'
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
adata_mouse_mesenchyme = sc.read_h5ad('data/for_mouse/adata_mouse.h5ad', chunk_size=100000)
adata_mouse_mesenchyme.obs['cluster'] = adata_mouse_mesenchyme.obs['cluster'].astype('str')


###########################
# 1A
###########################
# umap clusters
sc.pl.umap(adata_mouse_mesenchyme, color = 'cluster', save = '_mouse_clusters', title = '')


#############
## umap by genotype
adata_mouse_mesenchyme.obs.condition.value_counts()

# plot mutants
sc.pl.umap(
    adata_mouse_mesenchyme[adata_mouse_mesenchyme.obs["condition"] == "mutant"],
    color="key_new", size = 10, title = '', palette = 'viridis',
    save = '_mouse_models_mutants'
)

# plot Wildtype

sc.pl.umap(
    adata_mouse_mesenchyme[adata_mouse_mesenchyme.obs["condition"] == "wildtype"],
    color="key_new", size = 10, title = '',
    save = '_mouse_models_wildtype'
)


###########################
# 1B
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
# 1C
###########################
## heatmap of fraction of cells of each stromal cluster by mouse models
def distribution(
    adata,
    partition="cluster",
    labels="key_new",
    modelorder=['B6', 'B6.129', 'FVBN', 'WT for PN', 'WT for PRN', 'T-ERG', 'NP', 'Hi-MYC', 'PRN'],
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
# 1D
###########################
# parallel categories mouse models

# Arrange the cell annotation dataframe by the clusters
adata_mouse_mesenchyme.obs.sort_values(by = ['cluster'], inplace=True)

adata_mouse_mesenchyme.obs['cluster'] = adata_mouse_mesenchyme.obs['cluster'].astype('category')
color_map = dict(
    zip(adata_mouse_mesenchyme.obs['cluster'].cat.categories.tolist(), sc.pl.palettes.vega_20_scanpy)
)
adata_mouse_mesenchyme.obs["color"] = adata_mouse_mesenchyme.obs['cluster'].map(color_map)
fig = px.parallel_categories(
    adata_mouse_mesenchyme.obs,
    dimensions=['cluster', 'key_new'],
    color="color",
    labels={'key_new': "model", 'cluster': "cluster"},
)
fig.update_layout(autosize=False, width=400, height=600, font_size = 9, font_family="Arial Black")
fig.update_yaxes(tickfont_family="Arial Black")
fig.write_image("figures/figures_cell/parallel_categories_cluster.svg", scale = 2)

###############################
# 1E
###############################
# refer to the R script