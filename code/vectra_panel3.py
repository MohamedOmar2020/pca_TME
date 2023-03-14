
import os
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.pyplot import ion
from scanpy.plotting import _utils
from scanpy.plotting._anndata import *
from statannotations.Annotator import Annotator
import squidpy as sq
from statistics import mean
import warnings
warnings.filterwarnings("ignore")


sc.settings.figdir = 'figures/vectra_panel3'
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


###################################################################################
# define the functions to use downstream
def process_adata(adata, DAPI_thr, autofluor_thr, var_names, choices, POSTN, AR, Synaptophysin):
    # load the adata object
    import numpy as np
    import anndata as ad
    import pandas as pd

    adata = adata
    # Change the var names
    print('initial var names: ' + str(adata.var_names))
    adata.var_names = var_names
    print('given var names: ' + str(adata.var_names))

    # add raw attribute
    adata.raw = adata

    ## scanpy workflow
    import scanpy as sc

    # log
    print('log transformation')
    sc.pp.log1p(adata)

    # Batch correction
    #print('combat batch correction')
    #sc.pp.combat(adata, key='batch')

    # scale
    print('scaling and zero centering')
    sc.pp.scale(adata, max_value=10)

    # Filter by the DAPI and auto-fluor intensity
    #print('filtering by DAPI and autofluorsence')
    #print('number of cells before filtering: ' + str(len(adata.obs_names)))
    #adata = adata[(adata[:, 'DAPI'].X >= DAPI_thr) & (adata[:, 'autoflourscence'].X < autofluor_thr), :]
    #print('number of cells after filtering: ' + str(len(adata.obs_names)))

    # remove DAPI and autofluor
    #keep = var_names[1:(len(var_names)-1)]
    #adata = adata[:, keep]

    # PCA
    print('PCA')
    sc.tl.pca(adata, svd_solver='arpack')

    # Compute neighbors
    print('computing neighbors')
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=15)

    # UMAP
    print('umap')
    #from joblib import parallel_backend
    #with parallel_backend('threading', n_jobs=15):
    sc.tl.umap(adata)

    # Leiden clustering
    # print('leiden clustering')
    # with parallel_backend('threading', n_jobs=15):
    #    sc.tl.leiden(adata, resolution=0.5)

    ## Define conditions for cell types classification
    print('defining cell type conditions based on pre-defined thresholds')

    # Get the count matrix
    #adata_annot = adata.raw.to_adata()
    adata_annot = adata.copy()
    countMat = pd.concat([adata_annot.to_df(), adata_annot.obs], axis = 1)

    ## Annotate the cell types
    # initialize empty slot in adata.obs to store the new cell types
    #adata_annot.obs['cell_types'] = 'unknown'

    conditions = [
        (countMat[POSTN] >= countMat[POSTN].quantile(0.97)) & (countMat[AR] < countMat[AR].quantile(0.95)) & (countMat['halo_label'] == 'STROMA'),
        (countMat[POSTN] >= countMat[POSTN].quantile(0.97)) & (countMat[AR] < countMat[AR].quantile(0.95)) & (countMat['halo_label'] == 'EPITH'),
        (countMat[AR] >= countMat[AR].quantile(0.95)) & (countMat[POSTN] < countMat[POSTN].quantile(0.97)) & (countMat['halo_label'] == 'STROMA'),
        (countMat[AR] >= countMat[AR].quantile(0.95)) & (countMat[POSTN] < countMat[POSTN].quantile(0.97)) & (countMat['halo_label'] == 'EPITH'),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.95)) & (countMat['halo_label'] == 'STROMA'),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.95)) & (countMat['halo_label'] == 'EPITH')
    ]

    choices = choices

    countMat['cell_types'] = np.select(conditions, choices, default='unknown')

    adata_annot.obs['cell_types'] = countMat['cell_types']

#    for label in annot_dict.keys():
#        for key, value in annot_dict[label].items():
#            cond = np.logical_and.reduce([((countMat[k] >= countMat[k].quantile(list(v)[0])) & (countMat[k] <= countMat[k].quantile(list(v)[1]))) for k, v in annot_dict[label].items()])
#            adata_annot.obs.loc[cond, 'cell_types'] = label
            # replace nan with unknown
            # adata_annot.obs.cell_types.fillna('unknown', inplace=True)

            #################################################
        ## extract the coordinates
    #coords = {}
    #for i in adata.obs.cell_types:
    #    coords[i] = adata.obs[adata.obs.cell_types == i][['x', 'y']]

    ###############################
    ## Return the results
    print('done!')
    return adata_annot


##########################
def pipe(directory):
    alladata = None
    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            adata = halotoanndata(os.path.join(directory, filename))
            if alladata is None:
                alladata = adata
            else:
                alladata = alladata.concatenate(adata, join="outer")
    return alladata



def halotoanndata(path):
    """
    Convert Indica Labs HALO output to anndata format.
    """
    raw = pd.read_csv(path, decimal=',')
    cell = pd.DataFrame()
    for i in raw.columns:
        if i.endswith("Cell Intensity"):
            cell[i.replace(" Cell Intensity", "")] = raw[i]

    counts = anndata.AnnData(
        X=cell,
        obs=[
            tuple([np.mean([x2,x1]), np.mean([y2,y1])])
            for x1, x2, y1, y2 in zip(
                raw["XMin"], raw["XMax"], raw["YMin"], raw["YMax"]
            )
        ],
    )
    counts.obs = counts.obs.rename(columns={0: "x", 1: "y"})
    raw = raw.set_index(counts.obs.index)
    counts.obsm["spatial"] = np.array(counts.obs[["x", "y"]])
    counts.obs["cell_area"] = raw["Cell Area (µm²)"]
    counts.obs["cytoplasm_area"] = raw["Cytoplasm Area (µm²)"]
    counts.obs["nucleus_area"] = raw["Nucleus Area (µm²)"]
    counts.obs["nucleus_perimeter"] = raw["Nucleus Perimeter (µm)"]
    counts.obs["nucleus_roundness"] = raw["Nucleus Roundness"]
    counts.obs["halo_label"] = raw["Classifier Label"]
    #counts.obs["condition"] = raw["Algorithm Name"]
    counts.obs["name"] = os.path.basename(path).split("_P1")[0]
    return counts

#####################################################################

folder_path = '/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/from_Nicolo'

cells = pd.read_csv('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/from_Nicolo/MISI3542i_NB100_M2861_P1_Panel1_Scan1.qptiff_1389_job2590.object_results.csv', decimal=',')

#converters={'DAPI Nucleus Intensity': lambda x: float(x.replace(',','.'))}
#cell.apply(lambda x: x.str.replace(',','.'))
#df.stack().str.replace(',','.').unstack()
#pd.options.display.float_format = '{:,}'.format
#cell = cell.style.format('{:,}')

adata_3 = pipe(folder_path)

#adata_3.write_h5ad('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/panel3.h5ad')
# panel 3
#adata_3 = sc.read('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/panel3.h5ad')

#######################################
# modify the var names
adata_3.var_names = ['DAPI', 'Periostin', 'AR', 'autoflourscence', 'Chromogranin', 'PanCK', 'Synaptophysin']

# modify the slide ID
slide_dict = {
    'MISI3542i_NB100_M2861': 'MISI3542i_NB100_M2861',
    'MISI3542i_M3056_3_Panel1_Scan1.qptiff_1386_job2589.object_results.csv': 'MISI3542i_M3056_3'
}
adata_3.obs["slideID"] = adata_3.obs["name"].map(slide_dict)

adata_3.obs["slideID"].value_counts()
adata_3.obs["name"].value_counts()

#######################################
# check the range
#adata_3.X.min()
#adata_3.X.max()

# log-scaling and z score transformation
#sc.pp.log1p(adata_3)
#adata_3.raw = adata_3
#sc.pp.scale(adata_3, max_value=10)

############################################
# plot the spatial plots
adata_3.obs['slideID'].value_counts()
sq.pl.spatial_scatter(adata_3[adata_3.obs['slideID'] == 'MISI3542i_NB100_M2861'], shape=None, color="halo_label", size=1, cmap = 'inferno')
plt.savefig('figures/vectra_panel3/MISI3542i_NB100_M2861.tiff', dpi = 400)

sq.pl.spatial_scatter(adata_3[adata_3.obs['slideID'] == 'MISI3542i_M3056_3'], shape=None, color="halo_label", size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_M3056_3.tiff', dpi = 400)

sc.pl.spatial(adata_3[adata_3.obs['slideID'] == 'MISI3542i_NB100_M2861'], color="halo_label", spot_size=20,)
plt.show()


sq.pl.spatial_scatter(adata_3[adata_3.obs['slideID'] == 'MISI3542i_NB100_M2861'], shape=None, color=['halo_label', "Periostin", 'AR'], size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_NB100_M2861_AR_POSTN.tiff', dpi = 400)
sc.pl.spatial(adata_3[adata_3.obs['slideID'] == 'MISI3542i_NB100_M2861'], color = ['halo_label', "Periostin", 'AR'], spot_size=20, size=5, color_map = 'magma', save='MISI3542i_NB100_M2861_AR_POSTN2')


sq.pl.spatial_scatter(adata_3[adata_3.obs['slideID'] == 'MISI3542i_M3056_3'], shape=None, color=['halo_label', "Periostin", 'AR'], size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_M3056_3_AR_POSTN.tiff', dpi = 400)
sc.pl.spatial(adata_3[adata_3.obs['slideID'] == 'MISI3542i_M3056_3'], color = ['halo_label', "Periostin", 'AR'], spot_size=20, size = 5, color_map = 'magma', save='MISI3542i_M3056_3_AR_POSTN')


##################################################################################################
##################################################################################################
# cell type annotation
adata_MISI3542i_NB100_M2861_annot= process_adata(adata=adata_3[adata_3.obs['slideID'] == 'MISI3542i_NB100_M2861'],
                                 DAPI_thr=-2, autofluor_thr=3,
                                 var_names= ['DAPI', 'Periostin', 'AR', 'autoflourscence', 'Chromogranin', 'PanCK', 'Synaptophysin'],
                                 POSTN='Periostin',
                                 AR='AR',
                                Synaptophysin = 'Synaptophysin',
                                 choices=['Periostin+ stroma', 'Periostin+ epithelium', 'AR+ stroma', 'AR+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium']
                                 )

adata_MISI3542i_NB100_M2861_annot.obs.cell_types.value_counts()
# remove the unknown
adata_MISI3542i_NB100_M2861_annot = adata_MISI3542i_NB100_M2861_annot[~adata_MISI3542i_NB100_M2861_annot.obs['cell_types'].isin(['unknown'])]
# convert str to category
adata_MISI3542i_NB100_M2861_annot.obs['cell_types'] = adata_MISI3542i_NB100_M2861_annot.obs['cell_types'].astype('category')

# separate the stroma
#stroma_MISI3542i_NB100_M2861 = adata_MISI3542i_NB100_M2861_annot[adata_MISI3542i_NB100_M2861_annot.obs["halo_label"] == "STROMA"]
#stroma_MISI3542i_NB100_M2861.obs['cell_types'] = stroma_MISI3542i_NB100_M2861.obs['cell_types'].astype('category')

######
#adata_MISI3542i_NB100_M2861_annot.obs['x'].max()
#adata_MISI3542i_NB100_M2861_annot.obs['y'].max()
#cells['XMax'].max()
#cells['YMax'].max()

# save the coordinates to disk
for cell in adata_MISI3542i_NB100_M2861_annot.obs.cell_types.unique():
    coords = adata_MISI3542i_NB100_M2861_annot.obs[adata_MISI3542i_NB100_M2861_annot.obs.cell_types == cell][['x', 'y']]
    coords.to_csv('vectra/panel3/coords/MISI3542i_NB100_M2861/new/' + str(cell) + '.csv')


#############################

adata_MISI3542i_M3056_3_annot= process_adata(adata=adata_3[adata_3.obs['slideID'] == 'MISI3542i_M3056_3'],
                                 DAPI_thr=-2, autofluor_thr=3,
                                 var_names= ['DAPI', 'Periostin', 'AR', 'autoflourscence', 'Chromogranin', 'PanCK', 'Synaptophysin'],
                                 POSTN='Periostin',
                                 AR='AR',
                                 Synaptophysin = 'Synaptophysin',
                                 choices=['Periostin+ stroma', 'Periostin+ epithelium', 'AR+ stroma', 'AR+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium']
)

adata_MISI3542i_M3056_3_annot.obs.cell_types.value_counts()

# remove the unknown
adata_MISI3542i_M3056_3_annot = adata_MISI3542i_M3056_3_annot[~adata_MISI3542i_M3056_3_annot.obs['cell_types'].isin(['unknown'])]
# convert str to category
adata_MISI3542i_M3056_3_annot.obs['cell_types'] = adata_MISI3542i_M3056_3_annot.obs['cell_types'].astype('category')

# separate the stroma
#stroma_MISI3542i_M3056_3 = adata_MISI3542i_M3056_3_annot[adata_MISI3542i_M3056_3_annot.obs["halo_label"] == "STROMA"]
#stroma_MISI3542i_M3056_3.obs['cell_types'] = stroma_MISI3542i_M3056_3.obs['cell_types'].astype('category')

# save the coordinates to disk
for cell in adata_MISI3542i_M3056_3_annot.obs.cell_types.unique():
    coords = adata_MISI3542i_M3056_3_annot.obs[adata_MISI3542i_M3056_3_annot.obs.cell_types == cell][['x', 'y']]
    coords.to_csv('vectra/panel3/coords/MISI3542i_M3056_3/new/' + str(cell) + '.csv')

###############################################################
## plot spatial location of cell types:

# for MISI3542i_NB100_M2861:
sq.pl.spatial_scatter(adata_MISI3542i_NB100_M2861_annot, color = ['cell_types'], shape=None, size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_NB100_M2861_cellTypes.tiff', dpi = 400)
sc.pl.spatial(adata_MISI3542i_NB100_M2861_annot, color = ['cell_types'], groups = ['AR+ stroma', 'AR+ epithelium', 'Periostin+ stroma', 'Periostin+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium'], alpha_img = 1.2, spot_size=20, size = 5, color_map = 'magma', save='MISI3542i_NB100_M2861_cellTypes')

sq.pl.spatial_scatter(adata_MISI3542i_NB100_M2861_annot, color = ["Periostin", 'AR', 'cell_types'], groups = ['AR+ stroma', 'AR+ epithelium'], shape=None, size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_NB100_M2861_ARpositive_Cells.tiff', dpi = 400)
sc.pl.spatial(adata_MISI3542i_NB100_M2861_annot, color = ["Periostin", 'AR', 'cell_types'], groups = ['AR+ stroma', 'AR+ epithelium'], na_in_legend = False, alpha_img = 1.2, spot_size=20, size = 5, color_map = 'magma', save='MISI3542i_NB100_M2861_ARpositive_Cells')

sq.pl.spatial_scatter(adata_MISI3542i_NB100_M2861_annot, color = ["Periostin", 'AR', 'cell_types'], groups = ['Periostin+ stroma', 'Periostin+ epithelium'], shape=None, size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_NB100_M2861_POSTNpositive_Cells.tiff', dpi = 400)
sc.pl.spatial(adata_MISI3542i_NB100_M2861_annot, color = ["Periostin", 'AR', 'cell_types'], groups = ['Periostin+ stroma', 'Periostin+ epithelium'], na_in_legend = False, alpha_img = 1.2, spot_size=20, size = 5, color_map = 'magma', save='MISI3542i_NB100_M2861_POSTNpositive_Cells')

sq.pl.spatial_scatter(adata_MISI3542i_NB100_M2861_annot, color = ["Periostin", 'AR', 'cell_types'], groups = ['Periostin+ stroma', 'AR+ stroma'], shape=None, size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_NB100_M2861_AR_POSTN_stroma.tiff', dpi = 400)
sc.pl.spatial(adata_MISI3542i_NB100_M2861_annot, color = ["Periostin", 'AR'], groups = ['Periostin+ stroma', 'AR+ stroma'], ncols =1, na_in_legend = False, alpha_img = 1.2, spot_size=20, size = 5, color_map = 'magma', save='MISI3542i_NB100_M2861_AR_POSTN_stroma')
sc.pl.spatial(adata_MISI3542i_NB100_M2861_annot, color = ['cell_types'], groups = ['Periostin+ stroma', 'AR+ stroma'], palette = ['yellow', 'darkred', 'yellow', 'green', 'yellow', 'yellow'], na_color=None, na_in_legend = False, alpha_img = 1.2, spot_size=20, size = 5, color_map = 'magma', save='MISI3542i_NB100_M2861_AR_POSTN_stroma2')

sq.pl.spatial_scatter(adata_MISI3542i_NB100_M2861_annot, color = ["Synaptophysin", 'Chromogranin', 'cell_types'], groups = ['Synaptophysin+ stroma', 'Synaptophysin+ epithelium'], shape=None, size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_NB100_M2861_NEPC.tiff', dpi = 400)

# for MISI3542i_M3056_3:
sq.pl.spatial_scatter(adata_MISI3542i_M3056_3_annot, color = ["Periostin", 'AR', 'cell_types'], groups = ['AR+ stroma', 'AR+ epithelium'], shape=None, size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_M3056_3_ARpositive_Cells.tiff', dpi = 400)
sc.pl.spatial(adata_MISI3542i_M3056_3_annot, color = ['cell_types'], groups = ['AR+ stroma', 'AR+ epithelium', 'Periostin+ stroma', 'Periostin+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium'], alpha_img = 1.2, spot_size=20, size = 5, color_map = 'magma', save='MISI3542i_M3056_3_cellTypes')

sq.pl.spatial_scatter(adata_MISI3542i_M3056_3_annot, color = ["Periostin", 'AR', 'cell_types'], groups = ['Periostin+ stroma', 'Periostin+ epithelium'], shape=None, size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_M3056_3_POSTNpositive_Cells.tiff', dpi = 400)

sq.pl.spatial_scatter(adata_MISI3542i_M3056_3_annot, color = ["Periostin", 'AR', 'cell_types'], groups = ['Periostin+ stroma', 'AR+ stroma'], shape=None, size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_M3056_3_annot_AR_POSTN_stroma.tiff', dpi = 400)

sq.pl.spatial_scatter(adata_MISI3542i_M3056_3_annot, color = ["Synaptophysin", 'Chromogranin', 'cell_types'], groups = ['Synaptophysin+ stroma', 'Synaptophysin+ epithelium'], shape=None, size=1)
plt.savefig('figures/vectra_panel3/MISI3542i_M3056_3_NEPC.tiff', dpi = 400)

#############################
## neighborhoods enrichment analysis

# all
adata_MISI3542i_NB100_M2861_annot.uns['spatial'] = adata_MISI3542i_NB100_M2861_annot.obsm['spatial']
adata_MISI3542i_M3056_3_annot.uns['spatial'] = adata_MISI3542i_M3056_3_annot.obsm['spatial']

# stroma
#stroma_MISI3542i_NB100_M2861.uns['spatial'] = stroma_MISI3542i_NB100_M2861.obsm['spatial']
#stroma_MISI3542i_M3056_3.uns['spatial'] = stroma_MISI3542i_M3056_3.obsm['spatial']


# MISI3542i_NB100_M2861: All
sq.gr.spatial_neighbors(adata_MISI3542i_NB100_M2861_annot, spatial_key = 'spatial')
sq.gr.nhood_enrichment(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types")
sq.pl.nhood_enrichment(adata_MISI3542i_NB100_M2861_annot, mode = 'zscore', annotate=False, cluster_key="cell_types", figsize = [10,14], vmin=-100, vmax=100, dpi = 300, fontsize='x-large', save='MISI3542i_NB100_M2861_Neighborhoods.tiff')
# interactions
sq.gr.interaction_matrix(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types")
sq.pl.interaction_matrix(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types", figsize = [10,8], dpi = 400, save='MISI3542i_NB100_M2861_Interactions.tiff')

# cooccurence
#sq.gr.co_occurrence(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types")
#sq.pl.co_occurrence(
#    adata_MISI3542i_NB100_M2861_annot,
#    cluster_key="cell_types",
#    clusters='AR+ stroma',
#    palette = 'Accent',
#    #figsize=(10, 5),
#    save='MISI3542i_NB100_M2861_POSTN+stroma_cooccurence.tiff'
#)

# centrality
sq.gr.centrality_scores(
    adata_MISI3542i_NB100_M2861_annot,
    cluster_key="cell_types",
)
sq.pl.centrality_scores(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types", figsize=(20, 5), s=500, save='MISI3542i_NB100_M2861_centrality.tiff')

######
# MISI3542i_NB100_M2861: stroma
#sq.gr.spatial_neighbors(stroma_MISI3542i_NB100_M2861, spatial_key = 'spatial')
#sq.gr.nhood_enrichment(stroma_MISI3542i_NB100_M2861, cluster_key="cell_types")
#sq.pl.nhood_enrichment(stroma_MISI3542i_NB100_M2861, method='centroid', mode = 'zscore', annotate=False, cluster_key="cell_types", figsize = [10,12], vmin=-50, vmax=100, save='MISI3542i_NB100_M2861_Neighborhoods_stroma.tiff')
# interactions
#sq.gr.interaction_matrix(stroma_MISI3542i_NB100_M2861, cluster_key="cell_types")
#sq.pl.interaction_matrix(stroma_MISI3542i_NB100_M2861, cluster_key="cell_types", figsize = [10,12], save='MISI3542i_NB100_M2861_Interactions_stroma.tiff')

##############################
# MISI3542i_M3056_3: all
sq.gr.spatial_neighbors(adata_MISI3542i_M3056_3_annot, spatial_key = 'spatial')
sq.gr.nhood_enrichment(adata_MISI3542i_M3056_3_annot, cluster_key="cell_types")
sq.pl.nhood_enrichment(adata_MISI3542i_M3056_3_annot, mode = 'zscore', annotate=False, cluster_key="cell_types", figsize = [10,12], vmin=-100, vmax=100, dpi = 400,  save='MISI3542i_M3056_3_Neighborhoods.tiff')
# interactions
sq.gr.interaction_matrix(adata_MISI3542i_M3056_3_annot, cluster_key="cell_types")
sq.pl.interaction_matrix(adata_MISI3542i_M3056_3_annot, cluster_key="cell_types", figsize = [10,12], method="ward", vmax=20000, save='MISI3542i_M3056_3_Interactions.tiff')
# centrality
sq.gr.centrality_scores(
    adata_MISI3542i_M3056_3_annot,
    cluster_key="cell_types",
)
sq.pl.centrality_scores(adata_MISI3542i_M3056_3_annot, cluster_key="cell_types", figsize=(20, 5), s=500, save='MISI3542i_M3056_3_centrality.tiff')

###########
# MISI3542i_M3056_3: remove synapto+ cells
adata_MISI3542i_NB100_M2861_annot_noNEPC = adata_MISI3542i_NB100_M2861_annot[~adata_MISI3542i_NB100_M2861_annot.obs['cell_types'].isin(['Synaptophysin+ stroma', 'Synaptophysin+ epithelium'])]
sq.gr.spatial_neighbors(adata_MISI3542i_NB100_M2861_annot_noNEPC, spatial_key = 'spatial', coord_type = 'generic')
sq.gr.nhood_enrichment(adata_MISI3542i_NB100_M2861_annot_noNEPC, cluster_key="cell_types")
sq.pl.nhood_enrichment(adata_MISI3542i_NB100_M2861_annot_noNEPC, mode = 'zscore', annotate=False, cluster_key="cell_types", figsize = [6,6], dpi=500, vmin=-100, vmax=100,  palette = ['blue', 'darkred', 'orange', 'green'], kwargs={"fontsize": "xx-large"}, title='', save='MISI3542i_M3056_3_Neighborhoods_noNEPC.tiff')

# MISI3542i_M3056_3: stroma
#sq.gr.spatial_neighbors(stroma_MISI3542i_M3056_3, spatial_key = 'spatial')
#sq.gr.nhood_enrichment(stroma_MISI3542i_M3056_3, cluster_key="cell_types")
#sq.pl.nhood_enrichment(stroma_MISI3542i_M3056_3, method='centroid', mode = 'zscore', annotate=False, cluster_key="cell_types", figsize = [10,12], vmin=-50, vmax=100, save='MISI3542i_M3056_3_Neighborhoods_stroma.tiff')
# interactions
#sq.gr.interaction_matrix(stroma_MISI3542i_M3056_3, cluster_key="cell_types")
#sq.pl.interaction_matrix(stroma_MISI3542i_M3056_3, cluster_key="cell_types", figsize = [10,12], save='MISI3542i_M3056_3_Interactions_stroma.tiff')




#adata_3.obs["pheno"] = adata_3.obs["pheno"].astype("category")


for marker in [
        "DAPI",
        "Periostin",
        "AR",
        "autoflourscence",
        "Chromogranin",
        "PanCK",
        "Synaptophysin",
    ]:
    sanitize_anndata(stroma_MISI3542i_NB100_M2861)

    obs_df = get.obs_df(
        stroma_MISI3542i_NB100_M2861, keys=[marker] + ["cell_types"], layer=None, use_raw=False
    )
    obs_tidy = obs_df

    # compute p values using averages in each slide
    from scipy.stats import mannwhitneyu, ttest_ind, wilcoxon

    stats_df = get.obs_df(
        stroma_MISI3542i_NB100_M2861, keys=[marker] + ["name"] + ['cell_types'], layer=None, use_raw=False
    )
    stats_df = (
        stats_df.groupby(["name", 'cell_types'], as_index=False).mean().dropna()
    )
    pair = stats_df.cell_types.astype("category").cat.categories.tolist()
    assert len(pair) == 4
    data1 = stats_df.groupby("cell_types")[marker].get_group(pair[0])
    data2 = stats_df.groupby("cell_types")[marker].get_group(pair[2])
    stat, p = ttest_ind(data1, data2)
    pvalues = [p]
    print(pvalues)

    groupby = "cell_types"
    x = "cell_types"
    ylabel = [marker]
    ys = [marker]

    axs, _, _, _ = setup_axes(
        ax=None, panels=[marker], show_ticks=True, right_margin=0.3,
    )
    for ax, y, ylab in zip(axs, ys, ylabel):
        ax = sns.violinplot(
            x=x,
            y=y,
            data=obs_tidy,
            order=None,
            orient="vertical",
            ax=ax,
            cut=0,
        )
        ax = sns.stripplot(
            x=x,
            y=y,
            data=obs_tidy,
            order=None,
            jitter=True,
            color="black",
            size=0.005,
            ax=ax,
        )
        xlabel = groupby.replace("_", " ")
        ax.set_ylabel(ylab)
        ax.set_xlabel("condition")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        show = settings.autoshow
        pairs = [(pair[0], pair[2])]
        annotator = Annotator(ax, pairs, data=obs_tidy, x=x, y=y)
        annotator.configure(
            test=None, test_short_name="", text_format="star", loc="outside"
        )
        annotator.set_pvalues(pvalues=pvalues)
        annotator.annotate()
        # annotator.apply_and_annotate()
        _utils.savefig_or_show(
            "violin", show=show, save=f"_{marker}.png"
        )

###################################
from sklearn.neighbors import NearestNeighbors
import time
import sys
from sklearn.cluster import MiniBatchKMeans
import seaborn as sns

# Function for identifying the windows
def get_windows(job, n_neighbors):
    '''
    For each region and each individual cell in dataset, return the indices of the nearest neighbors.

    'job:  meta data containing the start time,index of region, region name, indices of region in original dataframe
    n_neighbors:  the number of neighbors to find for each cell
    '''
    start_time, idx, tissue_name, indices = job
    job_start = time.time()

    print("Starting:", str(idx + 1) + '/' + str(len(exps)), ': ' + exps[idx])

    # tissue_group: a grouped data frame with X and Y coordinates grouped by unique tissue regions
    tissue = compartment_groups.get_group(tissue_name)

    to_fit = tissue.loc[indices][['x', 'y']].values

    # Unsupervised learner for implementing neighbor searches.
    fit = NearestNeighbors(n_neighbors=n_neighbors).fit(tissue[['x', 'y']].values)

    # Find the nearest neighbors

    m = fit.kneighbors(to_fit)

    m = m[0], m[1]

    ## sort_neighbors
    args = m[0].argsort(axis=1)

    add = np.arange(m[1].shape[0]) * m[1].shape[1]

    sorted_indices = m[1].flatten()[args + add[:, None]]

    neighbors = tissue.index.values[sorted_indices]

    end_time = time.time()

    print("Finishing:", str(idx + 1) + "/" + str(len(exps)), ": " + exps[idx], end_time - job_start,
          end_time - start_time)
    return neighbors.astype(np.int32)

# Put in a dataframe for further analysis
countData_MISI3542i_NB100_M2861 = adata_MISI3542i_NB100_M2861_annot.copy().to_df()
#countData_MISI3542i_NB100_M2861.index = countData_MISI3542i_NB100_M2861.reset_index()
countData_MISI3542i_NB100_M2861.index = 'cell_' + countData_MISI3542i_NB100_M2861.index
obs_MISI3542i_NB100_M2861 = adata_MISI3542i_NB100_M2861_annot.obs.copy()
obs_MISI3542i_NB100_M2861.index = 'cell_' + obs_MISI3542i_NB100_M2861.index
#obs_MISI3542i_NB100_M2861.index = obs_MISI3542i_NB100_M2861.reset_index()
data_MISI3542i_NB100_M2861 = pd.concat([countData_MISI3542i_NB100_M2861, obs_MISI3542i_NB100_M2861], axis = 1)
data_MISI3542i_NB100_M2861.shape
#data_MISI3542i_NB100_M2861.index = data_MISI3542i_NB100_M2861.reset_index()
data_MISI3542i_NB100_M2861.index = range(len(data_MISI3542i_NB100_M2861))
data_MISI3542i_NB100_M2861['CellID'] = data_MISI3542i_NB100_M2861.index
data_MISI3542i_NB100_M2861['CellID'] = data_MISI3542i_NB100_M2861['CellID'].astype('str')
data_MISI3542i_NB100_M2861 = data_MISI3542i_NB100_M2861.loc[~data_MISI3542i_NB100_M2861['cell_types'].isin(['unknown'])]
data_MISI3542i_NB100_M2861 = data_MISI3542i_NB100_M2861.loc[~data_MISI3542i_NB100_M2861['cell_types'].isin(['artifact'])]
data_MISI3542i_NB100_M2861['cell_types'] = data_MISI3542i_NB100_M2861['cell_types'].cat.remove_unused_categories()
data_MISI3542i_NB100_M2861['cell_types'] = data_MISI3542i_NB100_M2861['cell_types'].astype('str')
data_MISI3542i_NB100_M2861.shape

# make dummy variables
data_MISI3542i_NB100_M2861 = pd.concat([data_MISI3542i_NB100_M2861,pd.get_dummies(data_MISI3542i_NB100_M2861['cell_types'])], axis = 1)
#data_MISI3542i_NB100_M2861 = data_MISI3542i_NB100_M2861.reset_index()
# Extract the cell types with dummy variables
sum_cols = data_MISI3542i_NB100_M2861['cell_types'].unique()
values = data_MISI3542i_NB100_M2861[sum_cols].values


# Keep the X and Y coordianates + the tissue regions >> then group by tissue regions (140 unique regions)
compartment_groups = data_MISI3542i_NB100_M2861[['x','y','halo_label']].groupby('halo_label')

# Create a list of unique tissue regions
exps = list(data_MISI3542i_NB100_M2861['halo_label'].unique())

# time.time(): current time is seconds
# indices: a list of indices (rownames) of each dataframe in tissue_group
# exps.index(t) : t represents the index of each one of the indices eg, exps.index("reg001_A") is 0 and exps.index("reg001_B") is 1 and so on
# t is the name of tissue regions eg, reg001_A
tissue_chunks = [(time.time(),exps.index(t),t,a) for t,indices in compartment_groups.groups.items() for a in np.array_split(indices,1)]

# Get the window (the 10 closest cells to each cell in each tissue region)
tissues = [get_windows(job,20) for job in tissue_chunks]

# for each cell and its nearest neighbors, reshape and count the number of each cell type in those neighbors.
ks = [20]
out_dict = {}
for k in ks:
    for neighbors, job in zip(tissues, tissue_chunks):
        chunk = np.arange(len(neighbors))  # indices
        tissue_name = job[2]
        indices = job[3]
        window = values[neighbors[chunk, :k].flatten()].reshape(len(chunk), k, len(sum_cols)).sum(axis=1)
        out_dict[(tissue_name, k)] = (window.astype(np.float16), indices)


## concatenate the summed windows and combine into one dataframe for each window size tested.
keep_cols = ['x','y','halo_label','cell_types']
windows = {}
for k in ks:
    window = pd.concat(
        [pd.DataFrame(out_dict[(exp, k)][0], index=out_dict[(exp, k)][1].astype(int), columns=sum_cols) for exp in
         exps], 0)
    window = window.loc[data_MISI3542i_NB100_M2861.index.values]
    window = pd.concat([data_MISI3542i_NB100_M2861[keep_cols], window], 1)
    windows[k] = window

neighborhood_name = "neighborhood"+str(k)
k_centroids = {}

windows2 = windows[10]


#####################
## Clustering the windows

km = MiniBatchKMeans(n_clusters = 5,random_state=0)

labelskm = km.fit_predict(windows2[sum_cols].values)
k_centroids[5] = km.cluster_centers_
data_MISI3542i_NB100_M2861['neighborhood20'] = labelskm
data_MISI3542i_NB100_M2861[neighborhood_name] = data_MISI3542i_NB100_M2861[neighborhood_name].astype('category')
data_MISI3542i_NB100_M2861[neighborhood_name].value_counts()

# This plot shows the cell types abundance in the different niches
cell_order = [
    'AR+ stroma',
    'AR+ epithelium',
    'Periostin+ stroma',
    'Periostin+ epithelium',
    'Synaptophysin+ stroma',
    'Synaptophysin+ epithelium'
    ]
niche_clusters = (k_centroids[5])
sns.set(font_scale=1)
tissue_avgs = values.mean(axis = 0)
fc = np.log2(((niche_clusters+tissue_avgs)/(niche_clusters+tissue_avgs).sum(axis = 1, keepdims = True))/tissue_avgs)
fc = pd.DataFrame(fc,columns = sum_cols)
s=sns.clustermap(fc.loc[[0,1,2, 3, 4],cell_order], vmin =-3,vmax = 3,cmap = 'bwr',row_cluster = False)
plt.savefig('figures/vectra_panel3/neighborhoodsHeatmap_MISI3542i_NB100_M2861.tiff')

# overlayed neighborhoods
sns.set(font_scale=3)
data_MISI3542i_NB100_M2861['neighborhood20'] = data_MISI3542i_NB100_M2861['neighborhood20'].astype('category')
g = sns.lmplot(data = data_MISI3542i_NB100_M2861,x = 'x',y='y',hue = 'neighborhood20',palette = 'bright',height = 20, col = 'halo_label', facet_kws={'legend_out': True}, col_wrap = 2,fit_reg = False)
# legend title
g._legend.set_title('Neighborhoods')
# replace labels
new_labels = ['AR+ epithelium', 'POSTN+ stroma', 'NE stroma', 'AR+ stroma', 'NE epithelium']
for t, l in zip(g._legend.texts, new_labels):
    t.set_text(l)
#plt.tight_layout()
plt.savefig('figures/vectra_panel3/neighborhoods_MISI3542i_NB100_M2861.tiff')







