
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


sc.settings.figdir = 'figures/Vectra_NEPC'
sc.set_figure_params(dpi_save = 400, transparent = False, fontsize =12, format='tiff')
plt.rcParams["font.family"] = "arial"
plt.rcParams['font.size'] = 14
plt.rcParams['font.style'] = 'italic'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 18


###################################################################################
# Thresholding with negative controls
###################################################################################

## PRN
def process_adata_PRN(adata, DAPI_thr, autofluor_thr, var_names, choices, POSTN, AR, Synaptophysin):
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

###################
## DKO
def process_adata_DKO(adata, DAPI_thr, autofluor_thr, var_names, choices, POSTN, AR, Synaptophysin):
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
        (countMat[POSTN] >= countMat[POSTN].quantile(0.99)) & (countMat[AR] < countMat[AR].quantile(0.45)) &  (countMat['PanCK'] < countMat['PanCK'].quantile(0.40)),
        (countMat[POSTN] >= countMat[POSTN].quantile(0.99)) & (countMat[AR] < countMat[AR].quantile(0.45)) & (countMat['PanCK'] > countMat['PanCK'].quantile(0.70)),
        (countMat[AR] >= countMat[AR].quantile(0.90)) & (countMat[POSTN] < countMat[POSTN].quantile(0.99)) & (countMat['PanCK'] < countMat['PanCK'].quantile(0.40)),
        (countMat[AR] >= countMat[AR].quantile(0.80)) & (countMat[POSTN] < countMat[POSTN].quantile(0.99)) & (countMat['PanCK'] > countMat['PanCK'].quantile(0.58)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.80)) & (countMat['PanCK'] < countMat['PanCK'].quantile(0.37)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.60)) &  (countMat['PanCK'] > countMat['PanCK'].quantile(0.43))
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

###################
## TKO
def process_adata_TKO(adata, DAPI_thr, autofluor_thr, var_names, choices, POSTN, AR, Synaptophysin):
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
        (countMat[POSTN] >= countMat[POSTN].quantile(0.93)) & (countMat[AR] < countMat[AR].quantile(0.60)) &  (countMat[PanCK] < countMat[PanCK].quantile(0.50)),
        (countMat[POSTN] >= countMat[POSTN].quantile(0.94)) & (countMat[AR] < countMat[AR].quantile(0.60)) & (countMat[PanCK] > countMat[PanCK].quantile(0.80)),
        (countMat[AR] >= countMat[AR].quantile(0.72)) & (countMat[POSTN] < countMat[POSTN].quantile(0.94)) & (countMat[PanCK] < countMat[PanCK].quantile(0.36)),
        (countMat[AR] >= countMat[AR].quantile(0.72)) & (countMat[POSTN] < countMat[POSTN].quantile(0.94)) & (countMat[PanCK] > countMat[PanCK].quantile(0.60)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.75)) & (countMat[PanCK] < countMat[PanCK].quantile(0.40)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.75)) &  (countMat[PanCK] > countMat[PanCK].quantile(0.70))
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


###################################################################################
# thresholding without negative controls
###################################################################################

## PRN
def process_adata_PRN2(adata, DAPI_thr, autofluor_thr, var_names, choices, POSTN, AR, Synaptophysin):
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
        (countMat[POSTN] >= countMat[POSTN].quantile(0.97)) & (countMat['halo_label'] == 'STROMA'),
        (countMat[POSTN] >= countMat[POSTN].quantile(0.97)) & (countMat['halo_label'] == 'EPITH'),
        (countMat[AR] >= countMat[AR].quantile(0.9)) &  (countMat['halo_label'] == 'STROMA'),
        (countMat[AR] >= countMat[AR].quantile(0.9)) & (countMat['halo_label'] == 'EPITH'),
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

###################
## DKO
def process_adata_DKO2(adata, DAPI_thr, autofluor_thr, var_names, choices, POSTN, AR, Synaptophysin):
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
        (countMat[POSTN] >= countMat[POSTN].quantile(0.99)) &  (countMat['PanCK'] < countMat['PanCK'].quantile(0.40)),
        (countMat[POSTN] >= countMat[POSTN].quantile(0.99)) & (countMat['PanCK'] > countMat['PanCK'].quantile(0.70)),
        (countMat[AR] >= countMat[AR].quantile(0.70))  & (countMat['PanCK'] < countMat['PanCK'].quantile(0.40)),
        (countMat[AR] >= countMat[AR].quantile(0.80)) & (countMat['PanCK'] > countMat['PanCK'].quantile(0.58)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.80)) & (countMat['PanCK'] < countMat['PanCK'].quantile(0.37)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.60)) &  (countMat['PanCK'] > countMat['PanCK'].quantile(0.43))
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

###################
## TKO
def process_adata_TKO2(adata, DAPI_thr, autofluor_thr, var_names, choices, POSTN, AR, Synaptophysin):
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
        (countMat[POSTN] >= countMat[POSTN].quantile(0.93)) & (countMat['PanCK'] < countMat['PanCK'].quantile(0.50)),
        (countMat[POSTN] >= countMat[POSTN].quantile(0.94)) & (countMat['PanCK'] > countMat['PanCK'].quantile(0.80)),
        (countMat[AR] >= countMat[AR].quantile(0.72)) & (countMat['PanCK'] < countMat['PanCK'].quantile(0.36)),
        (countMat[AR] >= countMat[AR].quantile(0.72)) & (countMat['PanCK'] > countMat['PanCK'].quantile(0.60)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.75)) & (countMat['PanCK'] < countMat['PanCK'].quantile(0.40)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.75)) &  (countMat['PanCK'] > countMat['PanCK'].quantile(0.70))
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


###################
## TKO
def process_adata_TRAMP2(adata, DAPI_thr, autofluor_thr, var_names, choices, POSTN, AR, Synaptophysin):
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
        (countMat[POSTN] >= countMat[POSTN].quantile(0.999)) & (countMat['PanCK'] < countMat['PanCK'].quantile(0.25)),
        (countMat[POSTN] >= countMat[POSTN].quantile(0.999)) & (countMat['PanCK'] > countMat['PanCK'].quantile(0.80)),
        (countMat[AR] >= countMat[AR].quantile(0.70)) & (countMat['PanCK'] < countMat['PanCK'].quantile(0.25)),
        (countMat[AR] >= countMat[AR].quantile(0.50)) & (countMat['PanCK'] > countMat['PanCK'].quantile(0.50)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.95)) & (countMat['PanCK'] < countMat['PanCK'].quantile(0.37)),
        (countMat[Synaptophysin] >= countMat[Synaptophysin].quantile(0.95)) &  (countMat['PanCK'] > countMat['PanCK'].quantile(0.43))
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

##############################################################################
##############################################################################
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


##############################################################################
##############################################################################
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

###############################################################################################
########################################################################################################
folder_path = '/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/from_Nicolo'

#cells = pd.read_csv('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/from_Nicolo/MISI3542i_NB100_M2861_P1_Panel1_Scan1.qptiff_1389_job2590.object_results.csv', decimal=',')

# converters={'DAPI Nucleus Intensity': lambda x: float(x.replace(',','.'))}
# cell.apply(lambda x: x.str.replace(',','.'))
# df.stack().str.replace(',','.').unstack()
# pd.options.display.float_format = '{:,}'.format
# cell = cell.style.format('{:,}')

adata_PRN = pipe(folder_path)

#adata_PRN.write_h5ad('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/panel3.h5ad')
# panel 3
adata_PRN = sc.read('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/PRN.h5ad')
adata_DKO = sc.read('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/adata_p1_DKO.h5ad')
adata_TKO = sc.read('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/adata_p1_TKO.h5ad')
adata_TRAMP = sc.read('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/panel3/adata_p1_TRAMP.h5ad')

#######################################
# modify the var names
adata_PRN.var_names = ['DAPI', 'Periostin', 'AR', 'autoflourscence', 'Chromogranin', 'PanCK', 'Synaptophysin']
adata_DKO.var_names = ['DAPI', 'Synaptophysin', 'Periostin', 'Chromogranin', 'AR', 'PanCK', 'autoflourscence']
adata_TKO.var_names = ['DAPI', 'Synaptophysin', 'Periostin', 'Chromogranin', 'AR', 'PanCK', 'autoflourscence']
adata_TRAMP.var_names = ['DAPI', 'Synaptophysin', 'Periostin', 'Chromogranin', 'AR', 'PanCK', 'autoflourscence']

# modify the slide ID
slide_dict = {
    'MISI3542i_NB100_M2861': 'MISI3542i_NB100_M2861',
    'MISI3542i_M3056_3_Panel1_Scan1.qptiff_1386_job2589.object_results.csv': 'MISI3542i_M3056_3'
}
adata_PRN.obs["slideID"] = adata_PRN.obs["name"].map(slide_dict)

adata_PRN.obs["slideID"].value_counts()
adata_PRN.obs["name"].value_counts()

##################################################################################################
# cell type annotation: PRN
##################################################################################################
adata_MISI3542i_NB100_M2861_annot= process_adata_PRN(adata=adata_PRN[adata_PRN.obs['slideID'] == 'MISI3542i_NB100_M2861'],
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

#save the coordinates to disk
for cell in adata_MISI3542i_NB100_M2861_annot.obs.cell_types.unique():
    coords = adata_MISI3542i_NB100_M2861_annot.obs[adata_MISI3542i_NB100_M2861_annot.obs.cell_types == cell][['x', 'y']]
    coords.to_csv('vectra/panel3/coords/MISI3542i_NB100_M2861/new/' + str(cell) + '.csv')

################################################
# re-annotate without negative controls
adata_MISI3542i_NB100_M2861_annot2= process_adata_PRN2(adata=adata_PRN[adata_PRN.obs['slideID'] == 'MISI3542i_NB100_M2861'],
                                 DAPI_thr=-2, autofluor_thr=3,
                                 var_names= ['DAPI', 'Periostin', 'AR', 'autoflourscence', 'Chromogranin', 'PanCK', 'Synaptophysin'],
                                 POSTN='Periostin',
                                 AR='AR',
                                Synaptophysin = 'Synaptophysin',
                                 choices=['Periostin+ stroma', 'Periostin+ epithelium', 'AR+ stroma', 'AR+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium']
                                 )

adata_MISI3542i_NB100_M2861_annot2.obs.cell_types.value_counts()
# remove the unknown
adata_MISI3542i_NB100_M2861_annot2 = adata_MISI3542i_NB100_M2861_annot2[~adata_MISI3542i_NB100_M2861_annot2.obs['cell_types'].isin(['unknown'])]
# convert str to category
adata_MISI3542i_NB100_M2861_annot2.obs['cell_types'] = adata_MISI3542i_NB100_M2861_annot2.obs['cell_types'].astype('category')

# save the coordinates to disk
for cell in adata_MISI3542i_NB100_M2861_annot2.obs.cell_types.unique():
     coords = adata_MISI3542i_NB100_M2861_annot2.obs[adata_MISI3542i_NB100_M2861_annot2.obs.cell_types == cell][['x', 'y']]
     coords.to_csv('vectra/panel3/coords/MISI3542i_NB100_M2861/new/' + str(cell) + '.csv')

##################################################################################################
# cell type annotation: DKO
##################################################################################################
adata_DKO_annot= process_adata_DKO(adata=adata_DKO,
                                 DAPI_thr=-2, autofluor_thr=3,
                                 var_names= ['DAPI', 'Synaptophysin', 'Periostin', 'Chromogranin', 'AR', 'PanCK', 'autoflourscence'],
                                 POSTN='Periostin',
                                 AR='AR',
                                Synaptophysin = 'Synaptophysin',
                                 choices=['Periostin+ stroma', 'Periostin+ epithelium', 'AR+ stroma', 'AR+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium']
                                 )

adata_DKO_annot.obs.cell_types.value_counts()
# remove the unknown
adata_DKO_annot = adata_DKO_annot[~adata_DKO_annot.obs['cell_types'].isin(['unknown'])]
# convert str to category
adata_DKO_annot.obs['cell_types'] = adata_DKO_annot.obs['cell_types'].astype('category')

################################################
# re-annotate without negative controls
adata_DKO_annot2= process_adata_DKO2(adata=adata_DKO,
                                 DAPI_thr=-2, autofluor_thr=3,
                                 var_names= ['DAPI', 'Synaptophysin', 'Periostin', 'Chromogranin', 'AR', 'PanCK', 'autoflourscence'],
                                 POSTN='Periostin',
                                 AR='AR',
                                Synaptophysin = 'Synaptophysin',
                                 choices=['Periostin+ stroma', 'Periostin+ epithelium', 'AR+ stroma', 'AR+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium']
                                 )

adata_DKO_annot2.obs.cell_types.value_counts()
# remove the unknown
adata_DKO_annot2 = adata_DKO_annot2[~adata_DKO_annot2.obs['cell_types'].isin(['unknown'])]
# convert str to category
adata_DKO_annot2.obs['cell_types'] = adata_DKO_annot2.obs['cell_types'].astype('category')

##################################################################################################
# cell type annotation: TKO
##################################################################################################
adata_TKO_annot= process_adata(adata=adata_TKO,
                                 DAPI_thr=-2, autofluor_thr=3,
                                 var_names= ['DAPI', 'Synaptophysin', 'Periostin', 'Chromogranin', 'AR', 'PanCK', 'autoflourscence'],
                                 POSTN='Periostin',
                                 AR='AR',
                                Synaptophysin = 'Synaptophysin',
                                 choices=['Periostin+ stroma', 'Periostin+ epithelium', 'AR+ stroma', 'AR+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium']
                                 )

adata_TKO_annot.obs.cell_types.value_counts()
# remove the unknown
adata_TKO_annot = adata_TKO_annot[~adata_TKO_annot.obs['cell_types'].isin(['unknown'])]
# convert str to category
adata_TKO_annot.obs['cell_types'] = adata_TKO_annot.obs['cell_types'].astype('category')

################################################
# re-annotate without negative controls
adata_TKO_annot2= process_adata_TKO2(adata=adata_TKO,
                                 DAPI_thr=-2, autofluor_thr=3,
                                 var_names= ['DAPI', 'Synaptophysin', 'Periostin', 'Chromogranin', 'AR', 'PanCK', 'autoflourscence'],
                                 POSTN='Periostin',
                                 AR='AR',
                                Synaptophysin = 'Synaptophysin',
                                 choices=['Periostin+ stroma', 'Periostin+ epithelium', 'AR+ stroma', 'AR+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium']
                                 )

adata_TKO_annot2.obs.cell_types.value_counts()
# remove the unknown
adata_TKO_annot2 = adata_TKO_annot2[~adata_TKO_annot2.obs['cell_types'].isin(['unknown'])]
# convert str to category
adata_TKO_annot2.obs['cell_types'] = adata_TKO_annot2.obs['cell_types'].astype('category')

################################################
# TRAMP: annotate without negative controls
adata_TRAMP_annot2= process_adata_TRAMP2(adata=adata_TRAMP,
                                 DAPI_thr=-2, autofluor_thr=3,
                                 var_names= ['DAPI', 'Synaptophysin', 'Periostin', 'Chromogranin', 'AR', 'PanCK', 'autoflourscence'],
                                 POSTN='Periostin',
                                 AR='AR',
                                Synaptophysin = 'Synaptophysin',
                                 choices=['Periostin+ stroma', 'Periostin+ epithelium', 'AR+ stroma', 'AR+ epithelium', 'Synaptophysin+ stroma', 'Synaptophysin+ epithelium']
                                 )

adata_TRAMP_annot2.obs.cell_types.value_counts()
# remove the unknown
adata_TRAMP_annot2 = adata_TRAMP_annot2[~adata_TRAMP_annot2.obs['cell_types'].isin(['unknown'])]
# convert str to category
adata_TRAMP_annot2.obs['cell_types'] = adata_TRAMP_annot2.obs['cell_types'].astype('category')

###############################################################
## plot spatial location of cell types:
cellTypes_colors = ['blue', 'deeppink', 'green', 'cyan', 'white', 'black']

sc.set_figure_params(dpi_save = 400, transparent = False, figsize = [5,5], fontsize =12, format='tiff')
plt.rcParams["font.family"] = "arial"
plt.rcParams['font.size'] = 14
plt.rcParams['font.style'] = 'italic'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['legend.fontsize'] = 14
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 18

# PRN
sc.pl.spatial(adata_MISI3542i_NB100_M2861_annot2, color = ['cell_types'], groups = ['AR+ stroma', 'Periostin+ stroma'], alpha_img = 1.2, spot_size=20, size = 4, palette =cellTypes_colors, title = '', na_in_legend = False, save='_PRN_spatial_celltypes')

# DKO
sc.pl.spatial(adata_DKO_annot2, color = ['cell_types'], groups = ['AR+ stroma', 'Periostin+ stroma'], alpha_img = 1.2, spot_size=20, size = 7, palette =cellTypes_colors, title = '', na_in_legend = False, save='_DKO_spatial_celltypes')

# TKO
sc.pl.spatial(adata_TKO_annot2, color = ['cell_types'], groups = ['AR+ stroma', 'Periostin+ stroma'], alpha_img = 1.2, spot_size=20, size = 3, palette =cellTypes_colors, title = '', na_in_legend = False, save='_TKO_spatial_celltypes')

# TRAMP
sc.pl.spatial(adata_TRAMP_annot2, color = ['cell_types'], groups = ['AR+ stroma', 'Periostin+ stroma'], alpha_img = 1.2, spot_size=20, size = 8, palette =cellTypes_colors, title = '', na_in_legend = False, save='_TRAMP_spatial_celltypes')

##############################################################
# Violin plots
##############################################################
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

## PRN

# Subset the data
adata_MISI3542i_NB100_M2861_annot_stroma2 = adata_MISI3542i_NB100_M2861_annot2[adata_MISI3542i_NB100_M2861_annot2.obs['cell_types'].isin(['AR+ stroma', 'Periostin+ stroma'])]
adata_MISI3542i_NB100_M2861_annot_stroma2.obs['cell_types'].value_counts()

AR_stroma2 = adata_MISI3542i_NB100_M2861_annot_stroma2[:, 'AR'].X.flatten()
Postn_stroma2 = adata_MISI3542i_NB100_M2861_annot_stroma2[:, 'Periostin'].X.flatten()


fig, axs = plt.subplots(1, 2, figsize=(6, 4), sharey=False)
# Plot the AR
sns.violinplot(x=adata_MISI3542i_NB100_M2861_annot_stroma2.obs['cell_types'], y=AR_stroma2, ax=axs[0], color='deeppink', scale='width')
sns.stripplot(x=adata_MISI3542i_NB100_M2861_annot_stroma2.obs['cell_types'], y=AR_stroma2, ax=axs[0], color='black', jitter=True, size=2)
axs[0].set_title('AR', fontsize=16, fontweight='bold')
axs[0].tick_params(axis='x', labelsize=12, rotation=35)  # Increase x-axis tick labels size
axs[0].tick_params(axis='y', labelsize=12)  # Increase y-axis tick labels size
axs[0].set_xlabel('', fontsize=14, fontweight='bold')  # Set x-axis title
axs[0].set_ylabel('normalized intensity', fontsize=12, fontweight='bold')  # Set y-axis label

# Plot the Periostin
sns.violinplot(x=adata_MISI3542i_NB100_M2861_annot_stroma2.obs['cell_types'], y=Postn_stroma2, ax=axs[1], color='orange', scale='width')
sns.stripplot(x=adata_MISI3542i_NB100_M2861_annot_stroma2.obs['cell_types'], y=Postn_stroma2, ax=axs[1], color='black', jitter=True, size=2)
axs[1].set_title('Periostin', fontsize=16, fontweight='bold')
axs[1].tick_params(axis='x', labelsize=12, rotation=35)  # Increase x-axis tick labels size
axs[1].tick_params(axis='y', labelsize=12)  # Increase y-axis tick labels size
axs[1].set_xlabel('', fontsize=14, fontweight='bold')  # Set x-axis title

plt.tight_layout()
plt.savefig('figures/vectra_panel3/violin_AR_Postn_stroma2.tiff')

##############
## compare the overall Ar and Postn intensity in the stroma of PRN
##########
# p-value

# Get protein expression
ar_expression_PRN = adata_MISI3542i_NB100_M2861_annot_stroma2.raw[:, 'AR'].X
postn_expression_PRN = adata_MISI3542i_NB100_M2861_annot_stroma2.raw[:, 'Periostin'].X

# Flatten the arrays if they are 2-dimensional
ar_expression_PRN = ar_expression_PRN.ravel()
postn_expression_PRN = postn_expression_PRN.ravel()

# Perform t-test
t_stat, p_value = stats.ttest_ind(ar_expression_PRN, postn_expression_PRN)

# Create DataFrame for plot
df = pd.DataFrame({
    'AR': ar_expression_PRN,
    'Periostin': postn_expression_PRN,
})

# Melt DataFrame
df_melted = df.melt(var_name='Protein', value_name='Expression')

# Create the violin plot
plt.figure(figsize=(5, 5))
sns.violinplot(x='Protein', y='Expression', palette =['deeppink', 'orange'], data=df_melted, scale='width')

# Add significance line
y_max = df_melted['Expression'].max()  # find maximum y value
plt.plot([0, 1], [y_max + 0.7, y_max + 0.7], lw=1.5, color='black')  # draw horizontal line

# Add downward-pointing edges
plt.plot([0, 0], [y_max + 0.40, y_max + 0.7], lw=1.5, color='black')  # left edge
plt.plot([1, 1], [y_max + 0.40, y_max + 0.7], lw=1.5, color='black')  # right edge

# Format p-value for display
#if p_value < 0.0001:
    #p_value_text = "p < 0.0001"
#else:
p_value_text = f"p = {p_value:.2e}"

plt.text(0.5, y_max + 0.8, p_value_text, ha='center', fontweight = 'bold', fontsize = 15)  # add p-value text

# Add title and labels with larger font size
plt.xlabel('', size=16)
plt.ylabel('Average Intensity', size=16)

# Increase the size of the tick labels
plt.xticks(fontsize=15, fontweight = 'bold')
plt.yticks(fontsize=15, fontweight = 'bold')

# Save the plot
plt.tight_layout()
plt.savefig('figures/vectra_panel3/violin_Ar_Postn_PRN_stroma.png')


###########################
## DKO

# Subset the data
adata_DKO_annot2_stroma = adata_DKO_annot2[adata_DKO_annot2.obs['cell_types'].isin(['AR+ stroma', 'Periostin+ stroma'])]
adata_DKO_annot2_stroma.obs['cell_types'].value_counts()

AR_stroma_DKO = adata_DKO_annot2_stroma[:, 'AR'].X.flatten()
Postn_stroma_DKO = adata_DKO_annot2_stroma[:, 'Periostin'].X.flatten()


fig, axs = plt.subplots(1, 2, figsize=(6, 4), sharey=False)
# Plot the AR
sns.violinplot(x=adata_DKO_annot2_stroma.obs['cell_types'], y=AR_stroma_DKO, order=['AR+ stroma', 'Periostin+ stroma'], cut=0.1, ax=axs[0], color='deeppink', scale='width')
sns.stripplot(x=adata_DKO_annot2_stroma.obs['cell_types'], y=AR_stroma_DKO, order=['AR+ stroma', 'Periostin+ stroma'], ax=axs[0], color='black', jitter=True, size=2)
axs[0].set_title('AR', fontsize=16, fontweight='bold')
axs[0].tick_params(axis='x', labelsize=12, rotation=35)  # Increase x-axis tick labels size
axs[0].tick_params(axis='y', labelsize=12)  # Increase y-axis tick labels size
axs[0].set_xlabel('', fontsize=14)  # Set x-axis title
axs[0].set_ylabel('normalized intensity', fontsize=14)  # Set y-axis label

# Plot the Periostin
sns.violinplot(x=adata_DKO_annot2_stroma.obs['cell_types'], y=Postn_stroma_DKO, order=['AR+ stroma', 'Periostin+ stroma'], cut = 0.1, ax=axs[1], color='cyan', scale='width')
sns.stripplot(x=adata_DKO_annot2_stroma.obs['cell_types'], y=Postn_stroma_DKO, order=['AR+ stroma', 'Periostin+ stroma'], ax=axs[1], color='black', jitter=True, size=2)
axs[1].set_title('Periostin', fontsize=16, fontweight='bold')
axs[1].tick_params(axis='x', labelsize=12, rotation=35)  # Increase x-axis tick labels size
axs[1].tick_params(axis='y', labelsize=12)  # Increase y-axis tick labels size
axs[1].set_xlabel('', fontsize=14)  # Set x-axis title

plt.tight_layout()
plt.savefig('figures/Vectra_NEPC/violin_AR_Postn_stroma_DKO.tiff')

##############
## compare the overall Ar and Postn intensity in the stroma of DKO
##########
# p-value

# Get protein expression
ar_expression_DKO = adata_DKO_annot2_stroma.raw[:, 'AR'].X
postn_expression_DKO = adata_DKO_annot2_stroma.raw[:, 'Periostin'].X

# Flatten the arrays if they are 2-dimensional
ar_expression_DKO = ar_expression_DKO.ravel()
postn_expression_DKO = postn_expression_DKO.ravel()

# Perform t-test
t_stat, p_value = stats.ttest_ind(postn_expression_DKO, ar_expression_DKO)

# Check if the mean of 'Periostin' is greater than the mean of 'AR'
# if postn_expression_DKO.mean() > ar_expression_DKO.mean():
#     p_value = p_value_two_tailed / 2
# else:
#     p_value = 1 - p_value_two_tailed / 2

# Create DataFrame for plot
df_DKO = pd.DataFrame({
    'AR': ar_expression_DKO,
    'Periostin': postn_expression_DKO,
})

# Melt DataFrame
df_DKO_melted = df_DKO.melt(var_name='Protein', value_name='Expression')

# Get the average expression of each protein
ar_avg = ar_expression_DKO.mean()
postn_avg = postn_expression_DKO.mean()

# Create the violin plot
plt.figure(figsize=(5, 5))
sns.violinplot(x='Protein', y='Expression', data=df_DKO_melted, palette =['deeppink', 'orange'], scale='width')
#sns.boxplot(x='Protein', y='Expression', data=df_DKO_melted)

# Add significance line
y_max = df_DKO_melted['Expression'].max()  # find maximum y value
plt.plot([0, 1], [y_max + 1.2, y_max + 1.2], lw=1.5, color='black')  # draw horizontal line

# Add downward-pointing edges
plt.plot([0, 0], [y_max + 0.90, y_max + 1.2], lw=1.5, color='black')  # left edge
plt.plot([1, 1], [y_max + 0.90, y_max + 1.2], lw=1.5, color='black')  # right edge

# Format p-value for display
#if p_value < 0.0001:
    #p_value_text = "p < 0.0001"
#else:
p_value_text = f"p = {p_value:.2e}"
plt.text(0.5, y_max + 1.3, p_value_text, ha='center', fontweight = 'bold', fontsize = 15)  # add p-value text

# Add title and labels with larger font size
plt.xlabel('', size=16)
plt.ylabel('Average Intensity', size=16)

# Increase the size of the tick labels
plt.xticks(fontsize=15, fontweight = 'bold')
plt.yticks(fontsize=15, fontweight = 'bold')

# Save the plot
plt.tight_layout()
plt.savefig('figures/Vectra_NEPC/violin_Ar_Postn_DKO_stroma.png')

###########################
## TKO

# Subset the data
adata_TKO_annot2_stroma = adata_TKO_annot2[adata_TKO_annot2.obs['cell_types'].isin(['AR+ stroma', 'Periostin+ stroma'])]
adata_TKO_annot2_stroma.obs['cell_types'].value_counts()

AR_stroma_TKO = adata_TKO_annot2_stroma[:, 'AR'].X.flatten()
Postn_stroma_TKO = adata_TKO_annot2_stroma[:, 'Periostin'].X.flatten()


fig, axs = plt.subplots(1, 2, figsize=(6, 4), sharey=False)
# Plot the AR
sns.violinplot(x=adata_TKO_annot2_stroma.obs['cell_types'], y=AR_stroma_TKO, ax=axs[0], color='deeppink', scale='width')
sns.stripplot(x=adata_TKO_annot2_stroma.obs['cell_types'], y=AR_stroma_TKO, ax=axs[0], color='black', jitter=True, size=2)
axs[0].set_title('AR', fontsize=16, fontweight='bold')
axs[0].tick_params(axis='x', labelsize=12, rotation=35)  # Increase x-axis tick labels size
axs[0].tick_params(axis='y', labelsize=12)  # Increase y-axis tick labels size
axs[0].set_xlabel('', fontsize=12)  # Set x-axis title
axs[0].set_ylabel('normalized intensity', fontsize=14)  # Set y-axis label

# Plot the Periostin
sns.violinplot(x=adata_TKO_annot2_stroma.obs['cell_types'], y=Postn_stroma_TKO, ax=axs[1], color='cyan', scale='width')
sns.stripplot(x=adata_TKO_annot2_stroma.obs['cell_types'], y=Postn_stroma_TKO, ax=axs[1], color='black', jitter=True, size=2)
axs[1].set_title('Periostin', fontsize=16, fontweight='bold')
axs[1].tick_params(axis='x', labelsize=12, rotation=35)  # Increase x-axis tick labels size
axs[1].tick_params(axis='y', labelsize=12)  # Increase y-axis tick labels size
axs[1].set_xlabel('', fontsize=14)  # Set x-axis title

plt.tight_layout()
plt.savefig('figures/Vectra_NEPC/violin_AR_Postn_stroma_TKO.tiff')

##############
## compare the overall Ar and Postn intensity in the stroma of TKO
##########
# p-value

# Get protein expression
ar_expression_TKO = adata_TKO_annot2_stroma.raw[:, 'AR'].X
postn_expression_TKO = adata_TKO_annot2_stroma.raw[:, 'Periostin'].X

# Flatten the arrays if they are 2-dimensional
ar_expression_TKO = ar_expression_TKO.ravel()
postn_expression_TKO = postn_expression_TKO.ravel()

# Perform t-test
t_stat, p_value = stats.ttest_ind(ar_expression_TKO, postn_expression_TKO)

# Create DataFrame for plot
df_TKO = pd.DataFrame({
    'AR': ar_expression_TKO,
    'Periostin': postn_expression_TKO,
})

# Melt DataFrame
df_TKO_melted = df_TKO.melt(var_name='Protein', value_name='Expression')

# Create the violin plot
plt.figure(figsize=(5, 5))
sns.violinplot(x='Protein', y='Expression', data=df_TKO_melted, palette =['deeppink', 'orange'], scale='width')

# Add significance line
y_max = df_TKO_melted['Expression'].max()  # find maximum y value
plt.plot([0, 1], [y_max + 0.7, y_max + 0.7], lw=1.5, color='black')  # draw horizontal line

# Add downward-pointing edges
plt.plot([0, 0], [y_max + 0.40, y_max + 0.7], lw=1.5, color='black')  # left edge
plt.plot([1, 1], [y_max + 0.40, y_max + 0.7], lw=1.5, color='black')  # right edge

# Format p-value for display
#if p_value < 0.0001:
    #p_value_text = "p < 0.0001"
#else:
p_value_text = f"p = {p_value:.2e}"

plt.text(0.5, y_max + 0.8, p_value_text, ha='center', fontsize = 16, fontweight = 'bold')  # add p-value text

# Add title and labels with larger font size
plt.xlabel('', size=15)
plt.ylabel('Average Intensity', size=15)

# Increase the size of the tick labels
plt.xticks(fontsize=15, fontweight = 'bold')
plt.yticks(fontsize=15, fontweight = 'bold')

# Save the plot
plt.tight_layout()
plt.savefig('figures/Vectra_NEPC/violin_Ar_Postn_TKO_stroma.png')



###########################
## TRAMP

# Subset the data
adata_TRAMP_annot2_stroma = adata_TRAMP_annot2[adata_TRAMP_annot2.obs['cell_types'].isin(['AR+ stroma', 'Periostin+ stroma'])]
adata_TRAMP_annot2_stroma.obs['cell_types'].value_counts()

AR_stroma_TRAMP = adata_TRAMP_annot2_stroma[:, 'AR'].X.flatten()
Postn_stroma_TRAMP = adata_TRAMP_annot2_stroma[:, 'Periostin'].X.flatten()


fig, axs = plt.subplots(1, 2, figsize=(6, 4), sharey=False)
# Plot the AR
sns.violinplot(x=adata_TRAMP_annot2_stroma.obs['cell_types'], y=AR_stroma_TRAMP, ax=axs[0], cut = 0, color='deeppink', scale='width')
sns.stripplot(x=adata_TRAMP_annot2_stroma.obs['cell_types'], y=AR_stroma_TRAMP, ax=axs[0], color='black', jitter=True, size=2)
axs[0].set_title('AR', fontsize=16, fontweight='bold')
axs[0].tick_params(axis='x', labelsize=12, rotation=35)  # Increase x-axis tick labels size
axs[0].tick_params(axis='y', labelsize=12)  # Increase y-axis tick labels size
axs[0].set_xlabel('', fontsize=12)  # Set x-axis title
axs[0].set_ylabel('normalized intensity', fontsize=14)  # Set y-axis label

# Plot the Periostin
sns.violinplot(x=adata_TRAMP_annot2_stroma.obs['cell_types'], y=Postn_stroma_TRAMP, ax=axs[1], color='orange', scale='width')
sns.stripplot(x=adata_TRAMP_annot2_stroma.obs['cell_types'], y=Postn_stroma_TRAMP, ax=axs[1], color='black', jitter=True, size=2)
axs[1].set_title('Periostin', fontsize=16, fontweight='bold')
axs[1].tick_params(axis='x', labelsize=12, rotation=35)  # Increase x-axis tick labels size
axs[1].tick_params(axis='y', labelsize=12)  # Increase y-axis tick labels size
axs[1].set_xlabel('', fontsize=14)  # Set x-axis title

plt.tight_layout()
plt.savefig('figures/Vectra_NEPC/violin_AR_Postn_stroma_TRAMP.tiff')

##############
## compare the overall Ar and Postn intensity in the stroma of TRAMP
##########
# p-value

# Get protein expression
ar_expression_TRAMP = adata_TRAMP_annot2_stroma.raw[:, 'AR'].X
postn_expression_TRAMP = adata_TRAMP_annot2_stroma.raw[:, 'Periostin'].X

# Flatten the arrays if they are 2-dimensional
ar_expression_TRAMP = ar_expression_TRAMP.ravel()
postn_expression_TRAMP = postn_expression_TRAMP.ravel()

# Perform t-test
t_stat, p_value = stats.ttest_ind(ar_expression_TRAMP, postn_expression_TRAMP)

# Create DataFrame for plot
df_TRAMP = pd.DataFrame({
    'AR': ar_expression_TRAMP,
    'Periostin': postn_expression_TRAMP,
})

# Melt DataFrame
df_TRAMP_melted = df_TRAMP.melt(var_name='Protein', value_name='Expression')

# Create the violin plot
plt.figure(figsize=(5, 5))
sns.violinplot(x='Protein', y='Expression', data=df_TRAMP_melted, palette =['deeppink', 'orange'], scale='width')

# Add significance line
y_max = df_TRAMP_melted['Expression'].max()  # find maximum y value
plt.plot([0, 1], [y_max + 0.7, y_max + 0.7], lw=1.5, color='black')  # draw horizontal line

# Add downward-pointing edges
plt.plot([0, 0], [y_max + 0.40, y_max + 0.7], lw=1.5, color='black')  # left edge
plt.plot([1, 1], [y_max + 0.40, y_max + 0.7], lw=1.5, color='black')  # right edge

# Format p-value for display
#if p_value < 0.0001:
    #p_value_text = "p < 0.0001"
#else:
p_value_text = f"p = {p_value:.2e}"

plt.text(0.5, y_max + 0.8, p_value_text, ha='center', fontsize = 16, fontweight = 'bold')  # add p-value text

# Add title and labels with larger font size
plt.xlabel('', size=15)
plt.ylabel('Average Intensity', size=15)

# Increase the size of the tick labels
plt.xticks(fontsize=15, fontweight = 'bold')
plt.yticks(fontsize=15, fontweight = 'bold')

# Save the plot
plt.tight_layout()
plt.savefig('figures/Vectra_NEPC/violin_Ar_Postn_TRAMP_stroma.png')


#############################
## neighborhoods enrichment analysis

# all
adata_MISI3542i_NB100_M2861_annot.uns['spatial'] = adata_MISI3542i_NB100_M2861_annot.obsm['spatial']
#adata_MISI3542i_M3056_3_annot.uns['spatial'] = adata_MISI3542i_M3056_3_annot.obsm['spatial']

# stroma
#stroma_MISI3542i_NB100_M2861.uns['spatial'] = stroma_MISI3542i_NB100_M2861.obsm['spatial']
#stroma_MISI3542i_M3056_3.uns['spatial'] = stroma_MISI3542i_M3056_3.obsm['spatial']


# MISI3542i_NB100_M2861: All
sq.gr.spatial_neighbors(adata_MISI3542i_NB100_M2861_annot, spatial_key = 'spatial')
sq.gr.nhood_enrichment(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types")
sq.pl.nhood_enrichment(adata_MISI3542i_NB100_M2861_annot, mode = 'zscore', annotate=False, title='', cluster_key="cell_types", figsize = [8,8], vmin=-100, vmax=100, dpi = 400, fontsize='x-large', save='MISI3542i_NB100_M2861_Neighborhoods.tiff')

# interactions
sq.gr.interaction_matrix(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types")
sq.pl.interaction_matrix(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types", figsize = [10,8], dpi = 400, save='MISI3542i_NB100_M2861_Interactions.tiff')

##################
# cooccurence
sq.gr.co_occurrence(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types", spatial_key='spatial', n_splits=2)

sq.pl.co_occurrence(
   adata_MISI3542i_NB100_M2861_annot,
   cluster_key="cell_types",
   clusters='AR+ stroma',
   palette = 'Accent',
   figsize=(10, 5),
   save='MISI3542i_NB100_M2861_POSTN+stroma_cooccurence.tiff'
)

#########
# centrality
sq.gr.centrality_scores(
    adata_MISI3542i_NB100_M2861_annot,
    cluster_key="cell_types",
)
sq.pl.centrality_scores(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types", figsize=(20, 5), s=500, save='MISI3542i_NB100_M2861_centrality.tiff')


########
# ripley
sq.gr.ripley(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types", mode='L')
sq.pl.ripley(adata_MISI3542i_NB100_M2861_annot, cluster_key="cell_types", mode='L', save='Ripley')

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



# MISI3542i_M3056_3: remove synapto+ cells no negative controls
adata_MISI3542i_NB100_M2861_annot_noNEPC2 = adata_MISI3542i_NB100_M2861_annot2[~adata_MISI3542i_NB100_M2861_annot2.obs['cell_types'].isin(['Synaptophysin+ stroma', 'Synaptophysin+ epithelium'])]
sq.gr.spatial_neighbors(adata_MISI3542i_NB100_M2861_annot_noNEPC2, spatial_key = 'spatial', coord_type = 'generic')
sq.gr.nhood_enrichment(adata_MISI3542i_NB100_M2861_annot_noNEPC2, cluster_key="cell_types")
sq.pl.nhood_enrichment(adata_MISI3542i_NB100_M2861_annot_noNEPC2, mode = 'zscore', annotate=False, cluster_key="cell_types", figsize = [6,6], dpi=500, vmin=-100, vmax=100,  palette = ['blue', 'darkred', 'orange', 'green'], kwargs={"fontsize": "xx-large"}, title='', save='MISI3542i_M3056_3_Neighborhoods_noNEPC2.tiff')


