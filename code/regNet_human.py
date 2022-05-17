# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import anndata as ad
from MulticoreTSNE import MulticoreTSNE as TSNE
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from cytoolz import compose
import loompy
from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss
import anndata as ad
import scipy.sparse as sp
import sys

#matplotlib.use('module://backend_interagg')

# read unfiltered data from a loom file
adata_human = sc.read_h5ad("data/scenic/erg_fibroblasts_scvi_v6_regulons_annot_new2.h5ad", chunk_size=100000)
adata_human = adata_human.raw.to_adata()
adata_human.X = sp.csr_matrix.todense(adata_human.X)
adata_human
adata_human.var_names_make_unique()

## Write to an unfiltered loom file
row_attrs = {"Gene": np.array(adata_human.var.index) ,}

col_attrs = {
    "CellID":  np.array(adata_human.obs.index) ,
    "nGene": np.array( np.sum(adata_human.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata_human.X.transpose() , axis=0)).flatten() ,
}

lp.create('data/scenic/SC_human_unfil.loom', adata_human.X.transpose(), row_attrs, col_attrs )

########################################################
## Initial/basic filtering:
nCountsPerGene = np.sum(adata_human.X, axis=0)
nCellsPerGene = np.sum(adata_human.X > 0, axis=0)
# Show info
print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - ", nCountsPerGene.max())
print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - ", nCellsPerGene.max())

nCells = adata_human.X.shape[0]

# pySCENIC thresholds
minCountsPerGene = 3 * .01 * nCells  # 3 counts in 1% of cells
print("minCountsPerGene: ", minCountsPerGene)

minSamples = .01 * nCells  # 1% of cells
print("minSamples: ", minSamples)

#################
# Basic filtering:

# All
sc.pp.filter_cells(adata_human, min_genes=300)
sc.pp.filter_genes(adata_human, min_cells=20)

# mito and genes/counts cuts
adata_human.var['mt'] = adata_human.var_names.str.startswith('mt-')

sc.pp.calculate_qc_metrics(adata_human, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata_human, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

##########
# Remove cells that have too many mitochondrial genes expressed or too many total counts:
sc.pl.scatter(adata_human, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata_human, x='total_counts', y='n_genes_by_counts')

# All
adata_human = adata_human[adata_human.obs.n_genes_by_counts < 4500, :]
adata_human = adata_human[adata_human.obs.pct_counts_mt < 20, :]

###############
## Diagnostic plots, post filtering:

sc.pl.violin(adata_human, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

## Save this loom file for Scenic
row_attrs = {"Gene": np.array(adata_human.var_names) ,}

col_attrs = {
    "CellID": np.array(adata_human.obs_names) ,
    "nGene": np.array( np.sum(adata_human.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata_human.X.transpose() , axis=0)).flatten() ,
}

lp.create('objs/Scenic/normal_tumor/SC_human_fil.loom', adata_human.X.transpose(), row_attrs, col_attrs)


#############################################################################
## Further pre-processing of expression data:
# << I guess we don't need this if we use the already filtred and annotated loom file >>

# save a copy of the raw data
adata_human.raw = adata_human # You can get back an AnnData of the object in .raw by calling .raw.to_adata().

# Total-count normalize (library-size correct) to 10,000 reads/cell
sc.pp.normalize_per_cell(adata_human, counts_per_cell_after=1e4)

# log transform the data.
sc.pp.log1p(adata_human)

# identify highly variable genes.
sc.pp.highly_variable_genes(adata_human, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_human)

# keep only highly variable genes:
adata_human = adata_human[:, adata_human.var['highly_variable']]

# regress out total counts per cell and the percentage of mitochondrial genes expressed
sc.pp.regress_out(adata_human, ['total_counts', 'pct_counts_mt'] , n_jobs=12)

# scale each gene to unit variance, clip values exceeding SD 10.
sc.pp.scale(adata_human, max_value=10)

#############
# principal component analysis
sc.tl.pca(adata_human, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_human, log=True)

# batch correction
sc.external.pp.harmony_integrate(adata_human, key='origin', max_iter_harmony = 20)

#adata_human.write_loom("objs/Scenic/all/all_SC_processed.loom", write_obsm_varm = True)
adata_human.write("objs/Scenic/normal_tumor/SC_human_harmony.h5ad")

# Read the saved adata object
adata_human = sc.read_h5ad('objs/Scenic/normal_tumor/SC_human_harmony.h5ad', chunk_size=100000)

##############
# Visualization of highly variable genes

# neighborhood graph of cells (determine optimal number of PCs here)
sc.pp.neighbors(adata_human, n_neighbors=15, n_pcs=30, use_rep='X_pca_harmony')

# compute UMAP
sc.tl.umap(adata_human)

# tSNE
#tsne = TSNE(n_jobs=12)
#adata_human.obsm['X_tsne'] = tsne.fit_transform(adata_human.X)

adata_human.write_loom("objs/Scenic/normal_tumor/SC_human_processed.loom", write_obsm_varm = True)
adata_human.write("objs/Scenic/normal_tumor/SC_human_processed.h5ad")

#######################
## Clustering:

# cluster the neighbourhood graph

sc.tl.louvain(adata_human,resolution=0.2)

#with rc_context({'figure.figsize': (10, 10)}):
#sc.pl.umap(adata_normal, color=['louvain'], save = "Normal_umap_louvain.png")

# find marker genes
sc.tl.rank_genes_groups(adata_human, 'louvain', method='t-test')
#sc.pl.rank_genes_groups(adata_human, n_genes=25, sharey=False)

# sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
#Markers = pd.DataFrame(adata_human.uns['rank_genes_groups']['names']).head(10)

adata_human.write_loom("objs/Scenic/normal_tumor/SC_human_processed.loom", write_obsm_varm=True)
adata_human.write("objs/Scenic/normal_tumor/SC_human_processed.h5ad")

####################################################################################
####################################################################################
## AFTER SCENIC steps:

#########################
## Visualization of SCENIC's AUC matrix
import json
import zlib
import base64

# Collect SCENIC AUCell output
lf = lp.connect('objs/Scenic/normal_tumor/normal_tumor_ScenicOutput.loom', mode='r+', validate=False )
auc_mtx = pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

## I guess we need to filter this to have just the cells in the processed loom
cellID_obs = pd.read_csv("objs/seurat_info/normal_tumor_cellID.csv")
cellID_obs.shape

# Filter the cells to keep only those in the seurat object
auc_mtx.shape
auc_mtx_fil = auc_mtx[np.isin(auc_mtx.index,cellID_obs["x"])]
auc_mtx_fil.shape



##################
## Run umap and tsne on the filted auc matrix
import umap

# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap(auc_mtx_fil)
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx_fil.index).to_csv("objs/Scenic/normal_tumor/normal_tumor_scenic_umap.txt", sep='\t')

# tSNE
tsne = TSNE(n_jobs = 14)
dr_tsne = tsne.fit_transform(auc_mtx_fil)
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx_fil.index).to_csv( "objs/Scenic/normal_tumor/normal_tumor_scenic_tsne.txt", sep='\t')

###############################
## Integrate the output:

# scenic output
# Filter the loom file to have just the cells in the processed loom (cellID) and also the cells in auc_mtx_fil
all_ScenicLoom = anndata.read_loom('objs/Scenic/normal_tumor/normal_tumor_ScenicOutput.loom', obs_names='obs_names', validate = False)
all_ScenicLoom = all_ScenicLoom[np.isin(all_ScenicLoom.obs.CellID, cellID_obs["x"])]
all_ScenicLoom

all_ScenicLoom.write_loom("objs/Scenic/normal_tumor/normal_tumor_ScenicLoom_fil.loom", write_obsm_varm=True)

lf = lp.connect('objs/Scenic/normal_tumor/normal_tumor_ScenicOutput.loom', mode='r+', validate=False )
lf_fil = lp.connect('objs/Scenic/normal_tumor/normal_tumor_ScenicLoom_fil.loom', mode='r+', validate=False )


meta = json.loads(zlib.decompress(base64.b64decode(lf.attrs.MetaData)))

#exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
#auc_mtx_fil = pd.DataFrame(lf_fil.ca.RegulonsAUC, index=lf_fil.ca.CellID)
regulons = lf.ra.Regulons
dr_umap = pd.read_csv('objs/Scenic/normal_tumor/normal_tumor_scenic_umap.txt', sep='\t', header=0, index_col=0 )
dr_tsne = pd.read_csv('objs/Scenic/normal_tumor/normal_tumor_scenic_tsne.txt', sep='\t', header=0, index_col=0 )

####################################################
## Fix regulon objects to display properly in SCope:
auc_mtx_fil.columns = auc_mtx_fil.columns.str.replace('\(','_(')
regulons.dtype.names = tuple([ x.replace("(","_(") for x in regulons.dtype.names ] )

# regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
    tmp = x.get('regulon').replace("(","_(")
    x.update( {'regulon': tmp} )


#################################
## Concatenate embeddings (tSNE, UMAP, etc.) >> ''From the processed loom/ Not the scenic ones '

## Read the annotated and filtered loom file
#loom_filtered_annot = anndata.read_loom('objs/loom/loom_filtered_annot.loom', obs_names="obs_names")
#loom_filtered_annot

procLoom = anndata.read_loom('objs/Scenic/normal_tumor/SC_human_processed.loom', obs_names="obs_names")
procLoom

## Filtering:

# Load the cell IDs from the seurat object
cellID_obs = pd.read_csv("objs/seurat_info/normal_tumor_cellID.csv")

# Read the umap coordinates from the seurat object metadata
umap_cord = pd.read_csv("objs/seurat_info/normal_tumor_umap.csv")

# Read the clusters information from the seurat object metadata
cell_clusters = pd.read_csv("objs/seurat_info/normal_tumor_clusters.csv")

# Filter the cells to keep only those in the auc_mtx_fil
procLoom_fil = procLoom[np.isin(procLoom.obs.index,cellID_obs["x"])]
procLoom_fil

## Filter the meta-data because they have some extra cells
cellID_obs = cellID_obs[np.isin(cellID_obs["x"], auc_mtx_fil.index)]
umap_cord = umap_cord[np.isin(umap_cord['Unnamed: 0'], auc_mtx_fil.index)]
cell_clusters = cell_clusters[np.isin(cell_clusters['CellID'], auc_mtx_fil.index)]

###### Attach the umap coords and cluster labels to procLoom_fil
## Change the name of the columns in loom and umap to Cell ID
procLoom_fil_index = pd.DataFrame(procLoom_fil.obs.index)
procLoom_fil_index = procLoom_fil_index.rename(columns = {'obs_names':'Cell ID'})
umap_cord = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered = procLoom_fil_index.merge(umap_cord, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]

## Add the umpa coordinate info to the loom file
procLoom_fil.obsm['X_umap'] = umap_ordered.values
procLoom_fil

# Change the column names of cell clusters dataframe
cell_clusters = cell_clusters.rename(columns = {'CellID':'Cell ID'})
cell_clusters = cell_clusters.rename(columns = {'Seurat_Normal.celltype':'cluster'})
cell_clusters_ordered = procLoom_fil_index.merge(cell_clusters, on = "Cell ID")
cell_clusters_ordered = cell_clusters_ordered.iloc[:,1:]

# Add the cluster info to the loom file
procLoom_fil.obs['clusters'] = cell_clusters_ordered.values
procLoom_fil

#tsneDF = pd.DataFrame(loom_filtered_annot.obsm['X_tsne'], columns=['_X', '_Y'])

Embeddings_X = pd.DataFrame( index=lf_fil.ca.CellID )

Embeddings_X = pd.concat( [
        pd.DataFrame(procLoom_fil.obsm['X_umap'],index=procLoom_fil.obs.index)[0] ,
        pd.DataFrame(procLoom_fil.obsm['X_pca'],index=procLoom_fil.obs.index)[0] ,
        dr_tsne['X'] ,
        dr_umap['X']
    ], sort=False, axis=1, join='outer' )

Embeddings_X.columns = ['1','2','3','4']

Embeddings_Y = pd.DataFrame(index=lf_fil.ca.CellID )

Embeddings_Y = pd.concat( [
        pd.DataFrame(procLoom_fil.obsm['X_umap'],index=procLoom_fil.obs.index)[1] ,
        pd.DataFrame(procLoom_fil.obsm['X_pca'],index=procLoom_fil.obs.index)[1] ,
        dr_tsne['Y'] ,
        dr_umap['Y']
    ], sort=False, axis=1, join='outer' )

Embeddings_Y.columns = ['1','2','3','4']

###############
### metadata
metaJson = {}

metaJson['embeddings'] = [
    {
        "id": -1,
        "name": f"Scanpy t-SNE (highly variable genes)"
    },
    {
        "id": 1,
        "name": f"Scanpy UMAP  (highly variable genes)"
    },
    {
        "id": 2,
        "name": "Scanpy PC1/PC2"
    },
    {
        "id": 3,
        "name": "SCENIC AUC t-SNE"
    },
    {
        "id": 4,
        "name": "SCENIC AUC UMAP"
    },
]

metaJson["clusterings"] = [{
    "id": 0,
    "group": "Scanpy",
    "name": "Scanpy louvain default resolution",
    "clusters": [],
}]

metaJson["metrics"] = [
    {
        "name": "nUMI"
    }, {
        "name": "nGene"
    }, {
        "name": "Percent_mito"
    }
]

metaJson["annotations"] = [
    {
        "name": "Louvain_clusters_Scanpy",
        "values": list(set(procLoom_fil.obs['clusters'].astype(np.str)))
    },
    # {
    #    "name": "Genotype",
    #    "values": list(set(adata.obs['Genotype'].values))
    # },
    # {
    #    "name": "Timepoint",
    #    "values": list(set(adata.obs['Timepoint'].values))
    # },
    # {
    #    "name": "Sample",
    #    "values": list(set(adata.obs['Sample'].values))
    # }
]

# SCENIC regulon thresholds:
metaJson["regulonThresholds"] = rt

for i in range(max(set([int(x) for x in procLoom_fil.obs['clusters']])) + 1):
    clustDict = {}
    clustDict['id'] = i
    clustDict['description'] = f'Unannotated Cluster {i + 1}'
    metaJson['clusterings'][0]['clusters'].append(clustDict)

clusterings = pd.DataFrame()
clusterings["0"] = procLoom_fil.obs['clusters'].values.astype(np.str)

########
## Assemble loom file row and column attributes
def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.values]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr

col_attrs = {
    "CellID": np.array(procLoom_fil.obs.index),
    "nUMI": np.array(procLoom_fil.obs['n_counts'].values),
    "nGene": np.array(procLoom_fil.obs['n_genes'].values),
    "Louvain_clusters_Scanpy": np.array(procLoom_fil.obs['clusters'].values ),
    #"Genotype": np.array(adata.obs['Genotype'].values),
    #"Timepoint": np.array(adata.obs['Timepoint'].values),
    #"Sample": np.array(adata.obs['Sample'].values),
    "Percent_mito": np.array(procLoom_fil.obs['pct_counts_mt'].values),
    #"Embedding": dfToNamedMatrix(tsneDF),
    "Embeddings_X": dfToNamedMatrix(Embeddings_X),
    "Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
    "RegulonsAUC": dfToNamedMatrix(auc_mtx_fil),
    #"Clusterings": dfToNamedMatrix(clusterings),
    "ClusterID": np.array(procLoom_fil.obs['clusters'].values)
}

row_attrs = {
    "Gene": lf_fil.ra.Gene,
    "Regulons": regulons,
}

attrs = {
    "title": "Scenic",
    "MetaData": json.dumps(metaJson),
    "Genome": 'mm9',
    "SCopeTreeL1": "",
    "SCopeTreeL2": "",
    "SCopeTreeL3": ""
}

# compress the metadata field:
attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

##########
## Create a new loom file, copying the expression matrix from the open loom connection:

lp.create(
    filename = 'objs/Scenic/all/all_IntegScenic.loom' ,
    layers=lf_fil[:,:],
    row_attrs=row_attrs,
    col_attrs=col_attrs,
    file_attrs=attrs
)

lf_fil.close() # close original pyscenic loom file
lf.close()


##############
## Check the saved loom file
IntegScenicLoom = anndata.read_loom('objs/Scenic/normal_tumor/normal_tumor_IntegScenic.loom', obs_names='CellID', validate = False, sparse=False)
IntegScenicLoom

