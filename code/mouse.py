
import sys
from pathlib import Path
import scanpy as sc
from matplotlib.pyplot import ion
from scutils.figures.base import basics
from scutils.figures.prostate import annotate_cell_types_prostate
from scutils.qc import PreprocessRNA


adata_mouse = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)
adata_mouse.obs['cluster'].value_counts()

adata_human = sc.read_h5ad('human/erg_fibroblasts_scvi_v6_regulons.h5ad', chunk_size=100000)

adata_human.var_names = [gene.lower() for gene in adata_human.var_names]
adata_mouse.var_names = [gene.lower() for gene in adata_mouse.var_names]


####################
## scanpy ingest

# find common genes
var_names = adata_mouse.var_names.intersection(adata_human.var_names)

# subset
adata_mouse = adata_mouse[:, var_names]
adata_human = adata_human[:, var_names]

# mapping the clusters from mouse to human using ingest
sc.pp.pca(adata_mouse)
sc.pp.neighbors(adata_mouse)
sc.tl.umap(adata_mouse)
sc.pl.umap(adata_mouse, color='cluster')

del adata_human.obsp['connectivities']
del adata_human.obsp['distances']
del adata_human.obsm['X_pca']
del adata_human.obsm['X_umap']
del adata_human.uns['neighbors']

sc.tl.ingest(adata_human, adata_mouse, obs='cluster')






















