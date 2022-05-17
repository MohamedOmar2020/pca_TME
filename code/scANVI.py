import os
import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import gdown
import scipy.sparse as sp
import scvi

#scvi.settings.seed = 94705
scvi.settings.batch_size = 10
#scvi.settings.progress_bar_style = "rich"
#scvi.settings.num_threads = 12

#sc.settings.set_figure_params(dpi=200, frameon=False)
#sc.set_figure_params(dpi=200)
#sc.set_figure_params(figsize=(4, 4))
#torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)

## Reference data: mouse
adata_mouse = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)

# filter the na
adata_mouse.obs['cluster'] = adata_mouse.obs['cluster'].astype('str')
adata_mouse = adata_mouse[adata_mouse.obs['cluster'] != 'nan', :]
adata_mouse.obs['cluster'].value_counts()

# raw
adata_mouse_raw = adata_mouse.copy()
adata_mouse_raw = adata_mouse_raw.raw.to_adata()
adata_mouse_raw.X = sp.csr_matrix.todense(adata_mouse_raw.X)
adata_mouse_raw.X = adata_mouse_raw.to_df()

adata_mouse_raw.layers['counts'] = adata_mouse.raw.X
adata_mouse_raw
adata_mouse_raw.obs.model.value_counts()

adata_mouse_raw.obs.drop( ['_scvi_batch','_scvi_labels','_scvi_local_l_mean','_scvi_local_l_var'], axis=1, inplace=True )
del adata_mouse_raw.uns['_scvi']
del adata_mouse_raw.obsm['X_scVI']



##########################
## Target data: human
# load the human data
adata_human = sc.read_h5ad('human/erg_fibroblasts_scvi_v6_regulons.h5ad', chunk_size=100000)
# raw
# save the normalized data (not z-scored) for cellchat
adata_human_raw = adata_human.copy()
adata_human_raw = adata_human_raw.raw.to_adata()
adata_human_raw.X = sp.csr_matrix.todense(adata_human_raw.X)
adata_human_raw.X = adata_human_raw.to_df()

# change the var_names to match mouse gene symbols
adata_human_raw.var_names = [gene.title() for gene in adata_human_raw.var_names]

## find common genes
var_names = adata_mouse_raw.var_names.intersection(adata_human_raw.var_names)
len(var_names)

# subset
adata_mouse_raw = adata_mouse_raw[:, var_names]
adata_human_raw = adata_human_raw[:, var_names]

###############################
# ## Create SCANVI model and train it on fully labelled reference dataset
adata_mouse_raw = remove_sparsity(adata_mouse_raw)
adata_mouse_raw2 = adata_mouse_raw.copy()

scvi.model.SCVI.setup_anndata(adata_mouse_raw2, batch_key="model", layer="counts")

arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=1,
    log_variational = 'False'
)

vae_ref = scvi.model.SCVI(
    adata_mouse_raw2,
    **arches_params
)
vae_ref.train()


# save the reference model
dirpath = 'scvi'
vae_ref.save('scvi/', overwrite=True)

adata_mouse_raw2.obsm["X_scVI"] = vae_ref.get_latent_representation()
sc.pp.neighbors(adata_mouse_raw2, use_rep="X_scVI")
sc.tl.leiden(adata_mouse_raw2)
sc.tl.umap(adata_mouse_raw2)


sc.pl.umap(
    adata_mouse_raw2,
    color=["cluster", "model"],
    frameon=False,
    ncols=2,
)

##############################################################
# Update with query
adata_human_raw.layers['counts'] = adata_human_raw.X
adata_query = adata_human_raw.copy()

adata_query.obs.drop( ['_scvi_batch','_scvi_labels'], axis=1, inplace=True )
del adata_query.uns['_scvi']
del adata_query.obsm['X_scVI']


# validate that our query data is ready to be loaded into the reference model
scvi.model.SCVI.prepare_query_anndata(adata_query, dirpath)
scvi.model.SCVI.prepare_query_anndata(adata_query, vae_ref)

# Now we create the new query model instance.

vae_q = scvi.model.SCVI.load_query_data(
    adata_query,
    dirpath,
)
vae_q = scvi.model.SCVI.load_query_data(
    adata_query,
    vae_ref,
)
