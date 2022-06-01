
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
sc.settings.figdir = ''


sc._settings.ScanpyConfig(figdir = 'figs_new', plot_suffix = '', n_jobs = 20, max_memory = 100)
sc.set_figure_params(dpi_save = 300)

#%%
adata_mouse_stroma = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)
adata_mouse_stroma.obs['cluster'] = adata_mouse_stroma.obs['cluster'].astype('str')
adata_mouse_stroma.obs['cluster'].value_counts()


adata_mouse_unfiltered = sc.read_h5ad('outs/h5ads/fapcm_unfiltered_v6.h5ad', chunk_size=100000)

adata_mouse_all = sc.read_h5ad('outs/h5ads/fapcm_all_v6.h5ad', chunk_size=100000)

adata_mouse_loda = sc.read_h5ad('outs/h5ads/prostate_mouse_loda.h5ad', chunk_size=100000)


# N of cells in table_3_alignment_qc.xlsx == 123652 cells
# N of cells in outs/h5ads/fapcm_unfiltered_v6.h5ad = 101853 cells == the N of cells in the manuscript





