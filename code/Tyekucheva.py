

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


## load the data
# mouse

sc.settings.figdir = 'figures'
sc.set_figure_params(dpi_save = 300)
#plt.rcParams.update({'xtick.labelsize' : '50'})
#plt.rcParams.update({'ytick.labelsize' : '50'})

adata_mouse = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)
adata_mouse.obs['cluster'] = adata_mouse.obs['cluster'].astype('str')
adata_mouse.obs['cluster'].value_counts()

# filter the na
adata_mouse = adata_mouse[adata_mouse.obs['cluster'] != 'nan', :]
adata_mouse.obs['cluster'].value_counts()

##############################
# load the annot human data

adata_BoneMet_stroma = sc.read_h5ad('data/kfoury/adata_BoneMet_stroma.h5ad')

# re-cap the gene symbols
adata_BoneMet_stroma.var_names = [gene.upper() for gene in adata_BoneMet_stroma.var_names]
# subset also the adata_human.raw
tempAdata = adata_BoneMet_stroma.raw.to_adata()
tempAdata.var_names = [gene.upper() for gene in tempAdata.var_names]
tempAdata.var_names_make_unique()
adata_BoneMet_stroma.raw = tempAdata

##########################################
# the stromal gene signature
stromaSignature_mouse = ['Aebp1', 'Antxr1', 'Bgn', 'C1qa', 'C1qb', 'C1qc', 'C1s1', 'C1ra', 'Cdh11', 'Col1a1', 'Col3a1', 'Fbln5', 'Fcgr2b', 'H2-Eb1', 'Lum', 'Moxd1', 'Prelp', 'Rnase1', 'Sfrp2', 'Sfrp4', 'Sulf1', 'Thbs2', 'Thy1', 'Tmsb4x']

stromaSignature_human = ['AEBP1', 'ANTXR1', 'BGN', 'C1QA', 'C1QB', 'C1QC', 'C1S', 'C1R', 'CDH11', 'COL1A1', 'COL3A1', 'FBLN5', 'FCGR2B', 'HLA-DRB1', 'LUM', 'MOXD1', 'PRELP', 'RNASE1', 'SFRP2', 'SFRP4', 'SULF1', 'THBS2', 'THY1', 'TMSB4X']

###############################
## dotplot mouse
dp_mouse = sc.pl.DotPlot(adata_mouse, var_names = stromaSignature_mouse,
                                #categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp_mouse.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/Tyekucheva_Signature_mouse.png')

sc.pl.umap(
    adata_mouse,
    color=stromaSignature_mouse,
    cmap="RdBu_r",
    vmax=5,
    save="_Tyekucheva_signature_mouse.png"
)

###############
## plot bone metastasis
dp_boneMets = sc.pl.DotPlot(adata_BoneMet_stroma, var_names = stromaSignature_human,
                                #categories_order = ['0','1','2','3','4','5','6','7'],
                                groupby='cluster', cmap = 'Reds'
                                )
dp_boneMets.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/Tyekucheva_Signature_BoneMets.png')

sc.pl.umap(
    adata_BoneMet_stroma,
    color=stromaSignature_human,
    cmap="RdBu_r",
    vmax=5,
    save="_Tyekucheva_signature_boneMets.png"
)



