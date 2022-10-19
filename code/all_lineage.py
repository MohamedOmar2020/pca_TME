
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
import seaborn as sb

sc.settings.figdir = 'figures'
sc.set_figure_params(dpi_save = 300, transparent = False)


# load the mouse data
adata_mouse_all = sc.read_h5ad('outs/h5ads/fapcm_all_v6.h5ad', chunk_size=100000)


adata_mouse_all.obs['leiden'] = adata_mouse_all.obs['leiden'].astype('str')
adata_mouse_all.obs['leiden'].value_counts()

adata_mouse_all.obs['key'] = adata_mouse_all.obs['key'].astype('str')
adata_mouse_all.obs['key'].value_counts()

# ##############################
# # annotation
# adata_mouse_all.obs["B_cells"] = (
#     adata_mouse_all.obs["leiden"].isin(["28"]).astype("category")
# )
#
# adata_mouse_all.obs["cDCs"] = (
#     adata_mouse_all.obs["leiden"].isin(["17"]).astype("category")
# )
#
# adata_mouse_all.obs["endothelial"] = (
#     adata_mouse_all.obs["leiden"].isin(["7", '36']).astype("category")
# )
#
# adata_mouse_all.obs["myeloid"] = (
#     adata_mouse_all.obs["leiden"].isin(["5", '30', '34']).astype("category")
# )
#
# adata_mouse_all.obs["pDCs"] = (
#     adata_mouse_all.obs["leiden"].isin(['3', '12']).astype("category")
# )
#
# adata_mouse_all.obs["t_NK"] = (
#     adata_mouse_all.obs["leiden"].isin(["20", "4", '11', '14']).astype("category")
#)

adata_mouse_all.obs["epithelium"] = (
    adata_mouse_all.obs["leiden"]
        .isin(["27", "24", "19", "18", "31", "11", "16", "21", "22", "6", "33", "23", "2", "29", '40', '38'])
        .astype("category")
)

adata_mouse_all.obs["stroma"] = (
    adata_mouse_all.obs["leiden"].isin(["8", "13", "15", "26", '32']).astype("category")
)

adata_mouse_all.obs["immune"] = (
    adata_mouse_all.obs["leiden"].isin(["28", '17', '5', '30', '34', '35', '3', '12', '20', '4', '14', '10', '37', '39', '41', '0', '9']).astype("category")
)

adata_mouse_all.obs["endothelium"] = (
     adata_mouse_all.obs["leiden"].isin(["7", '36']).astype("category")
 )
#adata.obs["basal"] = adata.obs["leiden"].isin(["2", "29"]).astype("category")


adata_mouse_all.obs["compartments"] = "Other"
adata_mouse_all.obs.loc[adata_mouse_all.obs.stroma == True, "compartments"] = "stroma"
adata_mouse_all.obs.loc[adata_mouse_all.obs.epithelium == True, "compartments"] = "epithelium"
adata_mouse_all.obs.loc[adata_mouse_all.obs.immune == True, "compartments"] = "immune"
adata_mouse_all.obs.loc[adata_mouse_all.obs.endothelium == True, "compartments"] = "endothelium"

adata_mouse_all.obs['compartments'].value_counts()
pd.crosstab( adata_mouse_all.obs['leiden'], adata_mouse_all.obs['compartments'])



adata_mouse_all.obs['b'].describe()
adata_mouse_all[adata_mouse_all.obs["compartments"] =='B_cells'].obs['b'].describe()

adata_mouse_all.obs['fibroblast'].describe()
adata_mouse_all[adata_mouse_all.obs["compartments"] =='fribroblasts'].obs['fibroblast'].describe()

adata_mouse_all[adata_mouse_all.obs['fibroblast'] > adata_mouse_all.obs['fibroblast'].quantile(0.90)]


#adata.obs.loc[adata.obs.basal == True, "phenotypes"] = "basal"
sc.pl.umap(
    adata_mouse_all,
    color="leiden",
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="leiden.png"
)

sc.pl.umap(
    adata_mouse_all,
    color="compartments",
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="compartments.png"
)

sc.pl.umap(
    adata_mouse_all,
    color=['compartments', 'dendritic', 'cCDs', 'langherhans_like', 'b', 't_nk', 'myeloid', 'mast', 'macrophages'],
    cmap="Reds",
    vmax=3,
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_immune.png"
)

sc.pl.umap(
    adata_mouse_all,
    color=['compartments', 'endothelial'],
    cmap="Reds",
    vmax=3,
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_endothelial.png"
)

sc.pl.umap(
    adata_mouse_all,
    color=['compartments', 'fibroblast', 'myofibroblast'],
    cmap="Reds",
    vmax=3,
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_stroma.png"
)

sc.pl.umap(
    adata_mouse_all,
    color=['compartments', 'luminal', 'basal', 'notluminal', 'neuroendocrine', 'seminal_vesicle_basal', 'seminal_vesicle_luminal', 'seminal_vesicle_ionocyte'],
    cmap="Reds",
    vmax=3,
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_epithelium.png"
)


# save
adata_mouse_all.write_loom('outs/h5ads/fapcm_all_v6_compartments_MO.loom', write_obsm_varm=True)

##########################################################################################
##########################################################################################
# load the human data
adata_human_all = sc.read_h5ad('outs_human/h5ads/erg_scvi_v6.h5ad', chunk_size=100000)
adata_human_all.obs['erg'].value_counts()
adata_human_all.obs['key'].value_counts()
adata_human_all.obs['tn'].value_counts()

######################
## annotation

sc.pl.umap(
    adata_human_all,
    color="leiden",
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_human_leiden.png"
)

sc.pl.umap(adata_human_all, color = ['endothelial'])
sc.pl.umap(adata_human_all, color = ['fibroblast', 'myofibroblast'])
sc.pl.umap(adata_human_all, color = ['dendritic', 'pDC', 'cCDs', 'langherhans_like', 'b', 't_nk', 'myeloid', 'mast', 'macrophages'])
sc.pl.umap(adata_human_all, color = ['luminal', 'basal', 'notbasal', 'neuroendocrine', 'seminal_vesicle_basal', 'seminal_vesicle_luminal', 'seminal_vesicle_ionocyte'])


adata_human_all.obs["epithelium"] = (
    adata_human_all.obs["leiden"].isin(['27', '5', '23', '20', '31', '28', '29', '40', '8', '36', '12', '39', '55', '26', '48', '61', '66']).astype("category")
)

adata_human_all.obs["stroma"] = (
    adata_human_all.obs["leiden"].isin(['16', '10', '9', '62']).astype("category")
)

adata_human_all.obs["immune"] = (
    adata_human_all.obs["leiden"].isin(['1', '2', '13', '3', '14', '22', '21', '18', '19', '24', '15', '7', '17', '32', '35', '37', '25']).astype("category")
)

adata_human_all.obs["endothelium"] = (
     adata_human_all.obs["leiden"].isin(['0', '44', '11', '30']).astype("category")
 )
#adata.obs["basal"] = adata.obs["leiden"].isin(["2", "29"]).astype("category")


adata_human_all.obs["compartments"] = "Other"
adata_human_all.obs.loc[adata_human_all.obs.stroma == True, "compartments"] = "stroma"
adata_human_all.obs.loc[adata_human_all.obs.epithelium == True, "compartments"] = "epithelium"
adata_human_all.obs.loc[adata_human_all.obs.immune == True, "compartments"] = "immune"
adata_human_all.obs.loc[adata_human_all.obs.endothelium == True, "compartments"] = "endothelium"


adata_human_all.obs['compartments'].value_counts()
pd.crosstab(adata_human_all.obs['leiden'], adata_human_all.obs['compartments'])
pd.crosstab(adata_human_all.obs['leiden'][adata_human_all.obs['leiden']=='Other'], adata_human_all.obs['compartments'])



#adata.obs.loc[adata.obs.basal == True, "phenotypes"] = "basal"
sc.pl.umap(
    adata_human_all,
    color="leiden",
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_human_leiden.png"
)

sc.pl.umap(
    adata_human_all,
    color="compartments",
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_human_compartments.png"
)

sc.pl.umap(
    adata_human_all,
    color=['compartments', 'dendritic', 'cCDs', 'langherhans_like', 'b', 't_nk', 'myeloid', 'mast', 'macrophages'],
    cmap="Reds",
    vmax=3,
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_human_immune.png"
)

sc.pl.umap(
    adata_human_all,
    color=['compartments', 'endothelial'],
    cmap="Reds",
    vmax=3,
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_human_endothelial.png"
)

sc.pl.umap(
    adata_human_all,
    color=['compartments', 'fibroblast', 'myofibroblast'],
    cmap="Reds",
    vmax=3,
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_human_stroma.png"
)

sc.pl.umap(
    adata_human_all,
    color=['compartments', 'luminal', 'basal', 'neuroendocrine', 'seminal_vesicle_basal', 'seminal_vesicle_luminal', 'seminal_vesicle_ionocyte'],
    cmap="Reds",
    vmax=3,
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save="_human_epithelium.png"
)



#################################################################
# adenosine related genes

#PPARG
#CYBB
#COL3A1
#FOXP3
#LAG3
#APP
#CD81
#GPI
#PTGS2
#CASP1
#FOS
#MAPK1
#MAPK3
#CREB1


# dotplot for adenosine related genes per different mouse models
dp = sc.pl.DotPlot(adata_mouse_all, var_names = ['Pparg', 'Cybb', 'Col3a1', 'Foxp3', 'Lag3', 'App', 'Cd81', 'Gpi1', 'Ptgs2', 'Casp1', 'Fos', 'Mapk1', 'Mapk3', 'Creb1'], groupby = ['key', 'compartments'], cmap = 'Reds', use_raw=True)
dp.legend(width=2.5).savefig('figures/dotplot_adenosine_mouse_all.png')

################
# dotplot for adenosine related genes per erg status
dp = sc.pl.DotPlot(adata_human_all, var_names = ['PPARG', 'CYBB', 'COL3A1', 'FOXP3', 'LAG3', 'APP', 'CD81', 'GPI', 'PTGS2', 'CASP1', 'FOS', 'MAPK1', 'MAPK3', 'CREB1'], groupby = ['erg', 'compartments'], cmap = 'Reds', use_raw=True)
dp.legend(width=2.5).savefig('figures/dotplot_adenosine_human_all_erg.png')

# dotplot for adenosine related genes per case
dp = sc.pl.DotPlot(adata_human_all, var_names = ['PPARG', 'CYBB', 'COL3A1', 'FOXP3', 'LAG3', 'APP', 'CD81', 'GPI', 'PTGS2', 'CASP1', 'FOS', 'MAPK1', 'MAPK3', 'CREB1'], groupby = ['case', 'compartments'], cmap = 'Reds', use_raw=True)
dp.legend(width=2.5).savefig('figures/dotplot_adenosine_human_all_case.png')

# dotplot for adenosine related genes per case
dp = sc.pl.DotPlot(adata_human_all, var_names = ['PPARG', 'CYBB', 'COL3A1', 'FOXP3', 'LAG3', 'APP', 'CD81', 'GPI', 'PTGS2', 'CASP1', 'FOS', 'MAPK1', 'MAPK3', 'CREB1'], groupby = ['condition', 'compartments'], cmap = 'Reds', use_raw=True)
dp.legend(width=2.5).savefig('figures/dotplot_adenosine_human_all_site.png')

## ############
# save
adata_human_all.write_loom('outs_human/h5ads/erg_scvi_v6_compartments.loom', write_obsm_varm=True)
