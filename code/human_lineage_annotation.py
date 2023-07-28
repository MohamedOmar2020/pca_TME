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
#import sc_toolbox.api as sct
import seaborn as sb
import plotly.express as px
import squidpy as sq
from sccoda.util import cell_composition_data as scc_dat
from sccoda.util import data_visualization as scc_viz

sc.settings.figdir = 'figures/human_compartments'
sc.set_figure_params(dpi_save = 300, transparent = False, fontsize =8, format='tiff')
plt.rcParams["font.family"] = "arial"
plt.rcParams['font.size'] = 8
plt.rcParams['font.style'] = 'italic'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['legend.fontsize'] = 8

############################################
# load the human data
adata_human_all = sc.read_h5ad("outs_human/h5ads/erg_scvi_v6.h5ad", chunk_size=100000)


adata_human_all.obs['leiden'] = adata_human_all.obs['leiden'].astype('str')
adata_human_all.obs['leiden'].value_counts()

adata_human_all.obs['key'] = adata_human_all.obs['key'].astype('str')
adata_human_all.obs['key'].value_counts()

###############################################
## annotate compartments

adata_human_all.obs["epithelium"] = (
    adata_human_all.obs["leiden"].isin(['27', '5', '23', '20', '31', '28', '29', '40', '8', '36', '12', '39', '55', '26', '48', '61', '66']).astype("category")
)

adata_human_all.obs["stroma"] = (
    adata_human_all.obs["leiden"].isin(['16', '10', '9', '62', '24']).astype("category")
)

adata_human_all.obs["immune"] = (
    adata_human_all.obs["leiden"].isin(['1', '2', '13', '3', '14', '22', '21', '18', '19', '15', '7', '17', '32', '35', '37', '25']).astype("category")
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

sc.pl.umap(adata_human_all, color = 'compartments', legend_loc = 'on data', save= '_human_compartments')


##########################################################
## annotate immune cells
human_immune = adata_human_all[adata_human_all.obs.compartments == 'immune']

human_immune_raw = human_immune.raw.to_adata()
human_immune_raw.write_h5ad('outs_human/h5ads/human_immune.h5ad')

# recluster c7
sc.tl.leiden(human_immune, restrict_to = ('leiden', ['7']), resolution=0.1)

human_immune.obs.leiden.value_counts()
human_immune.obs.leiden_R.value_counts()
sc.tl.rank_genes_groups(human_immune, 'leiden_R', pts=True, use_raw = True)
markers_immune = sc.get.rank_genes_groups_df(human_immune, group=None)

###############
sc.pl.umap(
    human_immune,
    color=["leiden_R", "CXCR3", "CLEC9A", "CD74", "XCR1", "CCR7", "CCL17", "CCL19", "CD207", "CCL17", 'dendritic', 'cCDs'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_dendritic'
)

sc.pl.umap(
    human_immune,
    color=["leiden_R", "CD79A", "CD79B", 'b'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_Bcells'
)

sc.pl.umap(
    human_immune,
    color=["leiden_R", "FOXP3", "CTLA4", 'TNFRSF4', 'IRF4', 'BATF', 'TNFRSF18', 'TOX2', 'PRDM1'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_Treg'
)

sc.pl.umap(
    human_immune,
    color=["leiden_R", "LEF1", "ATM", 'SELL', 'KLF2', 'ITGA6', 'IFNGR2', 'IL21R', 'FOXP1'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_naiveCD4Tcells'
)

sc.pl.umap(
    human_immune,
    color=["leiden_R", "CD68", "ADGRE1", "CD163", "MRC1", "CD14", "ITGAM", "CCR2", "FCGR1A", "MERTK", "ARG1", "NOS2"],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_macrophages'
)


sc.pl.umap(
    human_immune,
    color=["leiden_R", 't_nk', 'NCR1', 'PRF1'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6, use_raw = False,
    save = '_NKcells.png'
)


human_immune.obs["dendritic cells"] = (
    human_immune.obs["leiden_R"].isin(['7,0', '7,1']).astype("category")
)

human_immune.obs["B cells"] = (
    human_immune.obs["leiden_R"].isin(['25']).astype("category")
)

human_immune.obs["Treg"] = (
    human_immune.obs["leiden_R"].isin(['19']).astype("category")
)

human_immune.obs["CD4+ T lymphocytes"] = (
    human_immune.obs["leiden_R"].isin(['1']).astype("category")
)

human_immune.obs["NK/cytoxic T lymphocytes"] = (
    human_immune.obs["leiden_R"].isin(['2', '3', '13', '14', '18', '21', '22', '37']).astype("category")
)

human_immune.obs["mast cells"] = (
    human_immune.obs["leiden_R"].isin(['15']).astype("category")
)

human_immune.obs["monocytes/macrophages"] = (
    human_immune.obs["leiden_R"].isin(['7,2', '7,3', '7,4', '17', '35', '32']).astype("category")
)



human_immune.obs["immune"] = "Other"
human_immune.obs.loc[human_immune.obs['dendritic cells'] == True, "immune"] = "dendritic cells"
human_immune.obs.loc[human_immune.obs['B cells'] == True, "immune"] = "B cells"
human_immune.obs.loc[human_immune.obs['Treg'] == True, "immune"] = "Treg"
human_immune.obs.loc[human_immune.obs['NK/cytoxic T lymphocytes'] == True, "immune"] = "NK/cytoxic T lymphocytes"
human_immune.obs.loc[human_immune.obs['CD4+ T lymphocytes'] == True, "immune"] = "CD4+ T lymphocytes"
human_immune.obs.loc[human_immune.obs['mast cells'] == True, "immune"] = "mast cells"
human_immune.obs.loc[human_immune.obs['monocytes/macrophages'] == True, "immune"] = "monocytes/macrophages"

###
human_immune.obs['immune'].value_counts()
human_immune.obs['leiden_R'].value_counts()

# umap of annotated immune cell types
sc.pl.umap(human_immune, color = 'immune', save = "_human_immune")


################################################
# differential expression between the cell types
sc.tl.rank_genes_groups(human_immune, 'immune', pts=True, use_raw = True)

markers_immune = sc.get.rank_genes_groups_df(human_immune, group=None)


# dotplot
dp = sc.pl.DotPlot(human_immune,
                    var_names = ['CD4', 'CD69', 'CD40LG', 'LEF1',
                                 'FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'ICOS',
                                  'GZMK', 'KLRG1', 'NKG7', 'XCL1',
                                  'CD79A', 'CD19', 'MS4A1', 'CD74',
                                  'CCR7', 'CCL22',
                                 'CCR9', 'CYBB',
                                 'CCL5', 'CCL7', 'IL1B', 'IFIT1', 'APOE', 'C1QA', 'C1QB', 'C1QC',
                                 ],
                    categories_order = ['CD4+ T lymphocytes', 'Treg', 'NK/cytoxic T lymphocytes', 'B cells','dendritic cells','monocytes/macrophages', 'mast cells'],
                    groupby='immune', cmap = 'Reds', mean_only_expressed = True
)

dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_compartments/DotPlot_immune_human.png')


sc.pl.umap(
    human_immune,
    color=["leiden_R", "immune"],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_immuneCells'
)



sc.pl.violin(human_immune, ['C1QA', 'C1QB', 'C1QC'], groupby = 'immune', use_raw=True, rotation=90, save='_macrophages_human')
sc.pl.violin(human_immune, ['CD79A', 'CD19', 'MS4A1', 'CD74'], groupby = 'immune', use_raw=True, rotation=90, save='_Bcells_human')
sc.pl.violin(human_immune, ['GZMK', 'KLRG1', 'NKG7', 'XCL1'], groupby = 'immune', use_raw=True, rotation=90, save='_Tcells_human')
sc.pl.violin(human_immune, ['FOXP3', 'CTLA4', 'IL2RA', 'ICOS'], groupby = 'immune', use_raw=True, rotation=90, save='_Treg_human')
sc.pl.violin(human_immune, ['CD4', 'CD69', 'CD40LG', 'LEF1'], groupby = 'immune', use_raw=True, rotation=90, save='_Cd4_human')

##################################################
# stacked violin plot for immune cell types per ERG status

human_immune.obs['cluster'] = human_immune.obs['immune']
human_immune.obs['cluster'] = human_immune.obs['cluster'].astype('category')

sc.pl.umap(human_immune, color = ['cluster'])


relativeFrequency_immune = sct.tools.relative_frequency_per_cluster(human_immune, group_by='erg', xlabel='cluster')
sct.plot.cluster_composition_stacked_barplot(relativeFrequency_immune, xlabel='erg', figsize=(50, 40), width=0.7, margins=(0.02, 0.02), label_size=0, tick_size=50, colors=human_immune.uns['cluster_colors'], save = 'figures/human_compartments/immune_frequency_Per_ERG.png')

####################################################
# save immune
####################################################
human_immune.write('outs_human/h5ads/human_immune_annot.h5ad')

# for cell chat
adata_human_immune = human_immune.raw.to_adata()
adata_human_immune.X = sp.csr_matrix.todense(adata_human_immune.X)
adata_human_immune.X = adata_human_immune.to_df()
adata_human_immune.write('forCellChat/adata_human_immune_raw.h5ad')

########################################################
# read immune
########################################################
human_immune = sc.read_h5ad('outs_human/h5ads/human_immune_annot.h5ad', chunk_size=100000)

########################################################
########################################################


##################################################################################################################
# annotate human epithelium componenet
##################################################################################################################

# get the epithelium
human_epithelium = adata_human_all[adata_human_all.obs.compartments == 'epithelium']

sc.pl.umap(
    human_epithelium,
    color=["leiden", 'basal'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_basal'
)

sc.pl.umap(
    human_epithelium,
    color=["leiden", 'luminal'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_luminal'
)

sc.pl.umap(
    human_epithelium,
    color=["leiden", 'neuroendocrine'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_neuroendocrine'
)

sc.pl.umap(
    human_epithelium,
    color=["leiden", 'seminal_vesicle_basal'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_seminal_vesicle_basal'
)

sc.pl.umap(
    human_epithelium,
    color=["leiden", 'seminal_vesicle_luminal'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_seminal_vesicle_luminal'
)

sc.pl.umap(
    human_epithelium,
    color=["leiden", 'seminal_vesicle_ionocyte'],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_seminal_vesicle_ionocyte'
)

human_epithelium.obs["luminal"] = (
    human_epithelium.obs["leiden"].isin(['5', '8', '27', '20', '23', '31', '39', '36', '48', '66', '40', '28', '61', '26', '55']).astype("category")
)

human_epithelium.obs["basal"] = (
    human_epithelium.obs["leiden"].isin(['12', '29']).astype("category")
)


human_epithelium.obs["epithelium"] = "Other"
human_epithelium.obs.loc[human_epithelium.obs['luminal'] == True, "epithelium"] = "luminal"
human_epithelium.obs.loc[human_epithelium.obs['basal'] == True, "epithelium"] = "basal"

human_epithelium.obs['epithelium'].value_counts()
human_epithelium.obs['leiden'].value_counts()

# umap of annotated immune cell types
sc.pl.umap(human_epithelium, color = 'epithelium', save = "_human_epithelium")


################################################
# differential expression between the cell types
sc.tl.rank_genes_groups(human_epithelium, 'epithelium', pts=True, use_raw = True)

markers_immune = sc.get.rank_genes_groups_df(human_epithelium, group=None)


# dotplot
dp = sc.pl.DotPlot(human_epithelium,
                    var_names = [
                        "KLK3", "AR", "ACPP", "NKX3-1", "KLK2",
                        "KRT5", "KRT14", "TP63", "CD44", "KRT15"
                        ],
                    categories_order = ['luminal', 'basal'],
                    groupby='epithelium', cmap = 'Reds', mean_only_expressed = True
)

dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/human_compartments/DotPlot_epithelium_human.png')


sc.pl.umap(
    human_epithelium,
    color=["leiden", "epithelium"],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_epithCells'
)

sc.pl.violin(human_epithelium, ["KLK3", "AR", "ACPP", "NKX3-1", "KLK2"], groupby = 'epithelium', use_raw=True, rotation=90, save='_luminal_human')
sc.pl.violin(human_epithelium, ["KRT5", "KRT14", "TP63", "CD44", "KRT15"], groupby = 'epithelium', use_raw=True, rotation=90, save='_basal_human')

##################################################
# stacked violin plot for epithelial cell types per ERG status

human_epithelium.obs['cluster'] = human_epithelium.obs['epithelium']
human_epithelium.obs['cluster'] = human_epithelium.obs['cluster'].astype('category')

sc.pl.umap(human_epithelium, color = ['cluster'])


relativeFrequency_epithelium = sct.tools.relative_frequency_per_cluster(human_epithelium, group_by='erg', xlabel='cluster')
sct.plot.cluster_composition_stacked_barplot(relativeFrequency_epithelium, xlabel='erg', figsize=(50, 40), width=0.7, margins=(0.02, 0.02), label_size=0, tick_size=50, colors=human_epithelium.uns['cluster_colors'], save = 'figures/human_compartments/epithelium_frequency_Per_ERG.png')

####################################################
# save epithelium
####################################################
human_epithelium.write('outs_human/h5ads/human_epithelium_annot.h5ad')

# for cell chat
adata_human_epithelium = human_epithelium.raw.to_adata()
adata_human_epithelium.X = sp.csr_matrix.todense(adata_human_epithelium.X)
adata_human_epithelium.X = adata_human_epithelium.to_df()
adata_human_epithelium.write('forCellChat/adata_human_epithelium_raw.h5ad')

########################################################
# read epithelium
########################################################
human_epithelium = sc.read_h5ad('outs_human/h5ads/human_epithelium_annot.h5ad', chunk_size=100000)



##################################################################################
## integrate immune, stroma, and epithelium
##################################################################################

# add a key called cluster to human_immune with the cell type info
human_immune.obs['cluster'] = human_immune.obs['immune']

# add a key called cluster to human_epithelium with the cell type info
#human_epithelium.obs['cluster'] = 'epithelium'


#########
# load the human stroma
human_stroma = sc.read_h5ad('data/for_human/adata_human.h5ad', chunk_size=100000)
human_stroma.obs['cluster'] = human_stroma.obs['cluster'].astype('str')
human_stroma.obs['cluster'].value_counts()

human_stroma.obs['erg'].replace('ERG+', 'positive', inplace=True)
human_stroma.obs['erg'].replace('ERG-', 'negative', inplace=True)
human_stroma.obs['erg'].value_counts()
############################################
## combine

# Get cell identifiers for each AnnData object
human_epithelium_cells = human_epithelium.obs.index
human_stroma_cells = human_stroma.obs.index
human_immune_cells = human_immune.obs.index

# Find overlapping cells
overlapping_cells = human_epithelium_cells.intersection(human_stroma_cells)

# Print overlapping cells
print(overlapping_cells)
# Assume overlapping_cells contains the identifiers of the overlapping cells

# Get cluster assignments for the overlapping cells in the stromal object
epith_clusters = human_epithelium.obs.loc[overlapping_cells, 'cluster']

# Exclude overlapping cells from the epithelium object
human_epithelium = human_epithelium[~human_epithelium.obs.index.isin(overlapping_cells)]

# Now 'human_epithelium_filtered' is your new AnnData object for epithelium without the overlapping cells
human_epithelium


adata_list = [human_epithelium, human_immune, human_stroma]
human_all_annot = ad.concat(adata_list, join="outer", label = 'compartment', uns_merge = 'same')
human_all_annot.obs['compartment'].value_counts()
human_all_annot.obs['compartment'].replace('0', 'epthelium', inplace=True)
human_all_annot.obs['compartment'].replace('1', 'immune', inplace=True)
human_all_annot.obs['compartment'].replace('2', 'stroma', inplace=True)

human_all_annot.obs['cluster'].value_counts()


human_all_annot.obs['erg'].value_counts()

#############################################
# save for future use
# to remove
toRem =  ['seminal_vesicle_basal', 'seminal_vesicle_luminal', 'seminal_vesicle_ionocyte', 'epithelium', 'stroma', 'immune', 'endothelium', 'compartments', 'leiden_R', 'dendritic cells', 'B cells', 'Treg', 'CD4+ T lymphocytes', 'NK/cytoxic T lymphocytes', 'monocytes/macrophages', 'Regulon(Arid5a)', 'Regulon(Arid5b)', 'Regulon(Ascl1)', 'Regulon(Ascl2)', 'Regulon(Atf3)', 'Regulon(Bach1)', 'Regulon(Batf)', 'Regulon(Bcl3)', 'Regulon(Cebpa)', 'Regulon(Cebpb)', 'Regulon(Cebpd)', 'Regulon(Creb5)', 'Regulon(Crem)', 'Regulon(Dusp26)', 'Regulon(Egr1)', 'Regulon(Egr2)', 'Regulon(Egr3)', 'Regulon(Egr4)', 'Regulon(Eomes)', 'Regulon(Erg)', 'Regulon(Ets1)', 'Regulon(Fezf1)', 'Regulon(Fosb)', 'Regulon(Fosl2)', 'Regulon(Foxa1)', 'Regulon(Foxd3)', 'Regulon(Foxi1)', 'Regulon(Foxo1)', 'Regulon(Foxq1)', 'Regulon(Foxs1)', 'Regulon(Gabpb1)', 'Regulon(Gata2)', 'Regulon(Gata3)', 'Regulon(Gata6)', 'Regulon(Grhl3)', 'Regulon(Hnf4a)', 'Regulon(Hoxb6)', 'Regulon(Ikzf2)', 'Regulon(Irf1)', 'Regulon(Irf4)', 'Regulon(Irf5)', 'Regulon(Irf6)', 'Regulon(Irf7)', 'Regulon(Irf8)', 'Regulon(Junb)', 'Regulon(Jund)', 'Regulon(Klf2)', 'Regulon(Klf4)', 'Regulon(Klf5)', 'Regulon(Lhx6)', 'Regulon(Mafb)', 'Regulon(Maff)', 'Regulon(Mef2c)', 'Regulon(Myc)', 'Regulon(Myod1)', 'Regulon(Nfe2l2)', 'Regulon(Nfia)', 'Regulon(Nfil3)', 'Regulon(Nfix)', 'Regulon(Nfkb1)', 'Regulon(Nkx6-2)', 'Regulon(Onecut2)', 'Regulon(Pax3)', 'Regulon(Peg3)', 'Regulon(Pgr)', 'Regulon(Pou2f3)', 'Regulon(Pparg)', 'Regulon(Prrx2)', 'Regulon(Rel)', 'Regulon(Runx1)', 'Regulon(Runx3)', 'Regulon(Six2)', 'Regulon(Snai3)', 'Regulon(Sox10)', 'Regulon(Sox11)', 'Regulon(Sox18)', 'Regulon(Sox2)', 'Regulon(Sox4)', 'Regulon(Sox7)', 'Regulon(Sox9)', 'Regulon(Spi1)', 'Regulon(Spib)', 'Regulon(Spic)', 'Regulon(Srebf1)', 'Regulon(Stat3)', 'Regulon(Tagln2)', 'Regulon(Tal1)', 'Regulon(Tbx1)', 'Regulon(Tbx21)', 'Regulon(Tcf4)', 'Regulon(Tead1)', 'Regulon(Tff3)', 'Regulon(Trp63)', 'Regulon(Twist1)', 'name', 'name_nosuperscript']

for i in toRem:
    del human_all_annot.obs[i]

human_all_annot.write('data/human_all_annot.h5ad')


#############################################
# save raw for cellchat
human_all_annot_raw = human_all_annot.raw.to_adata()
#%%
human_all_annot_raw.X = sp.csr_matrix.todense(human_all_annot_raw.X)
human_all_annot_raw.X = human_all_annot_raw.to_df()
#%%
# to remove
#toRem =  ['_scvi_batch', '_scvi_labels', '_scvi_local_l_mean', '_scvi_local_l_var', 'endothelial', 'fibroblast', 'myofibroblast', 'dendritic', 'cCDs', 'langherhans_like', 'b', 't_nk', 'myeloid', 'mast', 'luminal', 'basal', 'notluminal', 'macrophages', 'neuroendocrine', 'seminal_vesicle_basal', 'seminal_vesicle_luminal', 'seminal_vesicle_ionocyte', 'epithelium', 'stroma', 'immune', 'endothelium', 'compartments', 'leiden_R', 'dendritic cells', 'B cells', 'Treg', 'CD4+ T lymphocytes', 'NK/cytoxic T lymphocytes', 'monocytes/macrophages', 'Regulon(Arid5a)', 'Regulon(Arid5b)', 'Regulon(Ascl1)', 'Regulon(Ascl2)', 'Regulon(Atf3)', 'Regulon(Bach1)', 'Regulon(Batf)', 'Regulon(Bcl3)', 'Regulon(Cebpa)', 'Regulon(Cebpb)', 'Regulon(Cebpd)', 'Regulon(Creb5)', 'Regulon(Crem)', 'Regulon(Dusp26)', 'Regulon(Egr1)', 'Regulon(Egr2)', 'Regulon(Egr3)', 'Regulon(Egr4)', 'Regulon(Eomes)', 'Regulon(Erg)', 'Regulon(Ets1)', 'Regulon(Fezf1)', 'Regulon(Fosb)', 'Regulon(Fosl2)', 'Regulon(Foxa1)', 'Regulon(Foxd3)', 'Regulon(Foxi1)', 'Regulon(Foxo1)', 'Regulon(Foxq1)', 'Regulon(Foxs1)', 'Regulon(Gabpb1)', 'Regulon(Gata2)', 'Regulon(Gata3)', 'Regulon(Gata6)', 'Regulon(Grhl3)', 'Regulon(Hnf4a)', 'Regulon(Hoxb6)', 'Regulon(Ikzf2)', 'Regulon(Irf1)', 'Regulon(Irf4)', 'Regulon(Irf5)', 'Regulon(Irf6)', 'Regulon(Irf7)', 'Regulon(Irf8)', 'Regulon(Junb)', 'Regulon(Jund)', 'Regulon(Klf2)', 'Regulon(Klf4)', 'Regulon(Klf5)', 'Regulon(Lhx6)', 'Regulon(Mafb)', 'Regulon(Maff)', 'Regulon(Mef2c)', 'Regulon(Myc)', 'Regulon(Myod1)', 'Regulon(Nfe2l2)', 'Regulon(Nfia)', 'Regulon(Nfil3)', 'Regulon(Nfix)', 'Regulon(Nfkb1)', 'Regulon(Nkx6-2)', 'Regulon(Onecut2)', 'Regulon(Pax3)', 'Regulon(Peg3)', 'Regulon(Pgr)', 'Regulon(Pou2f3)', 'Regulon(Pparg)', 'Regulon(Prrx2)', 'Regulon(Rel)', 'Regulon(Runx1)', 'Regulon(Runx3)', 'Regulon(Six2)', 'Regulon(Snai3)', 'Regulon(Sox10)', 'Regulon(Sox11)', 'Regulon(Sox18)', 'Regulon(Sox2)', 'Regulon(Sox4)', 'Regulon(Sox7)', 'Regulon(Sox9)', 'Regulon(Spi1)', 'Regulon(Spib)', 'Regulon(Spic)', 'Regulon(Srebf1)', 'Regulon(Stat3)', 'Regulon(Tagln2)', 'Regulon(Tal1)', 'Regulon(Tbx1)', 'Regulon(Tbx21)', 'Regulon(Tcf4)', 'Regulon(Tead1)', 'Regulon(Tff3)', 'Regulon(Trp63)', 'Regulon(Twist1)', 'name', 'name_nosuperscript']

#for i in toRem:
#    del human_all_annot_raw.obs[i]


human_all_annot_raw.obs_names_make_unique()
human_all_annot_raw.write('forCellChat/human_all_annot_raw.h5ad')

