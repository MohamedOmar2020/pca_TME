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
import sc_toolbox.api as sct
import seaborn as sb
import plotly.express as px
import squidpy as sq

############################################
# load the mouse data
adata_mouse_all = sc.read_h5ad('outs/h5ads/fapcm_all_v6.h5ad', chunk_size=100000)


adata_mouse_all.obs['leiden'] = adata_mouse_all.obs['leiden'].astype('str')
adata_mouse_all.obs['leiden'].value_counts()

adata_mouse_all.obs['key'] = adata_mouse_all.obs['key'].astype('str')
adata_mouse_all.obs['key'].value_counts()

###############################################
## annotate compartments

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

##########################################################
## annotate the different immune cells
mouse_immune = adata_mouse_all[adata_mouse_all.obs.compartments == 'immune']

# sc.pl.umap(
#     mouse_immune,
#     color=["leiden", "Cxcr3", "Clec9a", "Xcr1", "Ccr7", "Ccl17", "Ccl19", "Cd207", "Clec9a", "Ccl17", 'dendritic', 'cCDs'],
#     frameon=False,
#     legend_loc="on data",
#     size=5,
#     legend_fontsize=6,
#     save = '_dendritic_cells.png'
# )
#
# sc.pl.umap(
#     mouse_immune,
#     color=["leiden", "Cd79a", "Cd79b", 'b'],
#     frameon=False,
#     legend_loc="on data",
#     size=5,
#     legend_fontsize=6,
#     save = '_B_cells.png'
# )
#
# sc.pl.umap(
#     mouse_immune,
#     color=["leiden", "Foxp3", "Ctla4", 'Tnfrsf4', 'Irf4', 'Batf', 'Tnfrsf18', 'Tox2', 'Prdm1'],
#     frameon=False,
#     legend_loc="on data",
#     size=5,
#     legend_fontsize=6,
#     save = '_Treg.png'
# )
#
# sc.pl.umap(
#     mouse_immune,
#     color=["leiden", "Lef1", "Atm", 'Sell', 'Klf2', 'Itga6', 'Ifngr2', 'Il21r', 'Foxp1'],
#     frameon=False,
#     legend_loc="on data",
#     size=5,
#     legend_fontsize=6,
#     save = '_naiveCD4Tcells.png'
# )
#
# sc.pl.umap(
#     mouse_immune,
#     color=["leiden", "Ccl3", "Ifng", 'Ccl4', 'Xcl1', 'Csf2', 'Il10', 'Hopx', 'Lag3', 'Prf1', 'Tnfrsf9', 'Nkg7'],
#     frameon=False,
#     legend_loc="on data",
#     size=5,
#     legend_fontsize=6,
#     save = '_cytotoxicCD8Tcells.png'
# )
#
# sc.pl.umap(
#     mouse_immune,
#     color=["leiden", 't_nk', 'Ncr1', 'Prf1'],
#     frameon=False,
#     legend_loc="on data",
#     size=5,
#     legend_fontsize=6, use_raw = False,
#     save = '_NKcells.png'
# )
#
# mouse_immune.obs.leiden.value_counts()

# markers_mouse_c0 = sc.get.rank_genes_groups_df(mouse_immune, group = '0')
# markers_mouse_c3 = sc.get.rank_genes_groups_df(mouse_immune, group = '3')
# markers_mouse_c4 = sc.get.rank_genes_groups_df(mouse_immune, group = '4')
# markers_mouse_c5 = sc.get.rank_genes_groups_df(mouse_immune, group = '5')
# markers_mouse_c9 = sc.get.rank_genes_groups_df(mouse_immune, group = '9')# remove
# markers_mouse_c10 = sc.get.rank_genes_groups_df(mouse_immune, group = '10')
# markers_mouse_c12 = sc.get.rank_genes_groups_df(mouse_immune, group = '12')
# markers_mouse_c14 = sc.get.rank_genes_groups_df(mouse_immune, group = '14')
# markers_mouse_c17 = sc.get.rank_genes_groups_df(mouse_immune, group = '17')
# markers_mouse_c20 = sc.get.rank_genes_groups_df(mouse_immune, group = '20') # unknown
# markers_mouse_c28 = sc.get.rank_genes_groups_df(mouse_immune, group = '28')
# markers_mouse_c30 = sc.get.rank_genes_groups_df(mouse_immune, group = '30')
# markers_mouse_c34 = sc.get.rank_genes_groups_df(mouse_immune, group = '34')
# markers_mouse_c35 = sc.get.rank_genes_groups_df(mouse_immune, group = '35')
# markers_mouse_c37 = sc.get.rank_genes_groups_df(mouse_immune, group = '37')
# markers_mouse_c39 = sc.get.rank_genes_groups_df(mouse_immune, group = '39') # remove
# markers_mouse_c41 = sc.get.rank_genes_groups_df(mouse_immune, group = '41')

sc.tl.rank_genes_groups(mouse_immune, 'leiden', pts=True, use_raw = True)

sc.tl.leiden(mouse_immune, restrict_to = ('leiden', ['0']), resolution=0.9)
mouse_immune.obs.leiden.value_counts()
mouse_immune.obs.leiden_R.value_counts()
sc.tl.rank_genes_groups(mouse_immune, 'leiden_R', pts=True, use_raw = True)

# remove 9, 20, and 39, 41
mouse_immune = mouse_immune[~mouse_immune.obs["leiden_R"].isin(["9", "20", "2", '35', "39", "41"])]

markers_mouse_c0_0 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,0')
markers_mouse_c0_1 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,1')
markers_mouse_c0_2 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,2')
markers_mouse_c0_3 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,3')
markers_mouse_c0_4 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,4')
markers_mouse_c0_5 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,5')
markers_mouse_c0_6 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,6')
markers_mouse_c0_7 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,7')
markers_mouse_c0_8 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,8')
markers_mouse_c0_9 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,9')
markers_mouse_c0_10 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,10')
markers_mouse_c0_11 = sc.get.rank_genes_groups_df(mouse_immune, group = '0,11')
markers_mouse_c3 = sc.get.rank_genes_groups_df(mouse_immune, group = '3')
markers_mouse_c4 = sc.get.rank_genes_groups_df(mouse_immune, group = '4')
markers_mouse_c5 = sc.get.rank_genes_groups_df(mouse_immune, group = '5')
markers_mouse_c10 = sc.get.rank_genes_groups_df(mouse_immune, group = '10')
markers_mouse_c12 = sc.get.rank_genes_groups_df(mouse_immune, group = '12')
markers_mouse_c14 = sc.get.rank_genes_groups_df(mouse_immune, group = '14')
markers_mouse_c17 = sc.get.rank_genes_groups_df(mouse_immune, group = '17')
markers_mouse_c28 = sc.get.rank_genes_groups_df(mouse_immune, group = '28')
markers_mouse_c30 = sc.get.rank_genes_groups_df(mouse_immune, group = '30')
markers_mouse_c34 = sc.get.rank_genes_groups_df(mouse_immune, group = '34')
markers_mouse_c37 = sc.get.rank_genes_groups_df(mouse_immune, group = '37')




mouse_immune.obs["dendritic cells"] = (
    mouse_immune.obs["leiden_R"].isin(['17']).astype("category")
)

mouse_immune.obs["B cells"] = (
    mouse_immune.obs["leiden_R"].isin(['28']).astype("category")
)

mouse_immune.obs["Treg"] = (
    mouse_immune.obs["leiden_R"].isin(['14']).astype("category")
)

mouse_immune.obs["CD4+ T lymphocytes"] = (
    mouse_immune.obs["leiden_R"].isin(['10']).astype("category")
)

mouse_immune.obs["NK/cytoxic T lymphocytes"] = (
    mouse_immune.obs["leiden_R"].isin(['3', '4', '12', '0,11']).astype("category")
)

#mouse_immune.obs["cycling NK/cytoxic T lymphocytes"] = (
#    mouse_immune.obs["leiden"].isin(['35']).astype("category")
#)

#mouse_immune.obs["monocytes_Ccl5"] = (
#    mouse_immune.obs["leiden"].isin(['0']).astype("category")
#)
#mouse_immune.obs["monocytes_Lyz2"] = (
#    mouse_immune.obs["leiden"].isin(['5', '34']).astype("category")
#)


mouse_immune.obs["monocytes/macrophages"] = (
    mouse_immune.obs["leiden_R"].isin(['0,0', '0,1', '0,2', '0,3', '0,4', '0,5', '0,6', '0,7', '0,8', '0,9', '0,10', '5', '30', '34', '37']).astype("category")
)



mouse_immune.obs["immune"] = "Other"
mouse_immune.obs.loc[mouse_immune.obs['dendritic cells'] == True, "immune"] = "dendritic cells"
mouse_immune.obs.loc[mouse_immune.obs['B cells'] == True, "immune"] = "B cells"
mouse_immune.obs.loc[mouse_immune.obs['Treg'] == True, "immune"] = "Treg"
mouse_immune.obs.loc[mouse_immune.obs['NK/cytoxic T lymphocytes'] == True, "immune"] = "NK/cytoxic T lymphocytes"
mouse_immune.obs.loc[mouse_immune.obs['CD4+ T lymphocytes'] == True, "immune"] = "CD4+ T lymphocytes"
#mouse_immune.obs.loc[mouse_immune.obs['cycling NK/cytoxic T lymphocytes'] == True, "immune"] = "cycling NK/cytoxic T lymphocytes"
#mouse_immune.obs.loc[mouse_immune.obs['monocytes_Ccl5'] == True, "immune"] = "Ccl5+ monocytes"
#mouse_immune.obs.loc[mouse_immune.obs['monocytes_Lyz2'] == True, "immune"] = "Lyz2+ monocytes"
#mouse_immune.obs.loc[mouse_immune.obs['TAMs'] == True, "immune"] = "TAMs"
mouse_immune.obs.loc[mouse_immune.obs['monocytes/macrophages'] == True, "immune"] = "monocytes/macrophages"
#mouse_immune.obs.loc[mouse_immune.obs['pDCs'] == True, "immune"] = "pDCs"

###
mouse_immune.obs['immune'].value_counts()

################################################
# differential expression between the cell types
sc.tl.rank_genes_groups(mouse_immune, 'immune', pts=True, use_raw = True)

markers_immune = sc.get.rank_genes_groups_df(mouse_immune, group=None)


# dotplot
dp = sc.pl.DotPlot(mouse_immune,
                    var_names = ['Cd4', 'Cd69', 'Cd40lg', 'Lef1',
                                 'Foxp3', 'Il2ra', 'Ctla4', 'Ikzf2', 'Icos',
                                  'Gzmk', 'Klrg1', 'Nkg7', 'Xcl1',
                                  'Cd79a', 'Cd19', 'Ms4a1', 'Cd74',
                                  'Ccr7', 'Ccl22',
                                  'Siglech', 'Ccr9', 'Cybb',
                                 'Ccl5', 'Ccl7', 'Il1b', 'Ifit1','Lyz2', 'Apoe', 'C1qa',
                                 ],
                    categories_order = ['CD4+ T lymphocytes', 'Treg', 'NK/cytoxic T lymphocytes', 'B cells','dendritic cells','monocytes/macrophages'],
                    groupby='immune', cmap = 'Reds', mean_only_expressed = True
)

dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5).savefig('figures/DotPlot_immune_mouse.png')


sc.pl.umap(
    mouse_immune,
    color=["leiden", "immune"],
    frameon=False,
    legend_loc="on data",
    size=5,
    legend_fontsize=6,
    save = '_immune.png'
)



sc.pl.violin(mouse_immune, ['C1qa', 'C1qb', 'C1qc'], groupby = 'immune', use_raw=True, rotation=90, save='_macrophages_mouse.png')
sc.pl.violin(mouse_immune, ['Cd79a', 'Cd19', 'Ms4a1', 'Cd74'], groupby = 'immune', use_raw=True, rotation=90, save='_Bcells_mouse.png')
sc.pl.violin(mouse_immune, ['Gzmk', 'Klrg1', 'Nkg7', 'Xcl1'], groupby = 'immune', use_raw=True, rotation=90, save='_Tcells_mouse.png')
sc.pl.violin(mouse_immune, ['Foxp3', 'Ctla4', 'Il2ra', 'Icos'], groupby = 'immune', use_raw=True, rotation=90, save='_Treg_mouse.png')
sc.pl.violin(mouse_immune, ['Cd4', 'Cd69', 'Cd40lg', 'Lef1'], groupby = 'immune', use_raw=True, rotation=90, save='_Cd4_mouse.png')

sc.pl.matrixplot(
    mouse_immune,
    var_names=['C1qa', 'C1qb', 'C1qc', 'Cd79a', 'Cd19', 'Ms4a1', 'Cd74', 'Gzmk', 'Klrg1', 'Nkg7', 'Xcl1', 'Foxp3', 'Ctla4', 'Il2ra', 'Icos', 'Cd4', 'Cd69', 'Cd40lg', 'Lef1'],
    standard_scale="var",
    groupby="immune",
    cmap="RdBu_r",
    save='mouse.png'
)



##################################
## parallel categories
color_map = dict(
    zip(mouse_immune.obs['key'].cat.categories.tolist(), sc.pl.palettes.vega_20_scanpy)
)

# by GEMM
mouse_immune.obs["color"] = mouse_immune.obs['key'].map(color_map)
fig = px.parallel_categories(
    mouse_immune.obs,
    dimensions=['key', 'immune'],
    color="color",
    labels={'key': "model", 'cell type': "immune"},
)
fig.update_layout(autosize=True, width=500, height=1000, font_size = 8)
fig.write_image("figures/parallel_categories_mouse_immune_byMM.png")


################
# by cell type
mouse_immune.obs['immune'] = mouse_immune.obs['immune'].astype('category')
color_map = dict(
    zip(mouse_immune.obs['immune'].cat.categories.tolist(), sc.pl.palettes.vega_20_scanpy)
)
mouse_immune.obs["color"] = mouse_immune.obs['immune'].map(color_map)
fig = px.parallel_categories(
    mouse_immune.obs,
    dimensions=['immune', 'key'],
    color="color",
    labels={'key': "model", 'cell type': "immune"},
)
fig.update_layout(autosize=False, width=500, height=1000, font_size = 8)
fig.write_image("figures/parallel_categories_mouse_immune_byCellType.png")


########################################
# save
mouse_immune.write('outs/h5ads/mouse_immune.h5ad')
mouse_immune = sc.read_h5ad('outs/h5ads/mouse_immune.h5ad', chunk_size=100000)

# cell chat
adata_mouse_immune = mouse_immune.raw.to_adata()
adata_mouse_immune.X = sp.csr_matrix.todense(adata_mouse_immune.X)
adata_mouse_immune.X = adata_mouse_immune.to_df()
adata_mouse.write('forCellChat/adata_mouse_immune_raw.h5ad')


##################################################################################
## integrate with stroma

# add a key called cluster to mouse_immune with the cell type info
mouse_immune.obs['cluster'] = mouse_immune.obs['immune']

mouse_stroma = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)
mouse_stroma.obs['cluster'] = mouse_stroma.obs['cluster'].astype('str')
mouse_stroma.obs['cluster'].value_counts()

# filter the na
mouse_stroma = mouse_stroma[mouse_stroma.obs['cluster'] != 'nan', :]
mouse_stroma.obs['cluster'].value_counts()



adata_list = [mouse_immune, mouse_stroma]
mouse_immune_stroma = ad.concat(adata_list, join="outer", label = 'compartment', uns_merge = 'same')
mouse_immune_stroma.obs['compartment'].value_counts()
mouse_immune_stroma.obs['compartment'].replace('0', 'immune', inplace=True)
mouse_immune_stroma.obs['compartment'].replace('1', 'stroma', inplace=True)

mouse_immune_stroma.obs['cluster'].value_counts()


#############################################
# save for cellchat
mouse_immune_stroma_raw = mouse_immune_stroma.raw.to_adata()
#%%
mouse_immune_stroma_raw.X = sp.csr_matrix.todense(mouse_immune_stroma_raw.X)
mouse_immune_stroma_raw.X = mouse_immune_stroma_raw.to_df()
#%%
# to remove
toRem =  ['_scvi_batch', '_scvi_labels', '_scvi_local_l_mean', '_scvi_local_l_var', 'endothelial', 'fibroblast', 'myofibroblast', 'dendritic', 'cCDs', 'langherhans_like', 'b', 't_nk', 'myeloid', 'mast', 'luminal', 'basal', 'notluminal', 'macrophages', 'neuroendocrine', 'seminal_vesicle_basal', 'seminal_vesicle_luminal', 'seminal_vesicle_ionocyte', 'epithelium', 'stroma', 'immune', 'endothelium', 'compartments', 'leiden_R', 'dendritic cells', 'B cells', 'Treg', 'CD4+ T lymphocytes', 'NK/cytoxic T lymphocytes', 'monocytes/macrophages', 'Regulon(Arid5a)', 'Regulon(Arid5b)', 'Regulon(Ascl1)', 'Regulon(Ascl2)', 'Regulon(Atf3)', 'Regulon(Bach1)', 'Regulon(Batf)', 'Regulon(Bcl3)', 'Regulon(Cebpa)', 'Regulon(Cebpb)', 'Regulon(Cebpd)', 'Regulon(Creb5)', 'Regulon(Crem)', 'Regulon(Dusp26)', 'Regulon(Egr1)', 'Regulon(Egr2)', 'Regulon(Egr3)', 'Regulon(Egr4)', 'Regulon(Eomes)', 'Regulon(Erg)', 'Regulon(Ets1)', 'Regulon(Fezf1)', 'Regulon(Fosb)', 'Regulon(Fosl2)', 'Regulon(Foxa1)', 'Regulon(Foxd3)', 'Regulon(Foxi1)', 'Regulon(Foxo1)', 'Regulon(Foxq1)', 'Regulon(Foxs1)', 'Regulon(Gabpb1)', 'Regulon(Gata2)', 'Regulon(Gata3)', 'Regulon(Gata6)', 'Regulon(Grhl3)', 'Regulon(Hnf4a)', 'Regulon(Hoxb6)', 'Regulon(Ikzf2)', 'Regulon(Irf1)', 'Regulon(Irf4)', 'Regulon(Irf5)', 'Regulon(Irf6)', 'Regulon(Irf7)', 'Regulon(Irf8)', 'Regulon(Junb)', 'Regulon(Jund)', 'Regulon(Klf2)', 'Regulon(Klf4)', 'Regulon(Klf5)', 'Regulon(Lhx6)', 'Regulon(Mafb)', 'Regulon(Maff)', 'Regulon(Mef2c)', 'Regulon(Myc)', 'Regulon(Myod1)', 'Regulon(Nfe2l2)', 'Regulon(Nfia)', 'Regulon(Nfil3)', 'Regulon(Nfix)', 'Regulon(Nfkb1)', 'Regulon(Nkx6-2)', 'Regulon(Onecut2)', 'Regulon(Pax3)', 'Regulon(Peg3)', 'Regulon(Pgr)', 'Regulon(Pou2f3)', 'Regulon(Pparg)', 'Regulon(Prrx2)', 'Regulon(Rel)', 'Regulon(Runx1)', 'Regulon(Runx3)', 'Regulon(Six2)', 'Regulon(Snai3)', 'Regulon(Sox10)', 'Regulon(Sox11)', 'Regulon(Sox18)', 'Regulon(Sox2)', 'Regulon(Sox4)', 'Regulon(Sox7)', 'Regulon(Sox9)', 'Regulon(Spi1)', 'Regulon(Spib)', 'Regulon(Spic)', 'Regulon(Srebf1)', 'Regulon(Stat3)', 'Regulon(Tagln2)', 'Regulon(Tal1)', 'Regulon(Tbx1)', 'Regulon(Tbx21)', 'Regulon(Tcf4)', 'Regulon(Tead1)', 'Regulon(Tff3)', 'Regulon(Trp63)', 'Regulon(Twist1)', 'name', 'name_nosuperscript']

for i in toRem:
    del mouse_immune_stroma_raw.obs[i]


mouse_immune_stroma_raw.obs_names_make_unique()
mouse_immune_stroma_raw.write('forCellChat/mouse_immune_stroma_raw.h5ad')

###################################################################
## LR
for key in mouse_immune_stroma.obs.key.cat.categories.tolist():
    lradata = mouse_immune_stroma[mouse_immune_stroma.obs["key"] == key]
    lradata.raw = lradata
    sq.gr.ligrec(
        lradata,
        n_perms=1000,
        cluster_key="compartment",
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        seed=123, show_progress_bar=True, corr_method='fdr_bh', alpha=0.05,
    )

    for i in lradata.uns['compartment_ligrec'].keys():
        lradata.uns['compartment_ligrec'][i].to_csv(f"objs/lr_immune_stroma/{key}_" + str(i) + ".csv")

    sq.pl.ligrec(
        lradata,
        cluster_key="compartment",
        source_groups=["stroma", "immune"],
        target_groups=["stroma", "immune"],
        remove_nonsig_interactions=True,
        means_range=(0.3, np.inf),
        alpha=0.05, pvalue_threshold=0.05,
        swap_axes=False,
        save=f"_{key}_lr_interactions_immune_stroma.pdf",
    )

###################################################################
## LR
mouse_immune_stroma.obs['cluster'] = mouse_immune_stroma.obs['cluster'].astype('category')
mouse_immune_stroma.obs['cluster'].value_counts()

sq.gr.ligrec(
        mouse_immune_stroma,
        n_perms=1000,
        cluster_key="cluster",
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        seed=123, show_progress_bar=True, corr_method='fdr_bh', alpha=0.05,
)

sq.pl.ligrec(
        mouse_immune_stroma,
        cluster_key="cluster",
        source_groups=['5', '6', '7'],
        target_groups=['Treg'],
        remove_nonsig_interactions=True,
        means_range=(0.3, np.inf),
        alpha=0.05, pvalue_threshold=0.05,
        swap_axes=False,
        title = 'c5/c6/c7 to Treg',
        save="lr_interactions_immune_stroma_c5c6c7_Treg.pdf",
)

########################################################
sc.pl.umap(
     mouse_immune_stroma,
     color=["key"],
     frameon=False,
     legend_loc="on data",
     size=5,
     legend_fontsize=6)

ax = sc.pl.correlation_matrix(mouse_immune_stroma, groupby='cluster', figsize=(5,3.5), save="_corr.png")

Mki67_relativeFrequency_all = sct.tools.relative_frequency_per_cluster(mouse_immune_stroma, group_by='cluster', xlabel='Mki67_status', condition=None)
Mki67_relativeFrequency_all['cluster'] = 'c'+Mki67_relativeFrequency_all['cluster']
sct.plot.cluster_composition_stacked_barplot(Mki67_relativeFrequency_all, xlabel='cluster', figsize=(8, 10), width=0.8, label_size=20, tick_size=16, margins=(0.02, 0.04), colors=adata_mouse.uns['Mki67_status_colors'], save = 'figures/Mki67_status.png')


relativeFrequency_all = sct.tools.relative_frequency_per_cluster(mouse_immune_stroma, group_by='cluster', xlabel='key')
sct.plot.cluster_composition_stacked_barplot(relativeFrequency_all, xlabel='cluster', figsize=(25, 25), width=0.7, margins=(0.02, 0.04), label_size=30, tick_size=20, colors=mouse_immune_stroma.uns['cluster_colors'], save = 'figures/frequency_Per_GEMM.png')



sct.plot.relative_frequencies_boxplots(relativeFrequency_all, cols = mouse_immune_stroma.uns['cluster_colors'], xlabel='cluster', figsize=(15, 6), width=0.5, order = ['0', '1'], cluster = ['0', '1'])



x = pd.crosstab(mouse_immune_stroma.obs['key'], mouse_immune_stroma.obs['cluster'])
x2 = relativeFrequency_all.transpose()

