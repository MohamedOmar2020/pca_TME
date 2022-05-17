# preprocessing pipeline for functional atlas of prostate cancer mesenchyme
# dataset overview:
# fvbn
# b6
# 1029b6
# terg
# himyc
# pten
# mycn

import sys

sys.path.insert(0, "/home/ryc4001/Documents/scutils/")

from pathlib import Path

import scanpy as sc
import scvi
from matplotlib.pyplot import ion
from scutils.figures.base import basics
from scutils.figures.prostate import annotate_cell_types_prostate
from scutils.qc import PreprocessRNA

# set global paths
basepath = "/home/ryc4001/Documents/scproj/fapcm"
outspath = basepath + "/outs"
Path(outspath).mkdir(parents=True, exist_ok=True)
h5ads = outspath + "/h5ads"
Path(h5ads).mkdir(parents=True, exist_ok=True)


def load():
    ion()
    fvbn1 = sc.read_10x_h5(basepath + "/data/6mon_FVBN/filtered_feature_bc_matrix.h5")
    fvbn1.var_names_make_unique()
    fvbn1.obs["batch"] = 0
    fvbn2 = sc.read_10x_h5(basepath + "/data/6mon_FVBN_2/filtered_feature_bc_matrix.h5")
    fvbn2.var_names_make_unique()
    fvbn2.obs["batch"] = 1
    fvbn = fvbn1.concatenate(fvbn2, join="outer",)
    fvbn.obs["key"] = "fvbn"
    fvbn.obs["model"] = "Fvbn"
    fvbn.obs["condition"] = "wildtype"
    b6 = sc.read_10x_h5(basepath + "/data/BL6_1116/filtered_feature_bc_matrix.h5")
    b6.var_names_make_unique()
    b6.obs["batch"] = 2
    b6.obs["key"] = "b6"
    b6.obs["model"] = "B6"
    b6.obs["condition"] = "wildtype"

    b6129 = sc.read_10x_h5(basepath + "/data/129b6/filtered_feature_bc_matrix.h5")
    b6129.var_names_make_unique()
    b6129.obs["batch"] = 3
    b6129.obs["key"] = "129b6"
    b6129.obs["model"] = "B6.129"
    b6129.obs["condition"] = "wildtype"

    terg1 = sc.read_10x_h5(basepath + "/data/TRG-1/filtered_feature_bc_matrix.h5")
    terg1.var_names_make_unique()
    terg2 = sc.read_10x_h5(basepath + "/data/TRG-2/filtered_feature_bc_matrix.h5")
    terg2.var_names_make_unique()
    terg = terg1.concatenate(terg2, join="outer",)
    terg.obs["batch"] = 4
    terg.obs["key"] = "terg"
    terg.obs["model"] = "T-ERG"
    terg.obs["condition"] = "mutant"

    himyc1 = sc.read_10x_h5(basepath + "/data/Myc_1116/filtered_feature_bc_matrix.h5")
    himyc1.var_names_make_unique()
    himyc2 = sc.read_10x_h5(basepath + "/data/Myc_1210/filtered_feature_bc_matrix.h5")
    himyc2.var_names_make_unique()
    himyc = himyc1.concatenate(himyc2, join="outer",)
    himyc.obs["batch"] = 5
    himyc.obs["key"] = "himyc"
    himyc.obs["model"] = "Hi-MYC"
    himyc.obs["condition"] = "mutant"

    pten1 = sc.read_10x_h5(basepath + "/data/P10_1/filtered_feature_bc_matrix.h5")
    pten1.var_names_make_unique()
    pten2 = sc.read_10x_h5(basepath + "/data/P10_2/filtered_feature_bc_matrix.h5")
    pten2.var_names_make_unique()
    pten = pten1.concatenate(pten2, join="outer",)
    pten.obs["batch"] = 6
    pten.obs["key"] = "pten"
    pten.obs["model"] = "Pten$^{KO}$"
    pten.obs["condition"] = "mutant"

    # TODO: finalize naming
    ptenwt = sc.read_10x_h5(basepath + "/data/WT_Pten/filtered_feature_bc_matrix.h5")
    ptenwt.var_names_make_unique()
    ptenwt.obs["batch"] = 7
    ptenwt.obs["key"] = "129b6_pten"
    ptenwt.obs["model"] = "B6.129"
    ptenwt.obs["condition"] = "wildtype"

    mycn1 = sc.read_10x_h5(basepath + "/data/2861/filtered_feature_bc_matrix.h5")
    mycn1.var_names_make_unique()
    mycn1.obs["batch"] = 8
    mycn2 = sc.read_10x_h5(basepath + "/data/2872/filtered_feature_bc_matrix.h5")
    mycn2.var_names_make_unique()
    mycn2.obs["batch"] = 9
    mycn3 = sc.read_10x_h5(basepath + "/data/3056/filtered_feature_bc_matrix.h5")
    mycn3.var_names_make_unique()
    mycn3.obs["batch"] = 10
    mycn4 = sc.read_10x_h5(basepath + "/data/3489/filtered_feature_bc_matrix.h5")
    mycn4.var_names_make_unique()
    mycn4.obs["batch"] = 11
    mycn5 = sc.read_10x_h5(basepath + "/data/3490/filtered_feature_bc_matrix.h5")
    mycn5.obs["batch"] = 12
    mycn5.var_names_make_unique()
    mycn = mycn1.concatenate(mycn2, mycn3, mycn4, mycn5, join="outer",)
    mycn.var_names_make_unique()
    mycn.obs["key"] = "mycn"
    mycn.obs["model"] = "MNRP DKO"
    mycn.obs["condition"] = "mutant"

    mycnwt = sc.read_10x_h5(basepath + "/data/WT_8weeks/filtered_feature_bc_matrix.h5")
    mycnwt.var_names_make_unique()
    mycnwt.obs["batch"] = 13
    mycnwt.obs["key"] = "129b6_mycn"
    mycnwt.obs["model"] = "B6.129"
    mycnwt.obs["condition"] = "wildtype"

    adata = fvbn.concatenate(
        b6, b6129, terg, himyc, pten, ptenwt, mycn, mycnwt, join="outer",
    )
    adata.layers["counts"] = adata.X.copy()
    return adata


def preprocess(adata):
    ion()
    sc.set_figure_params(dpi=400, dpi_save=400, figsize=(4, 4))
    adata.write_h5ad(h5ads + "/fapcm_unfiltered_v6.h5ad")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save="qc1.png",
    )
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save="qc2.png")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save="qc3.png")
    """
    min_counts=400
    max_counts=500000
    min_genes=300
    max_genes=8000
    max_percent_mito=0.15
    """
    pre = PreprocessRNA(subset_hvg=False)
    adata = pre.preprocess(adata)
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=4000,
        layer="counts",
        batch_key="batch",
        subset=True,
    )
    scvi.data.setup_anndata(adata, layer="counts", batch_key="batch")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata, resolution=1.0)
    sc.tl.umap(adata)
    basics(adata)
    adata.write_h5ad(h5ads + "/fapcm_all_v6.h5ad")
    return adata


def fibroblasts(adata):
    ion()
    sc.set_figure_params(dpi=400, dpi_save=400, figsize=(4, 4))
    # subset fibroblasts
    adata = adata[
        adata.obs["leiden"].isin(
            ["1", "8", "13", "30", "32", "3", "4", "5", "15", "26", "35",]
        )
    ]
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.umap(adata)
    sc.pl.umap(
        adata,
        color=[
            "Col1a1",
            "Bgn",
            "Dcn",
            "Pdgfra",
            "Pdgfrb",
            "Vim",
            "Tnc",
            "Alcam",
            "Col1a2",
            "Fgfr2",
            "Fgfr1",
            "Fgfr3",
            "Fgfr4",
            "Igf1r",
            "Igf2r",
            "Insr",
            "Insrr",
            "Met",
            "Mst1r",
            "Csf1r",
            "Flt3",
            "Kit",
            "Myh11",
            "Myl9",
            "Acta2",
        ],
        cmap="Reds",
        vmax=4,
        vmin=0,
        save="_fibroblast_myofibroblast_lineage_markers.png",
    )
    scvi.data.setup_anndata(adata, layer="counts", batch_key="batch")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.umap(adata)
    sc.tl.dendrogram(adata, groupby="leiden")
    basics(adata)
    adata.write_h5ad(h5ads + "/fapcm_fibroblasts_v6.h5ad")
    return adata


def checkmarkers(adata):
    ion()
    sc.set_figure_params(dpi=400, dpi_save=400, figsize=(4, 4))
    sc.pl.umap(
        adata,
        color=["Lgr5", "Col3a1", "Ddr2", "Pdgfra", "Pdgfrb", "Col1a1", "Trp63"],
        cmap="Reds",
        vmax=3,
        save="_check_markers.png",
    )
    return adata


def cleanfibroblasts(adata):
    ion()
    sc.set_figure_params(dpi=400, dpi_save=400, figsize=(4, 4))
    adata = adata[~adata.obs["leiden"].isin(["0", "1", "2", "5", "6", "9"])]
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.umap(adata)
    sc.pl.umap(
        adata,
        color=[
            "Col1a1",
            "Bgn",
            "Dcn",
            "Pdgfra",
            "Pdgfrb",
            "Vim",
            "Tnc",
            "Alcam",
            "Col1a2",
            "Fgfr2",
            "Fgfr1",
            "Fgfr3",
            "Fgfr4",
            "Igf1r",
            "Igf2r",
            "Insr",
            "Insrr",
            "Met",
            "Mst1r",
            "Csf1r",
            "Flt3",
            "Kit",
            "Myh11",
            "Myl9",
            "Acta2",
        ],
        cmap="Reds",
        vmax=4,
        vmin=0,
        save="_fibroblast_myofibroblast_lineage_markers.png",
    )
    scvi.data.setup_anndata(adata, layer="counts", batch_key="batch")
    vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30)
    vae.train()
    adata.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.umap(adata)
    sc.tl.dendrogram(adata, groupby="leiden")
    basics(adata)
    adata.write_h5ad(h5ads + "/fapcm_fibroblasts_v6_clean.h5ad")
    return adata


if __name__ == "__main__":
    adata = sc.read_h5ad(h5ads + "/fapcm_fibroblasts_v6.h5ad")
    adata = cleanfibroblasts(adata)
