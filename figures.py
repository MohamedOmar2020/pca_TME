# figures for functional atlas of prostate cancer mesenchyme

import os
import sys

sys.path.insert(0, "/home/ryc4001/Documents/scutils/")
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import squidpy as sq
from matplotlib.pyplot import ion
from scipy.sparse import csr_matrix
from scutils.figures import (
    basics,
    canonicalfibroblastmarkers,
    correlationmatrix,
    distribution,
)
from scutils.qc import PreprocessRNA

import seaborn as sns

basepath = "/home/ryc4001/Documents/scproj/fapcm"
outspath = basepath + "/outs"
Path(outspath).mkdir(parents=True, exist_ok=True)
h5ads = outspath + "/h5ads"
Path(h5ads).mkdir(parents=True, exist_ok=True)


def main(adata):
    ion()
    sc.set_figure_params(
        dpi=800,
        dpi_save=800,
        frameon=True,
        figsize=(4, 4),
        format="png",
    )
    sns.set_style("ticks")
    # call analysis here

def figures_regulon(adata):
    ion()
    sc.set_figure_params(
        dpi=800, dpi_save=800, frameon=True, figsize=(4, 4), format="png"
    )
    regulons = [
        'Regulon(Arid5a)', 'Regulon(Arid5b)', 'Regulon(Ascl1)', 'Regulon(Ascl2)', 'Regulon(Atf3)', 'Regulon(Bach1)', 'Regulon(Batf)', 'Regulon(Bcl3)', 'Regulon(Cebpa)', 'Regulon(Cebpb)', 'Regulon(Cebpd)', 'Regulon(Creb5)', 'Regulon(Crem)', 'Regulon(Dusp26)', 'Regulon(Egr1)', 'Regulon(Egr2)', 'Regulon(Egr3)', 'Regulon(Egr4)', 'Regulon(Eomes)', 'Regulon(Erg)', 'Regulon(Ets1)', 'Regulon(Fezf1)', 'Regulon(Fosb)', 'Regulon(Fosl2)', 'Regulon(Foxa1)', 'Regulon(Foxd3)', 'Regulon(Foxi1)', 'Regulon(Foxo1)', 'Regulon(Foxq1)', 'Regulon(Foxs1)', 'Regulon(Gabpb1)', 'Regulon(Gata2)', 'Regulon(Gata3)', 'Regulon(Gata6)', 'Regulon(Grhl3)', 'Regulon(Hnf4a)', 'Regulon(Hoxb6)', 'Regulon(Ikzf2)', 'Regulon(Irf1)', 'Regulon(Irf4)', 'Regulon(Irf5)', 'Regulon(Irf6)', 'Regulon(Irf7)', 'Regulon(Irf8)', 'Regulon(Junb)', 'Regulon(Jund)', 'Regulon(Klf2)', 'Regulon(Klf4)', 'Regulon(Klf5)', 'Regulon(Lhx6)', 'Regulon(Mafb)', 'Regulon(Maff)', 'Regulon(Mef2c)', 'Regulon(Myc)', 'Regulon(Myod1)', 'Regulon(Nfe2l2)', 'Regulon(Nfia)', 'Regulon(Nfil3)', 'Regulon(Nfix)', 'Regulon(Nfkb1)', 'Regulon(Nkx6-2)', 'Regulon(Onecut2)', 'Regulon(Pax3)', 'Regulon(Peg3)', 'Regulon(Pgr)', 'Regulon(Pou2f3)', 'Regulon(Pparg)', 'Regulon(Prrx2)', 'Regulon(Rel)', 'Regulon(Runx1)', 'Regulon(Runx3)', 'Regulon(Six2)', 'Regulon(Snai3)', 'Regulon(Sox10)', 'Regulon(Sox11)', 'Regulon(Sox18)', 'Regulon(Sox2)', 'Regulon(Sox4)', 'Regulon(Sox7)', 'Regulon(Sox9)', 'Regulon(Spi1)', 'Regulon(Spib)', 'Regulon(Spic)', 'Regulon(Srebf1)', 'Regulon(Stat3)', 'Regulon(Tagln2)', 'Regulon(Tal1)', 'Regulon(Tbx1)', 'Regulon(Tbx21)', 'Regulon(Tcf4)', 'Regulon(Tead1)', 'Regulon(Tff3)', 'Regulon(Trp63)', 'Regulon(Twist1)',
    ]
    sc.tl.dendrogram(adata, groupby="cluster")
    sc.pl.heatmap(adata, regulons, dendrogram=True, show_gene_labels=True, groupby="cluster", swap_axes=True, cmap="Reds", save="_regulon_heatmap.png")


def figures_paper(adata):
    # final figures for paper
    genes = [
        "Ar",
        "Acta2",
        "Myl9",
        "Myh11",
        "Tagln",
        "Pdgfra",
        "Mustn1",
        "Angpt2",
        "Notch3",
        "Sfrp1",
        "Gpx3",
        "C3",
        "C7",
        "Cfh",
        "Ccl11",
        "Cd55",
        "Ptx3",
        "Thbd",
        "Ifi204",
        "Ifi205",
        "Ifi207",
        "Jun",
        "Junb",
        "Jund",
        "Fos",
        "Fosb",
        "Fosl2",
        "Atf3",
        "Mafb",
        "Maff",
        "Nek2",
        "Id1",
        "Id3",
        "Btg2",
        "Gadd45a",
        "Hes1",
        "Bcl3",
        "Socs1",
        "Socs3",
        "Il6",
        "Irf1",
        "Map3k8",
        "Gadd45b",
        "Gadd45g",
        "Dusp1",
        "Dusp6",
        "Klf4",
        "Sfrp2",
        "Wnt5a",
        "Lgr5",
        "Apc",
        "Wnt4",
        "Wnt6",
        "Notum",
        "Wif1",
        "Nkd1",
        "Fzd1",
        "Wnt2",
        "Wnt10a",
        "Dkk2",
        "Rorb",
        "Cxxc4",
        "Nfat5",
        "Apoe",
        "Dact1",
        "Ctnnb1",
        "Lef1",
        "Tcf4",
        "Myc",
        "Mki67",
        "H2afx",
        "Top2a",
        "Ccnb1",
        "Ccnb2",
        "Stmn1",
        "Ptn",
        "Mdk",
        "Tubb3",
        "Mrc2",
        "Fn1",
        "Tnc",
        "Col12a1",
        "Col14a1",
        "Col16a1",
        "Mmp19",
        "Cthrc1",
        "Wisp1",
        "Fzd1",
        "Fzd2",
        "Sfrp4",
        "Sfrp2",
        "Bmp1",
        "Tle3",
        "Tgfb1",
        "Tgfb2",
        "Tgfb3",
        "Postn",
    ]
    import matplotlib

    matplotlib.rcdefaults()
    sc.set_figure_params(
        dpi=800,
        dpi_save=800,
        frameon=True,
        figsize=(4, 4),
        fontsize=20,
        format="png",
    )
    # sns.set(font_scale=2)
    dp = sc.pl.dotplot(
        adata,
        var_names=genes,
        groupby="cluster",
        title="Marker Genes",
        cmap="Reds",
        dot_max=1.0,
        var_group_positions=[(0, 0), (1, 8), (9, 20), (21, 46), (47, 68), (69, 96)],
        var_group_labels=["misc", "0", "1", "2", "3, 4", "5, 6, 7"],
        var_group_rotation=0.0,
        expression_cutoff=0.5,
        return_fig=True,
    )
    dp.add_totals().savefig("markers.png")

    """
    regulons = [
        "Regulon(Srebf1)",
        "Regulon(Foxo1)",
        "Regulon(Arid5b)",
        "Regulon(Gata3)",
        "Regulon(Creb5)",
        "Regulon(Cebpa)",
        "Regulon(Gabpb1)",
        "Regulon(Onecut2)",
        "Regulon(Nfkb1)",
        "Regulon(Atf3)",
        "Regulon(Arid5a)",
        "Regulon(Arid5b)",
        "Regulon(Stat3)",
        "Regulon(Sox9)",
        "Regulon(Sox10)",
        "Regulon(Peg3)",
        "Regulon(Gata6)",
        "Regulon(Runx1)",
        "Regulon(Gata2)",
        "Regulon(Lhx6)",
        "Regulon(Snai3)",
    ]
    sc.pl.matrixplot(
        adata,
        var_names=regulons,
        standard_scale="var",
        swap_axes=True,
        groupby="cluster",
        cmap="RdBu_r",
        save="_regulons.png",
    )
    """


def figures_presentation(adata):
    # figures for presentation (results separated by category)
    ion()
    sc.set_figure_params(
        dpi=800,
        dpi_save=800,
        frameon=True,
        figsize=(4, 4),
        format="png",
    )
    sns.set_style("ticks")
    regulons = [
        "Regulon(Srebf1)",
        "Regulon(Foxo1)",
        "Regulon(Arid5b)",
        "Regulon(Gata3)",
        "Regulon(Creb5)",
        "Regulon(Cebpa)",
        "Regulon(Gabpb1)",
        "Regulon(Onecut2)",
        "Regulon(Nfkb1)",
        "Regulon(Atf3)",
        "Regulon(Arid5a)",
        "Regulon(Arid5b)",
        "Regulon(Stat3)",
        "Regulon(Sox9)",
        "Regulon(Sox10)",
        "Regulon(Foxo1)",
        "Regulon(Peg3)",
        "Regulon(Gata6)",
        "Regulon(Runx1)",
        "Regulon(Gata2)",
        "Regulon(Lhx6)",
        "Regulon(Snai3)",
    ]
    sc.pl.matrixplot(
        adata,
        var_names=regulons,
        standard_scale="var",
        groupby="leiden",
        cmap="RdBu_r",
        save="_regulons.png",
    )
    myo = [
        "Acta2",
        "Myl9",
        "Myh11",
        "Tagln",
        "Pdgfra",
        "Mustn1",
        "Angpt2",
        "Notch3",
    ]
    sc.pl.dotplot(
        adata,
        var_names=myo,
        groupby="leiden",
        cmap="Reds",
        vmax=6,
        save="_myofibroblasts.png",
    )
    myoadata = sc.read_h5ad(h5ads + "/fapcm_myo.h5ad")
    analyze_myofibroblasts(myoadata)
    c0 = [
        "Sfrp1",
        "Gpx3",
        "C3",
        "C7",
        "Cfh",
        "Ccl11",
        "Cd55",
        "Ptx3",
        "Thbd",
        "Ifi204",
        "Ifi205",
        "Ifi207",
    ]
    sc.pl.dotplot(
        adata,
        var_names=c0,
        groupby="leiden",
        cmap="Reds",
        vmax=6,
        save="_c0.png",
    )
    c11 = [
        "Jun",
        "Junb",
        "Jund",
        "Fos",
        "Fosb",
        "Fosl2",
        "Atf3",
        "Mafb",
        "Maff",
        "Nek2",
        "Id1",
        "Id3",
        "Btg2",
        "Gadd45a",
        "Hes1",
        "Bcl3",
        "Socs1",
        "Socs3",
        "Il6",
        "Irf1",
        "Map3k8",
        "Gadd45b",
        "Gadd45g",
        "Dusp1",
        "Dusp6",
        "Klf4",
    ]
    sc.pl.dotplot(
        adata,
        var_names=c11,
        groupby="leiden",
        cmap="Reds",
        vmax=6,
        save="_c11.png",
    )
    c3 = [
        "Sfrp2",
        "Wnt5a",
        "Lgr5",
        "Apc",
        "Wnt4",
        "Wnt6",
        "Notum",
        "Wif1",
        "Nkd1",
        "Fzd1",
        "Wnt2",
        "Wnt10a",
        "Dkk2",
        "Rorb",
        "Cxxc4",
        "Nfat5",
        "Apoe",
        "Dact1",
        "Ctnnb1",
        "Lef1",
        "Tcf4",
        "Myc",
    ]
    sc.pl.dotplot(
        adata,
        var_names=c3,
        groupby="leiden",
        cmap="Reds",
        vmax=6,
        save="_c3.png",
    )
    c1_4_5 = [
        "Mki67",
        "H2afx",
        "Top2a",
        "Ccnb1",
        "Ccnb2",
        "Stmn1",
        "Ptn",
        "Mdk",
        "Tubb3",
        "Mrc2",
        "Fn1",
        "Tnc",
        "Col12a1",
        "Col14a1",
        "Col16a1",
        "Mmp19",
        "Cthrc1",
        "Wisp1",
        "Fzd1",
        "Fzd2",
        "Lrp1b",
        "Sfrp4",
        "Sfrp2",
        "Bmp1",
        "Tle3",
        "Tgfb1",
        "Tgfb2",
        "Tgfb3",
    ]
    sc.pl.dotplot(
        adata,
        var_names=c1_4_5,
        groupby="leiden",
        cmap="Reds",
        vmax=6,
        save="_c1_4_5.png",
    )

    # ar
    sc.pl.umap(adata, color="Ar", cmap="RdBu_r", vmax=3, vmin=-1, save="ar.png")
    sc.pl.violin(adata, keys="Ar", groupby="cluster", save="ar.png")


def markers(adata, genedict, dotplot=True, umap=True):
    ion()
    if dotplot:
        for key, val in genedict.items():
            sc.pl.dotplot(
                adata,
                var_names=val,
                groupby="leiden",
                title=key,
                cmap="Reds",
                vmax=4,
                dendrogram=True,
                save=f"_{key}.png",
            )

        dp = sc.pl.dotplot(
            adata,
            var_names=genedict,
            groupby="leiden",
            title=None,
            cmap="Reds",
            vmax=4,
            dendrogram=True,
            return_fig=True,
        )
        dp.add_totals().show()
        dp.savefig("markers_dotplot.png")
        for genotype in adata.obs["key"].cat.categories.tolist():
            data = adata[adata.obs["key"] == genotype]
            sc.tl.dendrogram(data, groupby="leiden")
            dp = sc.pl.dotplot(
                data,
                var_names=genedict,
                groupby="leiden",
                title=None,
                cmap="Reds",
                vmax=4,
                dendrogram=True,
                return_fig=True,
            )
            dp.add_totals().show()
            dp.savefig(f"{genotype}_markers_dotplot.png")
    if umap:
        for key, val in genedict.items():
            sc.pl.umap(
                adata,
                color=val,
                title=key,
                cmap="Reds",
                vmax=4,
                save=f"_{key}_markers.png",
            )
            for genotype in adata.obs["key"].cat.categories.tolist():
                data = adata[adata.obs["key"] == genotype]
                sc.tl.dendrogram(data, groupby="leiden")
                sc.pl.umap(
                    data,
                    color=val,
                    title=key,
                    cmap="Reds",
                    vmax=4,
                    save=f"_{genotype}_{key}.png",
                )


def pipe_myofibroblasts(adata):
    """
    Generate myofibroblast subset, basic analysis.
    Adata from fapcm_fibroblasts_v6_clean_regulons_5.h5ad.
    """
    ion()
    sc.set_figure_params(
        dpi=800, dpi_save=800, frameon=True, figsize=(4, 4), format="png"
    )
    myo = adata[adata.obs["leiden"].isin(["2"])]
    sc.pl.umap(
        myo,
        color=["Acta2", "Myl9", "Myh11", "Pecam1", "Rgs5", "Mef2c"],
        cmap="Reds",
        vmax=5,
        save="_myofibroblast_subset_check.png",
    )
    marks = ["Acta2", "Myl9", "Myh11", "Pecam1", "Rgs5", "Mef2c"]
    filterbarcodes = None
    for mark in marks:
        if filterbarcodes is None:
            filterbarcodes = csr_matrix(myo[:, mark].X > 1).todense()
        else:
            filterbarcodes = filterbarcodes | csr_matrix(myo[:, mark].X > 1).todense()
    puremyo = myo[filterbarcodes].copy()
    sc.pl.umap(
        puremyo,
        color=[
            "Acta2",
            "Myl9",
            "Myh11",
            "Pecam1",
            "Rgs5",
            "Mef2c",
        ],
        vmax=5,
        cmap="Reds",
        save="_pure_myofibroblast_original_projection.png",
    )
    # compute new embedding
    scvi.data.setup_anndata(puremyo, layer="counts", batch_key="batch")
    vae = scvi.model.SCVI(puremyo, n_layers=2, n_latent=30)
    vae.train()
    puremyo.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(puremyo, use_rep="X_scVI")
    sc.tl.leiden(puremyo, resolution=0.3)
    sc.tl.umap(puremyo)
    sc.pl.umap(
        puremyo,
        color=[
            "Acta2",
            "Myl9",
            "Myh11",
            "Pecam1",
            "Rgs5",
            "Mef2c",
            "leiden",
            "key",
        ],
        cmap="Reds",
        vmax=8,
        save="_filtered_myofibroblasts_new_projection.png",
    )
    sc.pl.umap(
        puremyo, color="Ar", vmax=2, cmap="Reds", save="_filtered_myofibroblasts_ar.png"
    )
    # compute degs
    puremyo.obs["arp"] = csr_matrix(puremyo.raw[:, "Ar"].X > 0).todense()
    sc.pp.filter_cells(puremyo, min_genes=200)
    sc.pp.filter_genes(puremyo, min_cells=3)
    puremyo.obs["arp"] = puremyo.obs["arp"].astype("category")
    sc.tl.rank_genes_groups(
        puremyo, groupby="leiden", method="t-test_overestim_var", groups="all"
    )
    df = sc.get.rank_genes_groups_df(puremyo, group=None)
    df.to_csv("myo_differential_expression.csv")
    # regress against ar
    sc.tl.rank_genes_groups(
        puremyo, groupby="arp", method="t-test_overestim_var", groups="all"
    )
    sc.pl.rank_genes_groups_dotplot(
        puremyo, n_genes=20, vmax=14, save="_myofibroblasts_ar_markers.png"
    )
    # regress against ar, key as covariate
    sc.pp.regress_out(puremyo, keys=["key"])
    sc.tl.rank_genes_groups(
        puremyo, groupby="arp", method="t-test_overestim_var", groups="all"
    )
    sc.pl.rank_genes_groups_dotplot(
        puremyo, n_genes=20, vmax=14, save="_myofibroblasts_ar_markers_regressed.png"
    )
    sc.pp.filter_cells(puremyo, min_genes=200)
    sc.pp.filter_genes(puremyo, min_cells=3)
    sc.tl.rank_genes_groups(puremyo, method="t-test", groupby="arp", groups="all")
    sc.get.rank_genes_groups_df(puremyo, group=["True"])
    df = sc.get.rank_genes_groups_df(puremyo, group=["True"])
    df.to_csv("ar_positive_differential_expression.csv")
    del puremyo.obs["arp"]
    puremyo.write_h5ad("fapcm_myo.h5ad")


def analyze_myofibroblasts(adata):
    """
    Analyze myofibroblasts. Takes myofibroblast object from fapcm_myo.h5ad.
    """
    ion()
    sc.set_figure_params(
        dpi=300, dpi_save=300, frameon=True, figsize=(4, 4), format="png"
    )
    sns.set_style("ticks")
    mapdict = {
        "0": "myofibroblast",
        "1": "pericyte",
        "2": "myofibroblast",
        "3": "pericyte",
    }
    adata.obs["cluster"] = adata.obs["leiden"].map(mapdict)
    adata.obs["cluster"] = adata.obs["cluster"].astype("category")
    sc.tl.dendrogram(adata, groupby="cluster")
    sc.pl.dotplot(
        adata,
        var_names=[
            "Rgs5",
            "Mef2c",
            "Vtn",
            "Cygb",
            "Pdgfrb",
            "Rras",
            "Rasl12",
            "Rspo3",
            "Nrg2",
            "Hopx",
            "Actg2",
        ],
        groupby="cluster",
        cmap="Reds",
        vmax=6,
        save="_myo_markers.png",
    )
    regs = [
        "Regulon(Egr4)",
        "Regulon(Crem)",
        "Regulon(Sox4)",
        "Regulon(Egr2)",
        "Regulon(Cebpa)",
        "Regulon(Pparg)",
        "Regulon(Klf2)",
        "Regulon(Klf4)",
    ]
    sc.pl.matrixplot(
        adata,
        var_names=regs,
        groupby="cluster",
        standard_scale="var",
        cmap="RdBu_r",
        save="myo_regulons.png",
    )
    sc.tl.rank_genes_groups(
        adata, groupby="cluster", method="t-test_overestim_var", groups="all"
    )
    df = sc.get.rank_genes_groups_df(adata, group=None)
    df.to_csv("myo_t_test_overstim_var.csv")
    pd.DataFrame(adata.obs["cluster"]).to_csv("myo_phenotype.csv")


def ligandreceptor_permutation_test(adata):
    """
    Perform ligand-receptor analysis per CellPhoneDB (permutation test).
    Takes all cells, separated into fibroblast, luminal, basal compartments.
    """
    ion()
    adata.obs["fibroblasts"] = (
        adata.obs["leiden"].isin(["8", "13", "15", "26"]).astype("category")
    )
    adata.obs["luminal"] = (
        adata.obs["leiden"]
        .isin(["27", "24", "19", "18", "31", "11", "16", "21", "22", "6", "33", "23"])
        .astype("category")
    )
    adata.obs["basal"] = adata.obs["leiden"].isin(["2", "29"]).astype("category")
    adata.obs["phenotypes"] = "Other"
    adata.obs.loc[adata.obs.fibroblasts == True, "phenotypes"] = "fibroblasts"
    adata.obs.loc[adata.obs.luminal == True, "phenotypes"] = "luminal"
    adata.obs.loc[adata.obs.basal == True, "phenotypes"] = "basal"
    sc.pl.umap(adata, color="phenotypes", save="checkcompartments.png")

    # instead could look at unprocessed layer
    rawadata = sc.read_h5ad(
        "outs/h5ads/fapcm_unfiltered_v6.h5ad"
    )
    rawadata = rawadata[rawadata.obs_names.isin(adata.obs_names.tolist())]
    rawadata.obs["phenotypes"] = adata.obs["phenotypes"]
    adata = rawadata

    for key in adata.obs.key.cat.categories.tolist():
        lradata = adata[adata.obs["key"] == key]
        lradata.raw = lradata
        sq.gr.ligrec(
            lradata,
            n_perms=1000,
            cluster_key="phenotypes",
            transmitter_params={"categories": "ligand"},
            receiver_params={"categories": "receptor"},
            seed=123, show_progress_bar=True, corr_method='fdr_bh', alpha=0.05,
        )

        for i in lradata.uns['phenotypes_ligrec'].keys():
            lradata.uns['phenotypes_ligrec'][i].to_csv(f"objs/lr/{key}_"+str(i)+".csv")

        sq.pl.ligrec(
            lradata,
            cluster_key="phenotypes",
            source_groups=["fibroblasts", "luminal", "basal"],
            target_groups=["fibroblasts", "luminal", "basal"],
            means_range=(0.3, np.inf),
            alpha=0.05, pvalue_threshold=0.05,
            swap_axes=True,
            save=f"_{key}_lr_interactions.png",
        )


def regulons(adata):
    """
    Plot regulons. adata object must already have regulons added by pyscenic pipeline.
    """
    ion()
    sc.set_figure_params(
        dpi=800, dpi_save=800, frameon=True, figsize=(4, 4), format="png"
    )
    sc.pl.matrixplot(
        adata,
        var_names=[
            "Regulon(Egr4)",
            "Regulon(Gata2)",
            "Regulon(Sox11)",
            "Regulon(Ascl2)",
            "Regulon(Cebpa)",
            "Regulon(Creb5)",
            "Regulon(Nfkb1)",
            "Regulon(Gata6)",
            "Regulon(Runx1)",
            "Regulon(Sox4)",
            "Regulon(Pgr)",
            "Regulon(Arid5a)",
            "Regulon(Crem)",
            "Regulon(Mef2c)",
            "Regulon(Nfil3)",
            "Regulon(Pparg)",
            "Regulon(Arid5b)",
            "Regulon(Atf3)",
            "Regulon(Bach1)",
            "Regulon(Bcl3)",
            "Regulon(Cebpd)",
            "Regulon(Egr1)",
            "Regulon(Fosb)",
            "Regulon(Junb)",
            "Regulon(Jund)",
            "Regulon(Klf2)",
            "Regulon(Klf4)",
            "Regulon(Stat3)",
            "Regulon(Ascl1)",
            "Regulon(Pax3)",
            "Regulon(Snai3)",
            "Regulon(Srebf1)",
            "Regulon(Foxo1)",
            "Regulon(Gabpb1)",
            "Regulon(Gata3)",
            "Regulon(Onecut2)",
            "Regulon(Peg3)",
            "Regulon(Sox9)",
            "Regulon(Egr2)",
            "Regulon(Sox10)",
            "Regulon(Sox2)",
            "Regulon(Foxd3)",
            "Regulon(Lhx6)",
        ],
        cmap="RdBu_r",
        standard_scale="var",
        groupby="leiden",
        save="cluster_regulons.png",
    )


def ar(adata):
    ion()
    sc.set_figure_params(
        dpi=800, dpi_save=800, frameon=True, figsize=(4, 4), format="png"
    )
    sc.pl.violin(adata, keys="Ar", groupby="leiden", save="_ar.png")
    sc.pl.umap(adata, color="Ar", cmap="RdBu_r", vmin=-1, vmax=4, save="_ar.png")


def umapdegs(adata, path, bymodel=False, n_genes=200, method="wilcoxon"):
    """
    generate umap plots of degs
    """
    ion()
    sc.set_figure_params(dpi=100, dpi_save=100, figsize=(4, 4))
    for cluster in adata.obs.leiden.cat.categories.tolist():
        degs = pd.read_csv(
            os.path.join(
                path,
                f"c{cluster}.csv",
            )
        )
        clusterdegs = degs["Unnamed: 0"].tolist()[:n_genes]
        sc.pl.umap(
            adata,
            color=clusterdegs,
            cmap="Reds",
            vmin=0,
            vmax=4,
            size=7,
            frameon=True,
            save=f"_degs_cluster{cluster}.png",
        )
        if bymodel:
            for key in adata.obs.key.cat.categories.tolist():
                data = adata[adata.obs["key"] == key]
                sc.pl.umap(
                    data,
                    color=clusterdegs,
                    cmap="Reds",
                    vmin=0,
                    vmax=4,
                    size=7,
                    frameon=True,
                    save=f"_degs_cluster{cluster}_{key}.png",
                )


def adhoc(adata):
    genes = [
        "Col1a1",
        "Col1a2",
        "Col3a1",
        "Col4a1",
        "Col4a2",
        "Col5a1",
        "Col6a1",
        "Col8a1",
        "Col10a1",
        "Col11a1",
        "Col12a1",
        "Col13a1",
        "Col14a1",
        "Col15a1",
        "Col18a1",
        "Col23a1",
        "Bgn",
        "Dcn",
        "Lum",
        "Tagln",
        "Eln",
        "Fn1",
        "Mmp2",
        "Mmp3",
        "Mmp9",
        "Mmp10",
        "Mmp14",
        "Mmp19",
        "Serpine1",
        "Cthrc1",
        "Sulf1",
        "Tghfbi",
        "Inhba",
        "Egfl6",
        "Angpt2",
        "Pdgfa",
        "Acta2",
        "Myh11",
        "Pln",
        "Tpm1",
        "Tpm2",
        "Sorbs2",
    ]
    genes = tomouse(genes)
    genes = [gene for gene in genes if gene in adata.var_names.tolist()]
    sc.pl.dotplot(
        adata,
        groupby="leiden",
        var_names=genes,
        cmap="Reds",
        vmax=4,
        swap_axes=True,
        save="genes_cluster.png",
    )
    sc.pl.dotplot(
        adata,
        groupby="model",
        var_names=genes,
        cmap="Reds",
        vmax=4,
        swap_axes=True,
        save="genes_model.png",
    )


def tomouse(l):
    return [x.lower().capitalize() for x in l]


if __name__ == "__main__":
    adata = sc.read_h5ad(h5ads + "/fapcm_fibroblasts_v6_clean_regulons_5.h5ad")
    adata = adata[~adata.obs["leiden"].isin(["10"])]
    clusterdict = {
        "2": "0",
        "0": "1",
        "11": "2",
        "3": "3",
        "6": "4",
        "1": "5",
        "4": "6",
        "5": "7",
    }
    adata.obs["cluster"] = adata.obs["leiden"].map(clusterdict).astype("category")
    adata.obs["cluster"].cat.reorder_categories(
        ["0", "1", "2", "3", "4", "5", "6", "7"],
        inplace=True,
    )
    adata.obs["leiden"] = adata.obs["cluster"]
    main(adata)

