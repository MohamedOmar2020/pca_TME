
import os
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib.pyplot import ion
from scanpy.plotting import _utils
from scanpy.plotting._anndata import *
from statannotations.Annotator import Annotator

# panel 1
adata = sc.read('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/common/panel1.h5ad')

modeldict = {
    "38-1": "HiMyc",
    "4114": "MNRP$^{DKO}$ Wildtype",
    "4324": "MNRP$^{DKO}$ Wildtype",
    "4328": "MNRP$^{DKO}$ Wildtype",
    "3672650-164": "HiMyc",
    "3674698-12": "HiMyc",
    "3695824-1": "HiMyc Wildtype",
    "3695824-2": "HiMyc Wildtype",
    "3725965-171": "Pten Wildtype",
    "3778541-26": "HiMyc Wildtype",
    "3852587-54": "TRG",
    "3862892-201": "TRG",
    "3889446-217": "TRG",
    "3889446-220": "TRG",
    "3918561-44": "TRG Wildtype",
    "FVBN_6mo_2": "TRG Wildtype",
    "FVBN_6mo_3": "TRG Wildtype",
    "M2861": "MNRP$^{DKO}$",
    "M2872": "MNRP$^{DKO}$",
    "M3056": "MNRP$^{DKO}$",
    "Nkx_3.1_Pten_1": "PtenKO",
    "Nkx_3.1_Pten_2": "PtenKO",
    "Nkx_3.1_Pten_3": "PtenKO",
    "Nkx_3.1_WT_1": "PtenKO Wildtype",
    "Nkx_3.1_WT_2": "PtenKO Wildtype",
}
adata.obs["pheno"] = adata.obs["name"].map(modeldict)

sc.pp.log1p(adata)
adata.raw = adata
sc.pp.scale(adata, max_value=10)


stroma = adata[adata.obs["halo_label"] == "STROMA"]

adata.obs["pheno"] = adata.obs["pheno"].astype("category")





for model in [
    ["TRG", "TRG Wildtype", ],
    ["PtenKO", "PtenKO Wildtype", ],
    ["HiMyc", "HiMyc Wildtype", ],
    ["MNRP$^{DKO}$", "MNRP$^{DKO}$ Wildtype", ],
]:
    for marker in [
        "DAPI",
        "GPX3",
        "C3",
        "CD55",
        "NFKB1",
        "SFRP1",
        "CK",
    ]:
        stromax = stroma[stroma.obs["pheno"].isin(model)]

        sanitize_anndata(stromax)
        obs_df = get.obs_df(
            stromax, keys=[marker] + ["pheno"], layer=None, use_raw=False
        )
        obs_tidy = obs_df

        # compute p values using averages in each slide
        from scipy.stats import mannwhitneyu, ttest_ind, wilcoxon

        stats_df = get.obs_df(
            stromax, keys=[marker] + ["pheno"] + ["name"], layer=None, use_raw=False
        )
        stats_df = (
            stats_df.groupby(["name", "pheno"], as_index=False).mean().dropna()
        )
        pair = stats_df.pheno.astype("category").cat.categories.tolist()
        assert len(pair) == 2
        data1 = stats_df.groupby("pheno")[marker].get_group(pair[0])
        data2 = stats_df.groupby("pheno")[marker].get_group(pair[1])
        stat, p = ttest_ind(data1, data2)
        pvalues = [p]
        print(pvalues)

        groupby = "pheno"
        x = "pheno"
        ylabel = [marker]
        ys = [marker]

        axs, _, _, _ = setup_axes(
            ax=None, panels=[marker], show_ticks=True, right_margin=0.3,
        )
        for ax, y, ylab in zip(axs, ys, ylabel):
            ax = sns.violinplot(
                x=x,
                y=y,
                data=obs_tidy,
                order=None,
                orient="vertical",
                ax=ax,
                cut=0,
            )
            ax = sns.stripplot(
                x=x,
                y=y,
                data=obs_tidy,
                order=None,
                jitter=True,
                color="black",
                size=0.005,
                ax=ax,
            )
            xlabel = groupby.replace("_", " ")
            ax.set_ylabel(ylab)
            ax.set_xlabel("condition")
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            show = settings.autoshow
            pairs = [(model[0], model[1])]
            annotator = Annotator(ax, pairs, data=obs_tidy, x=x, y=y)
            annotator.configure(
                test=None, test_short_name="", text_format="star", loc="outside"
            )
            annotator.set_pvalues(pvalues=pvalues)
            annotator.annotate()
            # annotator.apply_and_annotate()
            _utils.savefig_or_show(
                "violin", show=show, save=f"_{model[0]}_{marker}.png"
            )


##################################

# panel 2
adata_panel2 = sc.read('/Users/mohamedomar/Documents/Research/Projects/Pca_TME/vectra/wnt/panel2.h5ad')

modeldict = {
    "38-1": "HiMyc",
    "4114": "MNRP$^{DKO}$ Wildtype",
    "4324": "MNRP$^{DKO}$ Wildtype",
    "4328": "MNRP$^{DKO}$ Wildtype",
    "3672650-164": "HiMyc",
    "3674698-12": "HiMyc",
    "3695824-1": "HiMyc Wildtype",
    "3695824-2": "HiMyc Wildtype",
    "3725965-171": "Pten Wildtype",
    "3778541-26": "HiMyc Wildtype",
    "3852587-54": "TRG",
    "3862892-201": "TRG",
    "3889446-217": "TRG",
    "3889446-220": "TRG",
    "3918561-44": "TRG Wildtype",
    "FVBN_6mo_2": "TRG Wildtype",
    "FVBN_6mo_3": "TRG Wildtype",
    "M2861": "MNRP$^{DKO}$",
    "M2872": "MNRP$^{DKO}$",
    "M3056": "MNRP$^{DKO}$",
    "Nkx_3.1_Pten_1": "PtenKO",
    "Nkx_3.1_Pten_2": "PtenKO",
    "Nkx_3.1_Pten_3": "PtenKO",
    "Nkx_3.1_WT_1": "PtenKO Wildtype",
    "Nkx_3.1_WT_2": "PtenKO Wildtype",
}

adata_panel2.obs["pheno"] = adata_panel2.obs["name"].map(modeldict)


sc.pp.log1p(adata_panel2)
adata_panel2.raw = adata_panel2
sc.pp.scale(adata_panel2, max_value=10)


stroma_panel2 = adata_panel2[adata_panel2.obs["halo_label"] == "STROMA"]

sc.pl.violin(
    stroma_panel2,
    groupby="pheno",
    rotation=90.0,
    keys=["DAPI", "LEF1", "WIF1", "WNT2", "SFRP2", "LGR5", "PanCK", ],
    order=[
        "TRG",
        "TRG Wildtype",
        "PtenKO",
        "PtenKO Wildtype",
        "HiMyc",
        "HiMyc Wildtype",
        "MNRP$^{DKO}$",
        "MNRP$^{DKO}$ Wildtype",
    ],
    save="all_violin_panel2.png",
)


sns.set(font_scale=0.8)
adata_panel2.obs["pheno"] = adata_panel2.obs["pheno"].astype("category")

for model in [
    ["TRG", "TRG Wildtype", ],
    ["PtenKO", "PtenKO Wildtype", ],
    ["HiMyc", "HiMyc Wildtype", ],
    ["MNRP$^{DKO}$", "MNRP$^{DKO}$ Wildtype", ],
]:
    for marker in [
        "DAPI",
        "LEF1",
        "WIF1",
        "WNT2",
        "SFRP2",
        "LGR5",
        "PanCK",
    ]:
        stromax = stroma_panel2[stroma_panel2.obs["pheno"].isin(model)]
        sanitize_anndata(stromax)
        obs_df = get.obs_df(
            stromax, keys=[marker] + ["pheno"], layer=None, use_raw=False
        )
        obs_tidy = obs_df
        print(obs_df)
        # compute p values using averages in each slide
        from scipy.stats import mannwhitneyu, ttest_ind, wilcoxon

        stats_df = get.obs_df(
            stromax, keys=[marker] + ["pheno"] + ["name"], layer=None, use_raw=False
        )
        stats_df = (
            stats_df.groupby(["name", "pheno"], as_index=False).mean().dropna()
        )
        pair = stats_df.pheno.astype("category").cat.categories.tolist()
        assert len(pair) == 2
        data1 = stats_df.groupby("pheno")[marker].get_group(pair[0])
        data2 = stats_df.groupby("pheno")[marker].get_group(pair[1])
        stat, p = ttest_ind(data1, data2)
        pvalues = [p]
        print(pvalues)

        groupby = "pheno"
        x = "pheno"
        ylabel = [marker]
        ys = [marker]

        axs, _, _, _ = setup_axes(
            ax=None, panels=[marker], show_ticks=True, right_margin=0.3,
        )
        for ax, y, ylab in zip(axs, ys, ylabel):
            ax = sns.violinplot(
                x=x,
                y=y,
                data=obs_tidy,
                order=None,
                orient="vertical",
                ax=ax,
                cut=0,
            )
            ax = sns.stripplot(
                x=x,
                y=y,
                data=obs_tidy,
                order=None,
                jitter=True,
                color="black",
                size=0.005,
                ax=ax,
            )
            xlabel = groupby.replace("_", " ")
            ax.set_ylabel(ylab)
            ax.set_xlabel("condition")
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            show = settings.autoshow
            pairs = [(model[0], model[1])]
            annotator = Annotator(ax, pairs, data=obs_tidy, x=x, y=y)
            annotator.configure(
                test=None, test_short_name="", text_format="star", loc="outside"
            )
            annotator.set_pvalues(pvalues=pvalues)
            annotator.annotate()
            # annotator.apply_and_annotate()
            _utils.savefig_or_show(
                "violin", show=show, save=f"_{model[0]}_{marker}.png"
            )