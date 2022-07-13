import pandas as pd
import squidpy as sq

adata = sc.read("outs/h5ads/fapcm_unfiltered_v6.h5ad")

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
sc.pl.umap(adata, color="phenotypes")


adata.obs['phenotypes'].value_counts()

for key in adata.obs.key.cat.categories.tolist():
    lradata = adata[adata.obs["key"] == key]
    lradata.raw = lradata
    sq.gr.ligrec(
        lradata,
        n_perms=1000,
        cluster_key="phenotypes",
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        seed = 123, show_progress_bar = True, corr_method='fdr_bh', alpha = 0.05,
    )

sq.pl.ligrec(
        lradata,
        cluster_key="phenotypes",
        source_groups=["fibroblasts", "luminal", "basal"],
        target_groups=["fibroblasts", "luminal", "basal"],
        #means_range=(0.3, np.inf),
        alpha=0.05,
        pvalue_threshold=0.05,
        swap_axes=True,
        save=f"_{key}_lr_interactions.png",
    )

ligandreceptor_permutation_test(adata)

sc.pl.umap(adata)
adata.obs['key'].value_counts()
adata.obs['model'].value_counts()


pd.crosstab(adata.obs['key'], adata.obs['model'])