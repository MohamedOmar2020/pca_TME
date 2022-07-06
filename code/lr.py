import pandas as pd
import squidpy

adata = sc.read("outs/h5ads/fapcm_unfiltered_v6.h5ad")
ligandreceptor_permutation_test(adata)


adata.obs['key'].value_counts()
adata.obs['model'].value_counts()


pd.crosstab(adata.obs['key'], adata.obs['model'])