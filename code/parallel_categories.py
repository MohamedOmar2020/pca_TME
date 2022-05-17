# parallel categories plot
import plotly.express as px
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

adata_mouse = sc.read_h5ad('outs/h5ads/fapcm_fibroblasts_v6_clean_regulons_5.h5ad', chunk_size=100000)
# filter the na
adata_mouse.obs['cluster'] = adata_mouse.obs['cluster'].astype('str')
adata_mouse = adata_mouse[adata_mouse.obs['cluster'] != 'nan', :]
adata_mouse.obs['cluster'].value_counts()
# Arrange the cell annotation dataframe by the clusters
adata_mouse.obs.sort_values(by = ['cluster'], inplace=True)

color_map = dict(
    zip(adata_mouse.obs['key'].cat.categories.tolist(), sc.pl.palettes.vega_20_scanpy)
)

adata_mouse.obs["color"] = adata_mouse.obs['key'].map(color_map)
fig = px.parallel_categories(
    adata_mouse.obs,
    dimensions=['key', 'cluster'],
    color="color",
    labels={'key': "model", 'cluster': "cluster"},
)
fig.update_layout(autosize=False, width=500, height=1000, font_size = 18)
fig.write_image("figs_new/parallel_categories_model.png")

##################
adata_mouse.obs['cluster'] = adata_mouse.obs['cluster'].astype('category')
color_map = dict(
    zip(adata_mouse.obs['cluster'].cat.categories.tolist(), sc.pl.palettes.vega_20_scanpy)
)
adata_mouse.obs["color"] = adata_mouse.obs['cluster'].map(color_map)
fig = px.parallel_categories(
    adata_mouse.obs,
    dimensions=['cluster', 'key'],
    color="color",
    labels={'key': "model", 'cluster': "cluster"},
)
fig.update_layout(autosize=False, width=500, height=1000, font_size = 18)
fig.write_image("figs_new/parallel_categories_cluster.png")

####################
# generate greyscaled parallel categories

def isgrey(key, current, expected):
    if key in expected:
        return current
    else:
        return "#808080"


## parallel categories for each cluster category
adata_mouse.obs['cluster'] = adata_mouse.obs['cluster'].astype('category')
clusters = ["0"]
color_map = dict(
    zip(adata_mouse.obs['cluster'].cat.categories.tolist(), sc.pl.palettes.vega_20_scanpy)
)
colormap = color_map
colormap = {key: isgrey(key, val, clusters) for key, val in colormap.items()}
adata_mouse.obs["color"] = adata_mouse.obs['cluster'].map(colormap)
fig = px.parallel_categories(
    adata_mouse.obs,
    dimensions=['cluster', 'key'],
    color="color",
    labels={'key': "model", 'cluster': "cluster"},
)
fig.update_layout(
    autosize=False, width=500, height=1000, font=dict(size=22), margin=dict(r=180)
)
fig.write_image(f"figs_new/parallel_categories_0.png")


######################

clusters = ["1"]
colormap = color_map
colormap = {key: isgrey(key, val, clusters) for key, val in colormap.items()}
adata_mouse.obs["color"] = adata_mouse.obs['cluster'].map(colormap)
fig = px.parallel_categories(
    adata_mouse.obs,
    dimensions=['cluster', 'key'],
    color="color",
    labels={'key': "model", 'cluster': "cluster"},
)
fig.update_layout(
    autosize=False, width=500, height=1000, font=dict(size=22), margin=dict(r=180)
)
fig.write_image(f"figs_new/parallel_categories_1.png")

########################

clusters = ["2"]
colormap = color_map
colormap = {key: isgrey(key, val, clusters) for key, val in colormap.items()}
adata_mouse.obs["color"] = adata_mouse.obs['cluster'].map(colormap)
fig = px.parallel_categories(
    adata_mouse.obs,
    dimensions=['cluster', 'key'],
    color="color",
    labels={'key': "model", 'cluster': "cluster"},
)
fig.update_layout(
    autosize=False, width=500, height=1000, font=dict(size=22), margin=dict(r=180)
)
fig.write_image(f"figs_new/parallel_categories_2.png")


######################

clusters = ["3", "4"]
colormap = color_map
colormap = {key: isgrey(key, val, clusters) for key, val in colormap.items()}
adata_mouse.obs["color"] = adata_mouse.obs['cluster'].map(colormap)
fig = px.parallel_categories(
    adata_mouse.obs,
    dimensions=['cluster', 'key'],
    color="color",
    labels={'key': "model", 'cluster': "cluster"},
)
fig.update_layout(
    autosize=False, width=500, height=1000, font=dict(size=22), margin=dict(r=180)
)
fig.write_image(f"figs_new/parallel_categories_3_4.png")

#########################

clusters = ["5", "6", "7"]
colormap = color_map
colormap = {key: isgrey(key, val, clusters) for key, val in colormap.items()}
adata_mouse.obs["color"] = adata_mouse.obs['cluster'].map(colormap)
fig = px.parallel_categories(
    adata_mouse.obs,
    dimensions=['cluster', 'key'],
    color="color",
    labels={'key': "model", 'cluster': "cluster"},
)
fig.update_layout(
    autosize=False, width=500, height=1000, font=dict(size=22), margin=dict(r=180)
)
fig.write_image(f"figs_new/parallel_categories_5_6_7.png")
