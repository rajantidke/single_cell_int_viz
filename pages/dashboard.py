import dash
from dash import dcc, html, Input, Output
import scanpy as sc
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from sklearn.preprocessing import StandardScaler
import logging
import dash_bootstrap_components as dbc
from pathlib import Path

dash.register_page(__name__, path='/dashboard')

# Get the directory where the current script is located
current_dir = Path(__file__).parent

# Construct the path to the data file
data_path = (current_dir / ".." / "data" / "processed_data.h5ad").resolve()

# Load processed data
adata = sc.read_h5ad(data_path)

# Precompute dimensionality reductions if not present
if 'X_umap' not in adata.obsm:
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
if 'X_tsne' not in adata.obsm:
    sc.tl.tsne(adata)

# Prepare data
if 'leiden' in adata.obs:
    clusters = adata.obs['leiden']
elif 'louvain' in adata.obs:
    clusters = adata.obs['louvain']
else:
    clusters = pd.Series(['0'] * adata.n_obs, index=adata.obs_names)

data = {
    'umap': pd.DataFrame(
        adata.obsm['X_umap_2d'],
        columns=['UMAP1', 'UMAP2'],
        index=adata.obs_names
    ),
    'tsne': pd.DataFrame(
        adata.obsm['X_tsne_2d'],
        columns=['tSNE1', 'tSNE2'],
        index=adata.obs_names
    ),
    'trimap': pd.DataFrame(
        adata.obsm['X_trimap'],
        columns=['TriMap1', 'TriMap2'],
        index=adata.obs_names
    ),
    'diffmap': pd.DataFrame(
        adata.obsm['X_diffmap'][:, :2],
        columns=['DiffMap1', 'DiffMap2'],
        index=adata.obs_names
    ),
    'phate': pd.DataFrame(
        adata.obsm['X_phate'],
        columns=['PHATE1', 'PHATE2'],
        index=adata.obs_names
    ),
    'genes': adata.var_names.tolist(),
    'expression': pd.DataFrame(
        adata.X.todense() if hasattr(adata.X, "todense") else adata.X,
        columns=adata.var_names,
        index=adata.obs_names
    ),
    'pca_variance': adata.uns['pca']['variance_ratio'] if 'pca' in adata.uns else np.zeros(50),
    'clusters': clusters
}

data['umap_3d'] = pd.DataFrame(
    adata.obsm['X_umap_3d'],
    columns=['UMAP1', 'UMAP2', 'UMAP3'],
    index=adata.obs_names
)
data['umap_3d']['Cluster'] = data['clusters'].values  # Add cluster info

# Add cluster information to all dimensionality reduction dataframes
for key in ['umap', 'tsne', 'trimap', 'diffmap', 'phate']:
    data[key]['Cluster'] = data['clusters'].values

# Theme styles
tab_style = {
    'backgroundColor': '#222',
    'color': 'white',
    'padding': '6px',
    'fontFamily': "'DM Mono', monospace",
    'fontWeight': 400,
    'letterSpacing': '0.03em'
}
tab_selected_style = {
    'backgroundColor': '#333',
    'color': 'white',
    'padding': '6px',
    'fontFamily': "'DM Mono', monospace",
    'fontWeight': 500,
    'letterSpacing': '0.03em'
}

def dark_plotly_layout(fig, is3d=False):
    font = dict(family="DM Mono, monospace", color="white")
    if is3d:
        fig.update_layout(
            scene=dict(
                xaxis=dict(
                    backgroundcolor='#111', color='white', gridcolor='grey', zerolinecolor='grey',
                    title_font=font, tickfont=font
                ),
                yaxis=dict(
                    backgroundcolor='#111', color='white', gridcolor='grey', zerolinecolor='grey',
                    title_font=font, tickfont=font
                ),
                zaxis=dict(
                    backgroundcolor='#111', color='white', gridcolor='grey', zerolinecolor='grey',
                    title_font=font, tickfont=font
                ),
            ),
            plot_bgcolor='#111',
            paper_bgcolor='#111',
            font=font,
            legend=dict(font=font),
            title_font=font
        )
    else:
        fig.update_layout(
            plot_bgcolor='#111',
            paper_bgcolor='#111',
            font=font,
            legend=dict(font=font),
            title_font=font,
            xaxis=dict(
                color='white', title_font=font, tickfont=font,
                gridcolor='grey', zerolinecolor='grey'
            ),
            yaxis=dict(
                color='white', title_font=font, tickfont=font,
                gridcolor='grey', zerolinecolor='grey'
            ),
        )
    return fig

layout = html.Div(
    style={
        'backgroundColor': '#111',
        'minHeight': '100vh',
        'color': 'white',
        'fontFamily': "'DM Mono', monospace",
        'padding': '0',
        'margin': '0'
    },
    children=[
        html.H1(
            "SINGLE CELL VISUALIZATION SUITE",
            style={
                'fontFamily': "'Montserrat', sans-serif",
                'fontWeight': 200,
                ...  # continues into next chunk