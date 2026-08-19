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
                'fontSize': '3vw',
                'letterSpacing': '0.14em',
                'marginBottom': '2vw',
                'lineHeight': '1.05',
                'textTransform': 'uppercase',
                'color': 'white',
                'marginLeft': '2vw',
                'marginTop': '2vw'
            }
        ),
        dcc.Tabs(
            className="custom-tabs",
            children=[
                dcc.Tab(label='QUALITY CONTROL', className="custom-tab", selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("QC Metrics", style={'color': 'White', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "Quality control (QC) metrics assess the quality of individual cells and genes in your dataset. "
                            "These plots help identify low-quality cells (e.g., with high mitochondrial content or low gene counts) "
                            "that may need to be filtered before downstream analysis. The plots below are quality control metrtics of your data post filtering",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        dbc.Row([
                            dbc.Col(dcc.Graph(id='violin-counts'), width=4),
                            dbc.Col(dcc.Graph(id='violin-genes'), width=4),
                            dbc.Col(dcc.Graph(id='violin-mito'), width=4),
                        ], className="mb-4"),
                        dbc.Row([
                            dbc.Col(dcc.Graph(id='ridge-counts'), width=4),
                            dbc.Col(dcc.Graph(id='ridge-genes'), width=4),
                            dbc.Col(dcc.Graph(id='ridge-mito'), width=4),
                        ])
                    ])
                ]),
                dcc.Tab(label='PCA', className="custom-tab", selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("PCA Variance Explained", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "View the proportion of variance explained by each principal component, providing insight into the main sources of variation in your data.",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        dcc.Graph(id='pca-variance-plot', style={'width': '80vw', 'height': '80vw'})
                    ])
                ]),
                dcc.Tab(label='UMAP', className="custom-tab", selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("2D UMAP colored by clusters", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "Visualize cell clusters in a 2D space using UMAP, revealing relationships and groupings among cells based on gene expression patterns.",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        dcc.Graph(id='umap-plot', style={'width': '80vw', 'height': '80vw'})
                    ])
                ]),
                dcc.Tab(label='UMAP 3D', className="custom-tab", 
                        selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("3D UMAP colored by clusters", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "Explore clusters in three dimensions with 3D UMAP, offering an enhanced perspective on cellular relationships and heterogeneity.",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        dcc.Graph(id='umap-3d-plot', style={'width': '80vw', 'height': '80vw'})
                    ])
                ]),
                dcc.Tab(label='t-SNE', className="custom-tab", selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("2D t-SNE colored by clusters", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "Examine cell populations in 2D using t-SNE, an alternative dimensionality reduction technique for visualizing complex single-cell data.",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        dcc.Graph(id='tsne-plot', style={'width': '80vw', 'height': '80vw'})
                    ])
                ]),
                dcc.Tab(label='TriMap', className="custom-tab", selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("2D TriMap colored by clusters", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "Investigate cell clusters using TriMap, which preserves global data structure in a 2D embedding for single-cell visualization.",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        dcc.Graph(id='trimap-plot', style={'width': '80vw', 'height': '80vw'})
                    ])
                ]),
                dcc.Tab(label='DIFFUSION MAP', className="custom-tab", selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("2D Diffusion Map colored by clusters", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "Analyze cellular trajectories and transitions with Diffusion Map, highlighting continuous processes such as differentiation.",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        dcc.Graph(id='diffmap-plot', style={'width': '80vw', 'height': '80vw'})
                    ])
                ]),
                dcc.Tab(label='PHATE', className="custom-tab", selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("2D PHATE colored by clusters", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "PHATE visualization captures both local and global structures, making it ideal for exploring developmental trajectories in single-cell data.",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        dcc.Graph(id='phate-plot', style={'width': '80vw', 'height': '80vw'})
                    ])
                ]),
                dcc.Tab(label='GENE EXPRESSION', className="custom-tab", selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("Gene Expression Visualization", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "Interactively visualize the expression of selected genes across cells and clusters using your chosen dimensionality reduction method.",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        html.Label("Select Reduction Method:", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        dcc.Dropdown(
                            id='reduction-selector',
                            options=[
                                {'label': 'UMAP', 'value': 'umap'},
                                {'label': 't-SNE', 'value': 'tsne'},
                                {'label': 'TriMap', 'value': 'trimap'},
                                {'label': 'Diffusion Map', 'value': 'diffmap'},
                                {'label': 'PHATE', 'value': 'phate'}
                            ],
                            value='umap',
                            style={'backgroundColor': '#222', 'color': 'black', 'fontFamily': "'DM Mono', monospace"}
                        ),
                        html.Label("Select Gene:", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        dcc.Dropdown(
                            id='gene-selector',
                            options=[{'label': gene, 'value': gene} for gene in data['genes']],
                            value='CD19' if 'CD19' in data['genes'] else data['genes'][0] if data['genes'] else None,
                            style={'backgroundColor': '#222', 'color': 'black', 'fontFamily': "'DM Mono', monospace"}
                        ),
                        dcc.Graph(id='gene-expression-plot', style={'width': '80vw', 'height': '80vw'})
                    ])
                ]),
                dcc.Tab(label='HEATMAP', className="custom-tab", selected_className="custom-tab--selected",
                        style=tab_style, selected_style=tab_selected_style, children=[
                    html.Div([
                        html.H3("Marker Gene Expression Heatmap", style={'color': 'white', 'fontFamily': "'DM Mono', monospace"}),
                        html.P(
                            "Compare marker gene expression across clusters in a heatmap, facilitating the identification of cluster-specific gene signatures.",
                            style={'color': '#ccc', 'fontFamily': "'DM Mono', monospace", 'fontSize': '1.1vw', 'marginBottom': '1.5vw'}
                            ),
                        dcc.Dropdown(
                            id='heatmap-cluster-selector',
                            options=[{'label': f'Cluster {i}', 'value': i} for i in sorted(data['clusters'].unique())],
                            value=[sorted(data['clusters'].unique())[0]],
                            multi=True,
                            style={'backgroundColor': '#222', 'color': 'black', 'width': '80%', 'fontFamily': "'DM Mono', monospace"}
                        ),
                        dcc.Graph(id='marker-heatmap', style={'width': '90vw', 'height': '80vh'})
                    ])
                ]),
            ]
        )
    ]
)

@dash.callback(
    [Output('violin-counts', 'figure'),
     Output('violin-genes', 'figure'),
     Output('violin-mito', 'figure'),
     Output('ridge-counts', 'figure'),
     Output('ridge-genes', 'figure'),
     Output('ridge-mito', 'figure')],
    [Input('violin-counts', 'id')]  # Dummy input
)
def update_qc_metrics(_):
    # Create QC DataFrame with cluster information
    qc_data = pd.DataFrame({
        'n_counts': adata.obs['total_counts'] if 'total_counts' in adata.obs else np.zeros(adata.n_obs),
        'n_genes': adata.obs['n_genes'] if 'n_genes' in adata.obs else np.zeros(adata.n_obs),
        'percent_mito': adata.obs['pct_counts_mt'] if 'pct_counts_mt' in adata.obs else np.zeros(adata.n_obs),
        'Cluster': data['clusters'].astype(str)
    })

    # Create individual violin plots
    violin_counts = px.violin(qc_data, y='n_counts', box=True, points='all',
                             title='Total UMI Counts Distribution')
    violin_genes = px.violin(qc_data, y='n_genes', box=True, points='all',
                            title='Genes Detected Distribution')
    violin_mito = px.violin(qc_data, y='percent_mito', box=True, points='all',
                           title='Mitochondrial % Distribution')

    # Create ridge plots using horizontal violins
    def create_ridge_plot(metric, title):
        fig = px.violin(qc_data, x=metric, y='Cluster', color='Cluster',
                        orientation='h', box=True, points=False,
                        title=title, height=800)
        fig.update_traces(side='positive', width=1.5, line_color='pink',
                         meanline_visible=True, scalemode='count')
        fig.update_layout(showlegend=False, xaxis_showgrid=False)
        return fig

    ridge_counts = create_ridge_plot('n_counts', 'Total Counts per Cluster')
    ridge_genes = create_ridge_plot('n_genes', 'Genes Detected per Cluster')
    ridge_mito = create_ridge_plot('percent_mito', 'Mitochondrial % per Cluster')

    # Apply dark theme to all figures
    return [dark_plotly_layout(fig) for fig in [
        violin_counts, violin_genes, violin_mito,
        ridge_counts, ridge_genes, ridge_mito
    ]]

@dash.callback(
    Output('pca-variance-plot', 'figure'),
    Input('pca-variance-plot', 'id')
)
def update_pca_variance(_):
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=list(range(1, min(51, len(data['pca_variance'])+1))),
        y=data['pca_variance'][:50],
        mode='lines+markers',
        name='Variance Explained'
    ))
    fig.update_layout(
        title='PCA Variance Explained (up to 50 dimensions)',
        xaxis_title='Principal Component',
        yaxis_title='Variance Explained',
        showlegend=True
    )
    return dark_plotly_layout(fig)

@dash.callback(
    Output('umap-plot', 'figure'),
    Input('umap-plot', 'id')
)
def update_umap(_):
    fig = px.scatter(
        data['umap'],
        x='UMAP1', y='UMAP2',
        color='Cluster',
        title='2D UMAP colored by clusters',
        color_discrete_sequence=px.colors.qualitative.Dark24
    )
    return dark_plotly_layout(fig)

@dash.callback(
    Output('umap-3d-plot', 'figure'),
    Input('umap-3d-plot', 'id')
)
def update_umap_3d(_):
    if 'umap_3d' not in data:
        return go.Figure().update_layout(
            title='3D UMAP Not Available',
            paper_bgcolor='#222',
            font=dict(color='white')
        )
    
    fig = px.scatter_3d(
        data['umap_3d'],
        x='UMAP1', y='UMAP2', z='UMAP3',
        color='Cluster',
        title='3D UMAP colored by clusters',
        color_discrete_sequence=px.colors.qualitative.Dark24,
        height=800
    )
    return dark_plotly_layout(fig, is3d=True)

@dash.callback(
    Output('tsne-plot', 'figure'),
    Input('tsne-plot', 'id')
)
def update_tsne(_):
    fig = px.scatter(
        data['tsne'],
        x='tSNE1', y='tSNE2',
        color='Cluster',
        title='2D t-SNE colored by clusters',
        color_discrete_sequence=px.colors.qualitative.Vivid
    )
    return dark_plotly_layout(fig)

@dash.callback(
    Output('trimap-plot', 'figure'),
    Input('trimap-plot', 'id')
)
def update_trimap(_):
    fig = px.scatter(
        data['trimap'],
        x='TriMap1', y='TriMap2',
        color='Cluster',
        title='2D TriMap colored by clusters',
        color_discrete_sequence=px.colors.qualitative.Bold
    )
    return dark_plotly_layout(fig)

@dash.callback(
    Output('diffmap-plot', 'figure'),
    Input('diffmap-plot', 'id')
)
def update_diffmap(_):
    fig = px.scatter(
        data['diffmap'],
        x='DiffMap1', y='DiffMap2',
        color='Cluster',
        title='2D Diffusion Map colored by clusters',
        color_discrete_sequence=px.colors.qualitative.Pastel
    )
    return dark_plotly_layout(fig)

@dash.callback(
    Output('phate-plot', 'figure'),
    Input('phate-plot', 'id')
)
def update_phate(_):
    fig = px.scatter(
        data['phate'],
        x='PHATE1', y='PHATE2',
        color='Cluster',
        title='2D PHATE colored by clusters',
        color_discrete_sequence=px.colors.qualitative.Prism
    )
    return dark_plotly_layout(fig)

@dash.callback(
    Output('gene-expression-plot', 'figure'),
    [Input('gene-selector', 'value'),
     Input('reduction-selector', 'value')]
)
def update_gene_expression(selected_gene, reduction_method):
    # Handle None or invalid gene selection
    if selected_gene is None or selected_gene not in data['genes']:
        return dash.no_update  # Prevent unnecessary updates

    # Get the appropriate dimensionality reduction data based on the selection
    reduction_data = data.get(reduction_method)
    if reduction_data is None:
        return dash.no_update
    
    # Get column names for the selected reduction method
    if reduction_method == 'umap':
        x_col, y_col = 'UMAP1', 'UMAP2'
        title_prefix = 'UMAP'
    elif reduction_method == 'tsne':
        x_col, y_col = 'tSNE1', 'tSNE2'
        title_prefix = 't-SNE'
    elif reduction_method == 'trimap':
        x_col, y_col = 'TriMap1', 'TriMap2'
        title_prefix = 'TriMap'
    elif reduction_method == 'diffmap':
        x_col, y_col = 'DiffMap1', 'DiffMap2'
        title_prefix = 'Diffusion Map'
    elif reduction_method == 'phate':
        x_col, y_col = 'PHATE1', 'PHATE2'
        title_prefix = 'PHATE'
    else:
        return dash.no_update

    # Extract gene expression values
    expression_values = adata[:, selected_gene].X.flatten()
    
    # Create the scatter plot with a new color scheme
    fig = px.scatter(
        x=reduction_data[x_col],
        y=reduction_data[y_col],
        color=expression_values,
        title=f'{title_prefix} Feature Plot for {selected_gene}',
        labels={'x': x_col, 'y': y_col, 'color': f'{selected_gene} Expression'},
        hover_name=adata.obs_names,
        color_continuous_scale='Agsunset'  
    )

    # Add a colorbar title and improve layout
    fig.update_coloraxes(colorbar_title=f'{selected_gene} Expression')
    fig.update_layout(
        coloraxis_colorbar=dict(
            title=f'{selected_gene} Expression',
            thicknessmode="pixels", thickness=20,
            lenmode="pixels", len=300
        )
    )

    # Apply dark theme styling
    return dark_plotly_layout(fig)

@dash.callback(
    Output('marker-heatmap', 'figure'),
    [Input('heatmap-cluster-selector', 'value')]
)
def update_heatmap(selected_clusters):
    # Validate input
    if not selected_clusters or not data['clusters'].isin(selected_clusters).any():
        return go.Figure()
    
    # Filter data
    mask = data['clusters'].isin(selected_clusters)
    cluster_subset = data['clusters'][mask]
    expr_subset = data['expression'].loc[mask]

    # Get top 5 marker genes per cluster
    top_genes = []
    for cluster in selected_clusters:
        cluster_expr = expr_subset[cluster_subset == cluster]
        if cluster_expr.empty:
            continue
        mean_expr = cluster_expr.mean().sort_values(ascending=False)
        top_genes.extend(mean_expr.head(5).index.tolist())
    
    # Deduplicate while preserving order
    seen = set()
    unique_genes = [g for g in top_genes if not (g in seen or seen.add(g))]
    
    # Sort cells by cluster
    sorted_cells = []
    sorted_cluster_labels = []
    for cluster in selected_clusters:
        cluster_cells = cluster_subset[cluster_subset == cluster].index.tolist()
        sorted_cells.extend(cluster_cells)
        sorted_cluster_labels.extend([cluster] * len(cluster_cells))
    
    # Downsample for performance
    max_cells = 10000
    if len(sorted_cells) > max_cells:
        step = len(sorted_cells) // max_cells
        sorted_cells = sorted_cells[::step]
        sorted_cluster_labels = sorted_cluster_labels[::step]
    
    heatmap_data = expr_subset.loc[sorted_cells, unique_genes]

    # Z-score normalization
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(heatmap_data.T).T

    # Create figure
    fig = go.Figure(go.Heatmap(
        z=scaled_data,
        x=unique_genes,
        y=[f"Cluster {cl}" for cl in sorted_cluster_labels],  # Simplified label
        colorscale='rdpu',
        colorbar=dict(title='Z-score'),
        hoverinfo='x+y+z'
    ))
    
    return dark_plotly_layout(fig)