# Single Cell RNA-seq Visualization Suite

An interactive Dash dashboard for exploring processed single-cell RNA-seq data --
quality control metrics, multiple dimensionality reduction embeddings (PCA, UMAP
2D/3D, t-SNE, TriMap, Diffusion Map, PHATE), gene expression overlays, and a
marker gene heatmap, all without writing code. This is a singular place where a non coding biologist could derive insights about single cell scRNAseq data without delving into computational methods and implementations. 

The dashboard has been built, tested, and confirmed working
end-to-end against a real processed dataset. The preprocessing pipeline script
is included and based on an envisioned approach, but **has not yet been
independently re-run or verified by me from raw data** -- see
[Status](#status) below.

---

## Features

- **Overview** -- summary of the dataset and pipeline
- **Quality Control** -- violin and ridge plots of per-cell QC metrics (UMI counts, genes detected, mitochondrial %)
- **PCA** -- variance explained per principal component
- **UMAP (2D & 3D)**, **t-SNE**, **TriMap**, **Diffusion Map**, **PHATE** -- dimensionality reduction embeddings, each colored by Leiden cluster
- **Gene Expression** -- overlay expression of any gene onto any embedding
- **Marker Heatmap** -- top marker genes per cluster, z-scored
- Consistent cluster colors across every embedding tab
- Dark, minimal UI

---

## Project structure

```
single_cell_int_viz/
├── app.py                        # Dash multipage entrypoint
├── single_cell_processing.py     # Scanpy-based preprocessing pipeline (see Status)
├── execution.sh                  # CLI wrapper: process data and/or launch dashboard
├── config/
│   └── requirements.txt          # Pinned dependencies
├── data/
│   └── processed_data.h5ad       # Processed dataset (subsampled, ~2,453 cells)
├── pages/
│   ├── home.py                   # Landing page
│   ├── dashboard.py              # Main visualization dashboard
│   └── assets/
│       └── custom.css            # Dark theme styling
└── screenshots/                  # (reference images, if added)
```

---

## Installation

This project uses conda for environment management.

```bash
# Clone this repository
git clone https://github.com/rajantidke/single_cell_int_viz.git
cd single_cell_int_viz

# Create and activate the environment
conda create -n scviz python=3.12 -y
conda activate scviz

# Install dependencies
pip install -r config/requirements.txt
```

---

## Usage

Make sure the `scviz` conda environment is activated first (`conda activate scviz`).

### Launch the dashboard

A processed dataset is already included at `data/processed_data.h5ad`, so you
can launch directly:

```bash
bash execution.sh
```

or, without the wrapper script:

```bash
python3 app.py
```

Then open `http://127.0.0.1:8050` in your browser.

### Re-run preprocessing (optional, see Status)

If you have your own raw `data/sample_data.h5ad` (10x Genomics-format,
converted to AnnData) and want to regenerate `processed_data.h5ad`:

```bash
bash execution.sh --process
```

or directly:

```bash
python3 single_cell_processing.py
```

---

## Status

- **Dashboard:** fully built, tested, and confirmed working against the
  included processed dataset -- every tab and callback has been verified to
  render correctly.
- **Preprocessing pipeline (`single_cell_processing.py`):** present and based
  on the the following method: 
  - QC filtering, 
  - normalization, 
  - HVG selection, 
  - scaling, 
  - PCA, 
  - Leiden clustering, and 
  - six dimensionality reduction embeddings 
  
  I have not yet personally re-run it end-to-end
  from raw data on my own machine. The `processed_data.h5ad` currently
  bundled with this repo is a pre-generated subsample (~2,453 cells, 2,000
  HVGs) used to build and test the dashboard. Verifying the preprocessing
  script against raw data of my own is a planned next step.

---

## Data

The bundled `processed_data.h5ad` is a subsample derived from a public 10x
Genomics breast cancer dataset (GEO), after QC filtering, normalization, HVG
selection, and scaling. It includes precomputed PCA, UMAP (2D & 3D), t-SNE,
TriMap, PHATE, and diffusion map embeddings, along with Leiden cluster
assignments (25 clusters).

---

## Roadmap

- Independently verify `single_cell_processing.py` against raw data
- Support for additional input formats (Seurat objects, Loom, CSV matrices)
- Multi-modal data support (ATAC-seq, spatial transcriptomics, CITE-seq)
- Automated cell type annotation against public reference atlases
- Exportable, publication-ready figures

---

## Acknowledgments

Built with
[Scanpy](https://scanpy.readthedocs.io/en/stable/),
[Dash](https://dash.plotly.com/),
[TriMap](https://github.com/eamid/trimap), and
[PHATE](https://github.com/KrishnaswamyLab/PHATE).