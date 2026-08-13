#!/usr/bin/env python3

#Primary imports
import scanpy as sc
import trimap  
import phate
from pathlib import Path

# Use pathlib to construct the path to the input and output data files
current_dir = Path(__file__).parent
input_file = (current_dir / "data" / "sample_data.h5ad").resolve()
output_file = (current_dir / "data" / "processed_data.h5ad").resolve()

# Check if the file exists
if not input_file.exists():
    raise FileNotFoundError(f"The file '{input_file}' does not exist.")

# Load the data
adata = sc.read_h5ad(input_file)
print("data loaded successfully.")

# Perform Quality Control (QC)
adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')

# calculate_qc_metrics computes total_counts, n_genes_by_counts, and
# pct_counts_mt (via qc_vars=['mt']) in a single pass, so no need to
# manually compute mitochondrial percentage beforehand.
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
adata = adata[adata.obs['n_genes_by_counts'] < 2500, :]
adata = adata[adata.obs['pct_counts_mt'] < 5, :]
print("QC metrics calculated and cells filtered based on n_genes_by_counts and pct_counts_mt.")

# Filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
print("Cells and genes filtered based on min_genes and min_cells.")

# Normalize the data
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("Data normalized and log-transformed.")

# Find variable features
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable]
print("Highly variable genes identified and data subsetted.")

# Scale the data
sc.pp.scale(adata, max_value=10)
print("Data scaled.")

# Run PCA
sc.tl.pca(adata, svd_solver='arpack')
print("PCA computed.")

# Find neighbors
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
print("Neighbors found.")

# Find clusters
sc.tl.leiden(adata)
print("Clusters found using Leiden algorithm.")

# Run UMAP for 2 components
# IMPORTANT: sc.tl.umap() writes its result to adata.obsm['X_umap'] each time it
# runs. Since we call it twice (once for 2D, once for 3D), we immediately
# copy each result into its own dedicated key before the next call
# overwrites 'X_umap'. This is intentional -- 'X_umap' itself should not be
# relied upon after this point, use 'X_umap_2d' / 'X_umap_3d' instead.
sc.tl.umap(adata, n_components=2)
adata.obsm["X_umap_2d"] = adata.obsm["X_umap"]
print("UMAP computed for 2 components.")

# Run UMAP for 3 components
sc.tl.umap(adata, n_components=3)
adata.obsm["X_umap_3d"] = adata.obsm["X_umap"]
print("UMAP computed for 3 components.")

# Run t-SNE for 2 components
sc.tl.tsne(adata, n_pcs=40, n_jobs=4, use_fast_tsne=True)
adata.obsm["X_tsne_2d"] = adata.obsm["X_tsne"]
print("t-SNE computed for 2 components.")

# TriMap implementation
trimap_embedding = trimap.TRIMAP().fit_transform(adata.obsm["X_pca"][:, :40])
adata.obsm["X_trimap"] = trimap_embedding
print("TriMap computed.")

# Diffusion Map
sc.tl.diffmap(adata)
print("Diffusion map computed.")

# PHATE implementation
phate_operator = phate.PHATE(n_components=2)
phate_embedding = phate_operator.fit_transform(adata.X)
adata.obsm["X_phate"] = phate_embedding
print("PHATE computed.")

# Save the processed data
adata.write(output_file)
print(f"Processed data saved to '{output_file}'.")