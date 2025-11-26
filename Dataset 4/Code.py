import scanpy as sc
import anndata as ad
import pandas as pd
import os
import scipy.sparse as sp



# ------------------------------- Load the main .h5ad file
adata = ad.read_h5ad("C:/Users/z5551702/OneDrive - UNSW/Aim 1/Data/Dataset 4/Chow_CanDisc_2023_data.h5ad")

print(adata)

# ------------------------------- output folders
output_dir_h5ad = "C:/Users/z5551702/OneDrive - UNSW/Aim 1/Data/Dataset 4/Processed/H5ad/"
output_dir_excel = "C:/Users/z5551702/OneDrive - UNSW/Aim 1/Data/Dataset 4/Processed/Excel/"


# ------------------------------- 1. Gene expression matrix (sparse) -> H5AD
X_sparse = sp.csr_matrix(adata.X) if not sp.issparse(adata.X) else adata.X
adata_genes = ad.AnnData(X=X_sparse, var=adata.var.copy(), obs=adata.obs.copy())
adata_genes.write_h5ad(os.path.join(output_dir_h5ad, "gene_expression_sparse.h5ad"))

# ------------------------------- 2. Metadata splitting
tcr_cols = ['TRAV', 'TRAJ', 'TRAC', 'cdr3_TRA', 'TRBV', 'TRBD', 'TRBJ', 'TRBC', 'cdr3_TRB', 'cdr3_full']
all_cols = adata.obs.columns.tolist()
meta_cols = [c for c in all_cols if c not in tcr_cols]

# Metadata without TCR
metadata = adata.obs[meta_cols]
metadata.to_csv(os.path.join(output_dir_excel, "metadata.csv"))
metadata_adata = ad.AnnData(obs=metadata)
metadata_adata.write_h5ad(os.path.join(output_dir_h5ad, "metadata.h5ad"))

# TCR metadata
tcr_data = adata.obs[tcr_cols]
tcr_data.to_csv(os.path.join(output_dir_excel, "tcr_data.csv"))
tcr_adata = ad.AnnData(obs=tcr_data)
tcr_adata.write_h5ad(os.path.join(output_dir_h5ad, "tcr_data.h5ad"))

# ------------------------------- 3. Save a 20x20 subset of gene expression to Excel
# Convert sparse to dense temporarily
X_dense = adata_genes.X.toarray() if sp.issparse(adata_genes.X) else adata_genes.X

# Take first 20 cells and first 20 genes
subset_df = pd.DataFrame(
    X_dense[:20, :20],
    index=adata_genes.obs.index[:20],
    columns=adata_genes.var.index[:20]
)

subset_df.to_excel(os.path.join(output_dir_excel, "gene_expression_subset_20x20.xlsx"))
