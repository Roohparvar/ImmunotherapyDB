import scanpy as sc
import anndata as ad
import pandas as pd
import os
import scipy.sparse as sp


# ------------------------------- Load the main .h5ad file
adata = ad.read_h5ad("C:/Users/z5551702/OneDrive - UNSW/Aim 1/Data/Dataset 2/Bassez_NatMed_2021_data.h5ad")

# ------------------------------- Print the contents of this file to see what it contains
print(adata)

adata.obs['cell_id'] = 'Dataset 2/' + adata.obs['cell_id'].astype(str)


import scanpy as sc
import anndata as ad
import pandas as pd
import scipy.sparse as sp

# ------------------------------- Load file
adata = ad.read_h5ad("C:/Users/z5551702/OneDrive - UNSW/Aim 1/Data/Dataset 2/Bassez_NatMed_2021_data.h5ad")

print(adata)

# ------------------------------- Extract first 5 cells
adata_5 = adata[:5, :].copy()

# ------------------------------- Select first 20 genes
gene_list = adata_5.var_names[:20]
adata_5_20 = adata_5[:, gene_list]

# ------------------------------- Convert Gene Expression (X) to DataFrame
if sp.issparse(adata_5_20.X):
    gene_expr = pd.DataFrame(adata_5_20.X.toarray(),
                             index=adata_5_20.obs_names,
                             columns=adata_5_20.var_names)
else:
    gene_expr = pd.DataFrame(adata_5_20.X,
                             index=adata_5_20.obs_names,
                             columns=adata_5_20.var_names)

# ------------------------------- Metadata
metadata = adata_5.obs.copy()

# ------------------------------- Write to Excel
output_path = "C:/Users/z5551702/OneDrive - UNSW/Aim 1/Data/Dataset 2/Processed/first5cells.xlsx"
with pd.ExcelWriter(output_path) as writer:
    gene_expr.to_excel(writer, sheet_name="GeneExpression_20genes")
    metadata.to_excel(writer, sheet_name="Metadata")



















"""
# ------------------------------- Add prefix to barcodes (convert to string first)

adata.obs['barcode'] = (
    "Dataset 2/" +
    adata.obs['barcode'].astype(str).str.split('_').str[-1]
)

# ------------------------------- Output directory and Excel file
output_dir = "C:/Users/z5551702/OneDrive - UNSW/Aim 1/Data/Dataset 2/Processed"
excel_file = os.path.join(output_dir, "ObsColumns_Descriptions.xlsx")

# ------------------------------- Define columns and their descriptions
columns_info = [
    ("original_id", "Original ID of the cell or sample in the raw dataset"),
    ("dataset", "Dataset from which this cell/sample comes"),
    ("batch", "Batch number or name for batch effect correction"),
    ("dataset_patient_id", "Patient ID within the specific dataset"),
    ("patient_id", "Standardized patient ID"),
    ("original_patient_id", "Original patient ID in the raw dataset"),
    ("cancer_type", "Type of cancer (e.g., Melanoma, Lung Cancer)"),
    ("ici_infused", "Whether the patient received Immune Checkpoint Inhibitor (ICI) treatment"),
    ("drug_name", "Name of drug administered"),
    ("N_ici_doses_infused", "Number of ICI doses received"),
    ("treatment_type", "Type of treatment (e.g., monotherapy, combination)"),
    ("previous_treatments_type", "Type of previous treatments"),
    ("previous_treatments_category", "Category of previous treatments (e.g., chemotherapy, radiation)"),
    ("group_response", "Patient response to treatment (responder/non-responder)"),
    ("group_response2", "Alternative or supplementary response grouping"),
    ("case_control_response", "Response in case/control format"),
    ("healthy_disease", "Whether the patient is healthy or has disease"),
    ("developed_irae", "Whether the patient developed immune-related adverse events (irAE)"),
    ("irae_type", "Type of irAE (e.g., colitis, dermatitis)"),
    ("group_irae", "Categorization based on severity of irAE"),
    ("case_control_irae", "irAE information in case/control format"),
    ("sex", "Sex of the patient"),
    ("age", "Age of the patient"),
    ("hla_haplotype", "HLA haplotype of the patient"),
    ("barcode", "Unique barcode for each cell"),
    ("nCounts_RNA", "Number of RNA reads per cell"),
    ("nFeatures_RNA", "Number of genes detected per cell"),
    ("percent_mito", "Percentage of mitochondrial gene expression per cell"),
    ("tissue", "Tissue of origin (e.g., blood, tumor)"),
    ("timepoint_pre_post", "Indicates whether sample is pre- or post-treatment"),
    ("timepoint_moment", "Exact timepoint of sample collection"),
    ("cell_sorting", "Method of cell sorting (e.g., FACS)"),
    ("original_celltype", "Cell type annotation in original dataset"),
    ("original_annotation", "Any additional annotation for the cell"),
    ("has_gex", "Whether gene expression data is available"),
    ("has_tcr", "Whether TCR data is available"),
    ("has_tcr_matched", "Whether TCR is matched to this cell"),
    ("has_vdj", "Whether V(D)J sequencing data is available"),
    ("has_gex_tcr", "Whether both gene expression and TCR data exist for the cell"),
    ("TRAV", "TCR alpha chain V gene segment"),
    ("TRAJ", "TCR alpha chain J gene segment"),
    ("TRAC", "TCR alpha chain C gene segment"),
    ("cdr3_TRA", "CDR3 sequence of TCR alpha chain"),
    ("TRBV", "TCR beta chain V gene segment"),
    ("TRBD", "TCR beta chain D gene segment"),
    ("TRBJ", "TCR beta chain J gene segment"),
    ("TRBC", "TCR beta chain C gene segment"),
    ("cdr3_TRB", "CDR3 sequence of TCR beta chain"),
    ("cdr3_full", "Full CDR3 sequence, possibly combined α and β chains")
]

# ------------------------------- Create DataFrame and save to Excel
df = pd.DataFrame(columns_info, columns=["Column Name", "Description"])
df.to_excel(excel_file, index=False)



# ------------------------------- Count unique patients and cells
print("Unique patients:")
print(adata.obs['patient_id'].nunique())

print("Cells per patient:")
print(adata.obs['patient_id'].value_counts())

print("Cells per tissue type:")
print(adata.obs['tissue'].value_counts())

print("Cells per timepoint:")
print(adata.obs['timepoint_pre_post'].value_counts())


print("treatment_type:")
print(adata.obs['treatment_type'].value_counts())



# ------------------------------- Extract and save the main gene expression matrix

# Convert to sparse if not already
X_sparse = sp.csr_matrix(adata.X) if not sp.issparse(adata.X) else adata.X

# Create new AnnData with only gene expression and barcodes
adata_gene = ad.AnnData(
    X=X_sparse,
    obs=pd.DataFrame({'barcode': adata.obs['barcode'].values}),
    var=pd.DataFrame(index=adata.var_names)
)

# Define output path
gene_file = os.path.join(output_dir, "GeneExpression_matrix.h5ad")

# Save
adata_gene.write(gene_file)

print(f"Gene expression matrix saved as sparse at: {gene_file}")



# ------------------------------- Extract and save the TCR (T-cell receptor) matrix

# Define the columns related to TCR sequences
tcr_columns = ['TRAV','TRAJ','TRAC','cdr3_TRA','TRBV','TRBD','TRBJ','TRBC','cdr3_TRB','cdr3_full']

# Copy the TCR columns from the AnnData obs (cell metadata)
tcr_df = adata.obs[tcr_columns].copy()

# Add the barcode as the first column
tcr_df.insert(0, 'barcode', adata.obs['barcode'])

# Create a new AnnData object for TCR data
adata_tcr = ad.AnnData(
    X=None,                 # X is empty because these are strings, not numeric expression values
    obs=tcr_df,             # Cell metadata containing TCR information
    var=pd.DataFrame(index=tcr_columns)  # Column metadata
)

# Define the output path for saving the TCR matrix
tcr_file = os.path.join(output_dir, "TCR_matrix.h5ad")

# Save the AnnData object to an .h5ad file
adata_tcr.write(tcr_file)




# ------------------------------- Extract and save all metadata excluding TCR columns

# List of TCR columns to exclude
tcr_columns = ['TRAV','TRAJ','TRAC','cdr3_TRA','TRBV','TRBD','TRBJ','TRBC','cdr3_TRB','cdr3_full']

# Copy the obs (metadata) from the original AnnData
metadata_df = adata.obs.copy()

# Drop TCR-related columns
metadata_df = metadata_df.drop(columns=tcr_columns)

# Move 'barcode' column to the first position
cols = ['barcode'] + [c for c in metadata_df.columns if c != 'barcode']
metadata_df = metadata_df[cols]

# Create a new AnnData object for metadata
adata_meta = ad.AnnData(
    X=None,              # No expression data
    obs=metadata_df,     # Cell metadata without TCR columns
    var=pd.DataFrame(index=[])  # No column metadata
)

# Define the output path for saving the metadata matrix
meta_file = os.path.join(output_dir, "Metadata_matrix.h5ad")

# Save the AnnData object to an .h5ad file
adata_meta.write(meta_file)




# ------------------------------- Extract and save all metadata excluding TCR columns

# -------------------------------------- Select 20 cells that have TCR


gene_file = os.path.join(output_dir, "GeneExpression_matrix.h5ad")
tcr_file = os.path.join(output_dir, "TCR_matrix.h5ad")
meta_file = os.path.join(output_dir, "Metadata_matrix.h5ad")

adata_gene = ad.read_h5ad(gene_file)
adata_tcr = ad.read_h5ad(tcr_file)
adata_meta = ad.read_h5ad(meta_file)


print("Gene barcodes example:", adata_gene.obs['barcode'][:5].tolist())
print("TCR barcodes example:", adata_tcr.obs.index[:5].tolist())
print("Metadata barcodes example:", adata_meta.obs.index[:5].tolist())
"""