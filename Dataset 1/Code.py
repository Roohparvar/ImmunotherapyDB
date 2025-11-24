import scanpy as sc
import anndata as ad
import pandas as pd
import os



# ------------------------------- Load the main .h5ad file
adata = ad.read_h5ad("C:/Users/z5551702/OneDrive - UNSW/Aim 1/Data/Dataset 1/Ali_ClinCanRes_2024_data.h5ad")

# ------------------------------- Print the contents of this file to see what it contains
print(adata)

# ------------------------------- Add prefix to barcodes (convert to string first)
adata.obs['barcode'] = 'Dataset 1/' + adata.obs['barcode'].astype(str)


# ------------------------------- Define columns and their descriptions
output_dir = "C:/Users/z5551702/OneDrive - UNSW/Aim 1/Data/Dataset 1/Processed"
excel_file = os.path.join(output_dir, "ObsColumns_Descriptions.xlsx")

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



# ------------------------------- Extract and save the main gene expression matrix
gene_expr_df = pd.DataFrame(
    adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
    index=adata.obs['barcode'],   # Use cell barcodes as row indices
    columns=adata.var_names       # Use gene names as column headers
)

# Add the barcode as the first column in the DataFrame
gene_expr_df.insert(0, 'barcode', gene_expr_df.index)

# Create a new AnnData object containing only the gene expression matrix
adata_gene = ad.AnnData(
    X=gene_expr_df.iloc[:,1:].values,              # Expression values only (exclude barcode)
    obs=pd.DataFrame({'barcode': gene_expr_df['barcode'].values}),  # Cell metadata with barcode
    var=pd.DataFrame(index=adata.var_names)       # Gene metadata
)

# Define the output path for saving the gene expression matrix
gene_file = os.path.join(output_dir, "GeneExpression_matrix.h5ad")

# Save the AnnData object to an .h5ad file
adata_gene.write(gene_file)



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




# ------------------------------- Load the saved AnnData objects
gene_file = os.path.join(output_dir, "GeneExpression_matrix.h5ad")
tcr_file = os.path.join(output_dir, "TCR_matrix.h5ad")
meta_file = os.path.join(output_dir, "Metadata_matrix.h5ad")

adata_gene = ad.read_h5ad(gene_file)
adata_tcr = ad.read_h5ad(tcr_file)
adata_meta = ad.read_h5ad(meta_file)

# ------------------------------- Convert X to dense if sparse and create a DataFrame for gene expression
if hasattr(adata_gene.X, "toarray"):  
    X_dense = adata_gene.X.toarray()
else:
    X_dense = adata_gene.X

gene_df = pd.DataFrame(X_dense, index=adata_gene.obs['barcode'], columns=adata_gene.var_names)
gene_df = gene_df.iloc[:20, :20]  # First 20 cells and 20 genes
gene_df.insert(0, 'barcode', gene_df.index)

# Take first 20 cells for TCR and Metadata
tcr_df = adata_tcr.obs.iloc[:20, :]
meta_df = adata_meta.obs.iloc[:20, :]

# ------------------------------- Create a DataFrame for matrix dimensions
matrix_info = pd.DataFrame({
    "Matrix": ["GeneExpression", "TCR", "Metadata"],
    "Rows": [adata_gene.X.shape[0], adata_tcr.obs.shape[0], adata_meta.obs.shape[0]],
    "Columns": [adata_gene.X.shape[1], adata_tcr.obs.shape[1], adata_meta.obs.shape[1]]
})

# ------------------------------- Save sample and dimensions to Excel with multiple sheets
excel_path = os.path.join(output_dir, "Sample_20_cells.xlsx")
with pd.ExcelWriter(excel_path) as writer:
    matrix_info.to_excel(writer, sheet_name="Matrix_Dimensions", index=False)  # New sheet with dimensions
    gene_df.to_excel(writer, sheet_name="GeneExpression", index=False)
    tcr_df.to_excel(writer, sheet_name="TCR", index=False)
    meta_df.to_excel(writer, sheet_name="Metadata", index=False)