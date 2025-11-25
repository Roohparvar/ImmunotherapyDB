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






# تعداد بیماران منحصربه‌فرد
print("Unique patients:", adata.obs['patient_id'].nunique())

# تعداد سلول‌ها یا نمونه‌ها برای هر بیمار
print("Cells per patient:")
print(adata.obs['patient_id'].value_counts())

# تعداد سلول‌ها در هر نوع نمونه (Tissue)
print("Cells per tissue type:")
print(adata.obs['tissue'].value_counts())

# تعداد سلول‌ها بر اساس زمان نمونه‌گیری (pre/post treatment)
print("Cells per timepoint:")
print(adata.obs['timepoint_pre_post'].value_counts())


