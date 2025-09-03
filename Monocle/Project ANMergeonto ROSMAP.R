expression_matrix_tmp <- read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/ADNI/Processed Data/ROSMAP_blood_gene_expression_monocyte_ACTL_04172023_ExpMatrix.csv", header=FALSE)
cell_metadata <- read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/ADNI/Processed Data/ROSMAP_blood_gene_expression_monocyte_ACTL_04172023_SampleInfo.csv")
gene_annotation <- read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/ADNI/Processed Data/ROSMAP_blood_gene_expression_monocyte_ACTL_04172023_GeneInfo.csv")

expression_matrix = as.matrix(expression_matrix_tmp)
row.names(expression_matrix) <- gene_annotation$X
row.names(gene_annotation) <- row.names(expression_matrix)
colnames(expression_matrix) <- row.names(cell_metadata)
row.names(cell_metadata) <- colnames(expression_matrix)

# create cds for monocle
cds_ref <- new_cell_data_set(expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_annotation)


expression_matrix_tmp <- read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/ANMerge/Data/ANMerge_blood_rna_gene_expr_processed_XH_10112023_ExpMatrix.csv", header=FALSE)
cell_metadata <- read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/ANMerge/Data/ANMerge_blood_rna_gene_expr_processed_XH_10112023_SampleInfo.csv")
gene_annotation <- read.csv("/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/ANMerge/Data/ANMerge_blood_rna_gene_expr_processed_XH_10112023_GeneInfo.csv")

expression_matrix = t(as.matrix(expression_matrix_tmp))
row.names(expression_matrix) <- gene_annotation$X
row.names(gene_annotation) <- row.names(expression_matrix)
colnames(expression_matrix) <- row.names(cell_metadata)
row.names(cell_metadata) <- colnames(expression_matrix)

# create cds for monocle
cds_qry <- new_cell_data_set(expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_annotation)


# Remove not shared genes
# Genes in reference.
genes_ref <- row.names(cds_ref)

# Genes in query.
genes_qry <- row.names(cds_qry)

# Shared genes.
genes_shared <- intersect(genes_ref, genes_qry)

# Remove non-shared genes.
cds_ref <- cds_ref[genes_shared,]
cds_qry <- cds_qry[genes_shared,]

# Process reference data
cds_ref <- preprocess_cds(cds_ref, num_dim=100)
cds_ref <- reduce_dimension(cds_ref, build_nn_index=TRUE)
# Save the PCA and UMAP transform models for use with projection.
save_transform_models(cds_ref, 'cds_ref_test_models')

# Project the query data set into the reference space
# Load the reference transform models into the query cds.
cds_qry <- load_transform_models(cds_qry, 'cds_ref_test_models')
# Apply the reference transform models to the query cds.
cds_qry <- preprocess_transform(cds_qry)
cds_qry <- reduce_dimension_transform(cds_qry)

# Plot the combined data sets
plot_cells(cds_ref)
plot_cells(cds_qry)


# Plot the combined data sets
# Label the data sets.
colData(cds_ref)[['data_set']] <- 'reference'
colData(cds_qry)[['data_set']] <- 'query'
# Combine the reference and query data sets.
cds_combined <- combine_cds(list(cds_ref, cds_qry),  keep_all_genes=TRUE, cell_names_unique=FALSE, keep_reduced_dims=TRUE)
plot_cells(cds_combined, color_cells_by='data_set')

#Transfer the reference cell labels to the query data set
cds_qry_lab_xfr <- transfer_cell_labels(cds_qry, reduction_method='UMAP', ref_coldata=colData(cds_ref), ref_column_name='braaksc', query_column_name='braak_xfr', transform_models_dir='cds_ref_test_models')
cds_qry_lab_fix <- fix_missing_cell_labels(cds_qry_lab_xfr, reduction_method='UMAP', from_column_name='braak_xfr', to_column_name='braak_fix')

DF <- as.data.frame(colData(cds_qry_lab_fix))

write.csv(DF, file="/Users/xiaoqing/Library/CloudStorage/OneDrive-IndianaUniversity/R.Projects/ADNI/Processed Data/Monocle_LabelTransfer_ROSMAP_ANMerge_Braak.csv")
