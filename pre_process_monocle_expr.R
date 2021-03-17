#### NOTES #####################################################################
# This script continues the process of creating a CDS object, started by the 
# 'pre_process_monocle' script. This script is called when expression data
# is going to be used.

#### Custom Functions ##########################################################
# Custom functions created for this script
missing_genes <- function(i, df){
  return(ens_IDs[which(!(ens_IDs %in% as.character(df[which(df$sample == sample_IDs[i]),2])))])
}
data_frames_na <- function(i, missing_list){
  return(data.frame("sample" = rep(sample_IDs[i], length(missing_list[[i]])),
                    "ens_ID" = missing_list[[i]],
                    "expr" = rep(NA, length(missing_list[[i]]))))
}
one_of_each <- function(i, df){
  return(df[which(df$sample == sample_IDs[i]),]
         [-which(duplicated(df[which(df$sample == sample_IDs[i]),2])),]
         [order(df[which(df$sample == sample_IDs[i]),][-which(duplicated(df[which(df$sample == sample_IDs[i]),2])),][,2]),3])
}


##### Create Expression Matrix #################################################
data_pre_proc_trimmed_expr <- data_pre_proc_trimmed[,1:3]
ens_IDs <- as.character(unique(data_pre_proc_trimmed_expr$ens_ID))
sample_IDs <- unique(data_pre_proc_trimmed_expr$sample)


list_missing_ens_IDs <- lapply(X = 1:length(sample_IDs),
                               FUN = missing_genes,
                               df = data_pre_proc_trimmed_expr)
list_missing_ens_IDs_dfs <- lapply(X = 1:length(list_missing_ens_IDs),
                                   FUN = data_frames_na,
                                   missing_list = list_missing_ens_IDs)
rm(list_missing_ens_IDs)
data_and_missing_ens_IDs <- rbind(data_pre_proc_trimmed_expr,
                                  bind_rows(list_missing_ens_IDs_dfs))
rm(data_pre_proc_trimmed_expr, list_missing_ens_IDs_dfs)

list_all_genes <- lapply(X = 1:length(sample_IDs), 
                         FUN = one_of_each, data_and_missing_ens_IDs)

# Convert to a matrix, give the columns and rows names, and change all NA's
# to 0.
expr_matrix <- as.matrix(bind_cols(list_all_genes))
expr_matrix[which(is.na(expr_matrix))] <- 0
row.names(expr_matrix) <- sort(ens_IDs)
colnames(expr_matrix) <- sample_IDs
rm(list_all_genes)
##### Convert Expression Matrix ################################################
# 
expr_matrix <- Matrix(data = expr_matrix, sparse = TRUE)
cds_object <- monocle3::new_cell_data_set(expression_data = expr_matrix,
                                          cell_metadata = NULL, 
                                          gene_metadata = NULL)

saveRDS(object = cds_object, file = "cds_expr.rds")

if(do_analysis){
  source(file = "pseudotime_analysis.R", local = FALSE)
}
  