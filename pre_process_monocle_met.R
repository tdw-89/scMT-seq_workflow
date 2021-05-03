#### NOTES #####################################################################
# This script continues the process of creating a CDS object, started by the 
# 'pre_process_monocle' script. This script is called when methylation data
# is going to be used.

#### Custom Functions ##########################################################
# Custom functions created for this script
separate_by_context <- function(i, df){
  return(df[which(df$met_context == unique(df$met_context)[i]),])
}

missing_genes <- function(i, df){
  return(ens_IDs[which(!(ens_IDs %in% as.character(df[which(df$sample == sample_IDs[i]),2])))])
}

missing_genes_context <- function(i, df){
  return(ensID_context[which(!(ensID_context %in% as.character(df[which(df$sample == sample_IDs_met[i]),6])))])
}

missing_genes_context_list <- function(i, list_samples){
  return(ensID_context[which(!(ensID_context %in% as.character(list_samples[[i]][,6])))])
}

data_frames_na_met_weight_context <- function(i, missing_list){
  return(data.frame("sample" = rep(sample_IDs[i], length(missing_list[[i]])),
                    "ens_ID" = substr(missing_list[[i]], start = 1, stop = 18),
                    "met" = rep(NA, length(missing_list[[i]])),
                    "met_weight" = rep(NA, length(missing_list[[i]])),
                    "met_context" = substr(missing_list[[i]], start = 19, stop = 34),
                    "ens_and_context" = missing_list[[i]]))
}

data_frames_na_met_weight <- function(i, missing_list){
  return(data.frame("sample" = rep(sample_IDs[i], length(missing_list[[i]])),
                    "ens_ID" = missing_list[[i]],
                    "met" = rep(NA, length(missing_list[[i]])),
                    "met_weight" = rep(NA, length(missing_list[[i]]))))
}

data_frames_na_met_context <- function(i, missing_list){
  return(data.frame("sample" = rep(sample_IDs[i], length(missing_list[[i]])),
                    "ens_ID" = missing_list[[i]],
                    "met" = rep(NA, length(missing_list[[i]])),
                    "met_context" = rep(NA, length(missing_list[[i]]))))
}

data_frames_na_met <- function(i, missing_list){
  return(data.frame("sample" = rep(sample_IDs[i], length(missing_list[[i]])),
                    "ens_ID" = missing_list[[i]],
                    "met" = rep(NA, length(missing_list[[i]]))))
}

one_of_each <- function(i, df){
  return(df[which(df$sample == sample_IDs[i]),]
         [-which(duplicated(df[which(df$sample == sample_IDs[i]),2])),]
         [order(df[which(df$sample == sample_IDs[i]),][-which(duplicated(df[which(df$sample == sample_IDs[i]),2])),][,2]),3])
}

combine_lists <- function(i, list_1, list_2){
  return(bind_rows(list_1[[i]], list_2[[i]]))
}

merge_id_sample_context <- function(i, list){
  return(cbind(list[[i]], paste0(list[[i]]$sample, list[[i]]$ens_ID, list[[i]]$met_context)))
}

merge_id_sample_context_2 <- function(df){
  return(cbind(df, paste0(df$sample, df$ens_ID, df$met_context)))
}

remove_sixth <- function(i, list){
  return(list[[i]][,-6])
}

find_duplicates <- function(i, list){
  return(which(!isUnique(list[[i]][,6])))
}

return_duplicates <- function(i, list_1, list_2){
  return(list_2[[i]][list_1[[i]],])
}

find_means_weight <- function(i, list){
  met_means <- c()
  weight_means <- c()
  met_means <- sapply(X = 1:length(unique(list[[i]][,6])),
                      FUN = calc_mean,
                      df = list[[i]][,c(3, 6)],
                      unique_items = unique(list[[i]][,6]))
  weight_means <- sapply(X = 1:length(unique(list[[i]][,6])),
                         FUN = calc_mean,
                         df = list[[i]][,c(4, 6)],
                         unique_items = unique(list[[i]][,6]))
  return(cbind(met_means, weight_means))
}

aggregate_list_mean <- function(df, col_num){
  return(aggregate(x = df[,col_num], by = list(df$combined), FUN = mean))
}

split_cols_context <- function(df, ls){
  new_df <- cbind(
  substr(df[,1],
               start = 1 ,
               stop = ls),
  substr(df[,1],
         start = (ls + 1) ,
         stop = max(nchar(df[,1]))),
  df[,2])
  return(as.data.frame(new_df))
}

create_matrix_weight_context <- function(gdf){
  matrix_new <- lapply(X = 1:length(unique(gdf$V1)), fill_matrix_weight_context, gdf)
  return(matrix_new)
}

fill_matrix_weight_context <- function(i, gdf){
  gdf[group_rows(gdf)[[i]],][order(gdf[group_rows(gdf)[[i]],2]),3]
}

test_length_met_context <- function(gdf){
  lengths <- sapply(X = group_rows(gdf), FUN = length) 
  return(lengths)
}

add_labels <- function (i, mat_list, list_labels, sample_ids){
  colnames(mat_list[[i]]) <- sample_ids
  row.names(mat_list[[i]]) <- list_labels[[i]]
  return(mat_list[[i]])
}

list_labels <- function(df){
  return(unique(df$V2)[order(unique((df$V2)))])
}

apply_weights_met <- function(i, met_list, weight_list){
  weighted_matrices <- sapply(X = 1:nrow(met_list[[i]]), FUN = function(j, matrix, weights){
    return(as.numeric(matrix[j,]) * ((as.numeric(weights[j,]))/sum(as.numeric(weights[j,])) * ncol(matrix)))
    }, met_list[[i]], weight_list[[i]])
  return(weighted_matrices)
}

##### Pre-process ##############################################################
if(!is_context && !is_weighted){
  #### Create Methylation Matrix ###############################################
  length_ensID <- nchar(as.character(data_pre_proc_trimmed$ens_ID[1]))
  length_sample <- nchar(as.character(data_pre_proc_trimmed$sample[1]))
  data_pre_proc_trimmed_met <- data_pre_proc_trimmed[,c(1,2,4)]
  # Incomplete.....---REMOVE ONCE COMPLETE 
}else if(!is_context && is_weighted){
  #### Create Methylation Matrix ###############################################
  length_ensID <- nchar(as.character(data_pre_proc_trimmed$ens_ID[1]))
  length_sample <- nchar(as.character(data_pre_proc_trimmed$sample[1]))
  data_pre_proc_trimmed_met_weight <- data_pre_proc_trimmed[,c(1,2,4,5)]
  # Incomplete.....---REMOVE ONCE COMPLETE 
}else if(is_context && !is_weighted){
  #### Create Methylation Matrix ###############################################
  length_ensID <- nchar(as.character(data_pre_proc_trimmed$ens_ID[1]))
  length_sample <- nchar(as.character(data_pre_proc_trimmed$sample[1]))
  data_pre_proc_trimmed_context <- data_pre_proc_trimmed[,c(1,2,4,6)]
  # Incomplete.....---REMOVE ONCE COMPLETE                                                       
}else if(is_context && is_weighted){
  #### Create Methylation Matrices ###############################################
  length_ensID <- nchar(as.character(data_pre_proc_trimmed$ens_ID[1]))
  length_sample <- nchar(as.character(data_pre_proc_trimmed$sample[1]))
  data_pre_proc_trimmed_met_weight_context <- data_pre_proc_trimmed[,-3]
  ens_IDs_met <- as.character(unique(data_pre_proc_trimmed_met_weight_context$ens_ID))
  sample_IDs_met <- unique(data_pre_proc_trimmed_met_weight_context$sample)
  context_names <- unique(data_pre_proc_trimmed_met_weight_context$met_context)
  names_id_context_paste <- unique(paste0(data_pre_proc_trimmed_met_weight_context$ens_ID,
                                          data_pre_proc_trimmed_met_weight_context$met_context))
  
  if(is_complete != TRUE)
  {
    # Create an additional column that contains the ensembl ID for each data point
    # concatenated with its methylation context.
    data_pre_proc_trimmed_met_weight_context_zero <- data_pre_proc_trimmed_met_weight_context
    data_pre_proc_trimmed_met_weight_context_zero[which(is.na(data_pre_proc_trimmed_met_weight_context_zero$met)),3] <- 0
    data_pre_proc_trimmed_met_weight_context_zero[which(is.na(data_pre_proc_trimmed_met_weight_context_zero$met)),4] <- 0
    
    data_pre_proc_id_context_paste <- cbind(data_pre_proc_trimmed_met_weight_context_zero,
                              paste0(data_pre_proc_trimmed_met_weight_context_zero$ens_ID,
                                     data_pre_proc_trimmed_met_weight_context_zero$met_context))
    colnames(data_pre_proc_id_context_paste)[6] <- "ens_and_context"
    rm(data_pre_proc_trimmed_met_weight_context_zero)
    
    ensID_context <- unique(data_pre_proc_id_context_paste$ens_and_context)
    
    list_sample_context <- lapply(1:length(unique(data_pre_proc_id_context_paste$sample)), 
         function(i) 
         data_pre_proc_id_context_paste[which(data_pre_proc_id_context_paste$sample == data_pre_proc_id_context_paste$sample[i]),])
    names(list_sample_context) <- unique(data_pre_proc_id_context_paste$sample)
    
    #Find missing gene/methylation context combinations in the list of samples
    list_missing_context_genes <- lapply(1:length(list_sample_context),
                                         FUN = missing_genes_context_list,
                                         list_samples = list_sample_context)
    list_missing_context_genes_dfs <- lapply(1:length(list_missing_context_genes),
                                         FUN = data_frames_na_met_weight_context,
                                         missing_list = list_missing_context_genes)
    rm(list_missing_context_genes)
    list_sample_context_with_missing <- lapply(1:length(list_sample_context),
                                         FUN = combine_lists,
                                         list_1 = list_sample_context,
                                         list_2 = list_missing_context_genes_dfs)
    list_context_weight <- lapply(1:length(list_sample_context_with_missing),
                                               FUN = remove_sixth, 
                                               list_sample_context_with_missing) # untested --- REMOVE ONCE COMPLETE
  }else{
    # convert NaNs and NAs to 0
    data_pre_proc_trimmed_met_weight_context_zero <- data_pre_proc_trimmed_met_weight_context
    data_pre_proc_trimmed_met_weight_context_zero[which(is.na(data_pre_proc_trimmed_met_weight_context_zero$met)),3] <- 0
    data_pre_proc_trimmed_met_weight_context_zero[which(is.na(data_pre_proc_trimmed_met_weight_context_zero$met)),4] <- 0
    
    # create a list with data separated by methylation context
    list_context_weight <- lapply(1:length(unique(data_pre_proc_trimmed_met_weight_context_zero$met_context)),
                                  FUN = separate_by_context, 
                                  df = data_pre_proc_trimmed_met_weight_context_zero)
    names(list_context_weight) <- unique(data_pre_proc_trimmed_met_weight_context_zero$met_context)
  }
  if(duplicates == TRUE)
  {
    num_cores <- detectCores()
    cluster_1 <- makeCluster(num_cores)
    
    list_context_weight_paste <- parLapply(cl = cluster_1, 
                                             X = list_context_weight, 
                                             fun = merge_id_sample_context_2)
    list_context_weight_paste <- parLapply(cl = cluster_1, 
              X = list_context_weight_paste, 
              fun = `names<-`, ... = append(colnames(list_context_weight_paste[[1]][-6]),
                                        values = "combined"))
    
    list_means_met <- parLapply(cl = cluster_1,
                            X = list_context_weight_paste,
                            fun = aggregate_list_mean, col_num = 3)
    list_means_weights <- parLapply(cl = cluster_1,
                                X = list_context_weight_paste,
                                fun = aggregate_list_mean, col_num = 4)
    list_means_met_split <- parLapply(cl = cluster_1,
                                     X = list_means_met, 
                                     fun = split_cols_context,
                                     ls = length_sample)
    list_means_weights_split <- parLapply(cl = cluster_1,
                                      X = list_means_weights, 
                                      fun = split_cols_context,
                                      ls = length_sample)
    list_means_met_split <- lapply(X = list_means_met_split, FUN = group_by, V1)
    list_means_weights_split <- lapply(X = list_means_weights_split, FUN = group_by, V1)
  
    #Check if all samples have all gene/contexts for each context
    test_results <- unlist(lapply(list_means_met_split,
                                  FUN = test_length_met_context))
    
    if((length(unique(test_results))) != length(list_means_met_split)){
      cat(sprintf("\nERROR: Not every sample has a data point for every\n
                   gene/methylation context."))
      stop()
    }
    
    test_results_2 <- unlist(lapply(list_means_weights_split,
                                  FUN = test_length_met_context))
    
    if((length(unique(test_results_2))) != length(list_means_weights_split)){
      cat(sprintf("\nERROR: Not every sample has a weight for every\n
                   gene/methylation context."))
      stop()
    }
    
    labels_list_met <- parLapply(cl = cluster_1,
                                 X = list_means_met_split,
                                 fun = list_labels)
    labels_list_weights <- parLapply(cl = cluster_1,
                                     X = list_means_weights_split,
                                     fun = list_labels)
    
    met_matrix_list <- lapply(X = list_means_met_split,
                              FUN = create_matrix_weight_context)# parLapply does not seem to work
    met_matrix_weights_list <- lapply(X = list_means_weights_split,
                                      FUN = create_matrix_weight_context)# parLapply does not seem to work
    
    met_matrix <- parLapply(cl = cluster_1,
                            X = met_matrix_list,
                            fun = bind_cols) 
    met_weights_matrix <- parLapply(cl = cluster_1,
                                    X = met_matrix_weights_list,
                                    fun = bind_cols) 
    
    names(met_matrix) <- names(list_means_met_split)
    names(met_weights_matrix) <- names(list_means_weights_split)
    
    met_matrix <- parLapply(cl = cluster_1,
                      X = 1:length(met_matrix),
                      fun = add_labels,
                      mat_list = met_matrix,
                      list_labels = labels_list_met,
                      sample_ids = sample_IDs_met)
    met_weights_matrix <- parLapply(cl = cluster_1,
                            X = 1:length(met_weights_matrix),
                            fun = add_labels,
                            mat_list = met_weights_matrix,
                            list_labels = labels_list_weights,
                            sample_ids = sample_IDs_met)
    met_matrix_context_weighted <- parLapply(cl = cluster_1, 
                                             X = 1:length(met_weights_matrix), 
                                             fun = apply_weights_met,
                                             met_list = met_matrix, 
                                             weight_list = met_weights_matrix)
    #Each matrix needs to be transposed in order to read correctly
    met_matrix_context_weighted_transpose <- parLapply(cl = cluster_1, 
                                                  X = met_matrix_context_weighted, 
                                                  fun = t) 
    
    met_matrix_context_weighted_transpose_labeled <- parLapply(cl = cluster_1,
                            X = 1:length(met_matrix_context_weighted_transpose),
                            fun = add_labels,
                            mat_list = met_matrix_context_weighted_transpose,
                            list_labels = labels_list_met,
                            sample_ids = sample_IDs_met)
    names(met_matrix_context_weighted_transpose_labeled) <- names(labels_list_met)
    list_matrix_context <- met_matrix_context_weighted_transpose_labeled
    stopCluster(cluster_1)
  }else{
    #DO STUFF --- REMOVE ONCE COMLETED
  }
  # Convert the list of methylation matrices to a list of sparse matrices, then
  # to a list of monocle3 CDS objects.
  
  list_matrix_context_sparse <- lapply(X = list_matrix_context, 
            FUN = function(x){Matrix(data = x, sparse = TRUE)})
  
  
  list_context_CDS <- lapply(X = list_matrix_context_sparse, 
                             FUN = new_cell_data_set,
                             cell_metadata = NULL,
                             gene_metadata = NULL) # Is there a way to suppress the warnings? --- REMOVE ONCE FINISHED
  saveRDS(object = list_context_CDS, file = "cds_met_weight_context_list.rds")
  
  if(option == 3){
    met_context_analysis <- TRUE
    source(file = "pseudotime_analysis.R", local = FALSE)
  }
}
