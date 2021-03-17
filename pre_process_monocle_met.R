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
  return(cbind(unique(list[,6]), met_means, weight_means))
}

calc_mean <- function(i, df, unique_items){
  return(sum(df[(df[,2] == unique_items[i]),1]) / sum(df[,2] == unique_items[i]))
}

##### Pre-process ##############################################################
if(!is_context && !is_weighted){
  #### Create Methylation Matrix ###############################################
  data_pre_proc_trimmed_met <- data_pre_proc_trimmed[,c(1,2,4)]
  # Incomplete.....---REMOVE ONCE COMPLETE 
}else if(!is_context && is_weighted){
  #### Create Methylation Matrix ###############################################
  data_pre_proc_trimmed_met_weight <- data_pre_proc_trimmed[,c(1,2,4,5)]
  # Incomplete.....---REMOVE ONCE COMPLETE 
}else if(is_context && !is_weighted){
  #### Create Methylation Matrix ###############################################
  data_pre_proc_trimmed_context <- data_pre_proc_trimmed[,c(1,2,4,6)]
  # Incomplete.....---REMOVE ONCE COMPLETE                                                       
}else if(is_context && is_weighted){
  #### Create Methylation Matrices ###############################################
  data_pre_proc_trimmed_met_weight_context <- data_pre_proc_trimmed[,-3]
  ens_IDs_met <- as.character(unique(data_pre_proc_trimmed_met_weight_context$ens_ID))
  sample_IDs_met <- unique(data_pre_proc_trimmed_met_weight_context$sample)
  context_names <- unique(data_pre_proc_trimmed_met_weight_context$met_context)
  if(is_complete != TRUE)
  {
    # Create an additional column that contains the ensembl ID for each data point
    # concatenated with its methylation context.
    
    data_pre_proc_id_context_paste <- cbind(data_pre_proc_trimmed_met_weight_context,
                              paste0(data_pre_proc_trimmed_met_weight_context$ens_ID,
                              data_pre_proc_trimmed_met_weight_context$met_context))
    colnames(data_pre_proc_id_context_paste)[6] <- "ens_and_context"
    rm(data_pre_proc_trimmed_met_weight_context)
    
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
    
    # create a list with data separated by methylation context
    list_context_weight <- lapply(1:length(unique(data_pre_proc_trimmed_met_weight_context_zero$met_context)),
                                  FUN = separate_by_context, 
                                  df = data_pre_proc_trimmed_met_weight_context_zero)
    names(list_context_weight) <- unique(data_pre_proc_trimmed_met_weight_context_zero$met_context)
  }
  if(duplicates == TRUE)
  {
    list_context_weight_paste <- lapply(1:length(list_context_weight),
                                        FUN = merge_id_sample_context,
                                        list_context_weight)
    list_duplicates_index <- lapply(1:length(list_context_weight_paste),
                          FUN = find_duplicates,
                          list_context_weight_paste)
    list_duplicates <- lapply(1:length(list_duplicates_index),
                              FUN = return_duplicates,
                              list_1 = list_duplicates_index,
                              list_2 = list_context_weight_paste)
    list_means <- lapply(1:length(list_duplicates),
                         FUN = find_means_weight,
                         list_duplicates)
    
    test <- sapply(X = 1:length(unique(list_duplicates[[1]][,6])),
           FUN = calc_mean,
           df = list_duplicates[[1]][,c(4,6)],
           unique_items = unique(list_duplicates[[1]][,6]))
    
  }
}
