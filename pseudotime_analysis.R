#### Author Info: ##############################################################
# Name: Tom Wolfe
# Institution: University of Massachusetts Lowell (student)

#### NOTES #####################################################################
# 
# This is an R script used to perform clustering and pseudotime analysis.
# 
# #Input/Pre-condition: A cell data set object (CDS) containing exprssion or 
# methylation data.
#   
# 
# #Output/post-condition: Multiple plots showing cell clusters and trajectory.

#setwd("C:/Users/eccwo/Dropbox/UML/Spring 2021/Honors Project/R-scripts")# Used for testing --- REMOVE ONCE COMPLETE
#setwd("/home/tom/Dropbox/UML/Spring 2021/Honors Project/R-scripts")# Used for testing --- REMOVE ONCE COMPLETE

# List files in current directory, then read in cds objects.
if(option == 1){
  
  file_list <- list.files()
  
  if("cds_expr.rds" %in% file_list){
    cds_object_expr <- readRDS(file = grep(pattern = "cds_expr.rds", x = file_list, value = TRUE))
    expr_analysis <- TRUE
  }
  
  if("cds_met.rds" %in% file_list){
    cds_object_met <- readRDS(file = grep(pattern = "cds_met.rds", x = file_list, value = TRUE))
  }
  
  if("cds_met_weight_context_list.rds" %in% file_list){
    library(parallel)
    cds_object_met_weight_context_list <- readRDS(file = grep(pattern = "cds_met_weight_context_list.rds", x = file_list, value = TRUE))
    met_context_analysis <- TRUE
  }
}else{
  if(expr_analysis){
    cds_object_expr <- cds_expr
  }
  if(met_analysis){
    cds_object_met <- cds_met
  }
  if(met_context_analysis){
    cds_object_met_weight_context_list <- list_context_CDS # Should it be the same for weighted and un-weighted ? --- REMOVE ONCE COMPLETE
  }
}

if(expr_analysis){
  # NOTE: need to give option to set number of dimensions for PCA after seeing scree plot --- REMOVE ONCE COMPLETE
  
  # PCA analysis:
  cds_object_expr <- preprocess_cds(cds_object_expr, norm_method = "none", method = "PCA",
                 num_dim = 10,  verbose = TRUE)
  
  # Scree plot:
  print(plot_pc_variance_explained(cds_object_expr))
  
  # Reduce the dimensions (using UMAP by default):
  cds_object_expr <- reduce_dimension(cds_object_expr, preprocess_method = "PCA")
  
  #Cluster the cells and plot the clusters:
  cds_object_expr <- cluster_cells(cds_object_expr)
  print(plot_cells(cds_object_expr, 
                   group_label_size = 2,
                   cell_size = 1,
                   color_cells_by = "cell"))
  
  # Pseudotime plot:
  cds_object_expr <- learn_graph(cds_object_expr)
  print(plot_cells(cds_object_expr, group_label_size = 6, cell_size = 1))
}

if(met_analysis){
  # do stuff --- REMOVE ONCE COMPLETE
}

if(met_context_analysis){
  num_cores <- detectCores()
  cluster_1 <- makeCluster(num_cores)
  
  # NOTE: need to give option to set number of dimensions for PCA after seeing scree plot --- REMOVE ONCE COMPLETE
  
  # PCA analysis:
  cds_object_met_weight_context_list <- parLapply(cl = cluster_1, 
                                                  X = cds_object_met_weight_context_list, 
                                                  fun = preprocess_cds,
                                                  norm_method = "none", 
                                                  method = "PCA",
                                                  num_dim = 10,  
                                                  verbose = TRUE)
  
  # Print a scree plot for each context
  for(i in 1:length(cds_object_met_weight_context_list)){
    print(plot_pc_variance_explained(cds_object_met_weight_context_list[[i]]) +
            labs(title = names(cds_object_met_weight_context_list)[i]))
  }
  
  # Reduce the dimensions of the context matrices (using UMAP)
  cds_object_met_weight_context_list <- parLapply(cl = cluster_1, 
            X = cds_object_met_weight_context_list, 
            fun = reduce_dimension, 
            preprocess_method = "PCA")
  
  cds_object_met_weight_context_list <- parLapply(cl = cluster_1,
                                                  X = cds_object_met_weight_context_list, 
                                                  fun = cluster_cells)
  
  for(i in 1:length(cds_object_met_weight_context_list)){
    print(plot_cells(cds_object_met_weight_context_list[[i]],
                     group_label_size = 2,
                     cell_size = 1,
                     color_cells_by = "cell") +
            labs(title = names(cds_object_met_weight_context_list)[i]))
  }
  
  stopCluster(cluster_1)
}

