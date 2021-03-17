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

setwd("C:/Users/eccwo/Dropbox/UML/Spring 2021/Honors Project/R-scripts")# Used for testing --- REMOVE ONCE COMPLETE
#setwd("/home/tom/Dropbox/UML/Spring 2021/Honors Project/R-scripts")# Used for testing --- REMOVE ONCE COMPLETE

# List files in current directory, then read in cds objects.
file_list <- list.files()
analyze_expr <- FALSE
analyze_met <- FALSE
analyze_met_weight <- FALSE
analyze_met_context <- FALSE
analyze_met_weight_context <- FALSE

if(sum(grepl(pattern = "cds_expr.rds", x = file_list)) == 1){
  cds_object_expr <- readRDS(file = grep(pattern = "cds_expr.rds", x = file_list, value = TRUE))
  analyze_expr <- TRUE
}

if(sum(grepl(pattern = "cds_met.rds", x = file_list)) == 1){
  cds_object_met <- readRDS(file = grep(pattern = "cds_met.rds", x = file_list, value = TRUE))
  # Check here if methylation values are weighted, and given context --- REMOVE ONCE COMPLETE
}

if(analyze_expr){
  # Give option to set number of dimensions for PCA after seeing scree plot --- REMOVE ONCE COMPLETE
  cds_object_expr <- preprocess_cds(cds_object_expr, norm_method = "none", method = "PCA",
                 num_dim = 10,  verbose = TRUE)
  
  # Scree plot
  print(plot_pc_variance_explained(cds_object_expr))
  
  # Reduce the dimensions (using UMAP by default)
  cds_object_expr <- reduce_dimension(cds_object_expr, preprocess_method = "PCA")
  
  #Cluster the cells and plot the clusters:
  cds_object_expr <- cluster_cells(cds_object_expr)
  print(plot_cells(cds_object_expr))
  
  # Pseudotime plot:
  cds_object_expr <- learn_graph(cds_object_expr)
  print(plot_cells(cds_object_expr))
}
if(analyze_met){
  # do stuff --- REMOVE ONCE COMPLETE
}
if(analyze_met_weight){
  # do stuff --- REMOVE ONCE COMPLETE
}
if(analyze_met_context){
  # do stuff --- REMOVE ONCE COMPLETE
}
if(analyze_met_weight_context){
  # do stuff --- REMOVE ONCE COMPLETE
}
