setwd("C:/Users/eccwo/Dropbox/UML/Spring 2021/Honors Project/R-scripts")# Used for testing --- REMOVE ONCE COMPLETE
#setwd("/home/tom/Dropbox/UML/Spring 2021/Honors Project/R-scripts")# Used for testing --- REMOVE ONCE COMPLETE
#### Author Info: ##############################################################
# Name: Tom Wolfe
# Institution: University of Massachusetts Lowell (student)

#### NOTES #####################################################################
# 
# This is an R script used to generate a cell data set (CDS) object suitable for
# pseudotime and clustering analysis in monocle3. Put this script and whatever
# single-cell data set file to be processed in an otherwise empty folder.
#
# ###
# #Input/Pre-condition: A scRNA-seq/scM&T-seq based table that contains a column
# with ensembl gene ID's, sample/cell ID's, and (a) column(s) with at least
# one of the following:
#   
# For Methylation Data:  
# 1: Methylation rate + methylation context + methylation weights 
# (3 columns total)
# -OR-
# 2: Weighted methylation rates + methylation context (2 columns total)
# -OR-
# 3: Weighted methylation rates normalized across contexts (1 column total)
# 
# For Expression Data:
# 1: Normalized expression rate as TPM, FPKM, or UMI (1 column total)
# -OR-
# 2: Raw read counts + respective gene lengths (2 columns total)
# 
# Also cannot have 'cds_expr.rds' or 'cds_met.rds' in the directory at the start
# of this script
#
# ###
# #Output/post-condition: A single-cell data object or multiple single cell
# data objects.

expr_analysis <- FALSE
met_context_analysis <- FALSE
met_analysis <- FALSE
option <- 0
cat(sprintf("\n
             For pseudotime analysis enter: 1\n
             For data pre-processing enter: 2\n
             For both enter: 3\n"))
option <- readline(prompt = ">>:")
while((option != 1) && (option != 2) && (option != 3)){
  cat(sprintf("Enter 1, 2, or 3\n"))
  option <- readline(prompt = ">>:")
}

if(option == 1){
  packages_installed <- installed.packages()
  if(!("monocle3" %in% installed.packages())){
    cat(sprintf("ERROR: Monocle 3 must be installed before\n"))
    cat(sprintf("psuedotime analysis can be performed.\n"))
    cat(sprintf("Aborting script.\n"))
    stop()
  }
  
  library(monocle3)
  library(ggplot2)
  source(file = "pseudotime_analysis.R", local = FALSE)
}else{
  packages_installed <- installed.packages()
  required_packages <- c("monocle3", "Matrix", "dplyr", "parallel")
  if(sum(required_packages %in% packages_installed) != 4){
    missing_packages <- required_packages[which(!(required_packages %in% packages_installed))]
    for(i in 1:length(missing_packages)){
      cat(sprintf("ERROR: the package '%s' is missing..\n", missing_packages[i]))
    }
    
    cat(sprintf("You must install the missing packages before pre-processing.\n"))
    if(!("monocle3" %in% packages_installed)){
      cat(sprintf("ERROR: Monocle 3 must be installed before\n"))
      cat(sprintf("psuedotime analysis can be performed.\n"))
    }
    stop()
  }
  
  ##### Load Libraries ###########################################################
  library(monocle3, quietly = TRUE, verbose = FALSE)
  library(Matrix, quietly = TRUE, verbose = FALSE)
  library(dplyr, quietly = TRUE, verbose = FALSE)
  library(parallel, quietly = TRUE, verbose = FALSE)
  ##### Check Directory/Files ################################################
  
  # Create a file list, and check to make sure there are the correct files in the 
  # current directory.
  file_list <- list.files()
  if((sum(file_list == "monocle_script.R") +
      sum(file_list == "pre_process_monocle_expr.R") +
      sum(file_list == "pseudotime_analysis.R") + 
      sum(file_list == "pre_process_monocle_met.R")) != 4){ # Not sure if markdown will be used for this --- REMOVE ONCE COMPLETE
    cat(sprintf("Wrong files in the working directory!\n"))
    cat(sprintf("Make sure that only this script or markdown file,\n")) # Not sure if markdown will be used for this --- REMOVE ONCE COMPLETE
    cat(sprintf("along with the single data file to be processed are\n"))
    cat(sprintf("in the current working directory.\n"))
    stop()
  }
  
  
  
  ##### Read in the File #######################################################
  # Figure out what the file extension is, then read it into a data frame.
  if(sum(grepl(pattern = ".rds", x = file_list)) == 1) {
    data_pre_proc <- readRDS(file = file_list[which(grepl(pattern = ".rds",
                                                          x = file_list))])
  }else if (sum(grepl(pattern = ".xlsx", x = file_list)) == 1){
    library(readxl)
    data_pre_proc <- read_xlsx(path = file_list[which(grepl(pattern = ".xlsx",
                                                            x = file_list))])
    detach("package:readxl")
  }else if (sum(grepl(pattern = ".csv", x = file_list)) == 1){
    data_pre_proc <- read.csv(file = file_list[which(grepl(pattern = ".csv",
                                                           x = file_list))])
  }else if (sum(grepl(pattern = ".txt", x = file_list)) == 1){
    data_pre_proc <- read.delim(file = file_list[which(grepl(pattern = ".txt",
                                                             x = file_list))])
  }else{
    cat(sprintf("Wrong files in working directory!\n"))
    stop()
  }
  
  ##### Get Info From User #######################################################
  # Set everything to false initially
  col_name_sample <- 0
  col_name_ensID <- 0
  col_name_expr <- 0
  col_name_met <- 0
  col_name_context <- 0
  col_name_met_weight <- 0
  norm_method <- 0
  is_weighted <- FALSE
  is_met <- FALSE
  is_expr <- FALSE
  is_context <- FALSE
  is_complete <- FALSE
  duplicates <- FALSE
  entry <- 'a'
  col_names <- colnames(data_pre_proc)
  col_names_std <- c("sample", "ens_ID", "expr", "met", "met_weight", "met_context")
  
  # Prompt the user for information on the data type
  cat(sprintf("\nAfter each of the following prompts, either enter the\n"))
  cat(sprintf("requested column name, or the number 0 if there is no such\n"))
  cat(sprintf("column in your data set.\n"))
  cat(sprintf("\nThe column names of the provided data set are:\n"))
  for(i in 1:length(col_names)){
    cat(sprintf("'%s'  ", col_names[i]))
    if(i %% 4 == 0){
      cat(sprintf("\n"))
    }
  }
  cat(sprintf("\n"))
  while(col_name_sample == 0){
    col_name_sample <- readline(prompt = "Sample/cell ID column >>: ")
  }
  while(col_name_ensID == 0){
    col_name_ensID <- readline(prompt = "Ensembl ID column >>: ")
  }
  
  col_name_expr <- readline(prompt = "Expression rate/count column >>: ")
  if(col_name_expr != 0){
    is_expr <- TRUE
  }
  col_name_met <- readline(prompt = "Methylation rate column >>: ")
  if(col_name_met != 0){
    is_met <- TRUE
  }
  col_name_met_weight <- readline(prompt = "Methylation rate weights column >>: ")
  if(col_name_met_weight != 0){
    is_weighted <- TRUE
  }
  col_name_context <- readline(prompt = "Methylation context column >>: ")
  if(col_name_context != 0){
    is_context <- TRUE
  }
  
  while(col_name_expr == 0 && col_name_met == 0)
  {
    cat(sprintf("\nNeed to set a name for expression column or\n"))
    cat(sprintf("methylationin column! At least one of the following must\n"))
    cat(sprintf("have a name..\n"))
    col_name_expr <- readline(prompt = "Expression rate/count column >>: ")
    if(col_name_expr != 0){
      is_expr <- TRUE
    }
    col_name_met <- readline(prompt = "Methylation rate column >>: ")
    if(col_name_met != 0){
      is_met <- TRUE
    }
  }
  
  
  while(is_expr && ((norm_method != "TPM") && 
                               (norm_method != "RPKM") && 
                               (norm_method != "FPKM") && 
                               (norm_method != "UMI") && 
                               (norm_method != "raw")))
  {
    cat(sprintf("\nPlease enter the normalization method used for expression.\n"))
    cat(sprintf("NOTE:\n"))
    cat(sprintf("Enter 'TPM' for transcripts per million\n"))
    cat(sprintf("      'RPKM' for reads per kilobase million\n"))
    cat(sprintf("      'FPKM' for fragments per kilobase million\n"))
    cat(sprintf("      'UMI' for unique molecular identifiers\n"))
    cat(sprintf("      'raw' for raw reads\n"))
    norm_method <- readline(prompt = ">>: ")
  }
  
  cat(sprintf("\nAre there any duplicates in the data, or gene's with\n"))
  cat(sprintf("contexts that were measured more than once per gene per sample?"))
  user_choice <- 0
  
  while((user_choice != 'y') && 
        (user_choice != 'n') && 
        (user_choice != 'no') && 
        (user_choice != 'yes')){
    user_choice <- readline(prompt = "(y/n) >>: ")
    user_choice <- tolower(user_choice)
  }
  if(user_choice == 'y' || user_choice == 'yes'){
    duplicates <- TRUE
  }
  
  cat(sprintf("\nDo all samples contain all genes/methylation contexts used in the data?"))
  user_choice <- 0
  
  while((user_choice != 'y') && 
        (user_choice != 'n') && 
        (user_choice != 'no') && 
        (user_choice != 'yes')){
    user_choice <- readline(prompt = "(y/n) >>: ")
    user_choice <- tolower(user_choice)
  }
  if(user_choice == 'y' || user_choice == 'yes'){
    is_complete <- TRUE
  }
  
  user_choice <- 0
  
  
  
  
  
  ##### Create a New Data Frame ##################################################
  # Enumerate the columns that are present in the data, and create a new data
  # frame ('data_pre_process_new') from the original.
  included_columns <- c(col_name_expr,
                        col_name_met,
                        col_name_met_weight,
                        col_name_context)[
                          which(c(col_name_expr,
                                  col_name_met,
                                  col_name_met_weight,
                                  col_name_context) != 0)]
  data_pre_proc_new <- cbind(data_pre_proc[,which(col_names == col_name_sample)],
                             data_pre_proc[,which(col_names == col_name_ensID)])
  for(i in 1:length(included_columns)){
    data_pre_proc_new <- cbind(data_pre_proc_new, data_pre_proc[,which(col_names == included_columns[i])])
  }
  data_pre_proc_trimmed <- data_pre_proc_new
  colnames(data_pre_proc_trimmed) <- c(col_names_std[1:2], col_names_std[which(c(is_expr, is_met, is_weighted, is_context)) + 2])
  rm(data_pre_proc)
  
  # #Trim the obs. that don't contain all elements needed ---REMOVE ONCE COMPLETE??
  # data_pre_proc_trimmed <- 
  #   data_pre_proc_new[which(rowSums(is.na(data_pre_proc_new[,1:ncol(data_pre_proc_new)])) == 0),]
  # colnames(data_pre_proc_trimmed) <- c(col_names_std[1:2], col_names_std[which(c(is_expr, is_met, is_weighted, is_context)) + 2]) 
  # rm(data_pre_proc_new)
  
  #### Call Next Script ##########################################################
  # Call (source()) the next script depending on what type of data is to be 
  # included in the pseudotime analysis.
  if(is_expr && is_met){
    source(file = "pre_process_monocle_expr.R", local = FALSE)
    source(file = "pre_process_monocle_met.R", local = FALSE)
  }
  if(is_expr && !is_met){
    source(file = "pre_process_monocle_expr.R", local = FALSE)
  }
  if(!is_expr && is_met){
    source(file = "pre_process_monocle_met.R", local = FALSE)
  }
}







