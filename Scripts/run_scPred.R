# load libraries
library(scPred)
library(tidyverse)
library(magrittr)
library(Seurat)
library(data.table)

run_scPred<-function(RefPath, LabelsPath, QueryPaths, OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: scPred 
	Wrapper script to run an scPred classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	RefPath : Training reference data file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv) .
	QueryPaths : Test dataset paths : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDirs : Output directory defining the path of the exported files for each query.
  "
  
  # Read the Reference expression matrix
  Ref <- fread(RefPath,data.table=FALSE)
  row.names(Ref) <- Ref$V1
  Ref <-  Ref[, 2:ncol(Ref)]

  # Read reference labels
  Labels <- as.data.frame(read.csv(LabelsPath))
  
  #############################################################################
  #scPred #
  #############################################################################

  Ref <- t(as.matrix(Ref))

  # create seurat object for the reference
  reference <- CreateSeuratObject(counts = Ref)
  reference$cell_type <- Labels$label

  # Normalize and find dim red
  reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)
  
  reference <- getFeatureSpace(reference, "cell_type")

  # scPred Training
  start_time <- Sys.time()
  set.seed(1234)
  reference <- trainModel(reference)
  end_time <- Sys.time()

  # get training time
  Training_Time_scPred <- as.numeric(difftime(end_time,start_time,units = 'secs'))

  # Loop over query samples
  message("@reading query")
  i = 1

  # get the trained model
  scpred <- get_scpred(reference)
  
  # Loop over query datasets
  message("@reading query and predicting")

  # Set a counter
  i = 1

  # Loop
  for(query in QueryPaths){
      # Get current output ir for current query
      OutputDir <- OutputDirs[[i]]

      # Read Query
      query <- fread(query,data.table=FALSE)
      row.names(query) <- query$V1
      query <-  query[, 2:ncol(query)]

      cellnames <- row.names(query)

      # create seurat obeject
      query <- t(query)
      query <- CreateSeuratObject(counts = query)      
      
      # scPred Prediction
      start_time <- Sys.time()
      query <- scPredict(query, scpred,recompute_alignment = FALSE)
      end_time <- Sys.time()
      Query_time_scPred <- as.numeric(difftime(end_time,start_time,units = 'secs'))

     # tidy up
      predictions <- data.frame(scPred = query$scpred_prediction, row.names = cellnames)

      # Create scPred subdir in target dir
      dir.create(file.path(OutputDir, "scPred"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "scPred"))

      # write down and save the output
      write.csv(predictions,paste('scPred','_pred.csv', sep = ''))
      write.csv(Training_Time_scPred,paste('scPred','_training_time.csv', sep = ''),row.names = FALSE)
      write.csv(Query_time_scPred,paste('scPred','_query_time.csv', sep = ''),row.names = FALSE)
      i = i+1
  }
}

  
# Get Command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Split the arguments to form lists
arguments <- paste(unlist(args),collapse=' ')
listoptions <- unlist(strsplit(arguments,'--'))[-1]
# Get individual argument names
options.args <- sapply(listoptions,function(x){
         unlist(strsplit(x, ' '))[-1]
        })
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})
# Set variables containing command line argument values
names(options.args) <- unlist(options.names)
ref <- unlist(options.args['ref'])
labs <- unlist(options.args['labs'])

query <- unlist(options.args['query'])
output_dir <- unlist(options.args['output_dir' ])

run_scPred(ref,labs, query, output_dir)


sessionInfo()



