# load libraries
library(data.table)
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))

run_SciBet <- function(RefPath,LabelsPath,QueryPaths,OutputDirs){

  "
	Author: Hussein Lakkis
	Date: 2022-06-24
	run  classifier: SciBet 
	Wrapper script to run an SciBet classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	RefPath : Training Reference file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv) .
	QueryPaths : query Refset paths : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDirs : Output directory defining the path of the exported files for each query.
  "
  
  # Read the reference expression matrix
  Ref <- fread(RefPath,data.table=FALSE)
  # Read reference labels
  Labels <- as.data.frame(read.csv(LabelsPath))

  #############################################################################
                                     #SciBet #
  #############################################################################

  # prepare Ref for training 
  Ref = as.data.frame(Ref)
  colnames(Ref) <- toupper(colnames(Ref))
  Ref[] <- lapply(Ref, as.numeric)
  Ref$label = Labels$label

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

      # SCIBET  Training 
      start_time <- Sys.time()
      pred <- SciBet(Ref, query)
      end_time <- Sys.time()

      # get the training and querying time
      Training_Time_SciBet <- as.numeric(difftime(end_time,start_time,units = 'secs'))
      predicted = data.frame(SciBet = pred, row.names = row.names(query))

      # Set prediction dataframe index
      row.names(predicted) <- row.names(query)
      print(head(predicted))

      # Create SciBet subdir in target dir
      dir.create(file.path(OutputDir, "SciBet"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "SciBet"))

      # write down and save the output
      write.csv(predicted,paste('SciBet','_pred.csv', sep = ''))
      write.csv(Training_Time_SciBet,paste('SciBet','_training_time.csv', sep = ''),row.names = FALSE)
      write.csv(Training_Time_SciBet,paste('SciBet','_query_time.csv', sep = ''),row.names = FALSE)
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

# Run SciBet
run_SciBet(ref,labs, query, output_dir)
# Output session info
sessionInfo()

