# load libraries
library(data.table)
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))

run_SciBet <- function(RefPath,LabelsPath,TestPaths,OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: SciBet 
	Wrapper script to run an SciBet classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	RefPath : Ref Data file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv) .
	TestPaths : Test dataset paths : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDirs : Output directory defining the path of the exported files for each query.
  "
  
  # Read the reference expression matrix
  Data <- fread(RefPath,data.table=FALSE)
  # Read reference labels
  Labels <- as.data.frame(read.csv(LabelsPath))
  #############################################################################
                                     #SciBet #
  #############################################################################
  # prepare data for training 
  Data = as.data.frame(Data)
  colnames(Data) <- toupper(colnames(Data))
  Data[] <- lapply(Data, as.numeric)
  Data$label = Labels$label
  # Loop over test datasets
  message("@reading test and predicting")
  # Set a counter
  i = 1
  # Loop
  for(Test in TestPaths){
      # Get current output ir for current query
      OutputDir <- OutputDirs[[i]]
      # Read Query
      test <- fread(Test,data.table=FALSE)
      row.names(test) <- test$V1
      test <-  test[, 2:ncol(test)]
      # SCIBET  Training 
      start_time <- Sys.time()
      pred <- SciBet(Data, test)
      end_time <- Sys.time()
      # get the training and testing time
      Training_Time_SciBet <- as.numeric(difftime(end_time,start_time,units = 'secs'))
      predicted = data.frame(SciBet_Pred = pred, row.names = row.names(test))
      # Set prediction Dataframe index
      row.names(predicted) <- row.names(test)
      print(head(predicted))
      # Create SciBet subdir in target dir
      dir.create(file.path(OutputDir, "SciBet"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "SciBet"))
      # write down and save the output
      write.csv(predicted,paste('SciBet','_pred.csv', sep = ''))
      write.csv(Training_Time_SciBet,paste('SciBet','_training_time.csv', sep = ''),row.names = FALSE)
      write.csv(Training_Time_SciBet,paste('SciBet','_test_time.csv', sep = ''),row.names = FALSE)
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

test <- unlist(options.args['test'])
output_dir <- unlist(options.args['output_dir' ])

# Run SciBet
run_SciBet(ref,labs, test, output_dir)
# Output session info
sessionInfo()

