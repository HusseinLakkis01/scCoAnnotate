args <- commandArgs(TRUE)

run_SciBet <- function(RefPath,LabelsPath,TestPath,OutputDir){
  "
	Author: Hussein Lakkis
	Date: 2021-07-24
	run  classifier: SciBet 
	Wrapper script to run an SciBet classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	RefPath : Ref Data file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv) .
	TestPath : Test dataset path : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDir : Output directory defining the path of the exported file.
  "
  # load libraries
  suppressMessages(library(tidyverse))
  suppressMessages(library(scibet))
  suppressMessages(library(viridis))
  suppressMessages(library(ggsci))
  
  # Read the ref data
  Data <- read.csv(RefPath,row.names = 1)
  # Read Ref labels
  Labels <- as.data.frame(read.csv(LabelsPath))
  # read test expression
  test <- read.csv(TestPath,row.names = 1)

  #############################################################################
  #SciBet #
  #############################################################################
  Training_Time_scibet <- list()
  # prepare data
  Data = as.data.frame(Data)
  test = as.data.frame(test)
  colnames(test) <- toupper(colnames(test))
  colnames(Data) <- toupper(colnames(Data))
  Data$label = Labels$label
 # SCIBET  Training 
  start_time <- Sys.time()
  pred <- SciBet(Data, test)
  end_time <- Sys.time()
  # get training and testing time
  Training_Time_SciBet <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  #tidy up the predictions
  predicted = data.frame(SciBet_Pred = pred, row.names = row.names(test))
  
  setwd(OutputDir)
  # write down and save output
  write.csv(predicted,paste('SciBet_pred.csv', sep = ''),row.names = TRUE)
  write.csv(Training_Time_SciBet,paste('SciBet_training_time.csv', sep = ''),row.names = FALSE)
  write.csv(Training_Time_SciBet,paste('SciBet_test_time.csv', sep = ''),row.names = FALSE)

}

run_SciBet(args[1], args[2], args[3], args[4])
sessionInfo()

