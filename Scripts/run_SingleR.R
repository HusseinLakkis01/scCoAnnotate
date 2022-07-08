# load libraries
library(data.table)
library(SingleR)
library(SingleCellExperiment)

run_SingleR <- function(RefPath,LabelsPath,TestPaths,OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: SingleR 
	Wrapper script to run an SingleR classifier 
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

  # Read the reference data
  Data <- as.data.frame(fread(RefPath,data.table=FALSE, header = T))
  row.names(Data) <- Data$V1
  Data <-  Data[, 2:ncol(Data)]
  # Read Reference labels
  print(Data[1:4,1:4])

  Labels <- as.data.frame(read.csv(LabelsPath,row.names = 1))
  print(head(Labels))
  #############################################################################
                                      #SingleR #
  #############################################################################
  # prepare data
  Data <-  t(Data)
  Labels <- as.factor(Labels$label)
  # Train SingleR
  start_time <- Sys.time()
  print(Data[1:4,1:4])
  singler <- trainSingleR(Data, labels=Labels)
  end_time <- Sys.time()
  Training_Time_SingleR <- as.numeric(difftime(end_time,start_time,units = 'secs'))

  # Loop over test datasets
  message("@reading test and predicting")
  # Set a counter
  i = 1
  # Loop
  for(Test in TestPaths){
      # Get current output for current query
      OutputDir <- OutputDirs[[i]]
      # Read Query
      test <- fread(Test,data.table=FALSE,header = TRUE)
      row.names(test) <- test$V1
      test <-  test[, 2:ncol(test)]
      # save cell names for prediction output
      cellnames <- row.names(test)
      test <- t(test)
      # train and predict
      start_time <- Sys.time()
      pred <- classifySingleR(test, singler, assay.type=1)
      end_time <- Sys.time()
      # get total time
      Test_Time_SingleR <- as.numeric(difftime(end_time,start_time,units = 'secs'))
     # tidy up prediction dataframe
      Pred_Labels_SingleR <- as.vector(pred$labels)
      Pred_Labels_SingleR <- data.frame(SingleR =Pred_Labels_SingleR,row.names = cellnames)
      Test_Time_SingleR <- as.numeric(difftime(end_time,start_time,units = 'secs'))
      # Create SingleR subdir in target dir
      dir.create(file.path(OutputDir, "SingleR"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "SingleR"))
      # write down and save the output
      write.csv(Pred_Labels_SingleR,paste('SingleR','_pred.csv', sep = ''))
      write.csv(Training_Time_SingleR,paste('SingleR','_training_time.csv', sep = ''),row.names = FALSE)
      write.csv(Test_Time_SingleR,paste('SingleR','_test_time.csv', sep = ''),row.names = FALSE)

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

run_SingleR(ref,labs, test, output_dir)



