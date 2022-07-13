# load libraries
library(CHETAH)
library(SingleCellExperiment)
library(data.table)

run_CHETAH <- function(RefPath,LabelsPath,TestPaths,OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: CHETAH 
	Wrapper script to run an CHETAH classifier 
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
  # read reference
  Data <- fread(RefPath,data.table=FALSE)
  row.names(Data) <- Data$V1
  Data <-  Data[, 2:ncol(Data)]
  # Read Ref labels
  Labels <- as.data.frame(read.csv(LabelsPath, row.names = 1))
  #############################################################################
                                      #CHETAH #
  #############################################################################
  # prepare data
  Data <- t(as.data.frame(Data))
  # create sce matrix as input
  message("@creating reference")
  sce <- SingleCellExperiment(list(counts = Data), 
                                  colData = data.frame(celltypes = Labels$label))

  # Loop over test datasets
  message("@reading test and predicting")
  # Set a counter
  i = 1
  # Loop
  for(Test in TestPaths){
      # Get current output for current query
      OutputDir <- OutputDirs[[i]]
      # read query 
      test <- fread(Test,data.table=FALSE)
      row.names(test) <- test$V1
      cellnames <- row.names(test)
      test <-  test[, 2:ncol(test)]
      # transpose it for SingleCellExperiment
      test <- t(test)
      sce_test <- SingleCellExperiment(assays = list(counts = test))
      # CHETAH  Training and testinh
      start_time <- Sys.time()
      sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
      end_time <- Sys.time()
      # get the training and testing time
      Training_Time_CHETAH <- as.numeric(difftime(end_time,start_time,units = 'secs'))
      # create Data frame with predictions and cell names
      prediction <-  data.frame(CHETAH = c(sce_test$celltype_CHETAH), row.names = cellnames) 
      colnames(prediction) <- c ("CHETAH")
      # Create CHETAH subdir per query
      dir.create(file.path(OutputDir, "CHETAH"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "CHETAH"))
      # write down and save the output
      write.csv(prediction,paste('CHETAH','_pred.csv', sep = ''))
      write.csv(Training_Time_CHETAH,paste('CHETAH','_total_time.csv', sep = ''),row.names = FALSE)
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


run_CHETAH(ref,labs, test, output_dir)
sessionInfo()
