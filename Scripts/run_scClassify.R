# load libraries
suppressMessages(library(scClassify))
suppressMessages(library(Seurat))
library(data.table)

run_scClassify <- function(RefPath,LabelsPath,TestPaths,OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: scClassify 
	Wrapper script to run an scClassify classifier 
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
  row.names(Data) <- Data$V1
  Data <-  Data[, 2:ncol(Data)]
  # Read Ref labels
  Labels <- as.data.frame(read.csv(LabelsPath))
  #############################################################################
  #scClassify #
  #############################################################################
  # prepare data and create Seurat Object 
  Data = t(as.data.frame(Data))
  colnames(Data) <- toupper(colnames(Data))
  Data <- CreateSeuratObject(counts = Data)
  Data@meta.data$celltype = as.factor(Labels$label)
  Data <- NormalizeData(Data)
  Data <- GetAssayData(Data, slot = 'data')
  Data <- as(Data, "dgCMatrix")
  # Loop over test datasets
  message("@reading test and predicting")
  # Set a counter
  i = 1
  for(Test in TestPaths){
      # Get current output ir for current query
      OutputDir <- OutputDirs[[i]]
      # Read Query
      test <- fread(Test,data.table=FALSE)
      row.names(test) <- test$V1
      test <-  test[, 2:ncol(test)]
      test <- t(test)
      # save cellnames
      cellnames <- colnames(test)
      # Read and normalize test
      test <- CreateSeuratObject(counts = test)
      test <- NormalizeData(test)
      test <- GetAssayData(test, slot = 'data')
      test <- as(test, "dgCMatrix")

      # scClassify  Training 
      start_time <- Sys.time()
      scClassify_res <- scClassify(exprsMat_train = Data,
                             cellTypes_train = Labels$label,
                             exprsMat_test = list(test = test),
                             tree = "HOPACH",
                             algorithm = "WKNN",
                             selectFeatures = c("limma"),
                             similarity = c("pearson"),
                             returnList = FALSE,
                             verbose = FALSE)       
                             
      end_time <- Sys.time()
      # get the training and testing time
      Training_Time_scClassify <- as.numeric(difftime(end_time,start_time,units = 'secs'))
      # Get predictions
      pred <- scClassify_res$testRes$test$pearson_WKNN_limma$predRes
      Pred_Labels_scClassify <- list(pred)
      # Create prediction dataframe 
      prediction <-  data.frame(scClassify = c(unlist(Pred_Labels_scClassify)), row.names = cellnames) 
      colnames(prediction) <- c ("scClassify")
      # Create scClassify subdir
      dir.create(file.path(OutputDir, "scClassify"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "scClassify"))
      # write down and save the output
      write.csv(prediction,paste('scClassify','_pred.csv', sep = ''))
      write.csv(Training_Time_scClassify,paste('scClassify','_training_time.csv', sep = ''),row.names = FALSE)
      write.csv(Training_Time_scClassify,paste('scClassify','_test_time.csv', sep = ''),row.names = FALSE)
      # update counter for next query
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

# Run scClassify
run_scClassify(ref,labs, test, output_dir)
# output info session
sessionInfo()

