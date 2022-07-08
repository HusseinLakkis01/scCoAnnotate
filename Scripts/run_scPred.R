

run_scPred<-function(DataPath,LabelsPath,TestPaths,OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: scPred 
	Wrapper script to run an scPred classifier 
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
  
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.data.frame(read.csv(LabelsPath, row.names = 1))
  cell_type <- Labels$label
  
  #############################################################################
  #scPred #
  #############################################################################
  library(scPred)
  library(tidyverse)
  library(magrittr)
  library(Seurat)

  Pred_Labels_scPred <- list()
  Training_Time_scPred <- list()
  Testing_Time_scPred <- list()
  Data = t(as.matrix(Data))

  reference <- CreateSeuratObject(counts = Data, meta.data = Labels)

  reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)
  
  reference <- getFeatureSpace(reference, "label")
  
  # scPred Training
  start_time <- Sys.time()
  set.seed(1234)
  # plotEigen(scp, group = 'cell_type1')
  reference <- trainModel(reference)
  # plotTrainProbs(scp)
  end_time <- Sys.time()
  Training_Time_scPred <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  message("@reading test")
  i = 1
  scpred <- get_scpred(reference)
  for(Test in TestPaths){

      OutputDir <- OutputDirs[[i]]
      test <- t(as.matrix(read.csv(Test,row.names = 1)))
      cellnames <- colnames(test)
      test <- CreateSeuratObject(counts = test)      # train and predict
      # scPred Prediction
      start_time <- Sys.time()
      test <- scPredict(test, scpred,recompute_alignment = FALSE)
      end_time <- Sys.time()
      Testing_Time_scPred <- as.numeric(difftime(end_time,start_time,units = 'secs'))
     # tidy up
      predicted <- data.frame(scPred = test$scpred_prediction, row.names = cellnames)

      colnames(predicted) = c("scPred")
      dir.create(file.path(OutputDir, "scPred"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "scPred"))
      # write down and save the output
      write.csv(predicted,paste('scPred','_pred.csv', sep = ''))
      write.csv(Training_Time_scPred,paste0('scPred_training_time.csv'),row.names = FALSE)
      write.csv(Testing_Time_scPred,paste0('scPred_test_time.csv'),row.names = FALSE)      
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

run_scPred(ref,labs, test, output_dir)



sessionInfo()



