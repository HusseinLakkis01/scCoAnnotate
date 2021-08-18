args <- commandArgs(TRUE)

run_scPred<-function(DataPath,LabelsPath,TestPath,OutputDir){
  "
  run scPred
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
   : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.data.frame(read.csv(LabelsPath, row.names = 1))
  cell_type <- Labels$label
  Test <- read.csv(TestPath,row.names = 1)
  
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
  Test = t(as.matrix(Test))

  reference <- CreateSeuratObject(counts = Data, meta.data = Labels)

  reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

  Test <- CreateSeuratObject(counts = Test)
  
  reference <- getFeatureSpace(reference, "label")
  
  # scPred Training
  start_time <- Sys.time()
  set.seed(1234)
  # plotEigen(scp, group = 'cell_type1')
  reference <- trainModel(reference)
  # plotTrainProbs(scp)
  end_time <- Sys.time()
  Training_Time_scPred <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  # scPred Prediction
  start_time <- Sys.time()
  scp <- scPredict(Test, reference)
  end_time <- Sys.time()
  Testing_Time_scPred <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  Pred_Labels_scPred <- getPredictions(scp)$predClass

  Training_Time_scPred <- as.vector(Training_Time_scPred)
  Testing_Time_scPred <- as.vector(Testing_Time_scPred)
  
  setwd(OutputDir)
  predicted <- data.frame(scPred_Type = Pred_Labels_scPred, row.names = colnames(Test))
  write.csv(predicted,paste('scPred','_pred.csv', sep = ''),row.names = True)
  write.csv(Training_Time_corr,paste('scPred','_Training_Time.csv', sep = ''),row.names = FALSE)
  write.csv(Training_Time_corr,paste('scPred','_Testing_Time.csv', sep = ''),row.names = FALSE)

  
}

run_scPred(args[1], args[2], args[3], args[4])
