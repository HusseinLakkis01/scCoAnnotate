args <- commandArgs(TRUE)

run_scPred<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir){
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
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  library(data.table)
  # Read the reference expression matrix
  Data <- fread(DataPath,data.table=FALSE)
  row.names(Data) <- Data$V1
  Data <-  Data[, 2:ncol(Data)]  
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]

  
  #############################################################################
  #                                scPred                                     #
  #############################################################################
  library("scPred")
  library("Seurat")
  library("magrittr")
  True_Labels_scPred <- list()
  Pred_Labels_scPred <- list()
  Training_Time_scPred <- list()
  Testing_Time_scPred <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){
      print(dim(Data[,Train_Idx[[i]]]))
      meta = data.frame(cell_type = Labels[Train_Idx[[i]]])
  
      reference <- CreateSeuratObject(counts = Data[,Train_Idx[[i]]])
      reference$cell_type <- Labels[Train_Idx[[i]]]
      reference <- reference %>% 
      NormalizeData() %>% 
      FindVariableFeatures() %>% 
      ScaleData() %>% 
      RunPCA()
      reference <- getFeatureSpace(reference, "cell_type")
      query <- CreateSeuratObject(counts = Data[,Test_Idx[[i]]])
      query <- NormalizeData(query)
      library(doParallel)
      cl <- makePSOCKcluster(2)
      registerDoParallel(cl)
      start_time <- Sys.time()
      set.seed(1234)
      reference <- trainModel(reference)
      end_time <- Sys.time()
      stopCluster(cl)
      Training_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

      start_time <- Sys.time()
      scpred <- get_scpred(reference)
      query <- scPredict(query, scpred,recompute_alignment = FALSE)
      end_time <- Sys.time()
      Testing_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

      True_Labels_scPred[i] <- list(Labels[Test_Idx[[i]]])
     Pred_Labels_scPred[i] <- list(query$scpred_prediction)}
  
  True_Labels_scPred <- as.vector(unlist(True_Labels_scPred))
  Pred_Labels_scPred <- as.vector(unlist(Pred_Labels_scPred))
  Training_Time_scPred <- as.vector(unlist(Training_Time_scPred))
  Testing_Time_scPred <- as.vector(unlist(Testing_Time_scPred))
  
  setwd(OutputDir)
  


    write.csv(True_Labels_scPred,'scPred_true.csv',row.names = FALSE)
    write.csv(Pred_Labels_scPred,'scPred_pred.csv',row.names = FALSE)
    write.csv(Training_Time_scPred,'scPred_training_time.csv',row.names = FALSE)
    write.csv(Testing_Time_scPred,'scPred_test_time.csv',row.names = FALSE)
  
}

run_scPred(args[1], args[2], args[3], args[4])
