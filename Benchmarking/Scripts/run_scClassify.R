args <- commandArgs(TRUE)
library(data.table)

run_scClassify<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir){
  "
  run scClassify
  Wrapper script to run scClassify on a benchmark dataset with 5-fold cross validation,
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
  #                                scClassify                                     #
  #############################################################################
# load libraries
  suppressMessages(library(scClassify))
  suppressMessages(library(Seurat))

  True_Labels_scClassify <- list()
  Pred_Labels_scClassify <- list()
  Total_Time_scClassify <- list()
  Data = t(as.data.frame(Data))

    
  Data <- CreateSeuratObject(counts = Data)
  Data@meta.data$celltype = as.factor(Labels)
  Data <- NormalizeData(Data)
  Data <- GetAssayData(Data, slot = 'data')
  for (i in c(1:n_folds)){
    
      Train <- as(Data[, Train_Idx[[i]]], "dgCMatrix")
      Test <- as(Data[, Test_Idx[[i]]], "dgCMatrix")
      label <- Labels[Train_Idx[[i]]]
      print(length(Test_Idx[[i]]))
      print(dim(Data))
    
      start_time <- Sys.time()
      scClassify_res <- scClassify(exprsMat_train = Train,
                             cellTypes_train = label,
                             exprsMat_test = list(test = Test),
                             tree = "HOPACH",
                             algorithm = "WKNN",
                             selectFeatures = c("limma"),
                             similarity = c("pearson"),
                             returnList = FALSE,
                             verbose = FALSE)     
      end_time <- Sys.time()
    pred <- scClassify_res$testRes$test$pearson_WKNN_limma$predRes
    Total_Time_scClassify[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_scClassify[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scClassify[i] <- list(pred)
  }
  True_Labels_scClassify <- as.vector(unlist(True_Labels_scClassify))
  Pred_Labels_scClassify <- as.vector(unlist(Pred_Labels_scClassify))
  Total_Time_scClassify <- as.vector(unlist(Total_Time_scClassify))
  write.csv(True_Labels_scClassify,paste0(OutputDir,'/scClassify_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_scClassify,paste0(OutputDir,'/scClassify_pred.csv'),row.names = FALSE)
  write.csv(Total_Time_scClassify,paste0(OutputDir,'/scClassify_total_time.csv'),row.names = FALSE)
}


  run_scClassify(args[1], args[2], args[3], args[4])

  