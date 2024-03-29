args <- commandArgs(TRUE)

run_SingleR<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir){
  "
  run SingleR
  Wrapper script to run SingleR on a benchmark dataset with 5-fold cross validation,
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
  #                               SingleR                                     #
  #############################################################################
  library(SingleR)
  library(Seurat)
  True_Labels_SingleR <- list()
  Pred_Labels_SingleR <- list()
  Total_Time_SingleR <- list()
  Data = t(as.matrix(Data))

  for (i in c(1:n_folds)){

    start_time <- Sys.time()
    singler = SingleR(method = "single", Data[,Test_Idx[[i]]], Data[,Train_Idx[[i]]], Labels[Train_Idx[[i]]])
    end_time <- Sys.time()

    Total_Time_SingleR[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_SingleR[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_SingleR[i] <- list(as.vector(singler$labels))
  }
  True_Labels_SingleR <- as.vector(unlist(True_Labels_SingleR))
  Pred_Labels_SingleR <- as.vector(unlist(Pred_Labels_SingleR))
  Total_Time_SingleR <- as.vector(unlist(Total_Time_SingleR))

  write.csv(True_Labels_SingleR,paste0(OutputDir,'/SingleR_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_SingleR,paste0(OutputDir,'/SingleR_pred.csv'),row.names = FALSE)
  write.csv(Total_Time_SingleR,paste0(OutputDir,'/SingleR_total_time.csv'),row.names = FALSE)
}


run_SingleR(args[1], args[2], args[3], args[4])

