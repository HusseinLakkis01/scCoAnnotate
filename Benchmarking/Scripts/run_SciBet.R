args <- commandArgs(TRUE)
library(data.table)
run_SciBet<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir){
  "
  run SciBet
  Wrapper script to run SciBet on a benchmark dataset with 5-fold cross validation,
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
  #                                SciBet                                     #
  #############################################################################
# load libraries
  suppressMessages(library(tidyverse))
  suppressMessages(library(scibet))
  suppressMessages(library(viridis))
  suppressMessages(library(ggsci))
  True_Labels_SciBet <- list()
  Pred_Labels_SciBet <- list()
  Total_Time_SciBet <- list()
  Data = as.data.frame(Data)
  
  for (i in c(1:n_folds)){
    
      train <-  Data[Train_Idx[[i]],]
      train[] <- lapply(train, as.numeric)
      train$label <- Labels[Train_Idx[[i]]]
      print(length(Test_Idx[[i]]))
      print(dim(Data))
      test <- Data[Test_Idx[[i]],] 
    
      start_time <- Sys.time()
      pred <- scibet::SciBet(train, test)
      print(pred[1:5])
      end_time <- Sys.time()
    
    Total_Time_SciBet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_SciBet[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_SciBet[i] <- list(pred)
  }
  True_Labels_SciBet <- as.vector(unlist(True_Labels_SciBet))
  Pred_Labels_SciBet <- as.vector(unlist(Pred_Labels_SciBet))
  Total_Time_SciBet <- as.vector(unlist(Total_Time_SciBet))
  write.csv(True_Labels_SciBet,paste0(OutputDir,'/SciBet_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_SciBet,paste0(OutputDir,'/SciBet_pred.csv'),row.names = FALSE)
  write.csv(Total_Time_SciBet,paste0(OutputDir,'/SciBet_total_time.csv'),row.names = FALSE)
}


  run_SciBet(args[1], args[2], args[3], args[4])

  