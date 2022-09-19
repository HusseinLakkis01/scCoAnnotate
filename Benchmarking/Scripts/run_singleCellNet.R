args <- commandArgs(TRUE)

run_singleCellNet<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir){
  "
  run singleCellNet
  Wrapper script to run singleCellNet on a benchmark dataset with 5-fold cross validation,
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
  colnames(Data) <- gsub('_','.',colnames(Data), fixed = TRUE)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]


  #############################################################################
  #                              singleCellNet                                #
  #############################################################################
  library(singleCellNet)
  library(Seurat)
  library(dplyr)
  True_Labels_singleCellNet <- list()
  Pred_Labels_singleCellNet <- list()
  Training_Time_singleCellNet <- list()
  Testing_Time_singleCellNet <- list()
  Data = t(as.matrix(Data))              # deals also with sparse matrix

  for(i in c(1:n_folds)){

    DataTrain <- Data[,Train_Idx[[i]]]
    DataTest <- Data[,Test_Idx[[i]]]
    start_time <- Sys.time()
      # Create seurat object
    meta <- data.frame(Annotation = as.factor(Labels[Train_Idx[[i]]]), row.names = colnames(DataTrain))
    Seurat <- CreateSeuratObject(DataTrain, meta.data = meta)
    print(head(Seurat@meta.data))
    stList<-splitCommon(sampTab = Seurat@meta.data, ncells = 80, dLevel = "Annotation")
     # get the downsampled list
    meta <- stList[[1]]
    # start training
    start_time <- Sys.time()
    class_info<-scn_train(stTrain = meta, expTrain = as.matrix(GetAssayData(Seurat))[,row.names(meta)], nTopGenes = 12, nRand = 70, nTrees = 350, nTopGenePairs = 25, dLevel = "Annotation")
    end_time <- Sys.time()
    # get training time
    Training_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    # do the testing
    start_time <- Sys.time()
    crPBMC <- scn_predict(class_info[['cnProc']], DataTest, nrand = 0)
    stQuery <- assign_cate(classRes = crPBMC, sampTab = data.frame(row.names = colnames(DataTest)), cThresh = 0.5) 
    end_time <- Sys.time()

    Testing_Time_singleCellNet[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))

    True_Labels_singleCellNet[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_singleCellNet[i] <- list(stQuery)
  }
  True_Labels_singleCellNet <- as.vector(unlist(True_Labels_singleCellNet))
  Pred_Labels_singleCellNet <- as.vector(unlist(Pred_Labels_singleCellNet))
  Training_Time_singleCellNet <- as.vector(unlist(Training_Time_singleCellNet))
  Testing_Time_singleCellNet <- as.vector(unlist(Testing_Time_singleCellNet))
  write.csv(True_Labels_singleCellNet,paste0(OutputDir,'/singleCellNet_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_singleCellNet,paste0(OutputDir,'/singleCellNet_pred.csv'),row.names = FALSE)
  write.csv(Training_Time_singleCellNet,paste0(OutputDir,'/singleCellNet_training_time.csv'),row.names = FALSE)
  write.csv(Testing_Time_singleCellNet,paste0(OutputDir,'/singleCellNet_test_time.csv'),row.names = FALSE)
}

run_singleCellNet(args[1], args[2], args[3], args[4])

