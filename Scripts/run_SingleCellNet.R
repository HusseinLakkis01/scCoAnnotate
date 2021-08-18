args <- commandArgs(TRUE)
run_singleCellNet <- function(DataPath, LabelsPath, TestPath, OutputDir){
  "
	Author: Hussein Lakkis
	Date: 2021-08-13
	run  classifier: SingleCellNet 
	Wrapper script to run an SingleCellNet classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	RefPath : Ref Data file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv) .
	TestPath : Test dataset path : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDir : Output directory defining the path of the exported file.
  "
  # load libraries
  library(singleCellNet)
  library(dplyr)
  library(Seurat)
  
  message("reading reference")
  
  # read data and labels
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.data.frame(read.csv(LabelsPath))
  Labels <- Labels$label
  message("reading Test")
  # read test
  Test <- read.csv(TestPath,row.names = 1)
  # Map gene names to upper
  colnames(Data) <- toupper(  colnames(Data))
  colnames(Test) <- toupper(  colnames(Test))
  # transpose data due to SingleCellNet requirements
  Test <- t(as.matrix(Test))
  Data = t(as.matrix(Data))  
  # subset based on common genes
  commonGenes<-intersect(rownames(Data), rownames(Test))
  Data <- Data[commonGenes, ]
  Test <- Test[commonGenes, ]
  # singleCellNet                  
  # Create seurat object
  meta <- data.frame(Annotation = Labels, row.names = colnames(Data))
  Seurat <- CreateSeuratObject(Data, meta.data = meta)
  # downsample for computational reasons
  stList<-splitCommon(sampTab = Seurat@meta.data, ncells = 80, dLevel = "Annotation")
  # get the downsampled list
  meta <- stList[[1]]
  # start training
  start_time <- Sys.time()
  class_info<-scn_train(stTrain = meta, expTrain = as.matrix(GetAssayData(Seurat))[,row.names(meta)], nTopGenes = 12, nRand = 70, nTrees = 350, nTopGenePairs = 25, dLevel = "Annotation")
  end_time <- Sys.time()
  # get training time
  Training_Time_SCN <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  # do the testing
  start_time <- Sys.time()
  crPBMC <- scn_predict(class_info[['cnProc']], Test, nrand = 0)
  stQuery <- assign_cate(classRes = crPBMC, sampTab = data.frame(row.names = colnames(Test)), cThresh = 0.5) 
  end_time <- Sys.time()
  # get testing time
  Test_Time_SCN <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  # tidy up results
  Pred_Labels_singleCellNet <-data.frame(SCN_Prediction = stQuery, row.names = colnames(Test))
  colnames(Pred_Labels_singleCellNet) = c("SCN_Prediction")
  # write down and save output
  write.csv(as.data.frame(Pred_Labels_singleCellNet),paste0(OutputDir,'/SingleCellNet_pred.csv'),row.names = TRUE)
  write.csv(Training_Time_SCN,paste0(OutputDir,'/SingleCellNet_training_time.csv'),row.names = FALSE)
  write.csv(Test_Time_SCN,paste0(OutputDir,'/SingleCellNet_test_time.csv'),row.names = FALSE)
}

run_singleCellNet(args[1], args[2], args[3], args[4])
sessionInfo()
