args <- commandArgs(TRUE)

run_SingleR<-function(DataPath,LabelsPath,TestPath,OutputDir){
  "
	Author: Hussein Lakkis
	Date: 2021-07-22
	run  classifier: SingleR 
	Wrapper script to run an SingleR classifier 
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
  library(SingleR)
  library(Seurat)
  # read datasets
  Data <- read.csv(DataPath,row.names = 1)
  Test <- read.csv(TestPath,row.names = 1)
  Labels <- as.data.frame(read.csv(LabelsPath))
  Labels <- as.vector(Labels$label)
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
  # train and predict
  start_time <- Sys.time()
  singler = SingleR(method = "single", Test, Data, Labels)
  end_time <- Sys.time()
  # get total time
  Total_Time_SingleR <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  # tidy up
  Pred_Labels_SingleR <- as.vector(singler$labels)
  Pred_Labels_SingleR <- data.frame(singleR_Prediction =Pred_Labels_SingleR,row.names = colnames(Test))
  # write down predictions
  write.csv(Pred_Labels_SingleR,paste0(OutputDir,'/SingleR_pred.csv'),row.names = TRUE )
  write.csv(Total_Time_SingleR,paste0(OutputDir,'/SingleR_test_time.csv'),row.names = FALSE)
  write.csv(Total_Time_SingleR,paste0(OutputDir,'/SingleR_training_time.csv'),row.names = FALSE)
  write.csv(singler,paste0(OutputDir,'/SingleR_allpred.csv') )

}

run_SingleR(args[1], args[2], args[3], args[4])
sessionInfo()


