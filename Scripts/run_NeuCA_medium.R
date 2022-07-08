args <- commandArgs(TRUE)
run_NeuCA<- function(DataPath,LabelsPath,TestPath,OutputDir){
  "
  run NeuCA
  Wrapper script to run NeuCA on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  "
  
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
  # transpose data due to CHETAH's requirements
  Test <- t(as.matrix(Test))
  Data = t(as.matrix(Data))  
  # subset based on common genes
  commonGenes<-intersect(rownames(Data), rownames(Test))
  Data <- Data[commonGenes, ]
  Test <- Test[commonGenes, ]


  library(NeuCA)
  library(dplyr)
  
#############################################################################
  #                                NeuCA                                     #
  #############################################################################
  library(NeuCA)
  library(SingleCellExperiment)
  Pred_Labels_NeuCA <- list()
  Total_Time_NeuCA <- list()
  Data <- t(as.matrix(Data))
  Test <- t(as.matrix(Test))

      sce <- SingleCellExperiment(assays = list(counts = Data, 
                                  colData = data.frame(celltypes = Labels)
      
      sce_test <- SingleCellExperiment(assays = list(counts = Test)
      start_time <- Sys.time()
      predicted.label = NeuCA(train = sce, test = sce_test, 
                        model.size = "medium", verbose = T)
      end_time <- Sys.time()
    
    Total_Time_NeuCA <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    Pred_Labels_NeuCA <- list(predicted.label)


  }
  Pred_Labels_NeuCA <- as.vector(unlist(Pred_Labels_NeuCA))
  Total_Time_NeuCA <- as.vector(unlist(Total_Time_NeuCA))
  write.csv(Pred_Labels_NeuCA,paste0(OutputDir,'/NeuCA_medium_pred.csv'),row.names = FALSE)
  write.csv(Total_Time_NeuCA,paste0(OutputDir,'/NeuCA_medium_total_time.csv'),row.names = FALSE)
}


run_NeuCA(args[1], args[2], args[3], args[4])