args <- commandArgs(TRUE)

run_scmapcell <- function(DataPath,LabelsPath,CV_RDataPath,OutputDir){
  "
  run scmapcell
  Wrapper script to run scmap on a benchmark dataset with 5-fold cross validation,
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
  #                                 scmap                                     #
  #############################################################################
  library(scmap)
  library(SingleCellExperiment)
  True_Labels_scmapcell <- list()
  Pred_Labels_scmapcell <- list()
  Training_Time_scmapcell <- list()
  Testing_Time_scmapcell <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){

      sce <- SingleCellExperiment(list(normcounts = Data[,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
      logcounts(sce) <- log2(normcounts(sce) + 1)
      # use gene names as feature symbols
      rowData(sce)$feature_symbol <- rownames(sce)
      sce <- selectFeatures(sce, suppress_plot = TRUE)
      
      sce_test <- SingleCellExperiment(list(normcounts = Data[,Test_Idx[[i]]]), 
                                       colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
      rowData(sce_test)$feature_symbol <- rownames(sce_test)
      sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
    
    
    # scmap-cell
    start_time <- Sys.time()
    set.seed(1)
    sce <- indexCell(sce)
    end_time <- Sys.time()
    Training_Time_scmapcell[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    start_time <- Sys.time()
    scmapCell_results <- scmapCell(sce_test,list(metadata(sce)$scmap_cell_index))
    scmapCell_clusters <- scmapCell2Cluster(scmapCell_results,list(as.character(colData(sce)$cell_type1)))
    end_time <- Sys.time()
    Testing_Time_scmapcell[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    
    True_Labels_scmapcell[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scmapcell[i] <- list(scmapCell_clusters$combined_labs)
  }
  
  True_Labels_scmapcell <- as.vector(unlist(True_Labels_scmapcell))
  Pred_Labels_scmapcell <- as.vector(unlist(Pred_Labels_scmapcell))
  Training_Time_scmapcell <- as.vector(unlist(Training_Time_scmapcell))
  Testing_Time_scmapcell <- as.vector(unlist(Testing_Time_scmapcell))
  
  write.csv(True_Labels_scmapcell,paste0(OutputDir,'/scmapcell_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_scmapcell,paste0(OutputDir,'/scmapcell_pred.csv'),row.names = FALSE)
  write.csv(Training_Time_scmapcell,paste0(OutputDir,'/scmapcell_training_time.csv'),row.names = FALSE)
  write.csv(Testing_Time_scmapcell,paste0(OutputDir,'/scmapcell_test_time.csv'),row.names = FALSE)
}

run_scmapcell(args[1], args[2], args[3], args[4])


