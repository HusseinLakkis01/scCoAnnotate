run_scmap <- function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scmap
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
  
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
  GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #   scmap   #
  #############################################################################
  library(scmap)
  library(SingleCellExperiment)
  True_Labels_scmapcluster <- list()
  Pred_Labels_scmapcluster <- list()
  True_Labels_scmapcell <- list()
  Pred_Labels_scmapcell <- list()
  Training_Time_scmapcluster <- list()
  Testing_Time_scmapcluster <- list()
  Training_Time_scmapcell <- list()
  Testing_Time_scmapcell <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:n_folds)){
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    sce <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]), 
    colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))
    logcounts(sce) <- log2(normcounts(sce) + 1)
    # use gene names as feature symbols
    rowData(sce)$feature_symbol <- rownames(sce)
    sce <- selectFeatures(sce, n_features = NumGenes, suppress_plot = TRUE)
    
    sce_test <- SingleCellExperiment(list(normcounts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]), 
     colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
    logcounts(sce_test) <- log2(normcounts(sce_test) + 1)
    rowData(sce_test)$feature_symbol <- rownames(sce_test)
    sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
  }
  else{
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
  }
  
  # scmap-cluster
  start_time <- Sys.time()
  sce <- indexCluster(sce)
  end_time <- Sys.time()
  Training_Time_scmapcluster[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  
  start_time <- Sys.time()
  scmapCluster_results <- scmapCluster(projection = sce_test,index_list = list(metadata(sce)$scmap_cluster_index))
  end_time <- Sys.time()
  Testing_Time_scmapcluster[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  
  True_Labels_scmapcluster[i] <- list(Labels[Test_Idx[[i]]])
  Pred_Labels_scmapcluster[i] <- list(scmapCluster_results$combined_labs)
  
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
  
  True_Labels_scmapcluster <- as.vector(unlist(True_Labels_scmapcluster))
  Pred_Labels_scmapcluster <- as.vector(unlist(Pred_Labels_scmapcluster))
  True_Labels_scmapcell <- as.vector(unlist(True_Labels_scmapcell))
  Pred_Labels_scmapcell <- as.vector(unlist(Pred_Labels_scmapcell))
  Training_Time_scmapcluster <- as.vector(unlist(Training_Time_scmapcluster))
  Testing_Time_scmapcluster <- as.vector(unlist(Testing_Time_scmapcluster))
  Training_Time_scmapcell <- as.vector(unlist(Training_Time_scmapcell))
  Testing_Time_scmapcell <- as.vector(unlist(Testing_Time_scmapcell))
  
  setwd(OutputDir)
  


  write.csv(True_Labels_scmapcluster,'scmapcluster_true.csv',row.names = FALSE)
  write.csv(Pred_Labels_scmapcluster,'scmapcluster_pred_Labels.csv',row.names = FALSE)
  write.csv(True_Labels_scmapcell,'scmapcell_True_Labels.csv',row.names = FALSE)
  write.csv(Pred_Labels_scmapcell,'scmapcell_Pred_Labels.csv',row.names = FALSE)
  write.csv(Training_Time_scmapcluster,'scmapcluster_Training_Time.csv',row.names = FALSE)
  write.csv(Testing_Time_scmapcluster,'scmapcluster_Testing_Time.csv',row.names = FALSE)
  write.csv(Training_Time_scmapcell,'scmapcell_Training_Time.csv',row.names = FALSE)
  write.csv(Testing_Time_scmapcell,'scmapcell_Testing_Time.csv',row.names = FALSE)
  }
}
