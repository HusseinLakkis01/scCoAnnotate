args <- commandArgs(TRUE)
library(data.table)

# Function to label by Spearman correlation by (Selin Jessa and Marie Coutlier)
label_correlation <- function(test_expr_mat,
                              ref_expr_mat,
                              threshold_common_genes = 0.5) {
  
  rownames(test_expr_mat) <- toupper(rownames(test_expr_mat))
  rownames(ref_expr_mat) <- toupper(rownames(ref_expr_mat))
  
  # Testing how many genes are in common and stopping if not enough
  common_genes <- intersect(rownames(test_expr_mat), rownames(ref_expr_mat))
  prop_common <- length(common_genes) / length(rownames(test_expr_mat))
  message("@@ ", round(prop_common*100, digits = 2), "% of test dataset genes are in the reference dataset")
  
  if (prop_common < threshold_common_genes) stop("Proportion of common genes below threshold.")
  
  # Reducing matrices to common subset to. both test and reference
  mat1 <- as.matrix(test_expr_mat[common_genes, ])
  mat2 <- as.matrix(ref_expr_mat[common_genes, ])
  
  # sanity check
  nrow(mat1) == nrow(mat2)
  
  # Computing correlations
  cor_matrix <- cor(mat1, mat2, method = "spearman", use = "complete.obs")
  
  # Getting the best one
  cor_label <- as.data.frame(cor_matrix) %>%
    mutate("cell" = rownames(cor_matrix)) %>%
    gather("celltype", "correlation", -cell) %>%
    group_by(cell) %>%
    top_n(1, correlation) %>%
    dplyr::select(cell, celltype, correlation) %>% 
    arrange(cell)
  
  # Returning the results
  return(cor_label)
  
  }



run_Correlation<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir){
  "
  run Correlation
  Wrapper script to run Correlation on a benchmark dataset with 5-fold cross validation,
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
  data <-  Data[, 2:ncol(Data)]
  colnames(data) <- gsub('_','.',colnames(data), fixed = TRUE)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  data <- data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]

  #############################################################################
  #                              Correlation                                #
  #############################################################################
# load libraries
  # loads needed libraries
  library(here)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(Seurat)
  True_Labels_Correlation <- list()
  Pred_Labels_Correlation <- list()
  Total_Time_Correlation <- list()
  data = as.data.frame(data)
  
  for (i in c(1:n_folds)){
    
     train <-  t(data[Train_Idx[[i]],])
     print(length(Test_Idx[[i]]))
     print(dim(data))
     test <- t(as.matrix(data[Test_Idx[[i]],]))
    
     Data <- CreateSeuratObject(counts = train)
     test <- CreateSeuratObject(counts = test)
     Data@meta.data$celltype = as.factor(Labels[Train_Idx[[i]]])
     # generate mean expression matrix
     Idents(Data) <- "celltype"
     Data <- NormalizeData(Data)
     test <- NormalizeData(test)

     Data <- AverageExpression(Data, return.seurat = TRUE)
     # get mean expression matrix
     Ref_mean <- GetAssayData(Data)
     # transpose test Dataset
     test <- as.matrix(GetAssayData(test))
     row.names(test) <- toupper(row.names(test))
     test[1:3,1:3]
     # Correlation  Training time
     start_time <- Sys.time()
     # Predict the cells
     predicted = as.data.frame(label_correlation(test,Ref_mean,0.45))
     colnames(predicted) <- c("", "Correlation_Prediction","Correlation_score")
     True_Labels_Correlation[i] <- list(Labels[Test_Idx[[i]]])

     Pred_Labels_Correlation[i] <- list(predicted$Correlation_Prediction)
     end_time <- Sys.time()
     # get the training and testing time
    Training_Time_corr <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  }
  print(Pred_Labels_Correlation)
  True_Labels_Correlation <- as.vector(unlist(True_Labels_Correlation))
  Pred_Labels_Correlation <- as.vector(unlist(Pred_Labels_Correlation))
  Total_Time_Correlation <- as.vector(unlist(Total_Time_Correlation))
  write.csv(True_Labels_Correlation,paste0(OutputDir,'/Correlation_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_Correlation,paste0(OutputDir,'/Correlation_pred.csv'),row.names = FALSE)
  write.csv(Training_Time_corr,paste0(OutputDir,'/Correlation_total_time.csv'),row.names = FALSE)
}


  run_Correlation(args[1], args[2], args[3], args[4])
