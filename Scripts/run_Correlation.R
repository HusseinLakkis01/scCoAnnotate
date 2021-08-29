args <- commandArgs(TRUE)

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

run_correlation<-function(RefPath,LabelsPath,TestPath,OutputDir){
  "
	Author: Hussein Lakkis
	Date: 2021-07-22
	run  classifier: Correlation by Selin and Marie 
	Wrapper script to run an Correlation classifier 
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
  # loads needed libraries
  library(here)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(Seurat)
  # read Data and labels
  Data <- read.csv(RefPath,row.names = 1)
  Labels <- as.data.frame(read.csv(LabelsPath, row.names = 1))
  # read test path
  test <- read.csv(TestPath,row.names = 1)
  # transpose dataset to create seurat object
  Data <- t(as.data.frame(Data))
  row.names(Data) <- toupper(row.names(Data))
  # test if the rownames are correct
  row.names(Data)[1:4]
  # cr
  Data <- CreateSeuratObject(counts = Data, meta.data = Labels)
  Data@meta.data$Cell.Type = Data@meta.data$label
  # generate mean expression matrix
  Idents(Data) <- "Cell.Type"
  Data <- AverageExpression(Data, return.seurat = TRUE)
  # get mean expression matrix
  Ref_mean <- GetAssayData(Data)
  # transpose test Dataset
  test <- t(as.matrix(test))
  row.names(test) <- toupper(row.names(test))
  test[1:3,1:3]
  # Correlation  Training time
  start_time <- Sys.time()
  # Predict the cells
  predicted = as.data.frame(label_correlation(test,Ref_mean,0.45))
  colnames(predicted) <- c("", "Correlation_Prediction","Correlation_score")
  end_time <- Sys.time()
  # get the training and testing time
  Training_Time_corr <- as.numeric(difftime(end_time,start_time,units = 'secs'))
  
  setwd(OutputDir)
  # write down and save the output
  write.csv(predicted,paste('correlation','_pred.csv', sep = ''),row.names = FALSE)
  write.csv(Training_Time_corr,paste('correlation','_training_time.csv', sep = ''),row.names = FALSE)
  write.csv(Training_Time_corr,paste('correlation','_test_time.csv', sep = ''),row.names = FALSE)

}

run_correlation(args[1], args[2], args[3], args[4])
sessionInfo()


