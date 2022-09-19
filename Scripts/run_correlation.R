# loads needed libraries
library(tidyr)
library(dplyr)
library(readr)
library(Seurat)
library(data.table)

# Function to label by Spearman correlation by (Selin Jessa and Marie Coutlier)
label_correlation <- function(Query_expr_mat,
                              Ref_expr_mat,
                              threshold_common_genes = 0.3) {
  
  rownames(Query_expr_mat) <- toupper(rownames(Query_expr_mat))
  rownames(Ref_expr_mat) <- toupper(rownames(Ref_expr_mat))
  
  # Querying how many genes are in common and stopping if not enough
  common_genes <- intersect(rownames(Query_expr_mat), rownames(Ref_expr_mat))
  prop_common <- length(common_genes) / length(rownames(Query_expr_mat))
  message("@@ ", round(prop_common*100, digits = 2), "% of Query dataset genes are in the Reference dataset")
  
  if (prop_common < threshold_common_genes) stop("Proportion of common genes below threshold.")
  
  # Reducing matrices to common subset to. both Query and Reference
  mat1 <- as.matrix(Query_expr_mat[common_genes, ])
  mat2 <- as.matrix(Ref_expr_mat[common_genes, ])
  
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


run_correlation<-function(RefPath, LabelsPath, QueryPaths, OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-23
	run  classifier: Spearnman Correlation 
	Wrapper script to run an Spearnman Correlation classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
	----------
	RefPath : Training Reference file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv) .
	QueryPaths : Query dataset paths : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDirs : Output directory defining the path of the exported files for each query.
  "

  # Read the Reference expression matrix
  Ref <- fread(RefPath,data.table=FALSE)
  row.names(Ref) <- Ref$V1
  Ref <-  Ref[, 2:ncol(Ref)]
  Labels <- as.data.frame(read.csv(LabelsPath, row.names = 1))

  # transpose dataset to create seurat object
  Ref <- t(as.data.frame(Ref))
  row.names(Ref) <- toupper(row.names(Ref))
  # check if the rownames are correct
  row.names(Ref)[1:4]

  # create seurat object
  Ref <- CreateSeuratObject(counts = Ref, meta.data = Labels)
  Ref@meta.data$Cell.Type = Ref@meta.data$label

  # generate mean expression matrix
  Idents(Ref) <- "Cell.Type"
  Ref <- AverageExpression(Ref, return.seurat = TRUE)

  # get mean expression matrix
  Ref_mean <- GetAssayData(Ref)

  # Loop over Query datasets
  message("@reading Query and predicting")
  # Set a counter
  i = 1

  for(Query in QueryPaths){
      # Get current output ir for current query
      OutputDir <- OutputDirs[[i]]

      # Read Query
      Query <- fread(Query,data.table=FALSE)
      row.names(Query) <- Query$V1
      Query <-  Query[, 2:ncol(Query)]

      # transpose Query Dataset
      Query <- t(as.matrix(Query))
      row.names(Query) <- toupper(row.names(Query))

      # Correlation  Training time
      start_time <- Sys.time()

      # Predict the cells
      predicted = as.data.frame(label_correlation(Query,Ref_mean,0.3))
      colnames(predicted) <- c("", "correlation","correlation_score")
      end_time <- Sys.time()

      # get the training and Querying time
      Training_Time_corr <- as.numeric(difftime(end_time,start_time,units = 'secs'))
      dir.create(file.path(OutputDir, "correlation"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "correlation"))

      # write down and save the output
      write.csv(predicted,paste('correlation','_pred.csv', sep = ''),row.names = FALSE)
      write.csv(Training_Time_corr,paste('correlation','_training_time.csv', sep = ''),row.names = FALSE)
      write.csv(Training_Time_corr,paste('correlation','_query_time.csv', sep = ''),row.names = FALSE)
      i = i+1
  }
}
# Get Command line arguments
args <- commandArgs(trailingOnly = TRUE)
# Split the arguments to form lists
arguments <- paste(unlist(args),collapse=' ')
listoptions <- unlist(strsplit(arguments,'--'))[-1]
# Get individual argument names
options.args <- sapply(listoptions,function(x){
         unlist(strsplit(x, ' '))[-1]
        })
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})
# Set variables containing command line argument values
names(options.args) <- unlist(options.names)
Ref <- unlist(options.args['ref'])
labs <- unlist(options.args['labs'])

query <- unlist(options.args['query'])
output_dir <- unlist(options.args['output_dir' ])

# Run correlation
run_correlation(Ref,labs, query, output_dir)
# output info session
sessionInfo()


