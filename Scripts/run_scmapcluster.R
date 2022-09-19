library(scmap)
library(SingleCellExperiment)
library(data.table)

run_scmapcluster <- function(RefPath, LabelsPath, QueryPaths, OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: scmap-cluster 
	Wrapper script to run an scmap-cluster classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	RefPath :  Training reference file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv) .
	QueryPaths : Query datasets paths : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDirs : Output directory defining the path of the exported files for each query.
  "
  
  message("reading Ref")
  print(RefPath)
  # Read the Ref expression matrix
  Ref <- as.data.frame(fread(RefPath,data.table=FALSE, header = T))
  row.names(Ref) <- Ref$V1
  Ref <-  Ref[, 2:ncol(Ref)]
  Labels <- as.data.frame(read.csv(LabelsPath))
  Labels <- Labels$label
  message("reading Test")


  #############################################################################
  #                                 scmap-cluster                                #
  #############################################################################

  # Create singlecell experiment object as input and normalize
  Ref <- t(Ref)
  sce <- SingleCellExperiment(list(normcounts = Ref), 
                                  colData = data.frame(cell_type1 = Labels))
  logcounts(sce) <- log2(normcounts(sce) + 1)

  # use gene names as feature symbols
  rowData(sce)$feature_symbol <- rownames(sce)

  # scmap-cell
  start_time <- Sys.time()
  sce <- selectFeatures(sce, suppress_plot = TRUE)

  # train scmap-cell
  set.seed(1)
  sce <- indexCluster(sce)
  end_time <- Sys.time()
  Training_Time_scmapcluster <- as.numeric(difftime(end_time,start_time,units = 'secs'))

  # loop over query datasets
  # Use counter
  i=1

  for(query in QueryPaths){

      # Get current output ir for current query
      OutputDir <- OutputDirs[[i]]

      # Read Query
      query <- fread(query,data.table=FALSE)
      row.names(query) <- query$V1
      query <-  query[, 2:ncol(query)]
      cellnames <- row.names(query)

      # transpose and create sce object
      query <- t(query)
      sce_test <- SingleCellExperiment(list(normcounts = query))
      logcounts(sce_test) <- log2(normcounts(sce_test) + 1)

      # Use gene names as features
      rowData(sce_test)$feature_symbol <- rownames(sce_test)
      sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData

      # predict on each query
      start_time <- Sys.time()
      scmapCluster_results <- scmapCluster(projection = sce_test,index_list = list(metadata(sce)$scmap_cluster_index))
      end_time <- Sys.time()

      # get labels and train/query time
      Query_time_scmapcluster <- as.numeric(difftime(end_time,start_time,units = 'secs'))
      Pred_Labels_scmapcluster <- list(scmapCluster_results$combined_labs)
      Pred_Labels_scmapcluster <- as.vector(Pred_Labels_scmapcluster)

      # Set prediction Dataframe index
      predicted = data.frame(scmapcluster = Pred_Labels_scmapcluster, row.names = cellnames)
      colnames(predicted) <- c('scmapcluster')

      # Create scmapcluster subdir in target dir
      dir.create(file.path(OutputDir, "scmapcluster"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "scmapcluster"))

      # write down and save the output
      write.csv(predicted,paste('scmapcluster','_pred.csv', sep = ''))
      write.csv(Training_Time_scmapcluster, paste('scmapcluster','_training_time.csv', sep = ''),row.names = FALSE)
      write.csv(Query_time_scmapcluster, paste('scmapcluster','_query_time.csv', sep = ''),row.names = FALSE)
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
ref <- unlist(options.args['ref'])
labs <- unlist(options.args['labs'])
query <- unlist(options.args['query'])
output_dir <- unlist(options.args['output_dir' ])

# Run scmap-cluster
run_scmapcluster(ref,labs, query, output_dir)
# Output session info
sessionInfo()



