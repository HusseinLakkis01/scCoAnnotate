library(singleCellNet)
library(Seurat)
library(dplyr)
library(data.table)
  

run_SingleCellNet <- function(RefPath, LabelsPath, QueryPaths, OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: SingleCellNet 
	Wrapper script to run an SingleCellNet classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
	----------
	RefPath : Training reference data file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv) .
	QueryPaths : Test dataset paths : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDirs : Output directory defining the path of the exported files for each query.
  "
  # Read the reference expression matrix
  Ref <- fread(RefPath,data.table=FALSE)
  row.names(Ref) <- Ref$V1
  cellnames <- row.names(Ref)
  Ref <-  Ref[, 2:ncol(Ref)]

  # Read Ref labels
  Labels <- as.data.frame(read.csv(LabelsPath))

  print(dim(Labels))
  print(dim(Ref))

  #############################################################################
  #SingleCellNet #
  #############################################################################

  # prepare data
  Ref <- t(as.data.frame(Ref))
  colnames(Ref) <- cellnames
  Labels <- as.vector(Labels$label)

  # Create seurat object
  meta <- data.frame(Annotation = Labels, row.names = cellnames)
  Seurat <- CreateSeuratObject(counts = Ref, meta.data = meta)  
  stList <- splitCommon(sampTab = Seurat@meta.data, ncells = 80, dLevel = "Annotation")

  # get the downsampled list
  meta <- stList[[1]]
  # start training
  start_time <- Sys.time()
  class_info <- scn_train(stTrain = meta, expTrain = as.matrix(GetAssayData(Seurat))[,row.names(meta)], nTopGenes = 12, nRand = 70, nTrees = 350, nTopGenePairs = 25, dLevel = "Annotation")
  end_time <- Sys.time()

  # get training time
  Training_Time_SCN <- as.numeric(difftime(end_time,start_time,units = 'secs')) 

  message("@reading test")
  i = 1
  for(query in QueryPaths){

      OutputDir <- OutputDirs[[i]]

      # read query
      query <- fread(query,data.table=FALSE)
      row.names(query) <- query$V1
      query <-  query[, 2:ncol(query)]
      cellnames <- row.names(query)

      # Transpose as input to scn
      query <- t(query)

      # train and predict
      start_time <- Sys.time()
      crPBMC <- scn_predict(class_info[['cnProc']], query, nrand = 0)
      stQuery <- assign_cate(classRes = crPBMC, sampTab = data.frame(row.names = cellnames), cThresh = 0.5) 
      end_time <- Sys.time()

      # get total time
      query_Time_SCN <- as.numeric(difftime(end_time,start_time,units = 'secs'))

     # tidy up
      Pred_Labels_SingleCellNet <- data.frame(SCN_Prediction = stQuery)
      row.names(Pred_Labels_SingleCellNet) = cellnames
      colnames(Pred_Labels_SingleCellNet) = c("SingleCellNet")

      # create dir for scn
      dir.create(file.path(OutputDir, "SingleCellNet"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "SingleCellNet"))

      # write down and save the output
      write.csv(Pred_Labels_SingleCellNet,paste('SingleCellNet','_pred.csv', sep = ''))
      write.csv(Training_Time_SCN,paste0('SingleCellNet_training_time.csv'),row.names = FALSE)
      write.csv(query_Time_SCN,paste0('SingleCellNet_query_time.csv'),row.names = FALSE)      
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

run_SingleCellNet(ref,labs, query, output_dir)

sessionInfo()
