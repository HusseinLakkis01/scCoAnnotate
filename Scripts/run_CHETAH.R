# load libraries
library(CHETAH)
library(SingleCellExperiment)
library(data.table)

run_CHETAH <- function(RefPath,LabelsPath,QueryPaths,OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-06-22
	run  classifier: CHETAH 
	Wrapper script to run an CHETAH classifier 
	outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
	Parameters
	----------
	RefPath : Training reference file path (.csv), cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	LabelsPath : Cell population annotations file path (.csv) .
	QueryPaths : Query dataset paths : cells-genes matrix with cell unique barcodes 
	as row names and gene names as column names.
	OutputDirs : Output directory defining the path of the exported files for each query.
  "
  # read reference
  Ref <- fread(RefPath,data.table=FALSE)
  row.names(Ref) <- Ref$V1
  Ref <-  Ref[, 2:ncol(Ref)]

  # Read Ref labels
  Labels <- as.data.frame(read.csv(LabelsPath, row.names = 1))
  #############################################################################
                                      #CHETAH #
  #############################################################################
  # prepare data
  Ref <- t(as.data.frame(Ref))
  
  # create sce matrix as input
  message("@creating reference")
  sce <- SingleCellExperiment(list(counts = Ref), 
                                  colData = data.frame(celltypes = Labels$label))

  # Loop over Query datasets
  message("@reading Query and predicting")

  # Set a counter
  i = 1
  # Loop
  for(Query in QueryPaths){
      # Get current output for current query
      OutputDir <- OutputDirs[[i]]

      # read query 
      Query <- fread(Query,data.table=FALSE)
      row.names(Query) <- Query$V1
      cellnames <- row.names(Query)
      Query <-  Query[, 2:ncol(Query)]

      # transpose it for SingleCellExperiment
      Query <- t(Query)
      sce_Query <- SingleCellExperiment(assays = list(counts = Query))

      # CHETAH  Training and Queryinh
      start_time <- Sys.time()
      sce_Query <- CHETAHclassifier(input = sce_Query, ref_cells = sce)
      end_time <- Sys.time()

      # get the training and Querying time
      Training_Time_CHETAH <- as.numeric(difftime(end_time,start_time,units = 'secs'))

      # create Data frame with predictions and cell names
      prediction <-  data.frame(CHETAH = c(sce_Query$celltype_CHETAH), row.names = cellnames) 
      colnames(prediction) <- c ("CHETAH")
      
      # Create CHETAH subdir per query
      dir.create(file.path(OutputDir, "CHETAH"), showWarnings = FALSE)
      setwd(file.path(OutputDir, "CHETAH"))

      # write down and save the output
      write.csv(prediction,paste('CHETAH','_pred.csv', sep = ''))
      write.csv(Training_Time_CHETAH,paste('CHETAH','_total_time.csv', sep = ''),row.names = FALSE)
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


run_CHETAH(ref,labs, query, output_dir)
sessionInfo()
