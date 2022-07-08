library(glue)
library(data.table)

preprocess_data<-function(RefPath,TestPaths,OutputDir){
  "
  run preproccessing, normalization, and writing down of matrices withb common genes
  Wrapper script to run preproccessing, normalization, and writing down of matrices withb common genes on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  RefPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations (per cell) file path (.csv).
  TestPaths : List of Query datasets paths
  OutputDir : Output directory defining the path of the exported file.
  "
  print(OutputDir)
  print(RefPath)

  message("reading reference")
  # read data and labels
  # Read the reference expression matrix
  Data <- fread(RefPath,data.table=FALSE)
  row.names(Data) <- Data$V1
  Data <-  Data[, 2:ncol(Data)]
  colnames(Data) <- toupper(colnames(Data))
  # set the genes of the reference to common_genes
  common_genes <- colnames(Data)
  # loop over all query datasets and get common genes with all
  for(TestPath in TestPaths){
    message("reading Test")
    # Read Query
    Test <- fread(TestPath,data.table=FALSE)
    row.names(Test) <- Test$V1
    Test <-  Test[, 2:ncol(Test)]
    # Map gene names to upper
    colnames(Test) <- toupper(colnames(Test))
    # subset based on common genes
    common_genes<-intersect(common_genes, colnames(Test))
  }
  # loop over query datasets and subset them to common genes and write it down in each query output dir
  for(TestPath in TestPaths){
    # Read Query
    Test <- fread(TestPath,data.table=FALSE)
    row.names(Test) <- Test$V1
    Test <-  Test[, 2:ncol(Test)]
    colnames(Test) <- toupper(colnames(Test))
    Test <- Test[,common_genes]
    dir.create(file.path(OutputDir, basename(dirname(TestPath))), showWarnings = FALSE)
    setwd(paste(OutputDir, basename(dirname(TestPath)),sep= '/'))
    fwrite(Test, "expression.csv",row.names=T)
   }
  Data <- Data[,common_genes ]
  fwrite(Data, glue("{OutputDir}/expression.csv"),row.names=T)
  fwrite(as.data.frame(common_genes), glue("{OutputDir}/common_genes.csv"))
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
test <- unlist(options.args['test'])
output_dir <- unlist(options.args['output_dir' ])

preprocess_data(ref,test,output_dir)


