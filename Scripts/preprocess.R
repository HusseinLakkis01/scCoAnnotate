library(glue)
library(data.table)

preprocess_data<-function(RefPath, QueryPaths, OutputDir, check_genes = FALSE, genes_required = NULL){
  "
  run preproccessing, normalization, and writing down of matrices with common genes
  Wrapper script to run preproccessing, normalization, and writing down of matrices with common genes between query and reference datasets
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  RefPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations (per cell) file path (.csv).
  QueryPaths : List of Query datasets paths
  OutputDir : Output directory defining the path of the exported files.
  genes_required: path to file containing required genes in common feature space with column name genes_required
  "

  if(check_genes){
  genes_required <- as.vector(read.csv(genes_required$genes_required))}

  message("@Reading reference")
  # read data and labels
  # Read the reference expression matrix
  Data <- fread(RefPath,data.table=FALSE)
  message("@Done reading reference")

  row.names(Data) <- Data$V1
  Data <-  Data[, 2:ncol(Data)]
  colnames(Data) <- toupper(colnames(Data))

  # set the genes of the reference to common_genes
  common_genes <- colnames(Data)

  # loop over all query datasets and get common genes with all
  for(query in QueryPaths){
    message("reading Test")
    # Read Query
    Test <- fread(query,data.table=FALSE)
    row.names(Test) <- Test$V1
    Test <-  Test[, 2:ncol(Test)]
    # Map gene names to upper
    colnames(Test) <- toupper(colnames(Test))
    # subset based on common genes
    common_genes<- as.vector(intersect(common_genes, colnames(Test)))
  }

  # check if the common genes
    if(check_genes){
  if(!(genes_required %in% common_genes)){
    warning("Not all genes required in Gene list are found in the common genes between reference and query")
  }}

  # loop over query datasets and subset them to common genes and write it down in each query output dir
  for(query in QueryPaths){
    # Read Query
    Test <- fread(query,data.table=FALSE)
    row.names(Test) <- Test$V1
    Test <-  Test[, 2:ncol(Test)]
    colnames(Test) <- toupper(colnames(Test))
    Test <- Test[,common_genes]

    # write down query expression matrices
    dir.create(file.path(OutputDir, basename(dirname(query))), showWarnings = FALSE)
    setwd(paste(OutputDir, basename(dirname(query)),sep= '/'))
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
test <- unlist(options.args['query'])
output_dir <- unlist(options.args['output_dir' ])
check_genes <- unlist(options.args['check_genes' ])
genes_required <- unlist(options.args['genes_required' ])


preprocess_data(ref, test, output_dir, check_genes, genes_required)


