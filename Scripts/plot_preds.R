# load R libraries
library(Seurat)
library(data.table)
library(ggplot2)
library(glue)
library(dplyr)
library(tidyverse)

source("/project/kleinman/hussein.lakkis/from_hydra/scCoAnnotate/Scripts/style/theme_min_SelinJessa.R")

plot_figures <- function(TestPaths,OutputDirs){
  "
	Author: Hussein Lakkis
	Date: 2022-08-21
	Wrapper script to run an plot figures
	outputs pdfs of predictions and UMAP embeddings
  
  Parameters
	----------
  TestPaths: file path for each query dataset containing raw counts in the cell x gene format
	OutputDirs : Output directory defining the path of the exported files for each query.
  "

  # Set a counter
  i = 1
  for(Test in TestPaths){
      # Get current output ir for current query
      OutputDir <- OutputDirs[[i]]
      # read predictions outputted by all classifiers
      predictions <- as.data.frame(fread(glue("{OutputDir}/Prediction_Summary.tsv")))
      predictions <- dplyr::select(predictions, where(is.character))
      row.names(predictions) <- predictions$cellname
      predictions <- predictions[, -c(1)]
      # Read Query
      test <- fread(Test,data.table=FALSE)
      row.names(test) <- test$V1
      test <-  test[, 2:ncol(test)]

      # transpose test Dataset
      test <- t(as.matrix(test))
      row.names(test) <- toupper(row.names(test))

      # Create Seurat 
      seurat <- CreateSeuratObject(counts = test, meta.data = predictions)
      seurat <- NormalizeData(seurat)
      all.genes <- rownames(seurat)
      seurat <- ScaleData(seurat, features = all.genes)
      seurat <- FindVariableFeatures(seurat)

      # Run DimRed and cluster
      seurat <- RunPCA(seurat)
      seurat <- RunUMAP(seurat, dims=1:10)
      seurat <- FindNeighbors(seurat, dims = 1:10)
      seurat = FindClusters(seurat, resolution = 0.5)
      cell_embeddings <- as.data.frame(seurat[["umap"]]@cell.embeddings[, c(1,2)])
      rownames(cell_embeddings) <- colnames(seurat)
      cell_embeddings <- cell_embeddings %>% rownames_to_column(var = "Cell")
      metadata <- seurat@meta.data
      colnames(cell_embeddings) <- c("Cell", "UMAP_1", "UMAP_2")

      dir.create(file.path(OutputDir, "output"))
      dir.create(file.path(OutputDir, "figures"))

      # save output
      fwrite(cell_embeddings, file = glue("{OutputDir}/output/UMAP_embeddings.csv"), sep = ",", row.names = T, col.names = T)
      fwrite(metadata, file = glue("{OutputDir}/output/metadata.csv"), sep = ",", row.names = T, col.names = T)
      row.names(cell_embeddings) <- cell_embeddings$Cell

      # loop over tools and plot predictions per tool (UMAPs)
      for(tool in colnames(predictions)){
        pred <-  predictions[,tool]
        umap_df <- cbind(cell_embeddings, Cell_Type = pred)

        umap_df %>% 
            ggplot(aes(x = UMAP_1, y = UMAP_2)) +
            geom_point(aes(colour = Cell_Type), size = 0.4, alpha = 1) +
            ggtitle(paste(tool)) +
            theme_min() +
            theme(legend.position = "none", 
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank()) 
        
        ggsave(glue("{OutputDir}/figures/{tool}.png"), dpi = 100, width = 4, height =4)
      }

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
test <- unlist(options.args['query'])
output_dir <- unlist(options.args['output_dir' ])

# Run plot function
plot_figures(test, output_dir)
# output info session
sessionInfo()
