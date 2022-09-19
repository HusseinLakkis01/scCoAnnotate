# load libraries
library(Seurat)
library(glue)
library(ggplot2)
library(pvclust)
require(data.tree)
library(apTreeshape)
library(dendextend)
library(Wind)
library(tidyverse)
args <- commandArgs(TRUE)

# function to get pvclust to phylo
as.phylo.hclust.with.nodenames <- function (x, nodenames, ...) #We add a nodenames argument
{
    N <- dim(x$merge)[1]
    edge <- matrix(0L, 2 * N, 2)
    edge.length <- numeric(2 * N)
    node <- integer(N)
    node[N] <- N + 2L
    cur.nod <- N + 3L
    j <- 1L
    for (i in N:1) {
        edge[j:(j + 1), 1] <- node[i]
        for (l in 1:2) {
            k <- j + l - 1L
            y <- x$merge[i, l]
            if (y > 0) {
                edge[k, 2] <- node[y] <- cur.nod
                cur.nod <- cur.nod + 1L
                edge.length[k] <- x$height[i] - x$height[y]
            }
            else {
                edge[k, 2] <- -y
                edge.length[k] <- x$height[i]
            }
        }
        j <- j + 2L
    }
    if (is.null(x$labels)) 
        x$labels <- as.character(1:(N + 1))
    node.lab <- nodenames[order(node)] #Here we define our node labels
    obj <- list(edge = edge, edge.length = edge.length/2, tip.label = x$labels, 
        Nnode = N, node.label = node.lab) #And you put them in the final object
    class(obj) <- "phylo"
    reorder(obj)
}

Prepare <- function(DataPath, LabelsPath, col_Index = 1, OutputDir){
  "
  Prepare Seurat 
  Function returns prepares the pca variables and VST genes for complexity metrics
  it also outputs umap embeddings

  Parameters
  ----------
  LabelsPath : Cell population annotations file path (.csv).
  col_Index : column index (integer) defining which level of annotation to use,
  in case of multiple cell type annotations (default is 1)
  OutputDir : Output directory defining the path of the exported file.
  "
  # READ data and labels
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  Labels <- as.vector(Labels[,col_Index])

  # create seurat object
  meta.data <- data.frame(label = Labels, row.names = row.names(Data))
  seurat <- CreateSeuratObject(counts = t(as.matrix(Data)), meta.data = meta.data)

  # Normalize and cluster
  seurat <- NormalizeData(seurat)
  all.genes <- rownames(seurat)
  seurat <- ScaleData(seurat, features = all.genes)
  seurat <- FindVariableFeatures(seurat)
  seurat <- RunPCA(seurat)
  seurat <- RunUMAP(seurat, dims=1:10)
  seurat <- FindNeighbors(seurat, dims = 1:10)
  seurat = FindClusters(seurat, resolution = 0.5)
  
  # plot umaps
  DimPlot(seurat, group.by = "label", pt.size = 1.5 )+NoAxes()

  # calculate average expreesion to cluster cell types
  Idents(seurat) <- seurat$label
  seurat_av <- AverageExpression(seurat, return.seurat = TRUE)
  expression_data <- GetAssayData(seurat_av)[VariableFeatures(seurat),]

  # use pvclust
  result <- pvclust(expression_data, method.dist="cor", method.hclust="complete", nboot=100)

  # get bootstrap values
  bootstraps <- (round(result$edges,2)*100)[,c(1,3)]
  
  # get newick format of the tree
  yy<-as.phylo.hclust.with.nodenames(result$hclust, nodenames=bootstraps[,2])
  newick = write.tree(yy,tree.names=TRUE,digits=2)


  setwd(OutputDir)

  # save embeddings
  cell_embeddings <- as.data.frame(seurat[["umap"]]@cell.embeddings[, c(1,2)])
  rownames(cell_embeddings) <- colnames(seurat)
  cell_embeddings <- cell_embeddings %>% rownames_to_column(var = "Cell")
  colnames(cell_embeddings) <- c("Cell", "UMAP_1", "UMAP_2")
  save(cell_embeddings, file = "UMAP_embeddings.Rda")
  save(seurat,file = "seurat.Rda" )
  
  # save dendrogram
  pdf("dendrogram.pdf",   
    width = 10, 
    height = 10) 
  plot(result)
  dev.off()


  write.csv(newick,"tree.newick_format.csv")
  write.csv(seurat[["pca"]]@cell.embeddings, "pca.csv")
  mat <- t(as.matrix(Seurat::GetAssayData(seurat, assay = "RNA")[VariableFeatures(seurat),]))
  write.csv(mat, "VariableFeatures.csv")
  ggsave("UMAP.png", dpi = 600)

  }
  
# run
Prepare(args[1],args[2],as.numeric(args[3]),args[4])


