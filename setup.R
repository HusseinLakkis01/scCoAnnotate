list.of.packages <- c("ggplot2", "scibet", "Seurat", "singleCellNet", "viridis", "ggsci", "SingleR", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
