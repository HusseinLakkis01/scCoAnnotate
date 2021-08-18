
args <- commandArgs(TRUE)

select_genes <- function(DataPath, LabelsPath, num_genes, OutputDir){
        library(praznik)
        message("reading data")
        data <- read.csv(DataPath, row.names = 1)
        message("reading labels")
        labels <- read.csv(LabelsPath, row.names = 1)
        message("MRMR")
        selection.result = MRMR(data, as.factor(labels$label), k = num_genes, threads = 0)
        write.csv(as.data.frame(selection.result), file = paste0(OutputDir,'/MRMR.csv'),row.names = TRUE)

}
select_genes(args[1], args[2], args[3], args[4])

