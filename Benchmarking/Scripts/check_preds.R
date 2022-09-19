args <- commandArgs(TRUE)

output_dir <- args[1]

message <- "done"

write.csv(as.data.frame(message), paste(output_dir, "/evaluation/finished.csv", sep=""))