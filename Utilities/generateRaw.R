#!/usr/bin/env Rscript
args <-  commandArgs(trailingOnly = T)

# check for required argument
if (length(args)==0) {
  print("=======================================================")
  print(" Usage = Rscript test.R < PIN > < .Rdata file > ")  
  print("=======================================================")
  stop("Both arguments must be supplied!!! \n", call.=FALSE)
  
} 

PIN <- args[1]
load(args[2])

write.table(counts, paste0(PIN, "_rawCounts.txt"), sep = "\t", quote = F, col.names = NA)