#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("dplyr"))

args <-  commandArgs(trailingOnly = T)
gff <- args[1] # path for hash/dictionary file

names   <- list.files(getwd(), pattern = ".rawCounts" )
#gff    <- "/home/fa286/bin/scripts/gff.Hash.Aegl5.0.txt"

swapFunc <- function(gff,countFile,names){

        gff.mosquito <- read.table(gff, header = F)
        one <- read.table(countFile, header = F)
        swap <- left_join(one, gff.mosquito, by = c("V1" = "V1"))
        finalCount <- swap[,c(3,2)]
        write.table(finalCount, paste0(names,".swapCount"), quote = F, sep = "\t", col.names = F, row.names = F)

}

for (i in 1:length(names)){

        swapFunc(gff = gff, countFile = names[i], names = names[i] )

}
