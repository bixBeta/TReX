#!/usr/bin/env Rscript

suppressWarnings(library("dplyr"))

names 	<- list.files(getwd(), pattern = ".rawCounts" )
gff 	<- "/Users/fa286/bin/gff.Hash.Aegl5.0.txt"

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
