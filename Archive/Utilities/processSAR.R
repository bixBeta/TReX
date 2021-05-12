#!/usr/bin/env Rscript

#fileNames <- list.files(paste0(getwd()), pattern = ".complete")
#filePath <- paste0(getwd(), "/", fileNames)
library(progress)
pb <- progress_bar$new(total = 100)


for (i in 1:20) {
  pb$tick()
  Sys.sleep(1 / 10)
}

arg <-  commandArgs(trailingOnly = T)

suppressPackageStartupMessages(library(dplyr))


if (length(arg)==0) {
  print(" Usage = Rscript processSAR.R < PIN > ")  
  stop("Please provide the PIN !!! \n", call.=FALSE)
  
} 

dir = paste0(getwd(), "/")

fileNames <- list.files( dir, pattern = ".complete")
filePath <- paste0(dir, fileNames)





# import SAR tool files as objects 


for (i in 1:length(fileNames)) {
      
      assign(strsplit(fileNames[i], "\\.")[[1]][1],
      value = read.table(filePath[i], header = TRUE),
      envir = .GlobalEnv)
}

for (i in 1:30) {
  pb$tick()
  Sys.sleep(1 / 100)
}

for (i in 1:length(objects(pattern = "vs"))) {
  
  z <- c("dispGeneEst","dispFit","dispMAP","dispersion","betaConv","maxCooks" )
  pattern <- as.data.frame(mget(strsplit(fileNames[i], "\\.")[[1]][1]), col.names = NULL)
  assign(paste0(strsplit(fileNames[i], "\\.")[[1]][1], ".rmLast6"),
         value = select(pattern, -z),
         envir = .GlobalEnv)
  
}

for (i in 1:30) {
  pb$tick()
  Sys.sleep(1 / 100)
}

obj <- objects(pattern = ".rmLast6")

for (i in 1:length(obj)) {
  assign(x = paste0(obj[i],".Annotated"),
        value = as.data.frame(mget(obj[i]),col.names = NULL) %>% 
          rename(!!paste0(strsplit(obj[i], "\\.")[[1]][1],".FoldChange") := FoldChange,
                 !!paste0(strsplit(obj[i], "\\.")[[1]][1],".log2FoldChange") := log2FoldChange,
                 !!paste0(strsplit(obj[i], "\\.")[[1]][1],".stat") := stat,
                 !!paste0(strsplit(obj[i], "\\.")[[1]][1],".pvalue") := pvalue,
                 !!paste0(strsplit(obj[i], "\\.")[[1]][1],".padj") := padj),
           envir = .GlobalEnv
        
          )
}


obj2 <- objects(pattern = "Annotated")

new.data.frame <- as.data.frame(mget(obj2[1]),col.names = NULL)

for (i in 2:length(obj2)){
  
  new.data.frame <- full_join(new.data.frame, 
  select(as.data.frame(mget(obj2[i]),col.names = NULL), Id, matches("vs")), 
  by = "Id")
  
}

write.table(new.data.frame, paste0(arg[1],".final.txt"), sep = "\t", quote = F, row.names = F)


for (i in 1:20) {
  pb$tick()
  Sys.sleep(1 / 100)
}
