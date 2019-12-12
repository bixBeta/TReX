suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(dplyr))
load("1034.RData")

## -------------------------------------------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = target,
                              design = ~ group)

## -------------------------------------------------------------------------------------------------------------------
dds <- DESeq(dds)
#resultsNames(dds)

## -------------------------------------------------------------------------------------------------------------------
# add .custom to all contrast objects 

RF_vs_FF.custom <- results(dds, contrast=c("group", "RF", "FF"), alpha = 0.05)  
RF2_vs_FF2.custom <- results(dds, contrast=c("group", "RF", "FF"), alpha = 0.1)
# pHot_vs_pCold.custom <- results(dds, contrast=c("group", "pHOT", "pCOLD"), alpha = 0.05)
# cHot_vs_pHot.custom <- results(dds, contrast=c("group", "cHOT", "pHOT"), alpha = 0.05)
# cCold_vs_pCold.custom <- results(dds, contrast=c("group", "cCOLD", "pCOLD"), alpha = 0.05)

pattern1 <- objects(pattern = ".custom")

sink(file = "Contrasts_Summary.log")
for (i in 1:length(pattern1)) {
  
  print(paste0("Res Summary - ", pattern1[i], ":"))
  summary.DESeqResults(get(pattern1[i]))
  
}
sink()
## -------------------------------------------------------------------------------------------------------------------
norm.Counts <- as.data.frame(counts(dds, normalized = T))

norm.Counts <- rename_all(norm.Counts, .funs = ~ paste0("norm.", colnames(norm.Counts)))



# y <- as.data.frame(cColdHot.custom)
# y$Id <- row.names(y)
# joined <- full_join(norm.Counts, y, by= "Id")

## -------------------------------------------------------------------------------------------------------------------
# coercing DESeqResult object to data.frame
for (i in 1:length(pattern1)) {
  assign(x = paste0(pattern1[i], ".df"), 
         value = as.data.frame(mget(pattern1[i]), col.names = NULL),
         envir = .GlobalEnv
  )
}

pattern2 <- objects(pattern = ".df")

#Id <- as.data.frame(row.names(norm.Counts))
#colnames(Id) <- "Id"
# norm.counts$Id and row.names(*.df) must match
for (i in 1:length(pattern2)) {
  
  assign(x = pattern2[i],
         value = merge(as.data.frame(mget(pattern2[i]), col.names = NULL), norm.Counts, by = "row.names"), 
         envir = .GlobalEnv)
  
}

norm.Counts$Id <- rownames(norm.Counts)

for (i in 1:length(pattern2)) {
  
  assign(paste0(pattern2[i], ".join"), 
         value = full_join(as.data.frame(mget(pattern2)[i], col.names = NULL), 
                           norm.Counts, by = c("Row.names" = "Id"), envir = .GlobalEnv ))
  
}

pattern3 <- objects(pattern = ".join")

for (i in 1:length(pattern3)) {
  
  assign(paste0(pattern3[i], ".final"), 
         value = ((as.data.frame(mget(pattern3)[i], col.names = NULL)) %>%
                    select(Row.names, matches("norm"), everything())), 
         envir = .GlobalEnv )
}


obj <- objects(pattern = ".final")

for (i in 1:length(obj)) {
  assign(x = paste0(obj[i],".Annotated"),
         value = as.data.frame(mget(obj[i]),col.names = NULL) %>% 
           rename(!!paste0(strsplit(obj[i], "\\.")[[1]][1],".baseMean") := baseMean,
                  !!paste0(strsplit(obj[i], "\\.")[[1]][1],".log2FoldChange") := log2FoldChange,
                  !!paste0(strsplit(obj[i], "\\.")[[1]][1],".stat") := stat,
                  !!paste0(strsplit(obj[i], "\\.")[[1]][1],".pvalue") := pvalue,
                  !!paste0(strsplit(obj[i], "\\.")[[1]][1],".lfcSE") := lfcSE,
                  !!paste0(strsplit(obj[i], "\\.")[[1]][1],".padj") := padj),
         envir = .GlobalEnv
         
  )
}

obj2 <- objects(pattern = "Annotated")

new.data.frame <- as.data.frame(mget(obj2[1]),col.names = NULL)

for (i in 2:length(obj2)){
  
  new.data.frame <- full_join(new.data.frame, 
                              select(as.data.frame(mget(obj2[i]),col.names = NULL), Row.names, matches("_vs_")), 
                              by = "Row.names")
  
}

new.data.frame <- 
  new.data.frame %>%
  mutate_at(vars(matches("vs")), ~ round(.,3)) %>%
  mutate_at(vars(matches("norm")), ~ round(.,0))

write.table(new.data.frame, "final.txt", sep = "\t", quote = F, row.names = F)
