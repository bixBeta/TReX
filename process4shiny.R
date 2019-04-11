library(DESeq2)
library(dplyr)
## -------------------------------------------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = target,
                              design = ~ group)

## -------------------------------------------------------------------------------------------------------------------
dds <- DESeq(dds)
#resultsNames(dds)

## -------------------------------------------------------------------------------------------------------------------
# add .custom to all contrast objects 

cColdHot.custom <- results(dds, contrast=c("group", "cHOT", "cCOLD"), alpha = 0.05)      
pColdHot.custom <- results(dds, contrast=c("group", "pHOT", "pCOLD"), alpha = 0.05)
cpHot.custom <- results(dds, contrast=c("group", "cHOT", "pHOT"), alpha = 0.05)
cpCold.custom <- results(dds, contrast=c("group", "cCOLD", "pCOLD"), alpha = 0.05)

pattern1 <- objects(pattern = ".custom")

for (i in 1:length(pattern1)) {
  
  print(paste0("Res Summary - ", pattern1[i], ":"))
  summary.DESeqResults(get(pattern1[i]))
  
}

## -------------------------------------------------------------------------------------------------------------------
norm.Counts <- as.data.frame(counts(dds, normalized = T))

norm.Counts <- rename_all(norm.Counts, .funs = ~ paste0("norm.", colnames(norm.Counts)))
  
norm.Counts$Id <- rownames(norm.Counts)

      # y <- as.data.frame(cColdHot.custom)
      # y$Id <- row.names(y)
      # joined <- full_join(norm.Counts, y, by= "Id")

## -------------------------------------------------------------------------------------------------------------------
# coercing DESeqResult object to data.frame
for (i in 1:length(pattern1)) {
    assign(x = paste0(pattern1[i], ".df"), 
    value = as.data.frame(get(pattern1[i])),
    envir = .GlobalEnv
    )
  }

pattern2 <- objects(pattern = ".df")

Id <- as.data.frame(row.names(norm.Counts))
colnames(Id) <- "Id"

for (i in 1:length(pattern2)) {
  
  assign(x = pattern2[i],
        value = cbind(as.data.frame(mget(pattern2)[i], col.names = NULL),Id), 
        envir = .GlobalEnv)

 }

for (i in 1:length(pattern2)) {
  
  assign(paste0(pattern2[i], ".join"), 
         value = full_join(as.data.frame(mget(pattern2)[i], col.names = NULL), 
                           norm.Counts, by = "Id"), envir = .GlobalEnv )
  
 }

pattern3 <- objects(pattern = ".join")

for (i in 1:length(pattern3)) {
  
  assign(paste0(pattern3[i], ".final"), 
         value = ((as.data.frame(mget(pattern3)[i], col.names = NULL)) %>%
                          select(Id, matches("norm"), everything())), 
         envir = .GlobalEnv )
 }
