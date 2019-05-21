library(dplyr)

keys <- as.data.frame(t(phenoData))
colnames(keys) <- c("file", "group")


counts <- counts %>%
          mutate_at(1:24, ~ round(.,0))

for (i in 1:24) {
  
  assign(as.character(colnames(counts)[i]),
         value = counts %>% 
                  select(colnames(counts)[i]) ,
         envir = .GlobalEnv)
}

# custom 
hds <- objects(pattern = "HD")
ncs <- objects(pattern = "NC")

for (i in 1:length(hds)) {
  
  write.table(mget(hds)[i], file = paste0(hds[i], ".txt"), 
              sep = "\t", quote = F, col.names = F)
  
}

for (i in 1:length(ncs)) {
  
  write.table(mget(ncs)[i], file = paste0(ncs[i], ".txt"), 
              sep = "\t", quote = F, col.names = F)
  
}



