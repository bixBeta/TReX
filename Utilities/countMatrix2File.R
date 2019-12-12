library(dplyr)

# keys <- as.data.frame(t(phenoData))
# colnames(keys) <- c("file", "group")

# arg <-  commandArgs(trailingOnly = T)
# counts <- read.table(arg[1], header = T, sep = "\t")

counts <- read.table("CPM_run1.txt", header = T, sep = "\t")

counts <- counts %>%
           mutate_at(2:ncol(counts), ~ round(.,0))
 
row.names(counts) <- counts$X

counts <- counts %>% 
  select(-1)

for (i in 1:ncol(counts)) {

  assign(as.character(colnames(counts)[i]),
         value = counts %>%
                  select(colnames(counts)[i]) ,
         envir = .GlobalEnv)
}

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



