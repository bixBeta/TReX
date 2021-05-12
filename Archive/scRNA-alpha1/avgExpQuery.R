geneExpression <- function(gene){
  
  queryExpression <- norm_1.1[rownames(norm_1.1) == gene, ]
  query <- as.data.frame(t(queryExpression))
  query$UMI <- rownames(query)
  jnd <- left_join(myMeta, query, by = "UMI")
  jnd <- jnd %>% 
    group_by(Phenotype, CPET.Day, seurat_clusters) %>%
    summarise_each(funs = (mean), !!gene)
  #print(jnd)
}

hh <- as.data.frame(geneExpression("GRASP"))
hh$PhenoDay <- paste0(hh$Phenotype,"-",hh$CPET.Day)

jj <- melt(hh)

library(ggplot2)
library(ggpubr)
library(ggplot2)
library(EnvStats)
library(ggbeeswarm)
library(viridis)


ggplot(jj, aes(seurat_clusters, value, fill=Phenotype)) +
geom_boxplot()+
geom_quasirandom(aes(color= Phenotype),alpha = 0.75,dodge.width=0.75)+
  theme_bw() + scale_color_manual(values=c(viridis(4)))+
  scale_fill_manual(values=c(viridis(4)))



ggplot(jj, aes(seurat_clusters, value)) +
  geom_boxplot(aes(fill = Phenotype))+
  geom_quasirandom(aes(color= PhenoDay),alpha = 0.75,dodge.width=0.75)+
  theme_bw() + scale_color_manual(values=c(viridis(5)))+
  scale_fill_manual(values=c(viridis(5))) + ylab("Mean reads per cell")


