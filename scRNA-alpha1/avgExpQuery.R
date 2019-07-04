geneExpression <- function(gene){
  
  queryExpression <- norm_1.1[rownames(norm_1.1) == gene, ]
  query <- as.data.frame(t(queryExpression))
  query$UMI <- rownames(query)
  jnd <- left_join(myMeta, query, by = "UMI")
  jnd <- jnd %>% 
    group_by(Phenotype, CPET.Day) %>%
    summarise_each(funs = (mean), !!gene)
  print(jnd)
}
geneExpression("GRASP")
