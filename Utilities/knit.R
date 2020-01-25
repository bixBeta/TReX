#!/usr/bin/env Rscript

args 	 <-  commandArgs(trailingOnly = T)


title 	 <- args[1]
genome	 <- args[2]
annot 	 <- args [3]
markdown <- list.files(pattern = ".Rmd")

rmarkdown::render(markdown, params = list(genome = genome, dynamictitle=title, annot = annot),
                  output_file = paste0(title, ".html"), quiet = T)


print("HTML report created!")


