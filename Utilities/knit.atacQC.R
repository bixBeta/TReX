#!/usr/bin/env Rscript

args 	 <-  commandArgs(trailingOnly = T)


title 	 <- args[1]
markdown <- list.files(pattern=".Rmd")

rmarkdown::render(markdown, output_file = paste0(title, ".atacQC.html"), quiet = T)

print("HTML report created!")




