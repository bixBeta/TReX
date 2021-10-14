#!/usr/bin/env Rscript
args <-  commandArgs(trailingOnly = T)
pin <- args[1]
ref <- args[2]

if (length(args)<=1) {
  print(" Usage = Rscript test.R <pin> <base-line>")  
  stop("Both arguments must be supplied!!! \n", call.=FALSE)
  
} 

library(SARTools)

workDir <- getwd()      								# working directory for the R session

projectName <- pin                      				# name of the project
author <- "RSC"                                			# author of the statistical analysis/report

targetFile <- "targetFile.txt"                           # path to the design/target file
rawDir <- getwd()                                      # path to the directory containing raw counts files
#featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
#                     "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
#                    "not_aligned", "too_low_aQual")# NULL if no feature to remove

featuresToRemove <- NULL
varInt <- "group"                                    # factor of interest
condRef <- ref                                    # reference biological condition




batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default), "local" or "mean"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors



#colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
#            "MediumVioletRed","SpringGreen", c2)



# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

library("viridisLite")
#colors <- viridisLite::viridis(n = nrow(target), option = "inferno", begin = 0, end = 1 )

colors<-  c("#EF8A62",
	     "#1f78b4",
	     "#1b9e77",
	     "purple3", 	      
	     "khaki4",
	     "#E9A3C9",
	     "#A1D76A",
	     "#FFFF33",
	     "grey", 
	     "#b3e2cd",
	     "#67A9CF",
	     "peachpuff2",
	     "red",
	     "magenta3",
	     "blue",
	     "yellow"
)

#colors <- c(c1,c2)
forceCairoGraph <- FALSE

################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)
library(SARTools)
if (forceCairoGraph) options(bitmapType="cairo")

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
# target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)



################################################

system(paste0("/Users/fa286/bin/generateRaw.R ", projectName, " *.RData"))

system(paste("mkdir",projectName ))
system(paste("mv figures *.html *.RData tables *.txt", projectName))

setwd(paste0(projectName, "/tables"))
#system(paste("pwd"))
system(paste("/Users/fa286/bin/vs2_vs_.sh"))

#setwd(projectName)
#system(paste0("/Users/fa286/bin/generateRaw.R ../ ", projectName, " ../*.RData"))
