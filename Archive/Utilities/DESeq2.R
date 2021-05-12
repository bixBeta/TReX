#!/usr/bin/env Rscript
# written by Faraz Ahmed 
# initial commit 11/15/2018
## -------------------------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------------------------

		# USAGE SYNTAX:
        # $ Rscript DESeq2.R <Experiment.Name> <Numerator> <Denominator>

## -------------------------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------------------------

suppressWarnings(library("dplyr"))
suppressWarnings(library("DESeq2"))
## -------------------------------------------------------------------------------------------------------------------
args <-  commandArgs(trailingOnly = T)
outputPrefix <- args[1] # Experiment Name e.g. Morgan.Sammons.0002
numerator   <- args[2]
denominator <- args[3]


## -------------------------------------------------------------------------------------------------------------------

phenoData <- read.csv(file = "phenoData.csv", header = T)
countMatrix <- read.table(file = "countMatrix.txt", header = F)


## -------------------------------------------------------------------------------------------------------------------

# --- storing meta info -------------------------

sampleNames <- phenoData$sampleName
condition <- phenoData$condition
myDim <- ncol(countMatrix)/2
rawCounts <- countMatrix

n <- (0:myDim)*2
n[1] <- 1

countTable <- rawCounts[,n]
rownames(countTable) <- countTable[,1]
countTable <- select(countTable, -V1)
colnames(countTable) <- phenoData$sampleName

# --- final check ------------------------------
phenoData$sampleName == colnames(countTable)



# FOR HTSEQ  -- CHECK 
#grep("__",x = row.names(countTable))
countTable[grep("__",x = row.names(countTable)),] # uncomment line 56 if line 54 is True

#countTable <- countTable[1:(nrow(countTable) - 5),] # remove last five rows contatining htseq stats
#countTable[grep("__",x = row.names(countTable)),] # verify

print(" ", quote = F)
print(" ", quote = F)
print("                           RUNNING DESEQ2                              ", quote = F)
print(" ", quote = F)
print(" ", quote = F)

## -------------------------------------------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = phenoData,
                              design = ~ condition)

treatments = unique(phenoData$condition) # Treatments of interest
dds$condition <- factor(colData(dds)$condition,
                        levels = treatments)

## -------------------------------------------------------------------------------------------------------------------
dds <- DESeq(dds)
resultsNames(dds)

## -------------------------------------------------------------------------------------------------------------------
# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

# save normalized values
write.table(as.data.frame(assay(rld)),file = paste0(outputPrefix, ".rlog-transformed-counts.txt"), sep = '\t')
write.table(as.data.frame(assay(vsd)),file = paste0(outputPrefix, ".vst-transformed-counts.txt"), sep = '\t')
system("mkdir normCounts")
system("mv *rlog-transformed-counts.txt *vst-transformed-counts.txt normCounts")

## -------------------------------------------------------------------------------------------------------------------
## use contrast to drop deseq2 results for specific conditions; following are some expamples

# WD.WN <- results(dds, contrast=c("condition","WD","WN"), alpha = 0.05)
# pd.pn <- results(dds, contrast=c("condition","PD","PN"), alpha = 0.05)
# Ad.An <- results(dds, contrast=c("condition","AD","ANU"), alpha = 0.05)
# write.csv(WD.WN, file = paste0(outputPrefix, ".WD.WN.results.csv"))
# write.csv(pd.pn, file = paste0(outputPrefix, ".pd.pn.results.csv"))
# write.csv(Ad.An, file = paste0(outputPrefix, ".Ad.An.results.csv"))
# write.csv(WD.pd, file = paste0(outputPrefix, ".WD.pd.results.csv"))


custom <- results(dds, contrast=c("condition", numerator, denominator), alpha = 0.05)
mcols(custom, use.names=TRUE) # brief description of res headers
summary(custom)

save(countTable, phenoData, dds, rawCounts, rld, vsd, custom, file = paste0(outputPrefix, ".Rdata"))
write.csv(custom, file = paste0(outputPrefix, ".", numerator, ".", denominator, ".results.csv"))
system("mkdir DGE_results")
system("mv *.results.csv DGE_results")

system("mkdir DESeq2_output")
system("mv normCounts DGE_results *.Rdata DESeq2_output")


custom

print(" ", quote = F)
print(" ", quote = F)
print("                              DONE   :)                          ", quote = F)
print(" ", quote = F)
print(" ", quote = F)



