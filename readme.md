beta4.sh  

```
Usage: bash beta4.sh [-h] [-p arg] [-t arg] [-g arg] [-r arg] [-s arg] [-c arg] 
---------------------------------------------------------------------------------------------------------------
[-h] --> Display Help
[-p] --> Project Identifier Number
[-t] --> Small RNA Trimming <yes, no or paired>
[-g] --> Reference Genome < hg38, GRCh38, mm10, GRCm38, rat, cat, chicken, horse, ATCC_13047, grape, ercc >
[-r] --> <SE> or <PE> 
[-s] --> Library Strandedness < 0, 1, 2 > where 1 = first strand, 2 = reverse strand, 0 for unstranded counts 
[-c] --> GeneBody Coverage < yes, no > 
---------------------------------------------------------------------------------------------------------------
```



processSAR.R  -> Rscript processSAR.R <path/to/SAR/tables/> 
smRNA.beta2.sh -> Current TREx workflow for smRNA 
geneSwitch.R  -> For gff annotations use geneSwitch.R to replace the gene0, gene1, gene2 ... etc. 
naming schema with appropriate gene names
