### RNA-seq Differential Gene Expression Workflow 

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
### smRNA-seq Workflow
smRNA.beta3.sh 

```
Usage: bash smRNA.beta3.sh [-h arg] [-p arg] [-t arg] [-g arg]

---------------------------------------------------------------------------------------------------------------
[-h] --> Display Help 
[-p] --> Project Identifier Number 
[-t] --> NextSeq run < yes, no, na > 
[-g] --> Mapper Genome < hsa, mmu, cel > 
---------------------------------------------------------------------------------------------------------------
```


#### R script to generate mega DE-reults table
processSAR.R  
`Usage: Rscript processSAR.R <path/to/SAR/tables/>`

geneSwitch.R 
`For gff annotations use geneSwitch.R to replace the gene0, gene1, gene2 ... etc. naming schema with appropriate gene names`
