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

```
Usage: bash smRNA.beta3.sh [-h arg] [-p arg] [-t arg] [-g arg]

---------------------------------------------------------------------------------------------------------------
[-h] --> Display Help 
[-p] --> Project Identifier Number 
[-t] --> NextSeq run < yes, no, na > 
[-g] --> Mapper Genome < hsa, mmu, cel > 
---------------------------------------------------------------------------------------------------------------
```
### ATAC-seq Workflow

```
Usage: bash 3.1.sh [-h arg] [-p arg] [-t arg] [-g arg]

---------------------------------------------------------------------------------------------------------------
[-h] --> Display Help
[-p] --> Project Identifier Number
[-t] --> Trimming <yes>; only use it if trimming is required
[-g] --> Reference Genome <mm10 or hg38>
---------------------------------------------------------------------------------------------------------------
``` 
```Rscript atacQC.R < human, mouse  or 'path to gtf annotation' >```

> Dir Structure

    .
    ├── *.CLEAN.bam                     # clean bams
    ├── *.bam.bai                       # samtools index *.CLEAN.bam
    ├── peaks.OUT                       # macs2 peak files 
    └── README.md


#### R script to generate mega DE-reults table
processSAR.R  
`Usage: Rscript processSAR.R <path/to/SAR/tables/>`

geneSwitch.R <br>
`For gff annotations use geneSwitch.R to replace the gene0, gene1, gene2 ... etc. naming schema with appropriate gene names`
