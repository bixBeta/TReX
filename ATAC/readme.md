#### 3.4.sh <- ATAC-seq workflow
Usage:

`nohup bash 3.4.sh "[-h arg] [-p arg] [-d args] [-t arg] [-g arg] [-q arg]" > 3.4.log &`

Params Description:

```bash
"[-h] --> Display Help"
"[-p] --> Project Identifier Number"
"[-d] --> Comma Spearated Values for Delimiter and Field <delim,field or default> default: _,5 "
"[-t] --> Trimming <nextseq or nova>;"
"[-g] --> Reference Genome <mm10 or hg38>"
"[-q] --> Execute atacQC.R script <yes>"
```

- Create a folder named fastqs and add all the PE fastq files to this folder
- Run the [3.4.sh](https://raw.githubusercontent.com/bixBeta/TReX/master/ATAC/3.4.sh) script from the same dir as the fastqs directory

```bash
.
├── **3.4.sh**
├── dedup-BAMS
│   ├── PIN.FRIP.multiqc.report_data
│   ├── atacQC.out
│   ├── featureCounts
│   ├── peaks.OUT
│   └── tagDirs
├── fastQC
├── **fastqs**
├── primary-BAMS
├── trimmed_fastqs
└── TrimQC_stats
```

Tools needed:

1. trim_galore
2. fastqc
3. bwa
4. samtools
5. picard tools
6. homer suite (w/ human and mouse genome config)
7. macs2
8. featureCounts (rsubread)
9. multiqc
