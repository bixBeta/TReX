```Rscript atacQC.R < human, mouse  or 'path to gtf annotation' >```

### Folder Layout
    .
    ├── *.CLEAN.bam                     # clean bams
    ├── *.bam.bai                       # `samtools index` *.CLEAN.bam
    ├── peaks.OUT                       # peak output (`macs2`)
    └── README.md
