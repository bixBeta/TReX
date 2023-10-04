runmode = params.mode

process STARM {

    tag "$id"
    label 'process_high'
    
    // publishDir "$baseDir/STAR_OUT", mode: "copy", overwrite: false
    
    input:
        tuple val(id), path(trimmed)
        path genome

    output:
        path "*ReadsPerGene.out.tab"            , emit: read_per_gene_tab 
        path "*Log.final.out"                   , emit: log_final
        path "*Log.out"                         , emit: log_out
        path "*Log.progress.out"                , emit: log_progress
        path "*SJ.out.tab"                      , emit: sj_out_tab
        path "*.out.mate*"                      , emit: unmapped
        path "*bam"                             , emit: bam_sorted

    script:

    if (runmode == "SE")
        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome} \
            --readFilesIn ${trimmed} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

        """
    else if (runmode == "SES")

        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome} \
            --readFilesIn ${trimmed} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts
        
        gzip *.out.mate1 

        """

    else if (runmode == "SEBS")

        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir  ${genome} \
            --readFilesIn ${trimmed} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000 

        gzip *.out.mate1

        """


    else if (params.mode == "PE")
       
        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome} \
            --readFilesIn ${trimmed[0]} ${trimmed[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts

        """

    else if (params.mode == "PES")
   
        """
            STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome} \
            --readFilesIn ${trimmed[0]} ${trimmed[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx 

        gzip *.out.mate1
        gzip *.out.mate2

        """
    else if (params.mode == "PEBS")

        """
             STAR \
            --runThreadN ${task.cpus} \
            --genomeDir ${genome} \
            --readFilesIn ${trimmed[0]} ${trimmed[1]} \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix ${id}. \
            --limitBAMsortRAM 61675612266 \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx \
            --alignIntronMax 1 \
            --alignMatesGapMax 45000                   

        gzip *.out.mate1
        gzip *.out.mate2

        """

    else {
        error "Invalid alignment mode: ${runmode}"
        exit 0
    } 
}



