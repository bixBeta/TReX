#!/bin/sh

#  1007.sh
#  RSC-alpha
#
#  Created by Faraz Ahmed on 3/29/19.
#  sample workflow for secondary and tertiary QC

start=`date +%s`
source activate RSC
source ~/.bash_profile

PIN='1007' # example
BED12='/path/to/bed12'

mkdir fastQC

for i in *.gz

do

$TRIM --nextseq 20 --gzip --length 50  --fastqc --fastqc_args "--outdir ./fastQC" $i

done


mkdir trimmed-fastq TrimQC_Stats
mv *trimmed.fq.gz trimmed-fastq
mv *trimming_report.txt TrimQC_Stats

cd trimmed-fastq

    for i in *_trimmed.fq.gz

    do
        DIR='/path/to/starIndex/'
        iSUB=`echo $i | cut -d '.' -f1`

        STAR \
        --runThreadN 12 \
        --genomeDir $DIR \
        --readFilesIn $i \
        --readFilesCommand gunzip -c \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix $iSUB. \
        --limitBAMsortRAM 61675612266 \
        --quantMode GeneCounts

    done



    multiqc -f -n ${PIN}.star.multiqc.report .

    mkdir STAR.COUNTS STAR.BAMS STAR.LOGS
    mv *.ReadsPerGene.out.tab STAR.COUNTS
    mv *.bam STAR.BAMS
    mv *.out *.tab *_STARtmp *.list STAR.multiqc.report_data STAR.LOGS

        cd STAR.BAMS
            for i in *.bam

            do
            samtools index -b $i

            done
        cd ..


    geneBody_coverage.py -r $BED12 -i STAR.BAMS/ -o ${PIN}
    mkdir geneBodyCov
    mv *geneBodyCoverage.* log.txt geneBodyCov


cd ..

end=`date +%s`
runtime=$((end-start))
echo 'runtime seconds '"$runtime" > runtime.${PIN}.txt



# change PIN, DIR and BED12 for each run!
