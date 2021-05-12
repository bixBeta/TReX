#!/bin/sh
#SBATCH -J ATACseq
#SBATCH -o %x.out
#SBATCH -n 12
#SBATCH --mem-per-cpu=18000

source ~/.bash_profile

usage(){

    echo "A T A C - S E Q   W O R K F L O W - @bixBeta"
    echo
    echo

    echo "Usage: bash" $0 "[-h arg] [-p arg] [-d args] [-t arg] [-g arg] [-q arg]"
    echo
    echo "---------------------------------------------------------------------------------------------------------------"
    echo "[-h] --> Display Help"
    echo "[-p] --> Project Identifier Number"
    echo "[-d] --> Comma Spearated Values for Delimiter and Field <delim,field or default> default: _,5 "
    echo "[-t] --> Trimming <nextseq or nova>;"
    echo "[-g] --> Reference Genome <mm10 or hg38>"
    echo "[-q] --> Execute atacQC.R script <yes>"
    echo "---------------------------------------------------------------------------------------------------------------"
}



declare -A genomeDir

genomeDir=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/BWAIndex/genome.fa" \
["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/bwa.index/Homo_sapiens.GRCh38.dna.toplevel.fa"
)

declare -A gtfs

gtfs=(
["mm10"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.96.gtf" \
["hg38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.96.gtf"
)

declare -A gAlias # for compatibility with atacQC.R
gAlias=(
["mm10"]="mouse" \
["hg38"]="human"
)

declare -A gSize # for macs2
gSize=(
["mm10"]="mm" \
["hg38"]="hs"
)




trimPE(){

        cd fastqs
        ls -1 *_R1.fastq* > .R1
        ls -1 *_R2.fastq* > .R2
        paste -d " " .R1 .R2 > Reads.list

        readarray fastqs < Reads.list
        mkdir fastQC

        for i in "${fastqs[@]}"
        do
                trim_galore --nextseq 20 --length 50  -j 8 --paired --gzip --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
        done

        mkdir TrimQC_stats trimmed_fastqs
        mv *_trimming_report.txt TrimQC_stats
        mv *_val* trimmed_fastqs
        mv TrimQC_stats fastQC trimmed_fastqs ../

        cd ..
}


trimHiSeqPE(){

                cd fastqs
                ls -1 *_1.fq* > .R1
                ls -1 *_2.fq* > .R2
                paste -d " " .R1 .R2 > Reads.list

                readarray fastqs < Reads.list
                mkdir fastQC

                for i in "${fastqs[@]}"
                do
                        trim_galore --quality 20 --gzip -j 8 --length 50  --paired --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
                done

                mkdir TrimQC_stats trimmed_fastqs
                mv *_trimming_report.txt TrimQC_stats
                mv *_val* trimmed_fastqs
                mv TrimQC_stats fastQC trimmed_fastqs ..

                cd ..
}




alignPE(){

        cd trimmed_fastqs

        ls -1 *_val_1.fq.gz > .trR1
        ls -1 *_val_2.fq.gz > .trR2
        paste -d " " .trR1 .trR2 > Trimmed.list

        readarray trimmedFastqs < Trimmed.list

        for i in "${trimmedFastqs[@]}"

        do
                # INDEX="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/BWAIndex/genome.fa"
                iSUB=`echo $i | cut -d ${DELIMITER} -f${FIELD}`

                bwa mem -t 24 -M -R "@RG\tID:${iSUB}\tSM:${iSUB}\tPL:ILLUMINA\tLB:${iSUB}\tPU:1" ${genomeDir[${DIR}]} $i \
                | samtools view -@ 24 -b -h -F 0x0100 -O BAM -o ${iSUB}.bam
        done

                mkdir primary-BAMS
                mv *.bam primary-BAMS
                mv primary-BAMS ..
                cd ..
}

sort(){
                cd primary-BAMS
            for i in *.bam
            do
            samtools sort $i > `echo  $i | cut -d "." -f1`.sorted.bam
            done

            for i in *.sorted.bam
            do
                samtools index $i
            done

                    # alignment stats etc. on raw bams
                    for i in *.sorted.bam
                    do
                        iSUB=`echo $i | cut -d "." -f1`
                        samtools flagstat $i > ${iSUB}.primary.flagstat
                        samtools idxstats $i > ${iSUB}.primary.idxstats
                    done

                cd ..
                pwd
}





rmMT(){
                cd primary-BAMS
                for i in *.sorted.bam
                do

                    iSUB=`echo $i | cut -d "." -f1`

                    samtools view -H `ls -1 *.sorted.bam | head -1` | cut -f2 | grep "SN:" |  cut -d ":" -f2 | grep -v "MT\|_\|\." | xargs samtools view -b $i > ${iSUB}.noMT.bam

                done
                cd ..
}

markDups(){
                cd primary-BAMS
        for i in *.noMT.bam
        do
            iSUB=`echo $i | cut -d "." -f1`
            java -jar /programs/bin/picard-tools/picard.jar \
            MarkDuplicates \
            INPUT=$i \
            OUTPUT=${iSUB}.dupMarked.noMT.bam \
            ASSUME_SORTED=true \
            REMOVE_DUPLICATES=false \
            METRICS_FILE=${iSUB}.MarkDuplicates.metrics.txt \
            VALIDATION_STRINGENCY=LENIENT \
            TMP_DIR=tmp

        done
                cd ..
}

dedupBAM(){
                cd primary-BAMS
                # alignment stats etc. on dupMarked no MT bams
                for i in *.dupMarked.noMT.bam
                do
                    iSUB=`echo $i | cut -d "." -f1`
                    samtools index $i
                    samtools flagstat $i > ${iSUB}.noMT.flagstat
                    samtools idxstats $i > ${iSUB}.noMT.idxstats
                done

        for i in *.dupMarked.noMT.bam
        do
                iSUB=`echo $i | cut -d "." -f1`
                samtools view -b -h -F 0X400 $i > ${iSUB}.DEDUP.bam
        done

                for i in *.DEDUP.bam; do samtools index $i ; samtools idxstats $i > `echo $i | cut -d "." -f1`.DEDUP.idxstats; done
                for i in *.DEDUP.bam; do samtools flagstat $i > `echo $i | cut -d "." -f1`.DEDUP.flagstat; done

                multiqc -n ${PIN}.multiqc.report .

                mkdir dedup-BAMS
                mv *.DEDUP* dedup-BAMS/
                mv dedup-BAMS ..
                cd ..

}

tagDir(){
  cd dedup-BAMS
  for i in *.DEDUP.bam
  do
  iSUB=`echo "$i" | cut -d'.' -f1` # subset to rename
  /home/fa286/bin/HOMER/bin/makeTagDirectory "$iSUB".tag.dir "$i"
  done
  cd ..
}

callPeak(){

  cd dedup-BAMS
    echo "calling peaks on DEDUP bams"
    mkdir peaks.OUT
    for i  in *.DEDUP.bam
      do
        iSUB=`echo $i | cut -d "." -f1`
        macs2 callpeak -t $i \
        -f BAMPE \
        -n ${iSUB} \
        -g ${gSize[${DIR}]} \
        -q 0.05 \
        --outdir peaks.OUT \
        --nomodel --shift 37 --ext 73 \
        --keep-dup all
      done

  cd ..

}


mergedPeaks(){

  echo "running module mergedPeaks"
  cd dedup-BAMS

    allBams=`echo *.DEDUP.bam`

    macs2 callpeak -t ${allBams} \
    -f BAMPE \
    -n allSamplesMergedPeakset \
    -g ${gSize[${DIR}]} \
    -q 0.05 \
    --outdir peaks.OUT \
    --nomodel --shift 37 --ext 73 \
    --keep-dup all

  cd ..
}


saf(){
  # awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' ${sample}_peaks.narrowPeak > ${sample}_peaks.saf
  cd dedup-BAMS/peaks.OUT
  awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $4, $1, $2+1, $3, "."}' allSamplesMergedPeakset_peaks.narrowPeak > allSamplesMergedPeakset.saf
  cd ../../
}

#saf

frip(){
  # featureCounts -p -a ${sample}_peaks.saf -F SAF -o readCountInPeaks.txt ${sample}.sorted.marked.filtered.shifted.bam
  cd dedup-BAMS/
    for i  in *.DEDUP.bam
    do
      iSUB=`echo $i | cut -d "." -f1`
      featureCounts -p -a peaks.OUT/allSamplesMergedPeakset.saf -F SAF -o "${iSUB}".readCountInPeaks.txt $i
    done

  cd ..
}

annotatePeaks(){

    cd dedup-BAMS/peaks.OUT
    /home/fa286/bin/HOMER/bin/annotatePeaks.pl allSamplesMergedPeakset.saf ${DIR} -gtf ${gtfs[${DIR}]} > allSamplesMergedPeakset.Annotated.saf
    cd ../..
}

bedGraphs(){
  cd dedup-BAMS
    for i in *.tag.dir
    do
        makeUCSCfile ${i} -o auto -fsize 1e10 -res 1 -color 106,42,73 -style chipseq
    done

    mkdir tagDirs
        mv *.tag.dir tagDirs
        cd tagDirs
        mkdir bedGraphs
            for i in *.tag.dir
            do
                cd $i
                zcat *.ucsc.bedGraph.gz | awk '{if(NR>1) print "chr"$0; else print $0}' | gzip > `basename *.ucsc.bedGraph.gz .ucsc.bedGraph.gz`.ucsc.bg.gz
                mv *.ucsc.bg.gz ../bedGraphs
                cd ..
            done
        cd ..

    mkdir featureCounts
    mv *.txt featureCounts

    multiqc -n ${PIN}.FRIP.multiqc.report -b "Please note that the featureCounts M Assigned Column refers to Fragments and Not Reads" --ignore tagDirs --ignore peaks.OUT .

  cd ..
}

atacQC(){

    cd dedup-BAMS
    echo "genome alias" = ${gAlias[${DIR}]}
    /programs/R-3.6.3/bin/Rscript /home/fa286/bin/scripts/atacQC.R ${gAlias[${DIR}]}
    # ${gAlias[${DIR}]}
    ~/bin/scripts/html.atacQC.sh `echo ${PIN}_atacQC`

    cd ..

    /home/fa286/bin/tree-1.7.0/tree > folder.structure

}


while getopts "hp:t:g:q:d:" opt; do
    case ${opt} in

    h)
        echo
        echo
        echo
        usage
        echo
        echo
        exit 1

    ;;

    p )

        PIN=$OPTARG
        echo "Project Identifier = " $PIN
    ;;

    t )

        T=$OPTARG

    ;;

    g)

        DIR=$OPTARG

    ;;

    q)

        QC=$OPTARG

    ;;

    d)
        DELIM=$OPTARG

    ;;


    \? )
        echo
        echo
        echo
        usage

    ;;

    esac

done
# shift $((OPTIND -1))

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if PIN is provided

if [[ -z "${PIN+x}" ]]; then

    PIN="PIN_Null"
fi

# PARAMETER CHECKS

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if delimiter parameter exists
                    if [[ ! -z "${DELIM+x}" ]]; then
                        #statements
                        if [[ $DELIM == default ]]; then

                        DELIMITER="_"
                        FIELD="5"
                        echo "file naming will be done using the default delimiter settings"
                      else

                        DELIMITER=`echo $DELIM | cut -d , -f1`
                        FIELD=`echo $DELIM | cut -d , -f2-`
                        echo "file naming will be done using the delim = $DELIMITER and field = $FIELD settings"

                      fi

                    fi

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if trimming parameter exists

                    if [[ ! -z "${T+x}" ]]; then

                        if [[ $T == nextseq ]]; then
                            trimPE
                        elif [[ $T == nova ]]; then
                            trimHiSeqPE
                        else
                        echo "-t only accepts nextseq or nova as arguments"
                        exit 1

                        fi
                    fi

                    #-------------------------------------------------------------------------------------------------------------
                    #-------------------------------------------------------------------------------------------------------------
                    ## check if genomeDir provided

                    if [[ ! -z "${DIR+x}" ]]; then
                        if [ ${genomeDir[${DIR}]+_} ]; then
                            echo Reference genome selected = $DIR
                            echo
                            alignPE
                            sort
                            rmMT
                            markDups
                            dedupBAM
                            callPeak
                            mergedPeaks
                            saf
                            frip
                            tagDir
                            annotatePeaks
                            bedGraphs
                        else
                            echo "The reference genome provided '"$DIR"' is not available"
                            exit 1

                        fi
                    fi

                    if [[ ! -z "${QC+x}" ]]; then

                        if [[ $QC == yes ]]; then
                            atacQC
                        else
                            echo "-q option only accepts yes as an argument"
                            exit 1
                        fi
                    fi






#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------

if [[ -z $1 ]] || [[  $1 = "help"  ]] ; then
    #statements
    echo
    echo
    usage
    echo
    echo
    exit 1

else
    echo >> beta5.atac.log
    echo `date -u` >> beta5.atac.log
    echo "Project Identifier Specified = " $PIN >> beta5.atac.log
    echo "Reference Genome Specified   = " $DIR >> beta5.atac.log
    echo "Trimming                     = " $T >> beta5.atac.log
    echo >> beta5.atac.log

    echo "ENV INFO: " >> beta5.atac.log
    echo >> beta5.atac.log
    echo "STAR version:" `~/bin/STAR-2.7.0e/bin/Linux_x86_64/STAR --version` >> beta5.atac.log
    echo "multiqc version:" `~/miniconda2/envs/RSC/bin/multiqc --version` >> beta5.atac.log
    echo "samtools version:" `/programs/bin/samtools/samtools --version` >> beta5.atac.log
    echo "macs2 version: macs2 2.1.0.20150731 " >> beta5.atac.log
    echo "HOMER version: v4.10.4" >> beta5.atac.log
    echo -------------------------------------------------------------------------------------------------- >> beta5.atac.log

fi
