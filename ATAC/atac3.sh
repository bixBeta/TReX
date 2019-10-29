#!/bin/sh

source ~/.bash_profile

usage(){

	echo "A T A C - S E Q   W O R K F L O W - @bixBeta"
	echo
	echo

	echo "Usage: bash" $0 "[-h arg] [-p arg] [-t arg] [-g arg]"
	echo
	echo "---------------------------------------------------------------------------------------------------------------"
	echo "[-h] --> Display Help"
	echo "[-p] --> Project Identifier Number"
	echo "[-t] --> Trimming <yes>; only use it if trimming is required"
	echo "[-g] --> Reference Genome <GRCm38 or GRCh38>"
	echo "---------------------------------------------------------------------------------------------------------------"
}



declare -A genomeDir

genomeDir=(
["GRCm38"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/BWAIndex/genome.fa" \
["GRCh38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/bwa.index/Homo_sapiens.GRCh38.dna.toplevel.fa"
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
                $TRIM --nextseq 20 --length 20  --paired --gzip --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
        done

        mkdir TrimQC_stats trimmed_fastqs
        mv *_trimming_report.txt TrimQC_stats
        mv *_val* trimmed_fastqs
        mv TrimQC_stats fastQC trimmed_fastqs ..

        cd ..
}


alignPE(){

        cd trimmed_fastqs

        ls -1 *_R1_val_1.fq.gz > .trR1
      	ls -1 *_R2_val_2.fq.gz > .trR2
      	paste -d " " .trR1 .trR2 > Trimmed.list

      	readarray trimmedFastqs < Trimmed.list

        for i in "${trimmedFastqs[@]}"

        do
                # INDEX="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/BWAIndex/genome.fa"
                iSUB=`echo $i | cut -d "_" -f5`

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
					samtools idxstats $i | cut -f1 | grep -v MT | xargs samtools view -b $i > ${iSUB}.noMT.bam

				done
				cd ..
}

markDups(){
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
}

cleanBAM(){
				# alignment stats etc. on dupMarked no MT bams
				for i in *.dupMarked.noMT.bam
				do
					iSUB=`echo $i | cut -d "." -f1`
					samtools index $i
					samtools flagstat $i > ${iSUB}.secondary.flagstat
					samtools idxstats $i > ${iSUB}.secondary.idxstats
				done

        for i in *.dupMarked.noMT.bam
        do
                iSUB=`echo $i | cut -d "." -f1`
                samtools view -b -h -F 0X400 $i > ${iSUB}.CLEAN.bam
        done
}



while getopts "hp:t:g:" opt; do
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
					## check if trimming parameter exists

					if [[ ! -z "${T+x}" ]]; then

						if [[ $T == yes ]]; then
							trimPE
						else
							echo "-t option only accepts yes as an argument"
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
							#alignPE
							#sort
							rmMT
							#markDups
							#cleanBAM
						else
							echo "The reference genome provided '"$DIR"' is not available"
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
	echo >> beta3.atac.log
	echo `date -u` >> beta3.atac.log
	echo "Project Identifier Specified = " $PIN >> beta3.atac.log
	echo "Reference Genome Specified   = " $DIR >> beta3.atac.log
	echo "Trimming                     = " $T >> beta3.atac.log
	echo >> beta3.atac.log

	echo "ENV INFO: " >> beta3.atac.log
	echo >> beta3.atac.log
	echo "STAR version:" `~/bin/STAR-2.7.0e/bin/Linux_x86_64/STAR --version` >> beta3.atac.log
	echo "multiqc version:" `~/miniconda2/envs/RSC/bin/multiqc --version` >> beta3.atac.log
	echo "samtools version:" `/programs/bin/samtools/samtools --version` >> beta3.atac.log
	echo "rseqc version:" `/programs/MACS2-2.1.0/bin/macs2 --version` >> beta3.atac.log
	echo "HOMER version: v4.10.4" >> beta3.atac.log
	echo -------------------------------------------------------------------------------------------------- >> beta3.atac.log

fi
