#!/bin/sh

source ~/.bash_profile

usage(){

	echo "Usage: bash" $0 "[-h arg] [-p arg] [-t arg] [-g arg] [-s arg] "
	echo
	echo "---------------------------------------------------------------------------------------------------------------"
	echo "[-h] --> Display Help"
	echo "[-p] --> Project Identifier Number"
	echo "[-t] --> Small RNA Trimming <yes, no>"
	echo "[-g] --> Reference Genome < hg38, GRCh38, mm10, GRCm38, cat, chicken, horse >"
	echo "[-s] --> Library Strandedness < 0, 1, 2 > where 1 = first strand, 2 = reverse strand, 0 for unstranded counts "
	echo "---------------------------------------------------------------------------------------------------------------"
}



trim(){

		mkdir TrimQC_stats fastQC trimmed_fastqs
		for i in fastqs/*.gz
		do
			$TRIM --nextseq 20 --gzip --length 50  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i 
		done
		mv *_trimming_report.txt TrimQC_stats
		mv *trimmed.fq.gz trimmed_fastqs

}

trimSmall(){

		mkdir TrimQC_stats fastQC trimmed_fastqs
		for i in fastqs/*.gz
		do
			$TRIM --nextseq 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i 
		done
		mv *_trimming_report.txt TrimQC_stats
		mv *trimmed.fq.gz trimmed_fastqs

}

declare -A genomeDir

genomeDir=( ["hg38"]="/workdir/genomes/Homo_sapiens/hg38/UCSC/hg38.star" \
	["mm10"]="/workdir/genomes/Mus_musculus/mm10/UCSC/mm10.star" \
	["GRCh38"]="/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/GRCh38.star" \
	["GRCm38"]="/workdir/genomes/Mus_musculus/mm10/ENSEMBL/GRCm38.star" \
	["cat"]="/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/genomeDir" \
	["chicken"]="/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/galgal5.star" \
	["horse"]="/workdir/genomes/Equus_caballus/ENSEMBL/Equus_caballus.star"  )


declare -A bed12

bed12=( ["hg38"]="/workdir/genomes/Homo_sapiens/hg38/UCSC/genes.bed12" \
	["mm10"]="/workdir/genomes/Mus_musculus/mm10/UCSC/BED12/mm10.ucsc.bed12" \
	["GRCh38"]=" " \
	["GRCm38"]=" " \
	["cat"]="/workdir/genomes/Felis_catus/Felis_catus9.0/Ensembl/Felis_catus.Felis_catus_9.0.95.bed12" \
	["chicken"]="/workdir/genomes/Gallus_gallus/Galgal5/ENSEMBL/Gallus_gallus.Gallus_gallus-5.0.bed12" \
	["horse"]="/workdir/genomes/Equus_caballus/ENSEMBL/Equus_caballus.EquCab3.0.96.bed12"  )



align(){

	cd trimmed_fastqs

	for i in *_trimmed.fq.gz

	    do

	        iSUB=`echo $i | cut -d '-' -f3-`

	        STAR \
	        --runThreadN 12 \
	        --genomeDir ${genomeDir[${DIR}]} \
	        --readFilesIn $i \
	        --readFilesCommand gunzip -c \
	        --outSAMstrandField intronMotif \
	        --outFilterIntronMotifs RemoveNoncanonical \
	        --outSAMtype BAM SortedByCoordinate \
	        --outFileNamePrefix $iSUB. \
	        --limitBAMsortRAM 61675612266 \
	        --quantMode GeneCounts

	    done

	    source activate RSC 
	    multiqc -f -n ${PIN}.star.multiqc.report .
	    mkdir STAR.COUNTS STAR.BAMS STAR.LOGS
		mv *.ReadsPerGene.out.tab STAR.COUNTS
		mv *.bam STAR.BAMS
		mv *.out *.tab *_STARtmp *.list *star.multiqc.report_data STAR.LOGS
		mkdir STAR 
		mv STAR.* *.html STAR
		
		mv STAR ..
	cd ..


	cd STAR/STAR.BAMS
		for i in *.bam 
		do
			samtools index -b $i
		done
	cd ..
	echo
	echo
	pwd
	echo
	echo
	source activate RSeQC 
	geneBody_coverage.py -r ${bed12[${DIR}]} -i STAR.BAMS/ -o ${PIN}
	mkdir geneBodyCov
	mv *geneBodyCoverage.* log.txt geneBodyCov
	cd ..


}


while getopts "hp:t:g:s:" opt; do
	case ${opt} in
	
	h)  
		echo
		echo
		echo
		usage
		echo
		echo

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

	s)

		STRAND=$OPTARG

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

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if trimming parameter exists

if [[ ! -z "${T+x}" ]]; then
	#statements

	if [[ $T == yes ]]; then
		trimSmall
	else
		trim
	fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if genomeDir provided 

if [[ ! -z "${DIR+x}" ]]; then
	if [ ${genomeDir[${DIR}]+_} ]; then
		echo Reference genome selected = $DIR
		echo
		align

	else 
		echo "The reference genome provided '"$DIR"' is not available"
		exit 1

	fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if strandedness info provided

if [[ ! -z "${STRAND+x}" ]]; then
	if [[ $STRAND = "1" ]]; then
		echo first strand selected

		for i in STAR/STAR.COUNTS/*.ReadsPerGene.out.tab
	        do
	        awk 'NR > 4 {print $1 "\t" $3}' $i > $i.rawCounts
	        done
	        mkdir STAR/STAR.COUNTS/rawCounts
	        mv STAR/STAR.COUNTS/*.rawCounts STAR/STAR.COUNTS/rawCounts

	elif [[ $STRAND = "2" ]]; then
		echo reverse strand selected 

		for i in STAR/STAR.COUNTS/*.ReadsPerGene.out.tab
	        do
	        awk 'NR > 4 {print $1 "\t" $4}' $i > $i.rawCounts
	        done
	        mkdir STAR/STAR.COUNTS/rawCounts
	        mv STAR/STAR.COUNTS/*.rawCounts STAR/STAR.COUNTS/rawCounts

	else 
		echo unstranded selected

		for i in STAR/STAR.COUNTS/*.ReadsPerGene.out.tab
	        do
	        awk 'NR > 4 {print $1 "\t" $2}' $i > $i.rawCounts
	        done
	        mkdir STAR/STAR.COUNTS/rawCounts
	        mv STAR/STAR.COUNTS/*.rawCounts STAR/STAR.COUNTS/rawCounts
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
else
	echo `date -u` >> beta2.run.log
	echo "Project Identifier Specified = " $PIN >> beta2.run.log
	echo "Reference Genome Specified   = " $DIR >> beta2.run.log
	echo "Trimming for smRNA seq       = " $T >> beta2.run.log
	echo "Strandedness specified 	   = " $STRAND >> beta2.run.log 	
fi

