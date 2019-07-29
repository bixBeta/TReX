#!/bin/sh

# source /programs/bin/util/setup_mirdeep2.sh

usage(){

  echo "Usage: bash" $0 "[-h arg] [-p arg] [-t arg]"
	echo
	echo "---------------------------------------------------------------------------------------------------------------"
	echo "[-h] --> Display Help"
	echo "[-p] --> Project Identifier Number"
	echo "[-t] --> NextSeq run <yes, no>"
	echo "---------------------------------------------------------------------------------------------------------------"

}


trimSmall(){
    echo "trimSmall"
		mkdir TrimQC_stats fastQC trimmed_fastqs
		for i in fastqs/*.gz
		do
			/home/fa286/bin/TrimGalore-0.6.0/trim_galore --nextseq 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
		done
		mv *_trimming_report.txt TrimQC_stats
		mv *trimmed.fq.gz trimmed_fastqs

}

trimHiSeq(){
    echo "trimHiSeq"
		mkdir TrimQC_stats fastQC trimmed_fastqs
		for i in fastqs/*.gz
		do
			/home/fa286/bin/TrimGalore-0.6.0/trim_galore --quality 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
		done
		mv *_trimming_report.txt TrimQC_stats
		mv *trimmed.fq.gz trimmed_fastqs

}

fastq2fasta(){
  echo "fastq2fasta"
  cd trimmed_fastqs
  gunzip *.gz

  for i in *.fq
  do
    iSUB=`echo $i | cut -d "." -f1`
    fastq2fasta.pl $i > ${iSUB}.fasta
  done
  cd ..
}

config(){
  echo "config"
  cd trimmed_fastqs
  ls -1 *.fasta > f1
  readarray fastas < f1
  for i in "${fastas[@]}"
    do
      echo $i | cut -d "-" -f2 >> f2
    done

  paste f1 f2 > config.txt

  CONFIG="config.txt"
  cd ..

}

mapper(){
  echo "mapper"
  cd trimmed_fastqs
  DATE=`date +"%m_%d_%H-%M"`
  mapper.pl $CONFIG -d -c -m -s ${PIN}_${DATE}.collapsed.fa
  cd ..
}

quantHuman(){
  cd trimmed_fastqs
  quantifier.pl -p /workdir/RSC/referenceFiles/miRBase/v22_1/hairpin.fa \
  -m /workdir/RSC/referenceFiles/miRBase/v22_1/mature.fa \
  -t hsa -y ${PIN}_${DATE} -r ${PIN}_${DATE}.collapsed.fa -W -d
  cd ..
}

cleanUp(){
  cd trimmed_fastqs
  mkdir fastas
  mv *.fasta fastas

  mkdir logFiles
  mv f1 f2 mapper.log dir_mapper* config.txt logFiles
  cd ..

  mkdir mirDeep2
  mv expression_analyses expression*.html miRNAs_expressed_all_samples*.csv mirDeep2

}

# PIN="1058"
# fastq2fasta
# config
# mapper


while getopts "hp:t:" opt; do
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

	\? )
		echo
		echo
		echo
		usage

	;;

	esac

done
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if PIN is provided

if [[ -z "${PIN+x}" ]]; then

	PIN="PIN_Null"
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if trimming parameter exists and run on nextseq 500 series (or 2 color bias )

if [[ ! -z "${T+x}" ]]; then
	#statements

	if [[ $T == yes ]]; then
		trimSmall
    fastq2fasta
    config
    mapper
    quantHuman
    cleanUp
	else
		trimHiSeq
    fastq2fasta
    config
    mapper
    quantHuman
    cleanUp
	fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------

if [[ -z $1 ]] || [[  $1 = "--help"  ]] ; then
	#statements
	echo
	echo
	usage
	echo
	echo
	exit 1
else
	echo
	echo `date` >> beta2.small.run.log
	echo "Project Identifier Specified = " $PIN >> beta2.small.run.log
	echo "Trimming for NextSeq         = " $T >> beta2.small.run.log
	echo -------------------------------------------------------------------------------------------------- >> beta2.small.run.log
fi
