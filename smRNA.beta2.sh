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
			$TRIM --nextseq 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
		done
		mv *_trimming_report.txt TrimQC_stats
		mv *trimmed.fq.gz trimmed_fastqs
    cd trimmed_fastqs
    gunzip *
    cd ..

}

trimHiSeq(){
    echo "trimHiSeq"
		mkdir TrimQC_stats fastQC trimmed_fastqs
		for i in fastqs/*.gz
		do
			$TRIM -q 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
		done
		mv *_trimming_report.txt TrimQC_stats
		mv *trimmed.fq.gz trimmed_fastqs
    cd trimmed_fastqs
    gunzip *
    cd ..

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
	else
		trimHiSeq
    fastq2fasta
    config
    mapper
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
