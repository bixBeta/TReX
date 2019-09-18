#!/bin/sh

# source /programs/bin/util/setup_mirdeep2.sh

usage(){

  echo "Usage: bash" $0 "[-h arg] [-p arg] [-t arg] [-g arg]"
	echo
	echo "---------------------------------------------------------------------------------------------------------------"
	echo "[-h] --> Display Help "
	echo "[-p] --> Project Identifier Number "
	echo "[-t] --> NextSeq run < yes, no, na > "
  echo "[-g] --> Mapper Genome < hsa, mmu, cel > "
	echo "---------------------------------------------------------------------------------------------------------------"

}


trimSmall(){
    echo "trimSmall"
		mkdir TrimQC_stats fastQC mirDeep2_results
		for i in fastqs/*.gz
		do
			/home/fa286/bin/TrimGalore-0.6.0/trim_galore --nextseq 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
		done
		mv *_trimming_report.txt TrimQC_stats
		mv *trimmed.fq.gz mirDeep2_results
    cd mirDeep2_results
    gunzip *
    cd ..

}

trimHiSeq(){
    echo "trimHiSeq"
		mkdir TrimQC_stats fastQC mirDeep2_results
		for i in fastqs/*.gz
		do
			/home/fa286/bin/TrimGalore-0.6.0/trim_galore --quality 20 --gzip --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
		done
		mv *_trimming_report.txt TrimQC_stats
		mv *trimmed.fq.gz mirDeep2_results
    cd mirDeep2_results
    gunzip *
    cd ..

}

fastq2fasta(){
  echo "fastq2fasta"
  cd mirDeep2_results
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
  cd mirDeep2_results
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
  cd mirDeep2_results
  DATE=`date +"%m_%d_%H-%M"`
  mapper.pl $CONFIG -d -c -m -s ${PIN}_${DATE}.collapsed.fa
  cd ..
}

quant(){
  echo "quant"
  cd mirDeep2_results
  quantifier.pl -p /workdir/RSC/referenceFiles/miRBase/v22_1/hairpin.fa \
  -m /workdir/RSC/referenceFiles/miRBase/v22_1/mature.fa \
  -t $G -y ${PIN}_${DATE} -r ${PIN}_${DATE}.collapsed.fa -W -d
  cd ..
}


cleanUp(){
  echo "cleanUp"
  cd mirDeep2_results
  # rm -r dir_mapper* f1 f2 *_trimmed.fasta *_trimmed.fq
  cd expression_analyses/expression_analyses_${PIN}_${DATE}
      mv *.mrd *.arf ../../
  cd ../../

  mv miRBase.mrd ${PIN}_${DATE}_miRBase.mrd
  mv mature_mapped.arf ${PIN}_${DATE}_mature_mapped.arf
  rm -r expression_analyses
  mkdir  expression_analyses_${PIN}_${DATE}
  mv *.arf *.fa *.mrd *.html *.csv  expression_analyses_${PIN}_${DATE}

  cd expression_analyses_${PIN}_${DATE}
    gzip *.arf *.fa *.mrd
  cd ..

  echo ""
  echo "DONE =)"

}

# PIN="1058"
# fastq2fasta
# config
# mapper


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


  g )

    G=$OPTARG

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

  else
    trimHiSeq

	fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if genome parameter is  provided

if [[ ! -z "${G+x}" ]]; then
	#statements
  echo "Genome selected --> $G "
  fastq2fasta
  config
  mapper
  quant
  cleanUp

elif [[ -z "$G"  ]]; then
  echo "Genome info not provided "
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
	echo `date` >> beta3.small.run.log
	echo "Project Identifier Specified = " $PIN >> beta3.small.run.log
	echo "Trimming for NextSeq         = " $T >> beta3.small.run.log
	echo "Selected Genome              = " $G >> beta3.small.run.log
	echo -------------------------------------------------------------------------------------------------- >> beta3.small.run.log
fi
