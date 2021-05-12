#!/bin/sh

#SBATCH -J smRNA
#SBATCH -o %x.out
#SBATCH -n 6
#SBATCH --mem-per-cpu=18000


# source /programs/bin/util/setup_mirdeep2.sh

usage(){

    echo "sm R N A - S E Q   W O R K F L O W - @bixBeta"
    echo
    echo

    echo "Usage: bash" $0 "[-h arg] [-p arg] [-t arg] [-g arg]"
    echo
    echo "-------------------------------------------------------------------------------------------------------------------------------------------------------"
    echo "[-h] --> Display Help "
    echo "[-p] --> Project Identifier Number "
    echo "[-d] --> Comma Spearated Values for Delimiter and Field <delim,field or default> default: -,2 (complex field example: 2 | tail -c 4 or grep -o '...$')"
    echo "[-t] --> NextSeq run < yes, no, na > "
    echo "[-g] --> Mapper Genome < hsa, mmu, cel > "
    echo "[-c] --> CleanUP < yes or no > "
    echo "-------------------------------------------------------------------------------------------------------------------------------------------------------"
}


trimSmall(){
    echo "trimSmall"
        mkdir TrimQC_stats fastQC mirDeep2_results
        for i in fastqs/*.gz
        do
            /home/fa286/bin/TrimGalore-0.6.0/trim_galore --nextseq 20 --gzip -j 8 --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
        done
        mv *_trimming_report.txt TrimQC_stats
        mv *trimmed.fq.gz mirDeep2_results

}

trimHiSeq(){
    echo "trimHiSeq"
        mkdir TrimQC_stats fastQC mirDeep2_results
        for i in fastqs/*.gz
        do
            /home/fa286/bin/TrimGalore-0.6.0/trim_galore --quality 20 --gzip -j 8 --length 10  --fastqc --fastqc_args "-t 4 --outdir ./fastQC" $i
        done
        mv *_trimming_report.txt TrimQC_stats
        mv *trimmed.fq.gz mirDeep2_results

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



  if [ -f config.txt ]; then
      echo "config file exists"
      CONFIG="config.txt"
      cd ..

  else

    ls -1 *.fasta > f1

    # readarray fastas < f1
    #
    #
    # for i in "${fastas[@]}"
    #   do
    #
    #             if  echo $DELIM | grep -q "|"
    #
    #             then
    #             echo $i | cut -d ${DELIMITER} -f${FIELD} | ${CCOUNT} >> f2
    #
    #             else
    #             echo $i | cut -d ${DELIMITER} -f${FIELD} >> f2
    #
    #             fi
    #
    # done

    COUNTER=`wc -l f1 | cut -d " " -f1`
    COUNTERC=`expr $COUNTER + 100 `

    seq 101 1 $COUNTERC > f2

    paste f1 f2 > config.txt
    CONFIG="config.txt"
    cd ..

  fi

}

mapper(){
  echo "mapper"
  cd mirDeep2_results

  if [ -f *collapsed.fa ]; then

    echo "collapsed fasta exists"

    PIN=`echo *.collapsed.fa | cut -d '_' -f1`
    DATE=`echo *.collapsed.fa | cut -d '_' -f2- | cut -d '.' -f1`

  else

    DATE=`date +"%m_%d_%H-%M"`
    mapper.pl $CONFIG -d -c -m -s ${PIN}_${DATE}.collapsed.fa

  fi

  cd ..


}

quant(){

  echo "quant"
  cd mirDeep2_results

    if [ -d expression_analyses_${PIN}_${DATE} ]; then

      echo
      echo "previous expression analyses detected, re-setting the DATE variable"
      echo

      COLLAPSED="${PIN}_${DATE}.collapsed.fa"

      DATE=`date +"%m_%d_%H-%M"`

      echo
      echo "New DATE = $DATE"
      echo

      quantifier.pl -p /workdir/genomes/smRNA/hairpin.fa \
      -m /workdir/genomes/smRNA/mature.fa \
      -t $G -y ${PIN}_${DATE} -r ${COLLAPSED} -W -d


    else

      quantifier.pl -p /workdir/genomes/smRNA/hairpin.fa \
      -m /workdir/genomes/smRNA/mature.fa \
      -t $G -y ${PIN}_${DATE} -r ${PIN}_${DATE}.collapsed.fa -W -d

    fi


  cd ..

}


cleanUp(){
  echo "cleanUp"
  cd mirDeep2_results
  rm -r dir_mapper* f1 f2 *_trimmed.fasta *_trimmed.fq
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
  echo "DONE =) "

}

# awk 'BEGIN {OFS="\t"; m=0} FNR==NR {d[$1]=1; next} { if(FNR%2==1) {s=substr($1,2,3); c=substr($1,match($1,"_x")+2,15); if(substr($1,2,25) in d) m=1} else { l=length($1); f=substr($1,1,1); a[s,l,f,m]++; b[s,l,f,m]=b[s,l,f,m]+c;s=0;c=0;l=0;f=0; m=0}} END {print "library\treadlength\tbase1\tmiRBaseMatch\t#distinctReads\t#reads"; for (var in a) {split(var,q,SUBSEP); print q[1], q[2], q[3], q[4], a[var], b[var]} }' <(zcat *miRBase.mrd.gz) <(zcat *collapsed.fa.gz) > mirmap_firstbase_readlengthcounts.txt



# PIN="1058"
# fastq2fasta
# config
# mapper


while getopts "hp:t:g:d:c:" opt; do
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

    d )

    DELIM=$OPTARG

  ;;

    c )

        CLEAN=$OPTARG

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
## check if delimiter parameter exists

if [[ ! -z "${DELIM+x}" ]]; then
    #statements
    if [[ $DELIM == default ]]; then

    DELIMITER="-"
    FIELD="2"
    echo "file naming will be done using the default delimiter settings"

    else

        DELIMITER=`echo $DELIM | cut -d , -f1`
        FIELD=`echo $DELIM | cut -d , -f2- | cut -d "|" -f1`
        CCOUNT=`echo $DELIM | cut -d , -f2- | cut -d "|" -f2-`

    echo "file naming will be done using the delim = $DELIMITER and field = $FIELD initial params for $DELIM"

  fi

fi
#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if trimming parameter exists and run on nextseq 500 series (or 2 color bias )

if [[ ! -z "${T+x}" ]]; then
    #statements

    if [[ $T == yes ]]; then
        trimSmall
        fastq2fasta

  else
    trimHiSeq
        fastq2fasta

    fi
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if genome parameter is  provided

if [[ ! -z "${G+x}" ]]; then
    #statements
  echo "Genome selected --> $G "
  config
  mapper
  quant

elif [[ -z "$G"  ]]; then
  echo "Genome info not provided "
fi

#-------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------
## check if clean up parameter exists

if [[ ! -z "${CLEAN+x}" ]]; then
    #statements

    if [[ $CLEAN == yes ]]; then
        cleanUp

  else
    echo 'cleanUP not required at this time'

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
    echo `date` >> beta5.small.run.log
    echo "Project Identifier Specified = " $PIN >> beta5.small.run.log
    echo "Trimming for NextSeq         = " $T >> beta5.small.run.log
    echo "Selected Genome              = " $G >> beta5.small.run.log
    echo -------------------------------------------------------------------------------------------------- >> beta5.small.run.log
fi
