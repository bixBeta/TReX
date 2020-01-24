#!/bin/sh
if [ "$1" = "help" ] || [ "$1" = "--help" ] || [  -z $1  ]; then
    echo ""
    echo "--------------------------------------------------------------------------------------------------------"
    echo "  To run this script, use the following syntax:"
    echo "     bash" $0 "</path/to/outs/>"
    echo "--------------------------------------------------------------------------------------------------------"
    echo ""
    echo ""
    echo ""
    exit 1

else

# commit from gg02


FASTQS=$1

cut -d "," -f3 ${FASTQS}/input_samplesheet.csv | cut -d "-" -f2 | sed 1,2d > .ids
cut -d "," -f3 ${FASTQS}/input_samplesheet.csv | sed 1,2d > .names


readarray sampleIDs < .ids
readarray sampleNames < .names
#readarray -t sampleIDs <<< `cut -d "," -f3 ${FASTQS}/input_samplesheet.csv | cut -d "-" -f2 | sed 1,2d`
#readarray -t sampleNames <<< `cut -d "," -f3 ${FASTQS}/input_samplesheet.csv | cut -d "-" -f1 | sed 1,2d`

length=${#sampleIDs[@]}
for (( i = 0; i < length; i++ ))
do
  /programs/cellranger-3.0.2/cellranger count --id=`echo ${sampleIDs[i]}` \
  --transcriptome=/workdir/singleCellData/10x_reference_files/refdata-cellranger-GRCh38-3.0.0/ \
  --fastqs=${FASTQS}/fastq_path/ \
  --sample=`echo ${sampleNames[i]}` \
  --localcores 20 --localmem 250

done

rm .ids .names
fi
