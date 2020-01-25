for i in *.bam
do
	iSUB=`echo $i | cut -d "." -f1`
	samtools reheader in.header.sam $i > ${iSUB}.bam

done

