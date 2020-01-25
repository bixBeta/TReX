for i in *.txt
do
	iSUB=`echo $i | cut -d "." -f1`
	awk 'NR > 2{print $0 }' $i | cut -f1,7 > $iSUB.rawCounts

done

