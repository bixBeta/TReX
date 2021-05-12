#!/bin/sh

for i in *.pdf 
do 
	iSUB=`echo $i | cut -d "." -f3` 
	sips -s format png $i --out $iSUB.png 

done


mv *.png figures



