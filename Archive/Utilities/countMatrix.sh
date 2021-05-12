#!/bin/sh

cat phenoData.csv | cut -d ',' -f2 | awk 'NR > 1 {print $0}' > .tmp.Names
paste `cat .tmp.Names` > countMatrix.txt

