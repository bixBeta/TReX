#!/usr/bin/env bash


for i in *complete*; do mv $i `echo $i | sed 's/vs/_vs_/g'` ; done
