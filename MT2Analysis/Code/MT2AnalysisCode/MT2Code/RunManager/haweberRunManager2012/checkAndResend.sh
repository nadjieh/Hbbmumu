#!/bin/bash

mkdir recovery

for i in `../../../../check_runmanager.py | awk '/***SERIOUS/,/*** TIMING/{print}' | grep job_ | cut -d . -f 2 | cut -c 6- | uniq`
do 
srmrm `grep SERESULTDIR job_$i.out | cut -c 13-`"/output_$i.root"
qsub -q all.q -N myRecoJob_$i -o recovery/job_$i.out -e recovery/job_$i.err myScript_$i.sh $i
done
