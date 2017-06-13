#!/bin/bash

dir=$1
outfilebase=$2

ls $dir/*t65000* > process_t65000_$outfilebase
ls $dir/*t63000* > process_t63000_$outfilebase
ls $dir/*t61000* > process_t61000_$outfilebase

./process_sims.pl process_out_t65000_$outfilebase process_t65000_$outfilebase
./process_sims.pl process_out_t63000_$outfilebase process_t63000_$outfilebase
./process_sims.pl process_out_t61000_$outfilebase process_t61000_$outfilebase

Rscript plot.R process_out_t61000_$outfilebase process_out_t63000_$outfilebase process_out_t61000_$outfilebase $outfilebase
