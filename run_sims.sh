#!/bin/bash
for f in 0.5 2 10 20 50 70
do 
    python3 doddModel_subplots.py 60 10000 $f 0 1 n60_e5000_f$f\_init0_div1.mp4
    python3 doddModel_subplots.py 60 10000 $f -2 1 n60_e5000_f$f\_init-2_div1.mp4
    python3 doddModel_subplots.py 60 10000 $f 0 0 n60_e5000_f$f\_init0_div0.mp4
    python3 doddModel_subplots.py 60 10000 $f -2 0 n60_e5000_f$f\_init-2_div0.mp4
done
