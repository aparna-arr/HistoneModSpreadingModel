#!/bin/bash

#for rt in 5000 3000 1000
#do
#    for f in 0.5 1 2 5 8 10
#    do
#        total_time=$(($rt+60000))
#        python3 MainSim.py -t $total_time --divisions-num 1000 -i A -r --recruit-time-init 30000 --recruit-time $rt --recruit-n 5 -f $f -o bigsim_t$total_time\_d1000_iA_rinit30000_rt$rt\_rn5_f$f
#    done
#done

for rt in 5000 3000 1000
do
    for f in 0.5 1 2 5 8 10
    do
        total_time=$(($rt+60000))
        python3 MainSim.py -t $total_time --divisions-num 1000 -i A -r --recruit-time-init 30000 --recruit-time $rt --recruit-n 5 -f $f --prob-spread powerlaw -o bigsim_t$total_time\_d1000_iA_rinit30000_rt$rt\_rn5_f$f\_powerlaw
    done
done
