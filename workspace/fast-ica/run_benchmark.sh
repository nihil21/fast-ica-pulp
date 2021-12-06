#!/bin/bash 

# Iterate over strategies
for s in 0 1; do
    # Iterate over G functions
    for g in 0 1 2; do
        echo "Strategy: $s, G: $g";
        make clean all run STRATEGY=$s G_FUNC=$g VERB=1 >> report_$s\_$g.log
        echo "Done";
        sleep 10;
    done
done
