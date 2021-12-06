#!/bin/bash 

# Iterate over strategies
for s in 0 1; do
    # Iterate over G functions
    for g in 0 1 2; do
        echo "Strategy: $s, G: $g" | tee -a report_strategy.txt;
        make clean all run STRATEGY=$s G_FUNC=$g VERB=1 | tee -a report_strategy.txt
        echo "Done" | tee -a report_strategy.txt;
        sleep 10;
    done
done
