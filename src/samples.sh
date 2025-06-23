#!/bin/bash

# Settings
shotsSeed=1
numObjectives=100

# Data generation for the measurement project
> data/measure.dat

# Pre-compute the steady state
./run --2d 4 2 --precompute 

# First-order measurements, bounds on many quantities
echo "First-order measurements" | tee -a data/measure.dat
for M in 0 100 1000 2500 5000 7500 10000 20000
do
    for A in 0 10 25 50 100 200
    do
        for shots in 0 100 1000 2500 5000 7500 10000 20000
        do
            for i in seq 1 ${numObjectives}
            do
                ./run --2d 4 2 --precomputed -M ${M} -A ${A} -S ${i} --objRand 2 --shots ${shots} --all 1 -B | tee -a data/measure.dat
            done
        done
    done
done

# Mag without z vs number of measurements
echo "Mag without z vs number of measurements" | tee -a data/measure.dat
for shots in 0 100 1000 2500 5000 7500 10000 20000
do
    ./run --2d 4 2 --precomputed -S 1 -M 1000 -A 100 --objZ -B | tee -a data/measure.dat
    ./run --2d 4 2 --precomputed -S 1 --objZ --shots ${shots} --all 1 -B | tee -a data/measure.dat
    ./run --2d 4 2 --precomputed -S 1 -M 1000 -A 100 --objZ --shots ${shots} --all 1 -B | tee -a data/measure.dat
    ./run --2d 4 2 --precomputed -S 1 -M 1000 -A 100 --objZ --shots ${shots} --all 2 -B | tee -a data/measure.dat
done

# Purity versus number of measurements
echo "Purity measurements" | tee -a data/measure.dat
for shots in 0 100 1000 2500 5000 7500 10000 20000
do
    ./run --2d 4 2 --precomputed -S 1 -M 1000 -A 100 --objPurity -B | tee -a data/measure.dat
    ./run --2d 4 2 --precomputed -S 1 --objPurity --shots ${shots} --all 1 -B | tee -a data/measure.dat
    ./run --2d 4 2 --precomputed -S 1 -M 1000 -A 100 --objPurity --shots ${shots} --all 1 -B | tee -a data/measure.dat
    ./run --2d 4 2 --precomputed -S 1 -M 1000 -A 100 --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat
done

