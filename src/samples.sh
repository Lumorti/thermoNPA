#!/bin/bash

# Settings
shotsSeed=1
numObjectives=50
systemSize="2 2"

# Data generation for the measurement project
> data/measure.dat

# Pre-compute the steady state
./run --2d ${systemSize} --precompute 

# First-order measurements, bounds on many quantities
#echo "First-order measurements" | tee -a data/measure.dat
#for M in 100 500 1000 
#do
    #for A in 10 25 50
    #do
        #for shots in 100 1000 2500 5000 7500 10000 15000 20000
        #do
            #for i in $(seq 1 ${numObjectives})
            #do
                #./run --2d ${systemSize} --precomputed -M ${M} -A ${A} -S ${i} --objRand 2 -S ${shotsSeed} --shots ${shots} --all 1 -B | tee -a data/measure.dat
            #done
        #done
    #done
#done

# Mag without z vs number of measurements
echo "file & magnoz" | tee -a data/measure.dat
A=20
M=25
for shots in 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 50000000 -1
do

    # SDP plus level 1 shots
    ./run -N "sdp+shots1" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objZ --noz --shots ${shots} --all 1 -B | tee -a data/measure.dat

    # SDP plus level 2 shots
    ./run -N "sdp+shots2" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objZ --noz --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 3 shots
    ./run -N "sdp+shots3" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objZ --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 4 shots
    ./run -N "sdp+shots4" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objZ --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus 10 automatic
    ./run -N "sdp+auto10" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objZ --noz --shots ${shots} --auto 10 -B | tee -a data/measure.dat

    # SDP plus 50 automatic
    ./run -N "sdp+auto50" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objZ --noz --shots ${shots} --auto 50 -B | tee -a data/measure.dat

done

# Purity versus number of measurements
echo "file & purity" | tee -a data/measure.dat
A=25
M=50
for shots in 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 50000000 -1
do

    # only level 1 shots
    ./run -N "shots1" -s M --2d ${systemSize} --precomputed --millis -S 1 --objPurity --shots ${shots} --all 1 -B | tee -a data/measure.dat

    # only level 4 shots
    ./run -N "shots4" -s M --2d ${systemSize} --precomputed --millis -S 1 --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 1 shots
    ./run -N "sdp+shots1" -s M --2d ${systemSize} --precomputed --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 1 -B | tee -a data/measure.dat

    # SDP plus level 4 shots
    ./run -N "sdp+shots4" -s M --2d ${systemSize} --precomputed --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

done

