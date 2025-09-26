#!/bin/bash

# Data generation for the measurement project
> data/measure.dat

# System settings
systemSize="3 3"
A=150
M=2000
numRepeats=50
filename="data/2d_3x3.dat"
filenameEnergy="data/2d_3x3_H.dat"

# Pre-compute both the steady state and the ground state
> data/precomputes.log
./run --2dtfi ${systemSize} --precompute ${filename} | tee -a data/precomputes.log
./run --2dtfi ${systemSize} --precompute ${filenameEnergy} -H | tee -a data/precomputes.log

# Heat current without something vs number of measurements
echo "file & heat" | tee -a data/measure.dat
./run -B -N "sdp" -s M --2dtfi ${systemSize} -M ${M} -A ${A} --objHC H | tee -a data/measure.dat
letter=z
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # Only 100 first
    ./run -B -N "first100-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S ${ind} --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat

    # SDP plus 100 first
    ./run -B -N "sdp+first100-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S ${ind} -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat

done

# Purity versus number of measurements
echo "file & purity" | tee -a data/measure.dat
./run -B -N "sdp" -s M --2dtfi ${systemSize} -A ${A} --objPurity | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do
    for ind in $(seq 1 $numRepeats)
    do

        # Only level 2 shots
        ./run -B -N "all2, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S ${ind} --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat

        # SDP plus level 2 shots
        ./run -B -N "sdp+all2, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S ${ind} -A ${A} --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat

    done
done

# Ground state energy vs number of measurements
echo "file & energy" | tee -a data/measure.dat
./run -B -N "sdp" -s M --2dtfi ${systemSize} -A ${A} -H | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do
    for ind in $(seq 1 $numRepeats)
    do

        # SDP plus automatic 100
        ./run -B -N "sdp+auto100, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S ${ind} -A ${A} -H --shots ${shots} --auto 100 | tee -a data/measure.dat

        # SDP plus measure the objective
        ./run -B -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S ${ind} -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat

        # SDP plus level 2 shots
        ./run -B -N "sdp+all2, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S ${ind} -A ${A} -H --shots ${shots} --all 2 | tee -a data/measure.dat

    done
done

# File showing the different confidence levels
echo "file & confidence" | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000
do
    for conf in 68 95 99.7
    do

        # SDP plus level 2 shots
        ./run -B -N "sdp+all2, ${conf}%" -s M -p ${conf} --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 -A ${A} --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat

    done
done

# Ground state energy for large systems
echo "file & large" | tee -a data/measure.dat
systemSize="50"
A=200
./run -B -S ${ind} -N "sdp" -s M --mg ${systemSize} -A 70 -H | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do
    for ind in $(seq 1 $numRepeats)
    do

        # Only the objective
        ./run -B -S ${ind} -N "onlyobj, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -H --shots ${shots} --onlyobj | tee -a data/measure.dat

        # SDP plus only the objective
        ./run -B -S ${ind} -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat

    done
done

# Automatically commit once done
git add .
git commit -m "automatic data commit"
git push



