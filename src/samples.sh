#!/bin/bash

# Data generation for the measurement project
> data/measure.dat

# Pre-compute the steady state
systemSize="3 3"
filename="data/2d_3x3.dat"
#./run --2d ${systemSize} --precompute ${filename}

# Heat current without Z vs number of measurements
echo "file & heat" | tee -a data/measure.dat
A=50
M=1000
./run -N "sdp" -s M --2d ${systemSize} -M ${M} -A ${A} --objHC H -B | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # SDP plus level 2 shots
    ./run -N "sdp+all2-z, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 2 -B | tee -a data/measure.dat
    ./run -N "sdp+all2-z, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 2 -B | tee -a data/measure.dat
    ./run -N "sdp+all2-z, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 3 shots
    ./run -N "sdp+all3-z, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat
    ./run -N "sdp+all3-z, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat
    ./run -N "sdp+all3-z, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 4 shots
    #./run -N "sdp+all4-z, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat
    #./run -N "sdp+all4-z, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat
    #./run -N "sdp+all4-z, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus 20 automatic
    #./run -N "sdp+auto20, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 20 -B | tee -a data/measure.dat
    #./run -N "sdp+auto20, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 20 -B | tee -a data/measure.dat
    #./run -N "sdp+auto20, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 20 -B | tee -a data/measure.dat

    # SDP plus 50 automatic
    ./run -N "sdp+auto50-z, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 50 -B | tee -a data/measure.dat
    ./run -N "sdp+auto50-z, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 50 -B | tee -a data/measure.dat
    ./run -N "sdp+auto50-z, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus 100 automatic
    #./run -N "sdp+auto100-z, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 100 -B | tee -a data/measure.dat
    #./run -N "sdp+auto100-z, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 100 -B | tee -a data/measure.dat
    #./run -N "sdp+auto100-z, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus only objective
    ./run -N "sdp+onlyobj-z, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --onlyobj -B | tee -a data/measure.dat
    ./run -N "sdp+onlyobj-z, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --onlyobj -B | tee -a data/measure.dat
    ./run -N "sdp+onlyobj-z, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --onlyobj -B | tee -a data/measure.dat

done


# Purity versus number of measurements
echo "file & purity" | tee -a data/measure.dat
A=50
M=1000
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # Only level 2 shots
    #./run -N "all2, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat
    #./run -N "all2, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat
    #./run -N "all2, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # Only level 3 shots
    ./run -N "all3, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 3 -B | tee -a data/measure.dat
    ./run -N "all3, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 3 -B | tee -a data/measure.dat
    ./run -N "all3, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # Only level 4 shots
    #./run -N "all4, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat
    #./run -N "all4, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat
    #./run -N "all4, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 2 shots
    #./run -N "sdp+all2, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat
    #./run -N "sdp+all2, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat
    #./run -N "sdp+all2, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 3 shots
    ./run -N "sdp+all3, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 3 -B | tee -a data/measure.dat
    ./run -N "sdp+all3, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 3 -B | tee -a data/measure.dat
    ./run -N "sdp+all3, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 4 shots
    #./run -N "sdp+all4, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat
    #./run -N "sdp+all4, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat
    #./run -N "sdp+all4, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus auto 50
    #./run -N "sdp+auto50, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 50 -B | tee -a data/measure.dat
    #./run -N "sdp+auto50, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 50 -B | tee -a data/measure.dat
    #./run -N "sdp+auto50, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus auto 100
    ./run -N "sdp+auto100, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 100 -B | tee -a data/measure.dat
    ./run -N "sdp+auto100, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 100 -B | tee -a data/measure.dat
    ./run -N "sdp+auto100, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus auto 300
    ./run -N "sdp+auto300, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 300 -B | tee -a data/measure.dat
    ./run -N "sdp+auto300, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 300 -B | tee -a data/measure.dat
    ./run -N "sdp+auto300, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 300 -B | tee -a data/measure.dat

done

# Pre-compute the ground state
systemSize="3 3"
filename="data/2d_3x3_H.dat"
#./run --2d ${systemSize} --precompute ${filename} -H

# Ground state energy vs number of measurements
echo "file & energy" | tee -a data/measure.dat
A=50
./run -N "sdp" -s M --2d ${systemSize} -A ${A} -H -B | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # Just level 2 shots
    #./run -N "all2, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} -H -S 1 --shots ${shots} --all 2 -B | tee -a data/measure.dat
    #./run -N "all2, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} -H -S 1 --shots ${shots} --all 2 -B | tee -a data/measure.dat
    #./run -N "all2, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} -H -S 1 --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # Just level 3 shots
    #./run -N "all3, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} -H -S 1 --shots ${shots} --all 3 -B | tee -a data/measure.dat
    #./run -N "all3, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} -H -S 1 --shots ${shots} --all 3 -B | tee -a data/measure.dat
    #./run -N "all3, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} -H -S 1 --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # Only measure the objective
    ./run -N "sdp+onlyobj, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --onlyobj -B | tee -a data/measure.dat
    ./run -N "sdp+onlyobj, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --onlyobj -B | tee -a data/measure.dat
    ./run -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --onlyobj -B | tee -a data/measure.dat

    # SDP plus level 2 shots
    ./run -N "sdp+all2, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat
    ./run -N "sdp+all2, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat
    ./run -N "sdp+all2, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 3 shots
    #./run -N "sdp+all3, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 3 -B | tee -a data/measure.dat
    #./run -N "sdp+all3, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 3 -B | tee -a data/measure.dat
    #./run -N "sdp+all3, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 4 shots
    #./run -N "sdp+all4, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 4 -B | tee -a data/measure.dat
    #./run -N "sdp+all4, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 4 -B | tee -a data/measure.dat
    #./run -N "sdp+all4, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus 50 automatic
    #./run -N "sdp+auto50, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 50 -B | tee -a data/measure.dat
    #./run -N "sdp+auto50, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 50 -B | tee -a data/measure.dat
    #./run -N "sdp+auto50, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus 100 automatic
    #./run -N "sdp+auto100, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 100 -B | tee -a data/measure.dat
    #./run -N "sdp+auto100, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 100 -B | tee -a data/measure.dat
    #./run -N "sdp+auto100, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus 300 automatic
    ./run -N "sdp+auto300, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 300 -B | tee -a data/measure.dat
    ./run -N "sdp+auto300, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 300 -B | tee -a data/measure.dat
    ./run -N "sdp+auto300, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 300 -B | tee -a data/measure.dat

done


# Heat current for large systems using symmetries
echo "file & large" | tee -a data/measure.dat
systemSize="50 2"
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # SDP plus symmetry
    ./run -N "sdp+sym, 68%" -p 68 -s M --2d ${systemSize} -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat
    ./run -N "sdp+sym, 95%" -p 95 -s M --2d ${systemSize} -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat
    ./run -N "sdp+sym, 99.7%" -p 99.7 -s M --2d ${systemSize} -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat

done




