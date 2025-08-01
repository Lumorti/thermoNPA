#!/bin/bash

# Data generation for the measurement project
> data/measure.dat

# Pre-compute the steady state
systemSize="3 3"
filename="data/2d_3x3.dat"
#./run --2d ${systemSize} --precompute ${filename}

# Heat current without Z vs number of measurements
echo "file & heat" | tee -a data/measure.dat
A=20
M=25
./run -N "sdp" -s M --2d ${systemSize} -M ${M} -A ${A} --objHC H -B | tee -a data/measure.dat
for shots in 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 -1
do

    # Only 200 automatic
    #./run -N "auto200" -s M --2d ${systemSize} --precomputed ${filename} -S 1 --objHC H --noz --shots ${shots} --auto 200 -B | tee -a data/measure.dat

    # Only level 4 shots
    #./run -N "shots4" -s M --2d ${systemSize} --precomputed ${filename} -S 1 --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 1 shots
    #./run -N "sdp+shots1" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 1 -B | tee -a data/measure.dat

    # SDP plus level 2 shots, 68%
    ./run -N "sdp+shots2, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 2 shots, 95%
    ./run -N "sdp+shots2, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 2 shots, 99.7%
    ./run -N "sdp+shots2, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 3 shots
    #./run -N "sdp+shots3, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 3 shots
    #./run -N "sdp+shots3, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 3 shots
    #./run -N "sdp+shots3, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 4 shots 68%
    ./run -N "sdp+shots4, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 4 shots 95%
    ./run -N "sdp+shots4, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 4 shots 99.7%
    ./run -N "sdp+shots4, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus 10 automatic
    #./run -N "sdp+auto10" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 10 -B | tee -a data/measure.dat

    # SDP plus 20 automatic 68%
    #./run -N "sdp+auto20, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 20 -B | tee -a data/measure.dat

    # SDP plus 20 automatic 95%
    #./run -N "sdp+auto20, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 20 -B | tee -a data/measure.dat

    # SDP plus 20 automatic 99.7%
    #./run -N "sdp+auto20, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 20 -B | tee -a data/measure.dat

    # SDP plus 50 automatic, 68%
    ./run -N "sdp+auto50, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus 50 automatic, 95%
    ./run -N "sdp+auto50, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus 50 automatic, 99.7%
    ./run -N "sdp+auto50, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus only objective, 68%
    ./run -N "sdp+onlyobj, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --onlyobj -B | tee -a data/measure.dat

    # SDP plus only objective, 95%
    ./run -N "sdp+onlyobj, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --onlyobj -B | tee -a data/measure.dat

    # SDP plus only objective, 99.7%
    ./run -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --onlyobj -B | tee -a data/measure.dat

    # SDP plus 100 automatic
    #./run -N "sdp+auto100" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus 200 automatic
    #./run -N "sdp+auto200" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 200 -B | tee -a data/measure.dat

done


# Purity versus number of measurements
echo "file & purity" | tee -a data/measure.dat
A=25
M=50
for shots in 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 -1
do

    # Only level 2 shots, 68%
    ./run -N "shots2, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # Only level 2 shots, 95%
    ./run -N "shots2, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # Only level 2 shots, 99.7%
    ./run -N "shots2, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # Only level 4 shots, 68%
    ./run -N "shots4, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # Only level 4 shots, 95%
    ./run -N "shots4, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # Only level 4 shots, 99.7%
    ./run -N "shots4, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 2 shots, 68%
    ./run -N "sdp+shots1, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 2 shots, 95%
    ./run -N "sdp+shots1, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 2 shots, 99.7%
    ./run -N "sdp+shots1, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 4 shots, 68%
    ./run -N "sdp+shots4, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 4 shots, 95%
    ./run -N "sdp+shots4, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 4 shots, 99.7%
    ./run -N "sdp+shots4, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus auto 50, 68%
    ./run -N "sdp+auto50, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus auto 50, 95%
    ./run -N "sdp+auto50, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus auto 50, 99.7%
    ./run -N "sdp+auto50, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 50 -B | tee -a data/measure.dat

done

# Pre-compute the ground state
systemSize="3 3"
filename="data/2d_3x3_H.dat"
#./run --2d ${systemSize} --precompute ${filename} -H

# Ground state energy vs number of measurements
echo "file & energy" | tee -a data/measure.dat
A=20
./run -N "sdp" -s M --2d ${systemSize} -A ${A} -H -B | tee -a data/measure.dat
for shots in 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 -1
do

    # Only 200 automatic
    #./run -N "auto200" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -H --shots ${shots} --auto 200 -B | tee -a data/measure.dat

    # Only level 4 shots
    #./run -N "shots4" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -H --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # Only measure the objective, 68%
    ./run -N "sdp+onlyobj, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --onlyobj -B | tee -a data/measure.dat

    # Only measure the objective, 95%
    ./run -N "sdp+onlyobj, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --onlyobj -B | tee -a data/measure.dat

    # Only measure the objective, 99.7%
    ./run -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --onlyobj -B | tee -a data/measure.dat

    # SDP plus level 1 shots
    #./run -N "sdp+shots1" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 1 -B | tee -a data/measure.dat

    # SDP plus level 2 shots, 68%
    ./run -N "sdp+shots2, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 2 shots, 95%
    ./run -N "sdp+shots2, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 2 shots, 99.7%
    ./run -N "sdp+shots2, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 3 shots, 68%
    #./run -N "sdp+shots3, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 3 shots, 95%
    #./run -N "sdp+shots3, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 3 shots, 99.7%
    #./run -N "sdp+shots3, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 4 shots, 68%
    ./run -N "sdp+shots4, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 4 shots, 95%
    ./run -N "sdp+shots4, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus level 4 shots, 99.7%
    ./run -N "sdp+shots4, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus 10 automatic
    #./run -N "sdp+auto10" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 10 -B | tee -a data/measure.dat

    # SDP plus 20 automatic
    #./run -N "sdp+auto20" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 20 -B | tee -a data/measure.dat

    # SDP plus 50 automatic, 68%
    #./run -N "sdp+auto50, 68%" -p 68 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus 50 automatic, 95%
    #./run -N "sdp+auto50, 95%" -p 95 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus 50 automatic, 99.7%
    #./run -N "sdp+auto50, 99.7%" -p 99.7 -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus 100 automatic, 68%
    ./run -N "sdp+auto100, 68%" -s M -p 68 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus 100 automatic, 95%
    ./run -N "sdp+auto100, 95%" -s M -p 95 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus 100 automatic, 99.7%
    ./run -N "sdp+auto100, 99.7%" -s M -p 99.7 --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus 200 automatic
    #./run -N "sdp+auto200" -s M --2d ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 200 -B | tee -a data/measure.dat

done


