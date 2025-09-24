#!/bin/bash

# Data generation for the measurement project
> data/measure.dat
> data/precomputes.log

# System settings
systemSize="3 3"
A=150
M=2000

# Pre-compute the steady state
filename="data/2d_3x3.dat"
./run --2dtfi ${systemSize} --precompute ${filename} | tee -a data/precomputes.log

# Pre-compute the ground state
filenameEnergy="data/2d_3x3_H.dat"
./run --2dtfi ${systemSize} --precompute ${filenameEnergy} -H | tee -a data/precomputes.log

# Heat current without something vs number of measurements
echo "file & heat" | tee -a data/measure.dat
./run -B -N "sdp" -s M --2dtfi ${systemSize} -M ${M} -A ${A} --objHC H | tee -a data/measure.dat
letter=z
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # Only level 2 shots
    ./run -B -N "all2-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -N "all2-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --all 2 | tee -a data/measure.dat

    # SDP plus level 2 shots
    ./run -B -N "sdp+all2-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -N "sdp+all2-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --all 2 | tee -a data/measure.dat

    # Only level 3 shots
    ./run -B -N "all3-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --all 3 | tee -a data/measure.dat
    ./run -B -N "all3-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --all 3 | tee -a data/measure.dat

    # SDP plus level 3 shots
    ./run -B -N "sdp+all3-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --all 3 | tee -a data/measure.dat
    ./run -B -N "sdp+all3-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --all 3 | tee -a data/measure.dat

    # Only automatic 100
    ./run -B -N "auto100-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --auto 100 | tee -a data/measure.dat
    ./run -B -N "auto100-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --auto 100 | tee -a data/measure.dat

    # SDP plus 100 automatic
    ./run -B -N "sdp+auto100-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --auto 100 | tee -a data/measure.dat
    ./run -B -N "sdp+auto100-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --auto 100 | tee -a data/measure.dat

    # Only automatic 200
    ./run -B -N "auto200-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --auto 200 | tee -a data/measure.dat
    ./run -B -N "auto200-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --auto 200 | tee -a data/measure.dat

    # SDP plus 200 automatic
    ./run -B -N "sdp+auto200-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --auto 200 | tee -a data/measure.dat
    ./run -B -N "sdp+auto200-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --auto 200 | tee -a data/measure.dat

    # Only 50 first
    ./run -B -N "first50-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --first 50 | tee -a data/measure.dat
    ./run -B -N "first50-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --first 50 | tee -a data/measure.dat

    # SDP plus 50 first
    ./run -B -N "sdp+first50-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --first 50 | tee -a data/measure.dat
    ./run -B -N "sdp+first50-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --first 50 | tee -a data/measure.dat

    # Only 100 first
    ./run -B -N "first100-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat
    ./run -B -N "first100-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat

    # SDP plus 100 first
    ./run -B -N "sdp+first100-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat
    ./run -B -N "sdp+first100-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat

    # Only 200 first
    ./run -B -N "first200-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --first 200 | tee -a data/measure.dat
    ./run -B -N "first200-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --first 200 | tee -a data/measure.dat

    # SDP plus 200 first
    ./run -B -N "sdp+first200-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --first 200 | tee -a data/measure.dat
    ./run -B -N "sdp+first200-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --first 200 | tee -a data/measure.dat

    # Only the objective
    ./run -B -N "onlyobj-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --onlyobj | tee -a data/measure.dat
    ./run -B -N "onlyobj-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 --objHC H --no${letter} --shots ${shots} --onlyobj | tee -a data/measure.dat

    # SDP plus only objective
    ./run -B -N "sdp+onlyobj-${letter}, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --onlyobj | tee -a data/measure.dat
    ./run -B -N "sdp+onlyobj-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --onlyobj | tee -a data/measure.dat

done

# Purity versus number of measurements
echo "file & purity" | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # Only level 2 shots
    ./run -B -N "all2, 95%" -s M -p 95 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -N "all2, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat

    # SDP plus level 2 shots
    ./run -B -N "sdp+all2, 95%" -s M -p 95 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 -A ${A} --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -N "sdp+all2, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 -A ${A} --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat

    # Only level 3 shots
    ./run -B -N "all3, 95%" -s M -p 95 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 --objPurity --shots ${shots} --all 3 | tee -a data/measure.dat
    ./run -B -N "all3, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 --objPurity --shots ${shots} --all 3 | tee -a data/measure.dat

    # SDP plus level 3 shots
    ./run -B -N "sdp+all3, 95%" -s M -p 95 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 -A ${A} --objPurity --shots ${shots} --all 3 | tee -a data/measure.dat
    ./run -B -N "sdp+all3, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 -A ${A} --objPurity --shots ${shots} --all 3 | tee -a data/measure.dat

done

# Ground state energy vs number of measurements
echo "file & energy" | tee -a data/measure.dat
./run -B -N "sdp" -s M --2dtfi ${systemSize} -A ${A} -H | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # Only automatic 100
    ./run -B -N "auto100, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --auto 100 | tee -a data/measure.dat
    ./run -B -N "auto100, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --auto 100 | tee -a data/measure.dat

    # SDP plus automatic 100
    ./run -B -N "sdp+auto100, 95%"   -s M -p 95   --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --auto 100 | tee -a data/measure.dat
    ./run -B -N "sdp+auto100, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --auto 100 | tee -a data/measure.dat

    # Only automatic 200
    ./run -B -N "auto200, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --auto 200 | tee -a data/measure.dat
    ./run -B -N "auto200, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --auto 200 | tee -a data/measure.dat

    # SDP plus automatic 200
    ./run -B -N "sdp+auto200, 95%"   -s M -p 95   --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --auto 200 | tee -a data/measure.dat
    ./run -B -N "sdp+auto200, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --auto 200 | tee -a data/measure.dat

    # Only automatic 300
    ./run -B -N "auto300, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --auto 300 | tee -a data/measure.dat
    ./run -B -N "auto300, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --auto 300 | tee -a data/measure.dat

    # SDP plus automatic 300
    ./run -B -N "sdp+auto300, 95%"   -s M -p 95   --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --auto 300 | tee -a data/measure.dat
    ./run -B -N "sdp+auto300, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --auto 300 | tee -a data/measure.dat

    # Only first 100
    ./run -B -N "first100, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --first 100 | tee -a data/measure.dat
    ./run -B -N "first100, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --first 100 | tee -a data/measure.dat

    # SDP plus first 100
    ./run -B -N "sdp+first100, 95%"   -s M -p 95   --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --first 100 | tee -a data/measure.dat
    ./run -B -N "sdp+first100, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --first 100 | tee -a data/measure.dat

    # Only first 200
    ./run -B -N "first200, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --first 200 | tee -a data/measure.dat
    ./run -B -N "first200, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --first 200 | tee -a data/measure.dat

    # SDP plus first 200
    ./run -B -N "sdp+first200, 95%"   -s M -p 95   --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --first 200 | tee -a data/measure.dat
    ./run -B -N "sdp+first200, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --first 200 | tee -a data/measure.dat

    # Only first 300
    ./run -B -N "first300, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --first 300 | tee -a data/measure.dat
    ./run -B -N "first300, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --first 300 | tee -a data/measure.dat

    # SDP plus first 300
    ./run -B -N "sdp+first300, 95%"   -s M -p 95   --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --first 300 | tee -a data/measure.dat
    ./run -B -N "sdp+first300, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --first 300 | tee -a data/measure.dat

    # Only measure the objective
    ./run -B -N "onlyobj, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --onlyobj | tee -a data/measure.dat
    ./run -B -N "onlyobj, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --onlyobj | tee -a data/measure.dat

    # SDP plus measure the objective
    ./run -B -N "sdp+onlyobj, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat
    ./run -B -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat

    # Only level 2 shots
    ./run -B -N "all2, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -N "all2, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --all 2 | tee -a data/measure.dat
    
    # SDP plus level 2 shots
    ./run -B -N "sdp+all2, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -N "sdp+all2, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --all 2 | tee -a data/measure.dat

    # Only level 3 shots
    ./run -B -N "all3, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --all 3 | tee -a data/measure.dat
    ./run -B -N "all3, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -H --shots ${shots} --all 3 | tee -a data/measure.dat

    # SDP plus level 3 shots
    ./run -B -N "sdp+all3, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --all 3 | tee -a data/measure.dat
    ./run -B -N "sdp+all3, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S 1 -A ${A} -H --shots ${shots} --all 3 | tee -a data/measure.dat

done

# Ground state energy for large systems
echo "file & large" | tee -a data/measure.dat
systemSize="50"
A=200
./run -B -S 1 -N "sdp" -s M --mg ${systemSize} -A 70 -H | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # Only the objective
    ./run -B -S 1 -N "onlyobj, 95%"   -p 95   -s M --mg ${systemSize} --known -H --shots ${shots} --onlyobj | tee -a data/measure.dat
    ./run -B -S 1 -N "onlyobj, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -H --shots ${shots} --onlyobj | tee -a data/measure.dat

    # SDP plus only the objective
    ./run -B -S 1 -N "sdp+onlyobj, 95%"   -p 95   -s M --mg ${systemSize} --known -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat
    ./run -B -S 1 -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat

    # Only level 2 shots
    ./run -B -S 1 -N "all2, 95%"   -p 95   -s M --mg ${systemSize} --known -H --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -S 1 -N "all2, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -H --shots ${shots} --all 2 | tee -a data/measure.dat

    # SDP plus level 2 shots
    ./run -B -S 1 -N "sdp+all2, 95%" -p 95 -s M --mg ${systemSize} --known -A ${A} -H --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -S 1 -N "sdp+all2, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -A ${A} -H --shots ${shots} --all 2 | tee -a data/measure.dat

    # Only auto 300
    ./run -B -S 1 -N "auto300, 95%" -p 95 -s M --mg ${systemSize} --known -H --shots ${shots} --auto 300 | tee -a data/measure.dat
    ./run -B -S 1 -N "auto300, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -H --shots ${shots} --auto 300 | tee -a data/measure.dat

    # SDP plus auto 300
    ./run -B -S 1 -N "sdp+auto300, 95%" -p 95 -s M --mg ${systemSize} --known -A ${A} -H --shots ${shots} --auto 300 | tee -a data/measure.dat
    ./run -B -S 1 -N "sdp+auto300, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -A ${A} -H --shots ${shots} --auto 300 | tee -a data/measure.dat

    # Only first 300
    ./run -B -S 1 -N "first300, 95%" -p 95 -s M --mg ${systemSize} --known -H --shots ${shots} --first 300 | tee -a data/measure.dat
    ./run -B -S 1 -N "first300, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -H --shots ${shots} --first 300 | tee -a data/measure.dat

    # SDP plus first 300
    ./run -B -S 1 -N "sdp+first300, 95%" -p 95 -s M --mg ${systemSize} --known -A ${A} -H --shots ${shots} --first 300 | tee -a data/measure.dat
    ./run -B -S 1 -N "sdp+first300, 99.7%" -p 99.7 -s M --mg ${systemSize} --known -A ${A} -H --shots ${shots} --first 300 | tee -a data/measure.dat

done

# File showing the different confidence levels TODO

# Automatically commit once done
git add .
git commit -m "automatic data commit"
git push



