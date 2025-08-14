#!/bin/bash

# Data generation for the measurement project
> data/measure.dat

# System settings
systemSize="3 3"
A=50
M=1000

# Pre-compute the steady state
filename="data/2d_3x3.dat"
./run --2dtfi ${systemSize} --precompute ${filename}

# Heat current without Z vs number of measurements
echo "file & heat" | tee -a data/measure.dat
./run -B -N "sdp" -s M --2dtfi ${systemSize} -M ${M} -A ${A} --objHC H | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # SDP plus level 2 shots
    ./run -B -N "sdp+all2-z, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -N "sdp+all2-z, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 2 | tee -a data/measure.dat

    # SDP plus level 3 shots
    ./run -B -N "sdp+all3-z, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 3 | tee -a data/measure.dat
    ./run -B -N "sdp+all3-z, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 3 | tee -a data/measure.dat

    # SDP plus 50 automatic
    ./run -B -N "sdp+auto50-z, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 50 | tee -a data/measure.dat
    ./run -B -N "sdp+auto50-z, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 50 | tee -a data/measure.dat

    # SDP plus only objective
    ./run -B -N "sdp+onlyobj-z, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --onlyobj | tee -a data/measure.dat
    ./run -B -N "sdp+onlyobj-z, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --onlyobj | tee -a data/measure.dat

done


# Purity versus number of measurements
echo "file & purity" | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # Only level 3 shots
    ./run -B -N "all3, 95%" -s M -p 95 --2dtfi ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 3 | tee -a data/measure.dat
    ./run -B -N "all3, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filename} --millis -S 1 --objPurity --shots ${shots} --all 3 | tee -a data/measure.dat

    # SDP plus level 3 shots
    ./run -B -N "sdp+all3, 95%" -s M -p 95 --2dtfi ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 3 | tee -a data/measure.dat
    ./run -B -N "sdp+all3, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 3 | tee -a data/measure.dat

    # SDP plus auto 100
    ./run -B -N "sdp+auto100, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 100 | tee -a data/measure.dat
    ./run -B -N "sdp+auto100, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 100 | tee -a data/measure.dat

    # SDP plus auto 300
    ./run -B -N "sdp+auto300, 95%" -p 95 -s M --2dtfi ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 300 | tee -a data/measure.dat
    ./run -B -N "sdp+auto300, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --auto 300 | tee -a data/measure.dat

done

# Pre-compute the ground state
filename="data/2d_3x3_H.dat"
./run --2dtfi ${systemSize} --precompute ${filename} -H

# Ground state energy vs number of measurements
echo "file & energy" | tee -a data/measure.dat
./run -B -N "sdp" -s M --2dtfi ${systemSize} -A ${A} -H | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # Just level 2 shots (these are left as comments using // in the data file)
    val=$(./run -B -N "all2, 95%" -s M -p 95 --2dtfi ${systemSize} --precomputed ${filename} -H -S 1 --shots ${shots} --all 2)
    val=$(echo "//${val}")
    echo ${val} | tee -a data/measure.dat
    val=$(./run -B -N "all2, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filename} -H -S 1 --shots ${shots} --all 2)
    val=$(echo "//${val}")
    echo ${val} | tee -a data/measure.dat

    # Only measure the objective
    ./run -B -N "sdp+onlyobj, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat
    ./run -B -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat

    # SDP plus level 2 shots
    ./run -B -N "sdp+all2, 95%"   -p 95   -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 2 | tee -a data/measure.dat
    ./run -B -N "sdp+all2, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --all 2 | tee -a data/measure.dat

    # SDP plus 200 automatic
    ./run -B -N "sdp+auto200, 95%"   -s M -p 95   --2dtfi ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 200 | tee -a data/measure.dat
    ./run -B -N "sdp+auto200, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filename} -S 1 -A ${A} -H --shots ${shots} --auto 200 | tee -a data/measure.dat

done

# Ground state energy for large systems using symmetries
echo "file & large" | tee -a data/measure.dat
systemSize="3 30"
./run -B -N "sdp" -s M --2dtfiperiodic ${systemSize} -S 1 -m 3 -H -O "<Z1>+<Z2>+<Z3>+<X1>+<X2>+<X3>" -P --shots ${shots} --sym | tee -a data/measure.dat
for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
do

    # SDP plus symmetry
    ./run -B -N "sdp+sym, 68%" -p 68 -Y -s M --2dtfiperiodic ${systemSize} -S 1 -m 3 -H -O "<Z1>+<Z2>+<Z3>+<X1>+<X2>+<X3>" -P --shots ${shots} --sym | tee -a data/measure.dat
    ./run -B -N "sdp+sym, 95%" -p 95 -Y -s M --2dtfiperiodic ${systemSize} -S 1 -m 3 -H -O "<Z1>+<Z2>+<Z3>+<X1>+<X2>+<X3>" -P --shots ${shots} --sym | tee -a data/measure.dat
    ./run -B -N "sdp+sym, 99.7%" -p 99.7 -Y -s M --2dtfiperiodic ${systemSize} -S 1 -m 3 -H -O "<Z1>+<Z2>+<Z3>+<X1>+<X2>+<X3>" -P --shots ${shots} --sym | tee -a data/measure.dat

done

# Automatically commit once done
git add .
git commit -m "automatic data commit"
git push



