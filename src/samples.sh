#!/bin/bash

# System settings
systemSize="3 3"
A=150
M=2000
numRepeats=50
filename="data/2d_3x3.dat"
filenameEnergy="data/2d_3x3_H.dat"
shellCon=" --shell -9.897068167 0.09897068167 "

# Reset the files
# > data/precomputes.log
# > data/measure.dat

# Pre-compute both the steady state and the ground state
# ./run --2dtfi ${systemSize} --precompute ${filename} | tee -a data/precomputes.log
# ./run --2dtfi ${systemSize} --precompute ${filenameEnergy} -H | tee -a data/precomputes.log

# Heat current without something vs number of measurements
# echo "file & heat" | tee -a data/measure.dat
# ./run -B -N "sdp" -s M --2dtfi ${systemSize} -M ${M} -A ${A} --objHC H | tee -a data/measure.dat
# letter=z
# for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
# do
#     for ind in $(seq 1 $numRepeats)
#     do

#         # Only 100 first
#         ./run -B -N "first100-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S ${ind} --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat

#         # SDP plus 100 first
#         ./run -B -N "sdp+first100-${letter}, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filename} -S ${ind} -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat

#     done
# done

# Purity versus number of measurements
# echo "file & purity" | tee -a data/measure.dat
# ./run -B -N "sdp" -s M --2dtfi ${systemSize} -A ${A} --objPurity | tee -a data/measure.dat
# for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
# do
#     for ind in $(seq 1 $numRepeats)
#     do

#         # Only second-order shots
#         ./run -B -N "all2, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S ${ind} --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat

#         # SDP plus second-order shots
#         ./run -B -N "sdp+all2, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S ${ind} -A ${A} --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat

#     done
# done

# Ground state energy vs number of measurements
# echo "file & energy" | tee -a data/measure.dat
# ./run -B -N "sdp" -s M --2dtfi ${systemSize} -A ${A} -H | tee -a data/measure.dat
# for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
# do
#     for ind in $(seq 1 $numRepeats)
#     do

#         # SDP plus the 100 most common
#         ./run -B -N "sdp+auto100, 99.7%" -s M -p 99.7 --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S ${ind} -A ${A} -H --shots ${shots} --auto 100 | tee -a data/measure.dat

#         # SDP plus measure the objective
#         ./run -B -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S ${ind} -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat

#         # SDP plus level 2 shots
#         ./run -B -N "sdp+all2, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S ${ind} -A ${A} -H --shots ${shots} --all 2 | tee -a data/measure.dat

#     done
# done

# Local purity vs number of measurements
# echo "file & localpurity" | tee -a data/measure.dat
# ./run -B -N "sdp" -s M --2dtfi ${systemSize} -A ${A} -H --objLocalPurity 1 2 3 ${shellCon} | tee -a data/measure.dat
# for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
# do
#     for ind in $(seq 1 $numRepeats)
#     do

#         # Just measure the objective
#         ./run -B -N "onlyobj, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S ${ind} -H --objLocalPurity 1 2 3 --shots ${shots} --onlyobj ${shellCon} | tee -a data/measure.dat

#         # SDP plus measure the objective
#         ./run -B -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --2dtfi ${systemSize} --precomputed ${filenameEnergy} -S ${ind} -A ${A} -H --objLocalPurity 1 2 3 --shots ${shots} --onlyobj ${shellCon} | tee -a data/measure.dat

#     done
# done

# File showing the different confidence levels
# echo "file & confidence" | tee -a data/measure.dat
# for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000
# do
#     for conf in 68 95 99.7
#     do

#         # SDP plus level 2 shots
#         ./run -B -N "sdp+all2, ${conf}%" -s M -p ${conf} --2dtfi ${systemSize} --precomputed ${filenameEnergy} --millis -S 1 -A ${A} --objPurity --shots ${shots} --all 2 | tee -a data/measure.dat

#     done
# done

# Ground state energy for large systems
# echo "file & large" | tee -a data/measure.dat
# systemSizeLarge="50"
# A=200
# ./run -B -N "sdp" -s M --mg ${systemSizeLarge} -A ${A} -H | tee -a data/measure.dat
# for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
# do
#     for ind in $(seq 1 $numRepeats)
#     do

#         # Only the objective
#         ./run -B -S ${ind} -N "onlyobj, 99.7%" -p 99.7 -s M --mg ${systemSizeLarge} --known -H --shots ${shots} --onlyobj | tee -a data/measure.dat

#         # SDP plus only the objective
#         ./run -B -S ${ind} -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --mg ${systemSizeLarge} --known -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat

#     done
# done

# Ground state energy for the J1-J2 chain
# numRepeatsLarge2=1
# echo "file & large2" | tee -a data/measure.dat
# systemSize2="50"
# J2="0.35"
# chiMPS=128
# A=200
# filenameMPS="data/j1j2_${systemSize2}_${J2}.mps"
# python3 src/dmrg.py -n ${systemSize2} --j2 ${J2} --chi ${chiMPS} -o ${filenameMPS} | tee -a data/precomputes.log
# ./run -B -N "sdp" -s M --j1j2 ${systemSize2} ${J2} -A ${A} -H | tee -a data/measure.dat
# for shots in 10000 50000 100000 500000 1000000 5000000 10000000 50000000 100000000 -1
# do
#     for ind in $(seq 1 $numRepeatsLarge2)
#     do

#         # Only the objective
#         ./run -B -S ${ind} -N "onlyobj, 99.7%" -p 99.7 -s M --j1j2 ${systemSize2} ${J2} --mps ${filenameMPS} -H --shots ${shots} --onlyobj | tee -a data/measure.dat

#         # SDP plus only the objective
#         ./run -B -S ${ind} -N "sdp+onlyobj, 99.7%" -p 99.7 -s M --j1j2 ${systemSize2} ${J2} --mps ${filenameMPS} -A ${A} -H --shots ${shots} --onlyobj | tee -a data/measure.dat

#     done
# done

# Local violation of the second law for a 1D chain
echo "file & secondlaw" | tee -a data/measure.dat
systemSize3="9"
A=150
M=2000
numRepeatsSecondLaw=50
lkArgs="${systemSize3} 0.5 -0.99 0.05 0.05"
filenameLK="data/lkchain_${systemSize3}.dat"
[ -f ${filenameLK} ] || ./run --lkchain ${lkArgs} --precompute ${filenameLK} | tee -a data/precomputes.log
./run -B -N "sdp" -s M --lkchain ${lkArgs} -M ${M} -A ${A} --objHC H | tee -a data/measure.dat
letter=z
for shots in 10000 100000 1000000 10000000 100000000 500000000 1000000000 5000000000 10000000000 -1
do
    for ind in $(seq 1 $numRepeatsSecondLaw)
    do

        # Only 100 first
        ./run -B -N "first100-${letter}, 99.7%" -p 99.7 -s M --lkchain ${lkArgs} --precomputed ${filenameLK} -S ${ind} --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat

        # SDP plus 100 first
        ./run -B -N "sdp+first100-${letter}, 99.7%" -p 99.7 -s M --lkchain ${lkArgs} --precomputed ${filenameLK} -S ${ind} -M ${M} -A ${A} --objHC H --no${letter} --shots ${shots} --first 100 | tee -a data/measure.dat

    done
done

# Automatically commit once done
# git add .
# git commit -m "automatic data commit"
# git push



