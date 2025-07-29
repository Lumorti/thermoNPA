#!/bin/bash

# Settings
shotsSeed=1
numObjectives=50
systemSize="2 2"
A=20
M=25

# Data generation for the measurement project
> data/measure.dat

# Pre-compute the steady state
#./run --2d ${systemSize} --precompute 

# Heat current without Z vs number of measurements vs computational time
#echo "file & heat3d" | tee -a data/measure.dat
#declare -a arr=("-A 10 -M 15" "-A 20 -M 25" "-A 30 -M 35" "-A 40 -M 45" "-A 50 -M 55" "-M 50" "-A 10 -M 50" "-A 50" "-M 0" "-r 2" "-r 2 -M 50")
#for AM in "${arr[@]}"
#do

    ## Just the SDP
    #./run -N "sdp" -s M --2d ${systemSize} ${AM} --objHC H --noz -B | tee -a data/measure.dat

    #for shots in 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000
    #do

        ## SDP plus level 1 shots
        #./run -N "sdp+shots1" -s M --2d ${systemSize} --precomputed -S 1 ${AM} --objHC H --noz --shots ${shots} --all 1 -B | tee -a data/measure.dat

        ## SDP plus level 2 shots
        #./run -N "sdp+shots2" -s M --2d ${systemSize} --precomputed -S 1 ${AM} --objHC H --noz --shots ${shots} --all 2 -B | tee -a data/measure.dat

        ## SDP plus level 3 shots
        #./run -N "sdp+shots3" -s M --2d ${systemSize} --precomputed -S 1 ${AM} --objHC H --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat

        ## SDP plus level 4 shots
        #./run -N "sdp+shots4" -s M --2d ${systemSize} --precomputed -S 1 ${AM} --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat

        ## SDP plus 20 automatic
        #./run -N "sdp+auto20" -s M --2d ${systemSize} --precomputed -S 1 ${AM} --objHC H --noz --shots ${shots} --auto 20 -B | tee -a data/measure.dat

        ## SDP plus 50 automatic
        #./run -N "sdp+auto50" -s M --2d ${systemSize} --precomputed -S 1 ${AM} --objHC H --noz --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    #done

#done


# Heat current without Z vs number of measurements
#echo "file & heat" | tee -a data/measure.dat
#./run -N "sdp" -s M --2d ${systemSize} -M ${M} -A ${A} --objHC H -B | tee -a data/measure.dat
#for shots in 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 -1
#do

    ## Only 200 automatic
    ##./run -N "auto200" -s M --2d ${systemSize} --precomputed -S 1 --objHC H --noz --shots ${shots} --auto 200 -B | tee -a data/measure.dat

    ## Only level 4 shots
    ##./run -N "shots4" -s M --2d ${systemSize} --precomputed -S 1 --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat

    ## SDP plus level 1 shots
    ##./run -N "sdp+shots1" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 1 -B | tee -a data/measure.dat

    ## SDP plus level 2 shots
    ##./run -N "sdp+shots2" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 2 -B | tee -a data/measure.dat

    ## SDP plus level 3 shots
    #./run -N "sdp+shots3" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat

    ## SDP plus level 4 shots
    #./run -N "sdp+shots4" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --all 4 -B | tee -a data/measure.dat

    ## SDP plus 10 automatic
    ##./run -N "sdp+auto10" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 10 -B | tee -a data/measure.dat

    ## SDP plus 20 automatic
    #./run -N "sdp+auto20" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 20 -B | tee -a data/measure.dat

    ## SDP plus 50 automatic
    #./run -N "sdp+auto50" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    ## SDP plus 100 automatic
    ##./run -N "sdp+auto100" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    ## SDP plus 200 automatic
    ##./run -N "sdp+auto200" -s M --2d ${systemSize} --precomputed -S 1 -M ${M} -A ${A} --objHC H --noz --shots ${shots} --auto 200 -B | tee -a data/measure.dat

#done


# Purity versus number of measurements
#echo "file & purity" | tee -a data/measure.dat
#A=25
#M=50
#for shots in 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 50000000 -1
#do

    ## only level 1 shots
    #./run -N "shots1" -s M --2d ${systemSize} --precomputed --millis -S 1 --objPurity --shots ${shots} --all 1 -B | tee -a data/measure.dat

    ## only level 4 shots
    #./run -N "shots4" -s M --2d ${systemSize} --precomputed --millis -S 1 --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

    ## SDP plus level 1 shots
    #./run -N "sdp+shots1" -s M --2d ${systemSize} --precomputed --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 1 -B | tee -a data/measure.dat

    ## SDP plus level 4 shots
    #./run -N "sdp+shots4" -s M --2d ${systemSize} --precomputed --millis -S 1 -M ${M} -A ${A} --objPurity --shots ${shots} --all 4 -B | tee -a data/measure.dat

#done

# Pre-compute the ground state
systemSize="3 2"
A=20
./run --2d ${systemSize} --precompute -H

# Ground state energy vs number of measurements
echo "file & energy" | tee -a data/measure.dat
./run -N "sdp" -s M --2d ${systemSize} -A ${A} -H -B | tee -a data/measure.dat
for shots in 1000 5000 10000 50000 100000 500000 1000000 5000000 10000000 -1
do

    # Only 200 automatic
    #./run -N "auto200" -s M --2d ${systemSize} --precomputed -S 1 -H --shots ${shots} --auto 200 -B | tee -a data/measure.dat

    # Only level 4 shots
    #./run -N "shots4" -s M --2d ${systemSize} --precomputed -S 1 -H --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # Only measure the objective
    ./run -N "sdp+onlyobj" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --onlyobj -B | tee -a data/measure.dat

    # SDP plus level 1 shots
    #./run -N "sdp+shots1" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --all 1 -B | tee -a data/measure.dat

    # SDP plus level 2 shots
    ./run -N "sdp+shots2" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --all 2 -B | tee -a data/measure.dat

    # SDP plus level 3 shots
    ./run -N "sdp+shots3" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 3 shots, no x
    ./run -N "sdp+shots3-x" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --nox --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 3 shots, no y
    #./run -N "sdp+shots3-y" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --noy --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 3 shots, no z
    #./run -N "sdp+shots3-z" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --noz --shots ${shots} --all 3 -B | tee -a data/measure.dat

    # SDP plus level 4 shots
    #./run -N "sdp+shots4" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --all 4 -B | tee -a data/measure.dat

    # SDP plus 10 automatic
    #./run -N "sdp+auto10" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --auto 10 -B | tee -a data/measure.dat

    # SDP plus 20 automatic
    #./run -N "sdp+auto20" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --auto 20 -B | tee -a data/measure.dat

    # SDP plus 50 automatic
    ./run -N "sdp+auto50" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --auto 50 -B | tee -a data/measure.dat

    # SDP plus 100 automatic
    ./run -N "sdp+auto100" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus 100 automatic, no x
    ./run -N "sdp+auto100-x" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --nox --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus 100 automatic, no y
    #./run -N "sdp+auto100-y" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --noy --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus 100 automatic, no z
    #./run -N "sdp+auto100-z" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --noz --shots ${shots} --auto 100 -B | tee -a data/measure.dat

    # SDP plus 200 automatic
    #./run -N "sdp+auto200" -s M --2d ${systemSize} --precomputed -S 1 -A ${A} -H --shots ${shots} --auto 200 -B | tee -a data/measure.dat

done


