#!/bin/bash

#> data.dat
#for A in 10 20 30 40 50 60 70 80 90 100 200
##for A in 110 120 130 140 150 160 170 180 190
#do
    ##for M in 10 100 1000 2000 5000 7500 10000 20000
    #for M in 25000 30000 35000
    #do
        #./run --tensor 12 -M $M -A $A > temp.log
        #diff=$(cat temp.log | grep "Difference" | awk '{print $2}')
        #timeSolve=$(cat temp.log | grep "Time to solve" | awk '{print $4}')
        #timeGen=$(cat temp.log | grep "Time to generate" | awk '{print $4}')
        #timeSolve=${timeSolve::-1}
        #timeGen=${timeGen::-1}
        #echo $A $M $diff $timeSolve $timeGen | tee -a data.dat
    #done
#done

