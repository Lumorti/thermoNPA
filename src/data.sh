#!/bin/bash

# This script is used to generate the data for the example
mkdir -p data
> data/X1.dat
> data/Z1.dat
> data/X1X2.dat
> data/Z1Z2.dat
> data/X5X6.dat
> data/Z5Z6.dat
> data/X5.dat
> data/Y5.dat
> data/Z5.dat
for g in $(LANG=en_US seq 0.1 0.1 3)
do
    #valZ1Z2=$(./run --david $g  -M 20 -m 2 -O "<Z1Z2>" | awk '/Bounds/ {print $2, $4}')
    #valZ2=$(./run --david $g    -M 20 -m 2 -O "<Z2>" | awk '/Bounds/ {print $2, $4}')
    #valX2=$(./run --david $g    -M 20 -m 2 -O "<X2>" | awk '/Bounds/ {print $2, $4}')
    #valX1X2=$(./run --david $g  -M 20 -m 2 -O "<X1X2>" | awk '/Bounds/ {print $2, $4}')
    valX1=$(./run --david $g    -M 1000 -O "<X1>" | awk '/Bounds/ {print $2, $4}')
    valZ1=$(./run --david $g    -M 1000 -O "<Z1>" | awk '/Bounds/ {print $2, $4}')
    valX1X2=$(./run --david $g  -M 1000 -O "<X1X2>" | awk '/Bounds/ {print $2, $4}')
    valZ1Z2=$(./run --david $g  -M 1000 -O "<Z1Z2>" | awk '/Bounds/ {print $2, $4}')
    valX5X6=$(./run --david $g  -M 1000 -O "<X5X6>" | awk '/Bounds/ {print $2, $4}')
    valZ5Z6=$(./run --david $g  -M 1000 -O "<Z5Z6>" | awk '/Bounds/ {print $2, $4}')
    valX5=$(./run --david $g    -M 1000 -O "<X5>" | awk '/Bounds/ {print $2, $4}')
    valY5=$(./run --david $g    -M 1000 -O "<Y5>" | awk '/Bounds/ {print $2, $4}')
    valZ5=$(./run --david $g    -M 1000 -O "<Z5>" | awk '/Bounds/ {print $2, $4}')

    #echo $g $valZ1Z2 | tee -a data/Z1Z2.dat
    #echo $g $valZ2 | tee -a data/Z2.dat
    #echo $g $valX2 | tee -a data/X2.dat
    #echo $g $valX1X2 | tee -a data/X1X2.dat
    echo $g $valX1 | tee -a data/X1.dat
    echo $g $valZ1 | tee -a data/Z1.dat
    echo $g $valX1X2 | tee -a data/X1X2.dat
    echo $g $valZ1Z2 | tee -a data/Z1Z2.dat
    echo $g $valX5X6 | tee -a data/X5X6.dat
    echo $g $valZ5Z6 | tee -a data/Z5Z6.dat
    echo $g $valX5 | tee -a data/X5.dat
    echo $g $valY5 | tee -a data/Y5.dat
    echo $g $valZ5 | tee -a data/Z5.dat


done

# Plot the data (lower and upper bounds)
#gnuplot -p -e "plot 'data/X2.dat' u 1:2 w l t 'X2 lower bound', 
            #'data/Z2.dat' u 1:2 w l t 'Z2 lower bound',
            #'data/Z1Z2.dat' u 1:2 w l t 'Z1Z2 lower bound',
            #'data/X1X2.dat' u 1:2 w l t 'X1X2 lower bound',
            #'data/X2.dat' u 1:3 w l t 'X2 upper bound',
            #'data/Z2.dat' u 1:3 w l t 'Z2 upper bound',
            #'data/Z1Z2.dat' u 1:3 w l t 'Z1Z2 upper bound',
            #'data/X1X2.dat' u 1:3 w l t 'X1X2 upper bound'"
gnuplot -p -e "plot 'data/X1.dat' u 1:2 w l t 'X1 lower bound', 
            'data/Z1.dat' u 1:2 w l t 'Z1 lower bound',
            'data/X1X2.dat' u 1:2 w l t 'X1X2 lower bound',
            'data/Z1Z2.dat' u 1:2 w l t 'Z1Z2 lower bound',
            'data/X5X6.dat' u 1:2 w l t 'X5X6 lower bound',
            'data/Z5Z6.dat' u 1:2 w l t 'Z5Z6 lower bound',
            'data/X5.dat' u 1:2 w l t 'X5 lower bound',
            'data/Y5.dat' u 1:2 w l t 'Y5 lower bound',
            'data/Z5.dat' u 1:2 w l t 'Z5 lower bound',
            'data/X1.dat' u 1:3 w l t 'X1 upper bound',
            'data/Z1.dat' u 1:3 w l t 'Z1 upper bound',
            'data/X1X2.dat' u 1:3 w l t 'X1X2 upper bound',
            'data/Z1Z2.dat' u 1:3 w l t 'Z1Z2 upper bound',
            'data/X5X6.dat' u 1:3 w l t 'X5X6 upper bound',
            'data/Z5Z6.dat' u 1:3 w l t 'Z5Z6 upper bound',
            'data/X5.dat' u 1:3 w l t 'X5 upper bound',
            'data/Y5.dat' u 1:3 w l t 'Y5 upper bound',
            'data/Z5.dat' u 1:3 w l t 'Z5 upper bound'"


