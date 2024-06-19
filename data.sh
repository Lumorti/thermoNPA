#!/bin/bash


# This script is used to generate the data for the example
mkdir -p data
> data/X2.dat
> data/Z2.dat
> data/Z1Z2.dat
> data/X1X2.dat
for g in $(LANG=en_US seq 0.1 0.1 3)
do
    valZ1Z2=$(./run --david $g  -M 40 -O "<Z1Z2>" | awk '/Bounds/ {print $2, $4}')
    valZ2=$(./run --david $g    -M 40 -O "<Z2>" | awk '/Bounds/ {print $2, $4}')
    valX2=$(./run --david $g    -M 40 -O "<X2>" | awk '/Bounds/ {print $2, $4}')
    valX1X2=$(./run --david $g  -M 40 -O "<X1X2>" | awk '/Bounds/ {print $2, $4}')
    echo $g $valZ1Z2 | tee -a data/Z1Z2.dat
    echo $g $valZ2 | tee -a data/Z2.dat
    echo $g $valX2 | tee -a data/X2.dat
    echo $g $valX1X2 | tee -a data/X1X2.dat
done

# Plot the data (lower and upper bounds)
gnuplot -p -e "plot 'data/X2.dat' u 1:2 w l t 'X2 lower bound', 
            'data/Z2.dat' u 1:2 w l t 'Z2 lower bound',
            'data/Z1Z2.dat' u 1:2 w l t 'Z1Z2 lower bound',
            'data/X1X2.dat' u 1:2 w l t 'X1X2 lower bound',
            'data/X2.dat' u 1:3 w l t 'X2 upper bound',
            'data/Z2.dat' u 1:3 w l t 'Z2 upper bound',
            'data/Z1Z2.dat' u 1:3 w l t 'Z1Z2 upper bound',
            'data/X1X2.dat' u 1:3 w l t 'X1X2 upper bound'"

 



