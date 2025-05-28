#!/bin/bash


# If we just have 10000 first-order X measurements
# and we want to get bounds on many second-order quantities
objSeed=1
for i in {1..10}
do
    shotsSeed=$i
    ./run -S ${objSeed} --tensor 5 -M 100 -A 50 --objRand 2 --shots 10000 --all 1 --noy --noz -S ${shotsSeed} -B
done
./run -S ${objSeed} --tensor 5 -M 100 -A 50 --objPurity --shots 10000 --all 1 --noy --noz -S ${shotsSeed} -B
