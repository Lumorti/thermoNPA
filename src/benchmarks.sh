#!/bin/bash

# For the two-qubit system
#echo "" > data/twoqubit.tex
#./run -s M --two -l 1 -B | tee -a data/twoqubit.tex
#./run -s M --two -m 1 -l 1 -B | tee -a data/twoqubit.tex
#./run -s M --two -l 2 -B | tee -a data/twoqubit.tex
#./run -s M --two -M 0 -B | tee -a data/twoqubit.tex

# For the tensor chain
#echo "" > data/tensorchain.tex
#./run -s M --tensor 12 -M 1000 -B | tee -a data/tensorchain.tex
#./run -s M --tensor 12 -M 10000 -B | tee -a data/tensorchain.tex
#./run -s M --tensor 12 -M 10000 -m 1 -B | tee -a data/tensorchain.tex
#./run -s M --tensor 12 -M 10000 -A 50 -B | tee -a data/tensorchain.tex
#./run -s M --tensor 12 -M 10000 -r 4 -B | tee -a data/tensorchain.tex
#./run -s M --tensor 12 -M 10000 -r 4 -A 50 -B | tee -a data/tensorchain.tex
#./run -s M --tensor 12 -M 10000 -r 4 -Y -B | tee -a data/tensorchain.tex
#./run -s M --tensor 12 -M 15000 -r 4 -Y -B | tee -a data/tensorchain.tex
#./run -s M --tensor 12 -M 30000 -r 4 -Y -B | tee -a data/tensorchain.tex

# For the 2D case
#echo "" > data/david2d.tex
#./run -s M --david2d 0.5 2 5 -M 1000 -B | tee -a data/david2d.tex
#./run -s M --david2d 0.5 2 5 -M 10000 -B | tee -a data/david2d.tex
#./run -s M --david2d 0.5 2 5 -M 10000 -r 4 -B | tee -a data/david2d.tex
#./run -s M --david2d 0.5 2 5 -M 10000 -A 100 -B | tee -a data/david2d.tex
#./run -s M --david2d 0.5 5 2 -M 30000 -A 150 -B | tee -a data/david2d.tex
#./run -s M --david2d 0.5 5 2 -M 50000 -A 200 -B | tee -a data/david2d.tex
#./run -s M --david2d 0.5 5 2 -M 50000 -A 200 -y 1 6 -y 2 7 -y 3 8 -y 4 9 -y 5 10 -B | tee -a data/david2d.tex
#./run -s M --david2d 0.5 5 2 -M 50000 -A 200 -r 4 -y 1 6 -y 2 7 -y 3 8 -y 4 9 -y 5 10 -B | tee -a data/david2d.tex
#./run -s M --david2d 0.5 5 2 -M 70000 -A 250 -y 1 6 -y 2 7 -y 3 8 -y 4 9 -y 5 10 -B | tee -a data/david2d.tex
./run -s M --david2d 0.5 5 2 -M 100000 -A 250 -y 1 6 -y 2 7 -y 3 8 -y 4 9 -y 5 10 -B | tee -a data/david2d.tex
./run -s M --david2d 0.5 5 2 -M 100000 -A 300 -y 1 6 -y 2 7 -y 3 8 -y 4 9 -y 5 10 -B | tee -a data/david2d.tex
./run -s M --david2d 0.5 5 2 -M 100000 -r 4 -y 1 6 -y 2 7 -y 3 8 -y 4 9 -y 5 10 -B | tee -a data/david2d.tex
./run -s M --david2d 0.5 5 2 -M 100000 -r 5 -y 1 6 -y 2 7 -y 3 8 -y 4 9 -y 5 10 -B | tee -a data/david2d.tex

# Graph versus system size
#echo "" > data/graph1.dat
#echo "" > data/graph1.tex
#echo "" > data/graph2.dat
#echo "" > data/graph2.tex
#echo "" > data/graph3.dat
#echo "" > data/graph3.tex
#for i in 12 15 20 25 30 35 40 45 50 60 70 80 90 100 120 140 160 180 200
#do

    #results=$(./run -s M --tensor $i -M 5000 -A 25 -Y -B)
    #diff=$(echo "$results" | awk '{print $12}')
    #echo "$results" | tee -a data/graph1.tex
    #echo "$i $diff" | tee -a data/graph1.dat

    #results=$(./run -s M --tensor $i -M 10000 -A 50 -Y -B)
    #diff=$(echo "$results" | awk '{print $12}')
    #echo "$results" | tee -a data/graph2.tex
    #echo "$i $diff" | tee -a data/graph2.dat

    #results=$(./run -s M --tensor $i -M 20000 -A 100 -Y -B)
    #diff=$(echo "$results" | awk '{print $12}')
    #echo "$results" | tee -a data/graph3.tex
    #echo "$i $diff" | tee -a data/graph3.dat

#done

git add .
git commit -m "Updated data"
git push

