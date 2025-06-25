#!/bin/bash

> data.dat

val3=$(./run --tensor 3 -s M -l 3 -I 2 -Y -r 3 -B)
echo "val3: $val3"
val4=$(./run --tensor 4 -s M -l 3 -I 2,3 -Y -r 3 -B)
echo "val4: $val4"
val5=$(./run --tensor 5 -s M -l 3 -I 2,3,4 -Y -r 3 -B)
echo "val5: $val5"
val6=$(./run --tensor 6 -s M -l 3 -I 2,3,4,5 -Y -r 3 -B)
echo "val6: $val6"
val7=$(./run --tensor 7 -s M -l 3 -I 2,3,4,5,6 -Y -r 3 -B)
echo "val7: $val7"
val8=$(./run --tensor 8 -s M -l 3 -I 2,3,4,5,6,7 -Y -r 3 -B)
echo "val8: $val8"
val9=$(./run --tensor 9 -s M -l 3 -I 2,3,4,5,6,7,8 -Y -r 3 -B)
echo "val9: $val9"
val10=$(./run --tensor 10 -s M -l 3 -I 2,3,4,5,6,7,8,9 -Y -r 3 -B)
echo "val10: $val10"
