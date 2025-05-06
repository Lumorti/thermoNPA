#!/usr/bin/env python

# simple python script to run ./run --phase gamma for many gammas
import numpy as np
import os
import matplotlib.pyplot as plt

# Define the range of gamma values
Jxs = np.linspace(-4.0, 4.0, 15)
Jys = np.linspace(-4.0, 4.0, 15)

# Other parameters
unknownVal = -1000000

# The upper and lower bound matrices
uppers = np.zeros((len(Jxs), len(Jys)))
avgs = np.zeros((len(Jxs), len(Jys)))
lowers = np.zeros((len(Jxs), len(Jys)))

# Loop over the gamma values
# os.system(f'rm data/monoms_*.csv')
# for x, Jx in enumerate(Jxs):
    # for y, Jy in enumerate(Jys):
        # # os.system(f'./run -S 2 --phase2 2 2 {Jx} {Jy} --objRand 2 -M 0 > data.dat')
        # # os.system(f'./run -S 3 --phase2 3 2 {Jx} {Jy} --objRand 2 -M 2000 -r 2 > data.dat')
        # # os.system(f'./run -S 3 --phase2 3 2 {Jx} {Jy} --objRand 2 -M 0 -F > data.dat')
        # # os.system(f'./run -S 2 --phase2 3 2 {Jx} {Jy} --objRand 2 -M 1000 -A 30 -s M > data.dat')
        # os.system(f'./run -S 1 --phase2 4 4 {Jx} {Jy} --objRand 2 -M 5000 -r 2 -s M -F > data.dat')
        # os.system(f'mv monoms.csv data/monoms_{x}_{y}.csv')
        # res = open('data.dat').readlines()
        # bounds = [unknownVal, unknownVal]
        # for line in res:
            # if 'Bounds:' in line:
                # split_line = line.split()
                # bounds = [float(split_line[1]), float(split_line[3])]
        # print(f'Jx: {Jx}, Jy: {Jy}, Bounds: {bounds}')
        # uppers[x, y] = bounds[0]
        # lowers[x, y] = bounds[1]

# Instead, just lot from the monoms.csv files
monoms = np.genfromtxt(f'data/monoms_0_0.csv', delimiter=',')
numMonoms = monoms.shape[0]
obj = np.zeros(numMonoms)
for i in range(numMonoms):
    if monoms[i, 0] == 1:
        obj[i] = np.random.uniform(-1, 1)
for x, Jx in enumerate(Jxs):
    for y, Jy in enumerate(Jys):
        monoms = np.genfromtxt(f'data/monoms_{x}_{y}.csv', delimiter=',')
        if monoms.shape[0] != numMonoms:
            continue
        val = np.dot(monoms[:, 2], obj)
        uppers[x, y] = val
        lowers[x, y] = val

# Get the average over all elements
upperAvg = 0
lowerAvg = 0
for x in range(len(Jxs)):
    for y in range(len(Jys)):
        if uppers[x, y] != unknownVal:
            upperAvg += uppers[x, y]
        if lowers[x, y] != unknownVal:
            lowerAvg += lowers[x, y]
upperAvg /= len(Jxs) * len(Jys)
lowerAvg /= len(Jxs) * len(Jys)

# Set all the -1000000 values to the average
uppers[uppers == unknownVal] = upperAvg
lowers[lowers == unknownVal] = lowerAvg

# Set the average values
avgs = (uppers + lowers) / 2

# Plot the result
plt.figure(figsize=(8, 8))
plt.subplot(131)
plt.imshow(uppers, extent=[Jxs[0], Jxs[-1], Jys[0], Jys[-1]], origin='lower', aspect='equal')
plt.title('Upper Bound')
plt.xlabel('Jx')
plt.ylabel('Jy')
plt.subplot(132)
plt.imshow(avgs, extent=[Jxs[0], Jxs[-1], Jys[0], Jys[-1]], origin='lower', aspect='equal')
plt.title('Average Bound')
plt.xlabel('Jx')
plt.ylabel('Jy')
plt.subplot(133)
plt.imshow(lowers, extent=[Jxs[0], Jxs[-1], Jys[0], Jys[-1]], origin='lower', aspect='equal')
plt.title('Lower Bound')
plt.xlabel('Jx')
plt.ylabel('Jy')
plt.tight_layout()
plt.savefig('phases.png')
plt.show()

