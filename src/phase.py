#!/usr/bin/env python

# simple python script to run ./run --phase gamma for many gammas
import numpy as np
import os
import matplotlib.pyplot as plt

# Define the range of gamma values
gammas = np.linspace(0.1, 1.0, 30)

# Loop over the gamma values
boundsPerGamma = []
for gamma in gammas:
    # os.system(f'./run --phase {gamma} -M 1000 > data.dat')
    os.system(f'./run --phase2 1 {gamma} 1.5 1 8 -A 50 -M 1000 > data.dat')
    res = open('data.dat').readlines()
    bounds = [-1000000,-100000]
    for line in res:
        if 'Bounds:' in line:
            split_line = line.split()
            bounds = [float(split_line[1]), float(split_line[3])]
    print(f'Gamma: {gamma}, Bounds: {bounds[0]}, {bounds[1]}')
    boundsPerGamma.append(bounds)

# Plot the result
plt.plot(gammas, boundsPerGamma)
plt.xlabel('Gamma')
plt.ylabel('Bounds')
plt.show()


