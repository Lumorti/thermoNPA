import matplotlib.pyplot as plt
import numpy as np
import sys

# Set the font to the LaTeX font
plt.rcParams['font.family'] = 'serif'

# Load data from the first argument
data1 = np.loadtxt("data/graph1.dat")
data2 = np.loadtxt("data/graph2.dat")
data3 = np.loadtxt("data/graph3.dat")

# Increase font size
plt.rcParams.update({'font.size': 15})

# Create a figure and primary axis
fig, ax1 = plt.subplots()

# Bounds rather than difference as above
mids1 = (data1[:, 1] + data1[:, 2]) / 2
err1 = (data1[:, 2] - data1[:, 1]) / 2
mids2 = (data2[:, 1] + data2[:, 2]) / 2
err2 = (data2[:, 2] - data2[:, 1]) / 2
mids3 = (data3[:, 1] + data3[:, 2]) / 2
err3 = (data3[:, 2] - data3[:, 1]) / 2
markers1, caps1, bars1 = ax1.errorbar(data1[:, 0], mids1, yerr=err1, linestyle='-', color='black', fmt='none', label='Constraint set 1 (10s)', capsize=0, elinewidth=8, capthick=1)
markers2, caps2, bars2 = ax1.errorbar(data2[:, 0], mids2, yerr=err2, linestyle='-', color='red', fmt='none',   label='Constraint set 2 (1m)', capsize=0, elinewidth=8, capthick=1)
markers3, caps3, bars3 = ax1.errorbar(data3[:, 0], mids3, yerr=err3, linestyle='-', color='blue', fmt='none', label='Constraint set 3 (1h)', capsize=0, elinewidth=8, capthick=1)
markers1c, caps1c, bars1c = ax1.errorbar(data1[:, 0], mids1, yerr=err1, linestyle='-', color='black', fmt='none', capsize=0, elinewidth=2, capthick=1)
markers2c, caps2c, bars2c = ax1.errorbar(data2[:, 0], mids2, yerr=err2, linestyle='-', color='red', fmt='none', capsize=0, elinewidth=1, capthick=1)
markers3c, caps3c, bars3c = ax1.errorbar(data3[:, 0], mids3, yerr=err3, linestyle='-', color='blue', fmt='none',  capsize=0, elinewidth=1, capthick=1)
ax1.legend()
for bar in bars1:
    bar.set_alpha(0.2)
for bar in bars2:
    bar.set_alpha(0.1)
for bar in bars3:
    bar.set_alpha(0.1)
ax1.set_xlabel('Number of qubits, n', labelpad=15)
ax1.set_ylabel('Bounded interval, $[\\langle M \\rangle_{lb}, \\langle M \\rangle_{ub}]$', color='black', labelpad=15)
ax1.tick_params(axis='y', labelcolor='black')

# Also plot a point at each max and min
deltaMax = np.ones(len(data1)) * (-0.034)
deltaMin = np.ones(len(data1)) * (0.016)
ax1.plot(data1[:, 0], data1[:, 1]+deltaMin, marker="$\u23B5$", color='black', markersize=10, linestyle='None')
ax1.plot(data1[:, 0], data1[:, 2]+deltaMax, marker="$\u23B4$", color='black', markersize=10, linestyle='None')
ax1.plot(data2[:, 0], data2[:, 1]+deltaMin, marker="$\u23B5$", color='red',   markersize=10, linestyle='None')
ax1.plot(data2[:, 0], data2[:, 2]+deltaMax, marker="$\u23B4$", color='red',   markersize=10, linestyle='None')
ax1.plot(data3[:, 0], data3[:, 1]+deltaMin, marker="$\u23B5$", color='blue',  markersize=10, linestyle='None')
ax1.plot(data3[:, 0], data3[:, 2]+deltaMax, marker="$\u23B4$", color='blue',  markersize=10, linestyle='None')

# scale should go from -1 to 1 on the y
ax1.set_ylim(-1, 1)

# Add grid and show the plot
plt.tight_layout()
plt.savefig('data/scaling.pdf', bbox_inches='tight')
plt.show()
