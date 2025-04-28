import matplotlib.pyplot as plt
import numpy as np
import sys

# Load data from the first argument
data1 = np.loadtxt("data/graph1.dat")
data2 = np.loadtxt("data/graph2.dat")
data3 = np.loadtxt("data/graph3.dat")

# Increase font size
plt.rcParams.update({'font.size': 15})

# Create a figure and primary axis
fig, ax1 = plt.subplots()

# Plot data on the primary y-axis as points
# ax1.plot(data1[:, 0], data1[:, 1], linestyle='none', marker='o', color='black', label='Set 1')
# ax1.plot(data2[:, 0], data2[:, 1], linestyle='none', marker='x', color='red', label='Set 2')
# ax1.plot(data3[:, 0], data3[:, 1], linestyle='none', marker='^', color='blue', label='Set 3')
# ax1.set_xlabel('System size', labelpad=15)
# ax1.set_ylabel('Bound difference %', color='black', labelpad=15)
# ax1.tick_params(axis='y', labelcolor='black')

# Bounds rather than difference as above
mids1 = (data1[:, 1] + data1[:, 2]) / 2
err1 = (data1[:, 2] - data1[:, 1]) / 2
mids2 = (data2[:, 1] + data2[:, 2]) / 2
err2 = (data2[:, 2] - data2[:, 1]) / 2
mids3 = (data3[:, 1] + data3[:, 2]) / 2
err3 = (data3[:, 2] - data3[:, 1]) / 2
ax1.errorbar(data1[:, 0], mids1, yerr=err1, linestyle='-', color='black', fmt='none', label='Set 1 (10s)', capsize=5)
ax1.errorbar(data2[:, 0], mids2, yerr=err2, linestyle='-', color='red', fmt='none', label='Set 2 (1m)', capsize=5)
ax1.errorbar(data3[:, 0], mids3, yerr=err3, linestyle='-', color='blue', fmt='none', label='Set 3 (1hr)', capsize=5)
ax1.set_xlabel('System Size', labelpad=15)
ax1.set_ylabel('Bounds', color='black', labelpad=15)
ax1.tick_params(axis='y', labelcolor='black')

# Create a secondary y-axis
# ax2 = ax1.twinx()

# Apply a transformation or use a different scale for the secondary y-axis
# transformed_data1 = np.array([int(x * 50) for x in data1[:, 1]])
# transformed_data2 = np.array([int(x * 50) for x in data2[:, 1]])
# transformed_data3 = np.array([int(x * 50) for x in data3[:, 1]])

# Plot transformed data on the secondary y-axis
# ax2.plot(data1[:, 0], transformed_data1, linestyle='none')
# ax2.plot(data2[:, 0], transformed_data2, linestyle='none')
# ax2.plot(data3[:, 0], transformed_data3, linestyle='none')
# ax2.set_ylabel('Bound difference (%)', color='black', labelpad=15)
# ax2.tick_params(axis='y', labelcolor='black')

# ax2.spines['top'].set_visible(False)
# ax2.spines['right'].set_visible(False)
# ax2.spines['bottom'].set_visible(False)
# ax2.spines['left'].set_visible(False)

# Add grid and show the plot
ax1.grid()
ax1.legend()
plt.tight_layout()
plt.savefig('data/scaling.pdf', bbox_inches='tight')
