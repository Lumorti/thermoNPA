import matplotlib.pyplot as plt
import numpy as np
import sys

# Load data from the first argument
data1 = np.loadtxt("data/graph1.dat")
data2 = np.loadtxt("data/graph2.dat")
data3 = np.loadtxt("data/graph3.dat")

# Increase font size
plt.rcParams.update({'font.size': 16})

# Create a figure and primary axis
fig, ax1 = plt.subplots()

# Plot data on the primary y-axis as points
# ax1.plot(data[:, 0], data[:, 1], color='black')
ax1.plot(data1[:, 0], data1[:, 1], linestyle='none', marker='o', color='black', label='constraint set 1')
ax1.plot(data2[:, 0], data2[:, 1], linestyle='none', marker='x', color='red', label='constraint set 2')
ax1.plot(data3[:, 0], data3[:, 1], linestyle='none', marker='^', color='blue', label='constraint set 3')
ax1.set_xlabel('system size', labelpad=15)
ax1.set_ylabel('bound difference', color='black', labelpad=15)
ax1.tick_params(axis='y', labelcolor='black')

# Create a secondary y-axis
ax2 = ax1.twinx()

# Apply a transformation or use a different scale for the secondary y-axis
transformed_data1 = np.array([int(x * 50) for x in data1[:, 1]])
transformed_data2 = np.array([int(x * 50) for x in data2[:, 1]])
transformed_data3 = np.array([int(x * 50) for x in data3[:, 1]])

# Plot transformed data on the secondary y-axis
ax2.plot(data1[:, 0], transformed_data1, linestyle='none')
ax2.plot(data2[:, 0], transformed_data2, linestyle='none')
ax2.plot(data3[:, 0], transformed_data3, linestyle='none')
ax2.set_ylabel('bound difference (%)', color='black', labelpad=15)
ax2.tick_params(axis='y', labelcolor='black')

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['bottom'].set_visible(False)
ax2.spines['left'].set_visible(False)

# Add grid and show the plot
ax1.grid()
ax1.legend()
plt.tight_layout()
plt.savefig('data/scaling.pdf', bbox_inches='tight')
