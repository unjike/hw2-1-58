import matplotlib.pyplot as plt
import numpy as np

# Data to be plotted 
particle_num = [1000, 5000, 10000, 15000, 20000]
serial_On2 = [0.899033, 21.6698, 88.2551, 223.524,  367.433]
serial_On1 = [0.0974842, 0.583493, 1.24503, 1.95879, 2.63991]
openmp = [0.0698181, 0.411675, 0.442495, 0.802078, 0.964335]


# Plot
plt.figure(figsize=(8, 6))
plt.plot(particle_num, serial_On2, label='Serial On2', marker='o', color='r')
plt.plot(particle_num, serial_On1, label='Serial On1', marker='s', color='g')
plt.plot(particle_num, openmp, label='OpenMP', marker='^', color='b')

# Logarithmic scaling for both axes
plt.xscale('log')
plt.yscale('log')

# Labels and Title
plt.xlabel('Number of Particles (Log Transformed)')
plt.ylabel('Time (seconds) (Log Transformed)')
plt.title('Performance Comparison (Log Scale)')

# Grid and Legend
plt.grid(True, which='both', ls='--', linewidth=0.5)
plt.legend()

# Show plot
plt.tight_layout()
plt.show()
