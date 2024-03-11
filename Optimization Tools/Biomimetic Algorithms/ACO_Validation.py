# Jackson Morgan
# Biomimetic Algorithms
# Task 1 - Traveling Sales-Ant Problem

import ACO
import numpy as np
import statistics as stat
import matplotlib.pyplot as plt

# Algorithm Parameters:
N = 30 # number of sales-ants
n = 200 # number of trials
rho = 0.5 # evaporation rate


# Step 1 - Define Problem Attributes:
nodes = [0,1,2,3,4]
d = np.matrix([[0,10,12,11,14], [10,0,13,15,8], [12,13,0,9,14], [11,15,9,0,16], [14,8,14,16,0]])
tau = np.ones((5,5)) # initial pheromones are all 1


# Loop through iterations:
meanval = []
for i in range(n):

    # generate ants' paths
    paths = []
    for i in range(N):
        paths.append(ACO.path(nodes, d, tau))
    paths = np.matrix(paths)


    # Update Pheromones:
    tau = ACO.evaporate(rho, tau) # evaporate
    for i in range(N): # add pheromones from each ant's path
        path = paths[i,:]
        delta_tau = ACO.pheromone(nodes, path, d)
        tau = tau + delta_tau


    # Record Path Value:
    totval = []
    for i in range(N):
        val = ACO.value(nodes, paths[i,:], d)
        totval.append(sum(val))
    meanval.append(stat.mean(totval))


# Plot Results:
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(meanval)
plt.grid()
plt.xlabel('Trial Number')
plt.ylabel('Average Distance')
plt.title('N = ' +str(N) +'      Rho = ' +str(rho))
plt.show()

























