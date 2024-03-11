# Jackson Morgan
# Genetic Algorithm Test on Knapsack Problem

import GA
import random
import statistics as stat
import matplotlib.pyplot as plt

# define knapsack problem:
class Knapsack:
    def __init__(self,values,weights,max_weight):
        self.values = values
        self.weights = weights
        self.max_weight = max_weight
values = [23,21,8,1,3,7,18,19,17,15,24,22,6,28,4,2,27,20,5,10] # [$]
weights = [7,2,6,9,1,5,6,1,3,4,7,9,3,7,3,4,5,1,5,4] # [kg]
max_weight = 45 # [kg]
knapsack = Knapsack(values, weights, max_weight)

# Variable Inputs ('Knobs')
Gens = 100
populationSize = 1000
selection_ratio = .4
mutation_percent = .4

mutation_number = int(populationSize*mutation_percent)

# Step 1 - Initialize Population
myPop = GA.Population(populationSize, numGenes = len(values))

# record best and mean values of valid solutions for initial population:
fit = [GA.FitnessFunction(mem,knapsack) for mem in myPop.members]
values = []
for i in range(len(fit)):
    if fit[i] != 0:
        values.append(fit[i])
mean = [stat.mean(values)]
best = [max(values)]
generations = [1]

for i in range(Gens-1):

    # Step 2 - Selection
    for membs in myPop.members:
        membs.fitness = GA.FitnessFunction(membs, knapsack)
    myPop.selection(selection_ratio)

    # Step 3 - Crossover and Mutate
    # Crossover:
    num_parents = len(myPop.members)-1
    while len(myPop.members) < populationSize:
        p1 = random.randint(0, num_parents) # select parent 1
        p2 = random.randint(0, num_parents) # select parent 2
        if p1 != p2: # checks that we are not using the same parent
            myPop.members.append(myPop.members[p1] + myPop.members[p2])

    # Safegaurd Best Candidate Solution
    fit = [c.fitness for c in myPop.members] # calculate fitness for new population
    max_fitness = max(fit)
    max_index = fit.index(max_fitness)

    # Mutation:
    for j in range(mutation_number):
        mutator = random.randint(0, populationSize-1) # index for member to mutate

        # Check that mutator index doesn't corrospond to best value and re-generate index if necessary:
        if mutator == max_index:
            while mutator == max_index:
                mutator = random.randint(0, populationSize-1)
        
        # Mutate selected candidate:
        myPop.members[mutator].mutate()

    # save mean and best values of valid solutions:
    fit = [GA.FitnessFunction(mem, knapsack) for mem in myPop.members]
    values = []
    for k in range(len(fit)):
        if fit[k] != 0:
            values.append(fit[k])
    mean.append(stat.mean(values))
    best.append(max(values))
    generations.append(i+2)

    # Step 4 - repeat


# Plot Results:
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(generations, mean, label='Mean Value')
ax.plot(generations, best, label='Best Value')
plt.xlabel('Generation')
plt.ylabel('Value [$]')
plt.grid()
plt.show()
