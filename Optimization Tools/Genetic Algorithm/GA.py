# Jackson Morgan
# Genetic Algorithm Classes and Functions

import numpy as np
import random

class Chromosome:

    # Initialize Chromosome:
    def __init__(self, numGenes):
        self.genes = random.choices([0, 1], k = numGenes)
        self.fitness = 0

    # Mutation Function:
    def mutate(self):
        gene_loc = random.randint(0, len(self.genes)-1) # select random gene location
        self.genes[gene_loc] = 1 - self.genes[gene_loc] # flip value of that location
    
    # Breeding Method (syntax: child = Chromosome1 + Chromosome2)
    def __add__(self, o):
        crossover_point = random.randint(0, len(self.genes)-1) # random crossover point
        child_genes = self.genes[:crossover_point] + o.genes[crossover_point:] # mix parent genes to create child genes
        child = Chromosome(len(child_genes)) # initialize child
        child.genes = child_genes # insert mixed parent genes into child genes
        return child


class Population:

    def __init__(self, populationSize, numGenes):
        self.members = [Chromosome(numGenes) for i in range(populationSize)]

    def selection(self, ratio):
        # Step 1 - Sort members by fitness (ascending order)
        self.members.sort(key=lambda x: x.fitness) # NOTE - this requires you to determine member fitness in script

        # Step 2 - return some number of members based on the ratio provided
        num_mems = int(ratio*int(len(self.members))) # calculates the number (x) of candidates we will keep
        self.members = self.members[num_mems:] # reduces the population to the best x candidates

def FitnessFunction(chrom: Chromosome, knapsack):

    # calculate total weight and total value of candidate solution:
    tot_value = sum(chrom.genes[i]*knapsack.values[i] for i in range(len(chrom.genes)))
    tot_weight = sum(chrom.genes[i]*knapsack.weights[i] for i in range(len(chrom.genes)))

    # penalize fitness for candidate solutions that exceed the max weight:
    if tot_weight > knapsack.max_weight:
        tot_value = 0
    
    # return total value as the fitness level:
    return(tot_value)
    
