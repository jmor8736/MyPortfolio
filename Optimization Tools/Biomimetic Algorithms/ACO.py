# Jackson Morgan
# Ant Colony Optimization Functions

import numpy as np
import random


# function to build visibility matrix:
def visibility(d, nodes, Mk):    
    vis = np.zeros(np.shape(d)) # preallocating visibility matrix
    for i in range(len(nodes)): # iterate through rows
        for j in range(len(nodes)): # iterate through columns
            if i != j:
                if nodes[j] in Mk:
                    vis[i,j] = 1 / d[i,j] # change entry from 0 only if the column is an unvisited city
    return(vis)


# function to build cumulative probability array:
def cum(nodes, tau, eta, path, Mk):
    alpha = 1
    beta = 2
    p = [] # initialize probability matrix
    i = path[-1]
    for j in range(len(nodes)):
        if nodes[j] in Mk:
            p.append(tau[i,j]**alpha * eta[i,j]**beta) # calculate probability of going to unvisited cities

    sum_probs = sum(p)
    cum_p = []
    last = 0
    for i in range(len(p)):
        cum_p.append(p[i]/sum_probs + last)
        last = cum_p[-1]
    return(cum_p)


# function to generate ant path:
def path(nodes, d, tau):
    path = [random.choice(nodes)] # randomly get starting point for ant

    Mk = nodes.copy()
    Mk.remove(path[0]) # remove starting node from Mk

    # generate next node in path:
    while len(Mk) > 0:
        eta = visibility(d, nodes, Mk) # get visibility matrix at current state
        sum_probs = cum(nodes, tau, eta, path, Mk)
        rand = random.random() # random number between 1 and 0 for next node selection
        for j in range(len(sum_probs)):
            if rand < sum_probs[j]:
                break
        path.append(Mk[j])
        if len(Mk) != 0:
            Mk.remove(Mk[j])
    path.append(path[0])
    return(path)


# pheromone evaporation function:
def evaporate(rho, tau):
    for i in range(len(tau[:,0])): # iterate through rows
        for j in range(len(tau[0,:])): # iterate through columns
            tau[i,j] = (1-rho)*tau[i,j] # evaporate pheromones
    return(tau)

# pheremone update function:
def pheromone(nodes, path, d):
    # Initialize pheromone contribution matrix
    delta_tau = np.zeros(np.shape(d))

    val = []
    for k in range(len(nodes)):
        i = nodes.index(path[0,k]) # get i index for distance matrix
        j = nodes.index(path[0,k+1]) # get j index for distance matrix
        val.append(d[i,j]) # add up distance values
    tot = 1/sum(val) # get ant's total phermone contribution to path

    for k in range(len(nodes)):
        i = nodes.index(path[0,k]) # get i index for delta_tau matrix
        j = nodes.index(path[0,k+1]) # j index
        delta_tau[i,j] = tot # add ant's contribution

    return(delta_tau)


# path value function:
def value(nodes, path, d):
    val = []
    for k in range(len(nodes)):
        i = nodes.index(path[0,k])
        j = nodes.index(path[0,k+1])
        val.append(d[i,j])
    return val





