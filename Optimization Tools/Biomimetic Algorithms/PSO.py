# Jackson Morgan
# Particle Swarm Optimization - Classes and Functions

import random
import numpy as np


# create particle class to store state of each particle and update state, allowing each particle to maintain it's own agency:
class Particle:

    def __init__(self, x, v, w, phi_p, phi_g):
        # x - position 
        # v - velocity 
        # w - inertia
        # phi_p - cognitive coefficient
        # phi_g - social coefficient
        self.x = x
        self.v = v
        self.w = w
        self.phi_p = phi_p
        self.phi_g = phi_g
        self.p = None # start best known value as zero

    def updatePosition(self):
        # implement x = x + v on self
        for i in range(len(self.x)):
            self.x[i] = self.x[i] + self.v[i]

    def updateVelocity(self, g):
        # random r_p and r_g
        # update self.v
        r_p = random.random()
        r_g = random.random()
        for i in range(len(self.x)):
            self.v[i] = self.w*self.v[i] + self.phi_p*r_p*(self.p[i] - self.x[i]) + self.phi_g*r_g*(g[i] - self.x[i])

    def updateValue(self):
        # calculate self.value, bsaed on f(x)
        self.value = f(self.x)

        # update best known value (if appliccable)
        if self.p is None or self.value < f(self.p):
            self.p = self.x[:]
        


# Evaluate particle location with the ackley function:
def f(x):
    # INPUTS:
    # x - vector with coordinates for function to be evaluated at (length = number of dimensions)

    # Constants:
    a = 20
    b = 0.2
    c = 2*np.pi

    # get constant dependant on number of dimensions:
    d = len(x)

    # initialize summation terms
    sum1 = 0
    sum2 = 0
    
    # calculate summation terms, looping through number of dimensions:
    for i in range(d):
        sum1 = sum1 + x[i]**2
        sum2 = sum2 + np.cos(c*x[i])

    # calculate mathematical terms in function:
    term1 = -a * np.exp(-b*np.sqrt(sum1/d))
    term2 = -np.exp(sum2/d)

    # calculate total function and return as output:
    y = term1 + term2 + a + np.exp(1)
    return(y)
