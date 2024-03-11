# Jackson Morgan
# Biomimetic Task 2 - Particle Swarm Optimization

import PSO
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Initialize:
S = 20
particles = []
domain = [-32.768, 32.768] # for [x, y, z]
D = 3 # number of dimensions

# Adjustable Constants:
w = 0.8
phi_p = 0.1
phi_g = 0.1

# Create all the particles:
for s in range(S):
    x = [random.uniform(domain[0], domain[1]) for i in range(D)]
    v = [random.uniform(-5, 5) for i in range(D)]
    particles.append(PSO.Particle(x, v, w, phi_p, phi_g))

# let the for loop know what type we're dealing with:
particle = PSO.Particle

g = particles[random.randint(0,S)].x # initialize g as position of a random particle
for particle in particles:
    particle.updateValue()
    if particle.value < PSO.f(g):
        g = particle.x[:]
        print(PSO.f(g))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Swarm:
n = 0

G = [] # initialize array to keep track of g's at each step
f_G = [] # initialize array to keep track of f(g) at each step
while n < 50:
    for particle in particles:
        # let the object do the work:
        particle.updateVelocity(g)
        particle.updatePosition()
        particle.updateValue()
        if particle.value < PSO.f(g):
            g = particle.x[:]
            print(PSO.f(g))
            print(g)
        ax.scatter(particle.x[0], particle.x[1], particle.x[2], 'b*')
    f_G.append(PSO.f(g)) # add current value of g to array
    G.append(g) # add current g coords to array
    n += 1



plt.show()

fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
ax1.plot(f_G)
plt.xlabel('Trial Number')
plt.ylabel('f(g)')
plt.grid()
plt.show()











