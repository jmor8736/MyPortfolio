# Jackson Morgan
# Monte Carlo Definite Integral Estimate

import random
import math
import matplotlib.pyplot as plt
import numpy as np

n = 1000 # number of trials
a = -2
b = 3

# estimate count with mean value theorem
x1 = np.arange(a,b,(b-a)/n)
y1 = []
for i in range(len(x1)):
    y1.append((x1[i]**3)*math.sin(x1[i]))
En = ((b-a)/n)*sum(y1) # abosolute integral estimate

x = [random.random()*5-2 for i in range(n)] # random number between -2 and 3
y = [random.random()*10 for i in range(n)] # random number between 0 and 10

# Solve using Monte Carlo Method:
nUnder = 0 # initialize count for number of points under curve
integral = []
tot = 0
total = []
est = []
for i in range(len(x)):
    if y[i] <= x[i]**3*math.sin(x[i]):
        nUnder = nUnder + 1
    tot = tot + 1
    total.append(tot)
    integral.append((b-a)*10*nUnder/tot)
    est.append(En)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(total,integral)
ax.plot(total,est,'k--')
ax.set_title('Monte Carlo Integration')
ax.set_xlabel('Trial Number')
ax.set_ylabel('Integral Estimate')
plt.show()

print('Integral = ' + str(integral[-1]))
