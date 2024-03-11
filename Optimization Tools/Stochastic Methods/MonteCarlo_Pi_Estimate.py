# Jackson Morgan
# Monte Carlo Pi Estimate

import random
import math
import numpy as np
import matplotlib.pyplot as plt

## TASK 1)
numin = 0

x = [random.random()]
y = [random.random()]
if y[0] <= math.sqrt(1**2 - x[0]**2):
    numin = numin + 1

pi = [numin/1]
Pi = [math.pi]


while abs(pi[-1]-math.pi)/math.pi > .01:
    x.append(random.random())
    y.append(random.random())

    tot = len(x) + 1
    if y[-1] <= math.sqrt(1**2 - x[-1]**2):
        numin = numin + 1

    pi.append(4*numin/tot)
    Pi.append(math.pi)
    print(pi[-1])

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(pi)
ax.plot(Pi,'k--')
ax.set_title('Monte Carlo Pi Estimate')
ax.set_xlabel('Trial Number')
ax.set_ylabel('Estimate Value')
plt.show()
print(abs(pi[-1]-math.pi)/math.pi)


