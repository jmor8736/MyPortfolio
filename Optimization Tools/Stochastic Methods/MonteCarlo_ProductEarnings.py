# Jackson Morgan
# Monte Carlo Product Earnings Sensitivity Assessment

import statistics as stat
import numpy as np
import random 
import matplotlib.pyplot as plt
import math

n = 10000

# Given Triangular Distributions of inputs:
unit_price = [random.triangular(50, 70, 55) for i in range(n)]
unit_sales = [random.triangular(2000, 3000, 2440) for i in range(n)]
variable_costs = [random.triangular(50000, 65000, 55200) for i in range(n)]
fixed_costs = [random.triangular(10000, 20000, 14000) for i in range(n)]

# Calculate Earnings for Each of 10,000 combinations of inputs:
earnings = []
usup = []
for i in range(n):
    earnings.append((unit_price[i])*(unit_sales[i]) - fixed_costs[i] - variable_costs[i])
    usup.append(unit_sales[i]*unit_price[i])


# Calculate 95% confidence interval:
alpha = .05
sigma = stat.stdev(earnings)
theta_n = stat.mean(earnings)
z = 1.96
L = theta_n - z*sigma/math.sqrt(n)
U = theta_n + z*sigma/math.sqrt(n)
print('95% Confidence Interval = [' +str(L)+', '+str(U)+']')


# Create Histogram of Data and Label:
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.hist(earnings, color='black', ec='gray', bins=65)
ax.set_title('Earnings Distribution')
ax.set_xlabel('Earnings [$]')
ax.set_ylabel('Frequency')
plt.show()


# Input Sensitivities Plots:
fig1 = plt.figure()
ax1 = fig1.add_subplot(3,1,1)
ax2 = fig1.add_subplot(3,1,2)
ax3 = fig1.add_subplot(3,1,3)
ax1.scatter(variable_costs,earnings,marker='.')
ax1.set_ylabel('Earnings')
ax1.set_xlabel('Variabel Costs')
ax2.scatter(fixed_costs,earnings,marker='.')
ax2.set_ylabel('Earnings')
ax2.set_xlabel('Fixed Costs')
ax3.scatter(usup,earnings,marker='.')
ax3.set_ylabel('Earnings')
ax3.set_xlabel('Unit Sales x Unit Price')
plt.tight_layout()
plt.show()


# Sigma Normalized Derivitives:
vc_std = stat.stdev(variable_costs) # variable costs standard deviation
fc_std = stat.stdev(fixed_costs) # fixed costs standard deviation
usup_std = stat.stdev(usup) # unit price * units sold standard deviation
vc_snd = -vc_std/sigma
fc_snd = -fc_std/sigma
usup_snd = usup_std/sigma
print('Variable Costs Sigma-normalized Derivative = ' +str(vc_snd))
print('Fixed Costs Sigma-normalized Derivative = ' +str(fc_snd))
print('Unit Sales x Unit Price Sigma-normalized Derivative = ' +str(usup_snd))













