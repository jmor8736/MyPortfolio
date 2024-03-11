# Jackson Morgan
# MBFA Validation on Knapsack Problem

import MBFA
import matplotlib.pyplot as plt

# Define knapsack problem:
values = [23,21,8,1,3,7,18,19,17,15,24,22,6,28,4,2,27,20,5,10]
weights = [7,2,6,9,1,5,6,1,3,4,7,9,3,7,3,4,5,1,5,4] # [kg]
max_weight = 45 # kg
knapsack = MBFA.Knapsack(values, weights, max_weight)

# number of candidates to be created/tested:
n = 100000

# create dummy 1st, 2nd, and 3rd best candidates to test against and set total value to 0
mem1 = mem2 = mem3 = MBFA.candidate(knapsack)
mem1.tot_v = mem2.tot_v = mem3.tot_v = 0
all_vals = [0]
best_val = [0]
mean_val = [0]
trial = [0]

for j in range(n): # 

    mem = MBFA.candidate(knapsack) # generate candidate solution
    all_vals.append(mem.tot_v)

    # Check if current member is in the top 3 solutions and unique. Index as one of best if it qualifies:
    if mem.tot_v > mem3.tot_v and mem.tot_v != mem2.tot_v and mem.tot_v != mem1.tot_v:
        if mem.tot_v > mem2.tot_v and mem.tot_v != mem1.tot_v:
            if mem.tot_v > mem1.tot_v:
                mem1 = mem
            else:
                mem2 = mem
        else:
            mem3 = mem

    trial.append(trial[-1]+1)
    best_val.append(mem1.tot_v)
    mean_val.append(sum(all_vals)/trial[-1])

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(trial,best_val, label='Best Value')
ax.plot(trial,mean_val, label='Mean Value')
plt.xlabel('Trial Number')
plt.ylabel('Best Value [$]')
plt.grid()
plt.show()

# Third Best:
print(mem3.w)
print(mem3.v)
print('Third Weight: ' +str(mem3.tot_w))
print('Third Value: ' +str(mem3.tot_v))

# Second Best:
print(mem2.w)
print(mem2.v)
print('Second Weight: ' +str(mem2.tot_w))
print('Second Value: ' +str(mem2.tot_v))

# First Best:
print(mem1.w)
print(mem1.v)
print('First Weight: ' +str(mem1.tot_w))
print('First Value: ' +str(mem1.tot_v))








