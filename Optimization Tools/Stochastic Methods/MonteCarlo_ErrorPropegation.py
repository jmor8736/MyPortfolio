# Jackson Morgan
# Monte Carlo Error Propegation
import matplotlib.pyplot as plt
import math
import random
import statistics as stat
import numpy as np

# Ideal Gas Law: p = rho*R*T

# State 1:
R_air = 287.053 # [J/kgK]
T0 = 15 # [C]
T0 = 15 + 273.15 # [K]
rho0 = 1.225 # [kg/m^3]

# State 2:
T2 = 25 + 273.15 # [K]
T_std = .2 # standard deviation [K]
P2 = 104847 # [Pa]
P_std = 52 # standard deviation [Pa]

n = 100000

# Analytic Error Calculation
delta_rho = math.sqrt((P_std/P2)**2 + (T_std/T2)**2)
print(delta_rho)

# Monte Carlo Error Estimate
T_range = [random.random()*(2*T_std)+(T2-T_std) for i in range(n)]
P_range = [random.random()*(2*P_std)+(P2-P_std) for i in range(n)]
rho = []
for i in range(n):
    rho.append(P_range[i]/(R_air*T_range[i]))
delta_rho_MC = stat.stdev(rho)

print('Analytic Uncertainty = ' + str(delta_rho))
print('Monte Carlo Uncertainty = ' + str(delta_rho_MC))


# Varying Monte Carlo Uncertainty
P_std_range = np.arange(0,850,50)
print(P_std_range)

delta_rho_MCrange = []
for i in range(len(P_std_range)):
    P_rangeMC = [random.random()*(2*P_std_range[i])+(P2-P_std_range[i]) for j in range(n)]
    rhoMC = []
    for i in range(n):
        rhoMC.append(P_rangeMC[i]/(R_air*T_range[i]))
    delta_rho_MCrange.append(stat.stdev(rhoMC))

# Plot Results:
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(P_std_range,delta_rho_MCrange)
ax.set_title('Density Error vs Pressure Error')
ax.set_xlabel('Pressure Error [Pa]')
ax.set_ylabel('Density Error [kg/m^3]')
ax.grid()
plt.show()