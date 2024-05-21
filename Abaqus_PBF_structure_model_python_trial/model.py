# Microstructure evolution model that calculates phase fractions for a given thermal history:
import numpy as np
import matplotlib.pyplot as plt
import math as m
import sys

# read temperature history:
# data = np.genfromtxt('../input/Abaqus_one_elem_tests/amp_temp_amprint.txt', delimiter=',')
data = np.genfromtxt('../input/Abaqus_one_elem_tests/one_elem.csv', delimiter=',')
t = data[:,0] # time
T = data[:,1] # temperature

#Discretize the given temperature history with smaller time steps using interpolation:
time_step = 0.001
t = np.arange(t[0], t[-1], time_step)
T = np.interp(t, data[:,0], data[:,1])

# initialize arrays:
beta_f = np.zeros(len(t))
alpha_f = np.zeros(len(t))
mart_f = np.zeros(len(t))
RLS = np.zeros(len(t))

#Initial Phase fractions:
beta_f[0] = 0.0
alpha_f[0] = 0.0
mart_f[0] = 0.0
RLS[0] = -1

### Melting temperature:
T_melt = 1610

#### Beta to Alpha transformation data:
# JMAK parameters:
data = np.genfromtxt('./JMAK_params_BtoA.csv', delimiter=',', skip_header=1)
k_BtoA = data[:,0]
n_BtoA = data[:,1]
temp_jmak_BtoA = data[:,2]

# Alpha equilibrium phase fraction
data = np.genfromtxt('./Equilibrium_alpha_fraction.csv', delimiter=',', skip_header=1)
T_alpha_eq = data[:,1]
alpha_eq = data[:,0]

# KM parameters:
b_km = 0.005
T_mart = 851

# Mart to alpha + beta transformation:
T_mart_dis = 400
data = np.genfromtxt('./JMAK_params_mart_to_AandB.csv', delimiter=',', skip_header=1)
k_mart_to_AandB = data[:,0]
n_mart_to_AandB = data[:,1]
temp_jmak_mart_to_AandB = data[:,2]

for i in range(1, len(t)): #iterate over all time steps

    dt = t[i]-t[i-1] # time increment
    dT = T[i]-T[i-1] # temperature increment

    # print('Time:', t[i])
    # print('Temperature:', T[i])


  #First asses RLS state:  
    # if it is powder:
    if RLS[i-1] == -1:
        if T[i] <= T_melt: #if next temperature is below melting temperature:
            RLS[i] = -1 # keep powder state
            continue # still powder, nothing happens, go to next iteration
        else: # if temperature is above melting temperature:
            RLS[i] = 0 # change to liquid state
            continue # no transformation happens, go to next iteration
    
    # if it is liquid:
    if RLS[i-1] == 0:
        if T[i] >= T_melt: #if above T_melt:
            RLS[i] = 0 # keep liquid state
            continue # still liquid, nothing happens, go to next iteration
        else: # if below T_melt:
            RLS[i] = 1 # change to solid state
            beta_f[i] = 1.0 # set beta phase fraction to 1
            alpha_f[i] = 0.0 # set alpha phase fraction to 0
            mart_f[i] = 0.0 # set martensite phase fraction to 0
            continue
    
    if RLS[i-1] == 1: # if it is solid:
        if T[i] >= T_melt:
            RLS[i] = 0 # change to liquid state
            beta_f[i] = 0.0 # set phase fractions to zero
            alpha_f[i] = 0.0
            mart_f[i] = 0.0
            continue # no transformation happens, go to next iteration

        else:
            RLS[i] = 1 # keep solid state

            # if it is solid, we calculate phase transformations depending on heating:
            if dT/dt >= -20.0 and dT/dt <= 0: # JMAK model alpha to beta:
                # calculate alpha equilibrium phase fraction:
                # alpha_eq_i = 0 # alpha equilibrium phase fractio
                # P = [-31188.514, 170526.26, -388991.69, 471927.45, -315178.49, 99079.891, 1667.1991, -9726.8403, 1884.7280]
                # for n in range(0, 9):
                #     alpha_eq_i += P[n]*(T[i]/1000)**(8-n)

                alpha_eq_i = np.interp(T[i], T_alpha_eq, alpha_eq)

                # linearly interpolate JMAK parameters k_i and n_i according to current temperature:
                k_i = np.interp(T[i], temp_jmak_BtoA, k_BtoA)
                n_i = np.interp(T[i], temp_jmak_BtoA, n_BtoA)

                #Calculate Tau:
                tau = ( -1/k_i*m.log(1 - (alpha_f[i-1]/alpha_eq_i)/(beta_f[i-1]+alpha_f[i-1]) ))**(1/n_i)

                # print('tau:',tau)

                # caculate change in alpha phase fraction:
                alpha_f[i] = (1 - m.exp(-k_i * (tau+dt)**n_i)) * alpha_eq_i * (beta_f[i-1]+alpha_f[i-1])

                # calculate change in beta phase fraction:
                beta_f[i] = 1 - alpha_f[i]
                mart_f[i] = mart_f[i-1]

            if dT/dt <= -410.0 and mart_f[i-1] <= 1 and T[i] < T_mart: # KM model beta to martensite fast cooling:
                
                mart_f[i] = ( 1 - m.exp( -b_km*(T_mart - T[i]) ) ) * (beta_f[i-1]+mart_f[i-1])
                beta_f[i] = beta_f[i-1] - (mart_f[i] - mart_f[i-1])
                alpha_f[i] = alpha_f[i-1]

            elif (dT/dt > -410.0) and dT/dt <= -20.0 and mart_f[i-1] <= 1 and T[i] < T_mart: # KM model beta to martensite slow cooling:
                
                beta_eq_i = 1 - np.interp(T[i], T_alpha_eq, alpha_eq)

                mart_f[i] = ( 1 - m.exp( -b_km*(T_mart - T[i]) ) ) * (beta_f[i-1]+mart_f[i-1]-beta_eq_i)
                beta_f[i] = beta_f[i-1] - (mart_f[i] - mart_f[i-1])
                alpha_f[i] = alpha_f[i-1]

            # elif (dT/dt >= -0.0001 and T[i] >= T_mart_dis): # JMAK model of mart to alpha + beta transformation:
                
            #     # linearly interpolate JMAK parameters k_i and n_i according to current temperature:
            #     k_mart_i = np.interp(T[i], temp_jmak_mart_to_AandB, k_mart_to_AandB)
            #     n_mart_i = np.interp(T[i], temp_jmak_mart_to_AandB, n_mart_to_AandB)

            #     #Calculate Tau:
            #     tau_mart_i = ( -1/k_mart_i * m.log( (mart_f[i-1] - f_mart_eq) / (1 - f_mart_eq) )) ** (1/n_mart_i)

            #     f_mart_eq = 0.5 * (1 - m.tanh( (450 - T[i]) / 80 ))

            #     # mart_f[i] = 1 - (1 - m.exp(-k_mart_i*(tau_mart_i+dt)**n_mart_i)) * (1 - f_mart_eq)
            #     mart_f[i] = f_mart_eq - ( m.exp( -k_mart_i ( tau_mart_i +dt  )**n_mart_i) ) * (mart_f[i-1] + beta_f[i-1] - f_mart_eq)

            #     alpha_f[i] = alpha_f[i-1] + (mart_f[i-1] - mart_f[i])*f_mart_eq
            #     beta_f[i] = beta_f[i-1] + (mart_f[i-1] - mart_f[i])*(1 - f_mart_eq)

            else: # if no phase transfpormation happens, they stay the same

                alpha_f[i] = alpha_f[i-1]
                beta_f[i] = beta_f[i-1]
                mart_f[i] = mart_f[i-1]

        # if dT/dt <= -410.0:
        #     # KM model beta to martensite:
        #     mart_f[i] = (1 - m.exp(-b_km*(T_mart-T[i])))*(beta_f[i-1]+mart_f[i-1])
        #     beta_f[i] = beta_f[i-1] - mart_f[i]

        # if dT/dt < -20.0 and dT/dt >= -410.0:
        #     # KM model beta to martensite:
        #     mart_f[i] = (1 - m.exp(-b_km*(T_mart-T[i])))*(beta_f[i-1]+mart_f[i-1])
        #     beta_f[i] = beta_f[i-1] - mart_f[i]

# End of time-series iteration

# Plots:

# Print temperature:
plt.figure()
plt.plot(t, T)
plt.xlabel('Time (s)')
plt.ylabel('Temperature (C)')
plt.title('Temperature vs Time')
plt.savefig('temp_vs_time.png')


#plot RSL vs temperature:
plt.figure()
plt.plot(t, RLS)
plt.xlabel('Time (s)')
plt.ylabel('RLS')
plt.title('RLS vs Time')
plt.savefig('RLS_vs_time.png')

#plot graph with two lines, with one of them having their values represented on the right axis and the other on the left axis:
plt.figure()
fig, ax1 = plt.subplots()
ax1.plot(t, T, 'b-')
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Temperature (°C)', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(t, alpha_f, 'r-')
ax2.plot(t, beta_f, 'r--')
ax2.plot(t, mart_f, 'r:')
ax2.set_ylabel('Phase Fraction', color='r')
ax2.tick_params('y', colors='r')

ax1.legend(['Temperature'], loc='upper left')
ax2.legend(['Alpha', 'Beta', 'Mart'], loc='upper right')

fig.tight_layout()
plt.savefig('phase_fractions.png')


# # plot this function: alpha_eq_i += P[n]*(T[i]/1000)**(8-n)
# plt.figure()
# T_plot = np.linspace(600, 1100, 100)
# alpha_eq = np.zeros(len(T_plot))

# for i in range(0, len(T_plot)):
#     for n in range(0, 9):
#         alpha_eq[i] += P[n]*(T_plot[i]/1000)**(8-n)

# plt.plot(T_plot, alpha_eq)
# #save plot to file:
# plt.savefig('alpha_eq.png')

# 

# Print final phase fractions:
print('Final phase fractions:')
print('Alpha:', alpha_f[-1])
print('Beta:', beta_f[-1])
print('Martensite:', mart_f[-1])



