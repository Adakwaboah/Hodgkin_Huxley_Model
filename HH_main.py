'''
Written by: AKwasi Darkwah Akwaboah
Description: Forward Euler Implementation of the Hodgkin Huxley Model
Date: October 6, 2019
'''

import numpy as np
import matplotlib.pyplot as plt

Vrest = -60
ENa = 54.2 #mV
EK = -74.7 #mV
EL = -43.256 #mV

Cmem = 1 #uF/cm2
gK_bar = 36 #mS/cm2
gNa_bar = 120
gL_bar = 0.3

dt = 0.001
tStart = -1.000
tEnd = 15.000
nStep = np.int(np.ceil((tEnd-tStart)/dt)) #number of steps
print(nStep)
StimDur = 2
Is = 150 #stimulus strength

Vm = Vrest;

m = 0.05293
h = 0.59612; #initial value for h gate
n = 0.31768; #initial value of n gate

#create storage for state variables
plot_Vm = np.zeros((nStep, 1), dtype=np.float)
plot_m = np.zeros((nStep, 1), dtype=np.float)
plot_h = np.zeros((nStep, 1), dtype=np.float)
plot_n = np.zeros((nStep, 1), dtype=np.float)
plot_time = np.empty((nStep, 1), dtype=np.float)

print('The Hodgkin Huxley Model')
tNow = tStart
for iStep in range(nStep):
    JNa = gNa_bar*m*m*m*h*(Vm - ENa)
    JK = gK_bar*n*n*n*n*(Vm - EK)
    JL = gL_bar*(Vm - EL)

    Iss = (0.0 < tNow < 2.0) * Is + (tNow >= 2.0 or tNow < 0) * 0
    alpha_m = (0.1 * (25 - Vm)) / (np.exp((25 - Vm) / 10) - 1)
    beta_m = 4 * np.exp(-Vm / 18)
    m = m + (dt*((alpha_m*(1-m)) - (beta_m*m)))

    alpha_n = (0.01 * (10 - Vm)) / (np.exp((10 - Vm) / 10) - 1)
    beta_n = 0.125 * np.exp((- Vm) / 80)
    n = n + (dt * ((alpha_n * (1 - n)) - (beta_n * n)))

    alpha_h = 0.07 * np.exp((- Vm) / 20)
    beta_h = 1 / (1 + np.exp((30 - Vm) / 10))
    h = h + dt * ((alpha_h * (1 - h)) - (beta_h * h))
    dVm = -dt * (JNa + JK + JL - Iss) / Cmem
    plot_Vm[iStep] = Vm
    plot_m[iStep] = m
    plot_n[iStep] = n
    plot_h[iStep] = h
    plot_time[iStep] = tNow

    Vm = Vm + dVm
    tNow = tStart + iStep * dt

plt.plot(plot_time, plot_Vm)
plt.title('time course of Transmembrane Potential')
plt.xlabel('time(ms)')
plt.ylabel('Vm (mV)')
plt.grid()
plt.show()

plt.figure()
plt.plot(plot_time, plot_m)
plt.plot(plot_time, plot_h)
plt.plot(plot_time, plot_n)
plt.title('time course of Gating probabilties, (m, h, n)')
plt.xlabel('time (ms)')
plt.ylabel('gating (activation/ inactivation) probability')
plt.legend(['m', 'h', 'n'], loc='best')
plt.grid()
plt.show()


