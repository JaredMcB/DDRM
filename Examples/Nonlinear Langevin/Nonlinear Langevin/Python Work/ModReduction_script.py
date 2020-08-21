# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from ModReduction import modReduction
import scipy.optimize as so

start = 0
stop = 10**3
steps = 10**5 + 1

x_init = .1
h = (stop - start)/(steps - 1)
t = np.linspace(start,stop,steps)

X = np.zeros(steps)
X[0] = x_init

dVdx = lambda x : -x*(x**2 - 1)

for n in range(steps - 1):
    X[n+1] = X[n] + h*dVdx(X[n]) + np.sqrt(h)*np.random.normal(0,1)

Z = X[100:]
t = t[100:]
plt.plot(t[:1000],Z[:1000])

Psi = lambda x : np.array([dVdx(x), x])
A_init = [0,0,0,.01,1]

[Y,A_sol] = modReduction(Z,Psi,A_init)
plt.plot(t[:1000],Y[:1000])
print(A_sol)