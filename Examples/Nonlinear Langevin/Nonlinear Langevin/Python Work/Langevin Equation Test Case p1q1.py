# -*- coding: utf-8 -*-
"""
Created on Fri May 17 09:08:19 2019

@author: JaredMcBride
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sa
import scipy.optimize as so

start = 0
stop = 10**4
steps = stop*100 + 1

dVdx = lambda x : -x*(x**2 -1) 

x_init = .5

h = (stop - start)/(steps - 1)
t = np.linspace(start,stop,steps)

X = np.zeros(steps)
X[0] = x_init

for n in range(steps - 1):
    X[n+1] = X[n] + h*dVdx(X[n]) + np.sqrt(h)*np.random.normal(0,1)
    
def sim(X,Psi,A_init):
    N = len(X)

    A_sol = A_init
        
    # Then we run it
    a = A_sol[0]
    b = np.reshape(A_sol[1:],[2,2])

    Y = np.zeros(N)
    y = np.zeros(N-1)
    Y[:1] = X[:1]
    y[0] = X[1]
    for i in range(1,N-1):
        y[i] = -a*y[i-1] + b[1] @ Psi(Y[i]) + b[0] @ Psi(Y[i-1])
        Y[i+1] = y[i] + np.sqrt(abs(b[1][0]))*np.random.normal(0,1)
        
    return Y

Psi = lambda x : np.array([dVdx(x), x])
A_init = [0,0,0,h,1]
Y = sim(X,Psi,A_init)
plt.plot(t[:1000],X[:1000])
plt.plot(t[:1000],Y[:1000])

def modReduction(X,Psi,A_init):
    N = len(X)
    
    def aux_fun(A, X, Psi, N):
        a = A[0]
        b = np.reshape(A[1:],[2,2])
        
        y = np.zeros(N-1)
        y[0] = X[1]
        for i in range(1,N-1):
            y[i] = -a*y[i-1] + b[1] @ Psi(X[i]) + b[0] @ Psi(X[i-1])
        return X[1:] - y

    obj_fun = lambda A : aux_fun(A, X, Psi, N)
    
    A_sol = so.least_squares(obj_fun,A_init).x
        
    # Then we run it
    a = A_sol[0]
    b = np.reshape(A_sol[1:],[2,2])

    Y = np.zeros(N)
    y = np.zeros(N-1)
    
    Y[:1] = X[:1]
    y[0] = X[1]
    for i in range(1,N-1):
        y[i] = -a*y[i-1] + b[1] @ Psi(Y[i]) + b[0] @ Psi(Y[i-1])
        Y[i+1] = y[i] + np.sqrt(abs(b[1][0]))*np.random.normal(0,1)
        
    return [Y,A_sol]

n_disc = 60
Z = X[n_disc:]
A_init = [0,0,0,0,0]
[Y,A_sol] = modReduction(Z,Psi,A_init)
A_sol

plt.plot(t[:1000],Z[:1000])
plt.plot(t[:1000],Y[:1000])