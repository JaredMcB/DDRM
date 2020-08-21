# -*- coding: utf-8 -*-
"""
Created on Wed May 15 16:06:13 2019

@author: JaredMcBride
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 15:00:42 2019

@author: JaredMcBride
"""

import numpy as np
import scipy.optimize as so


def modReduction_1(X,Psi,A_init):
    N = len(X)
    
    def aux_fun(A, X, Psi, N):
        a = 0
        b = np.array([[0, 0],[A,1]])
        
        y = np.zeros(N-1)
        y[0] = X[1]
        for i in range(1,N-1):
            y[i] = -a*y[i-1] + b[1] @ Psi(X[i]) + b[0] @ Psi(X[i-1])
        return X[1:] - y

    obj_fun = lambda A : aux_fun(A, X, Psi, N)
    
    A_sol = so.least_squares(obj_fun,A_init).x
        
    # Then we run it
    A_sol = A_init
        
    # Then we run it
    a = 0
    b = np.array([[0, 0],[A,1]])

    Y = np.zeros(N)
    y = np.zeros(N-1)
    
    Y[:1] = X[:1]
    y[0] = X[1]
    for i in range(1,N-1):
        y[i] = -a*y[i-1] + b[1] @ Psi(Y[i]) + b[0] @ Psi(Y[i-1])
        Y[i+1] = y[i] + np.sqrt(abs(b[1][0]))*np.random.normal(0,1)
        
    return [Y,A_sol]





 


    
    