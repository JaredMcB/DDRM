# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:12:05 2020

@author: JaredMcBride
"""

import acor as ac
import numpy as np

t_start = 0
t_stop = 1000

dt = 2**(-8)
T = (t_stop - t_start)/dt + 1
T_disc = 10**4

Time = np.linspace(t_start,t_stop,num = T, retstep = True)
x = np.zeros(T+T_disc)
u = np.rand(T+T_disc)

x[0] = 1.5
for i in range(T+T_disc - 1):
    x[i+1] = r*x[i] + u[i]
