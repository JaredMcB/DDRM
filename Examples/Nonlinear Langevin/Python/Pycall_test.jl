# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:12:05 2020

@author: JaredMcBride
"""
using PyCall
np = pyimport("numpy")
using Plots
pyplot()

# import acor as ac
# import numpy as np

t_start = 0
t_stop = 1000

dt = 2^(-8)
T = Int((t_stop - t_start)/dt + 1)
T_disc = 10^4

np.sin(pi/2)

Time = range(t_start,t_stop,step = dt)
x = zeros(T+T_disc)
u = randn(T+T_disc)

r = .6

x[1] = 1.5
for i = 1 : T+T_disc - 1
    x[i+1] = r*x[i] + u[i]
end
lags = 0:T
A = autocov(x,lags)

end_ind = findall(A.<0)[1] - 1
plot(log.(A[1:end_ind]))
plot(A[1:end_ind])

A_mat = [ones(end_ind,1) reshape(1:end_ind,end_ind,1)]
b = inv(A_mat'*A_mat)*A_mat'*log.(A[1:end_ind])

b[2]
τ_exp = -1/b[2]
plot!(1:end_ind,A_mat*b)
N_disc = 20*τ_exp

## Integrated autocorrelation time

plot(A)
end_ind_int = findall(A.<0)[1] - 1
plot(A[1:end_ind])


A ./= A[1]

τ_int = .5 + sum(A[2:end_ind_int])
N_eff = T/τ_int

plot!(1:15,map(x -> .6^x/.64,0:14))

T = 2000
T_disc = 20

x = zeros(T+T_disc)
u = randn(T+T_disc)

x[1] = 1.5
for i = 1 : T+T_disc - 1
    x[i+1] = r*x[i] + u[i]
end

x = x[T_disc+1:end]

A = autocor(x,0:1999)

plot(1:20,A[1:20])
plot!(1:15,map(x -> .6^x,0:14))

end_ind = findall(A.<0)[1] - 2
end_ind = 15
plot(log.(A[1:end_ind]))
plot(A[1:end_ind])

A_mat = [ones(end_ind,1) reshape(1:end_ind,end_ind,1)]
b = inv(A_mat'*A_mat)*A_mat'*log.(A[1:end_ind])

b[2]
τ_exp = -1/b[2]
plot!(1:end_ind,A_mat*b)
N_disc = 20*τ_exp

plot(A)
end_ind_int = findall(A.<0)[1] - 1
plot(A[1:end_ind])


A ./= A[1]

τ_int = .5 + sum(A[2:end_ind_int])
N_eff = T/τ_int
