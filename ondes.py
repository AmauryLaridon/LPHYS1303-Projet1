# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi ,exp


def onde1(x):
    return 10*(1/sqrt(2*np.pi*0.1)) * exp(-5*(x-1)**2)
    
def onde2(x):
    return 10*(1/sqrt(2*np.pi*0.4)) * exp(-1.25*(x-5)**2)
    
def u0(x):
    return onde1(x) + onde2(x)

delta = 0.02
    

def time_ev(u_0, t_f, k, x_f, h):
    x_range = [i*h for i in range(int(x_f/h))]
    U = []
    U.append([u_0(x) for x in x_range])
    
    t_range = [j*k for j in range(int(t_f/k))]
    for t in t_range[1:]:
        u1 = U[-1]
        
        if t == k:
            u2 = u1
        else:
            u2 = U[-2]
            
        U.append([0,0] + [u2[i] - (k/(3*h)) *(u1[i+1] + u1[i] + u1[i-1]) * (u1[i+1] - u1[i-1]) + (delta**2) * (k/(h**3)) * (u1[i+2] - 2*u1[i+1] + 2*u1[i-1] - u1[i-2])  for i in range(2,len(u1)-2)] + [0,0])
    
    [xx,tt]=np.meshgrid(x_range,t_range)
    plt.contourf(xx,tt, np.array(U))
    plt.xlabel("x")
    plt.ylabel("t")
    plt.colorbar()
    plt.show()
    
    
    
if __name__ == "__main__":
    time_ev(u0, 1, 0.0001, 10, 0.01)