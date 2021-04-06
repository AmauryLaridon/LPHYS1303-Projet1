# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi ,exp, cos
from scipy.stats import hypsecant

delta = 0.02




# Gaussiennes
def gauss1(x):
    return (1/sqrt(2*np.pi*0.1)) * exp(-50*(x-2)**2)
    
def gauss2(x):
    return (1/sqrt(2*np.pi*0.2)) * exp(-25*(x-4)**2)
    
def gauss(x):
    return gauss1(x) + gauss2(x)






# Solitons
def soliton1(x,c,a):
    return 0.5*c*(pi*hypsecant.pdf((x-a)*sqrt(c)/2))**2

def soliton2(x,c,a):
    return 0.5*c*(pi*hypsecant.pdf((x-a)*sqrt(c)/2))**2

def solit(c1,a1,c2,a2):
    def sol_(x):
        return soliton1(x,c1,a1) + soliton2(x,c2,a2)
    return sol_




# Cosinus
def cos2(x):
    return 0.5*cos(0.5*x-0.5*pi)**2
    





# Résolution numérique
def time_ev(u_0, t_f, k, x_f, h, snaps = []):
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
            
        U.append([0,0] + [u2[i] - (k/(3*h)) *(u1[i+1] + u1[i] + u1[i-1]) * (u1[i+1] - u1[i-1]) - (delta**2) * (k/(h**3)) * (u1[i+2] - 2*u1[i+1] + 2*u1[i-1] - u1[i-2])  for i in range(2,len(u1)-2)] + [0,0])
        
        
    for s in snaps:
        plt.plot(x_range, U[int(s/k)], label = "t = {}".format(s))
    plt.legend()
    plt.show()
    plt.clf()
    
    [xx,tt]=np.meshgrid(x_range,t_range)
    plt.contourf(xx,tt, np.array(U), levels = 10)
    plt.xlabel("x")
    plt.ylabel("t")
    plt.colorbar()
    plt.show()
    
    
    
    
if __name__ == "__main__":
    #time_ev(gauss, 12, 0.0005, 10, 0.05, [0,4,8,11.5])
    #time_ev(solit(8,15,16,5), 13, 0.001, 45, 0.5, [0,4,8,12])
    time_ev(cos2, 8, 0.0001, 2*pi, 0.05, [0,1,4,5.5])