# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi ,exp
from scipy.stats import hypsecant
from matplotlib.animation import FuncAnimation
from matplotlib.pyplot import cm


def onde1(x):
    return (1/sqrt(2*np.pi*0.01)) * exp(-50*(x-1)**2)

def onde2(x):
    return (1/sqrt(2*np.pi*0.4)) * exp(-1.25*(x-5)**2)

def u0(x):
    return onde1(x) + onde2(x)

delta = 0.022

def time_ev(u_0, t_f, k, x_f, h):
    x_0 = 0
    L = x_f- x_0
    N = L/h
    t_0 = 0
    T = t_f - t_0
    M = T/k
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

    print("Résolution numérique avec une grille spatiale de {} points".format(N))
    print("Résolution numérique avec une grille temporelle de {} points".format(M))
    print("Paramètres numérique : L = {}, T = {}s, h = {}, k = {}, delta = {}".format(L, T, h, k, delta))

    #Plot semi-animé d'Instantanés
    t_span = np.arange(0,8,0.4)
    n = np.zeros(np.shape(t_span)[0])
    for t in t_span:
        n = [int(t/k)]
        #for j in range(np.shape(t_span)[0]):
        if n[0] < len(U):
            plt.plot(x_range, U[n[0]], label = "$t={:2.2f}\; s$".format(t), marker ='.')
            plt.title('Instantanés dépassement de solitons.\n $L = {}, h = {}, k = {}, T = {}, \delta = {}$'.format(x_f, h, k, t_f, delta))
            plt.xlabel('$x$')
            plt.ylabel('Amplitude')
            plt.legend()
            plt.show(block=False)
            plt.pause(0.5)
            plt.close()
    #Plot d'Instantanés
    t_span = [0, 0.1, 3, 5]
    n = np.zeros(np.shape(t_span)[0])
    n = [int(t_span[i]/k) for i in range(np.shape(t_span)[0])]
    for i in range(np.shape(t_span)[0]):
        plt.plot(x_range, U[n[i]], label="$t = {:2.2f}s$".format(t_span[i]), marker=',')
    plt.title('Instantanés dépassement de solitons.\n $L = {}, h = {}, k = {}, T = {}, \delta = {}$'.format(x_f, h, k, t_f, delta))
    plt.xlabel('$x$')
    plt.ylabel('Amplitude')
    plt.legend()
    plt.show()

    """
    #Animation 2D
    x_data = x_range
    y_data = U
    U0 = [u_0(x) for x in x_range]
    #class matplotlib.animation.FuncAnimation(fig, func,
    #frames = None, init_func=None, fargs=None, save_count = None, *, cache_frame_data = True, **kwargs)


    x_0 = 0
    fig, ax = plt.subplots()
    ax.set_xlim(x_0,x_f)
    ax.set_ylim(-5,20)
    line, = ax.plot(U0)

    def animation_frame(i):
        x_data.append(i)
        y_data.append(i)
        line.set_xdata(x_data)
        line.set_ydata(y_data)
        return line,

    animation = FuncAnimation(fig, func=animation_frame, frames = np.arange(0, 10, 0.01), interval = 1)
    plt.show()
    """
    #Plot 2D
    [xx,tt]=np.meshgrid(x_range,t_range)
    plt.contourf(xx,tt, np.array(U), levels = 10)
    plt.title('Dépassement de solitons.\n $L = {}, h = {}, k = {}, T = {}, \delta = {}$'.format(x_f, h, k, t_f, delta))
    plt.xlabel("x")
    plt.ylabel("t")
    plt.colorbar()
    plt.show()
    plt.clf()

def soliton1(x,c,a):
    return 0.5*c*(pi*hypsecant.pdf((x-a)*sqrt(c)/2))**2

def soliton2(x,c,a):
    return 0.5*c*(pi*hypsecant.pdf((x-a)*sqrt(c)/2))**2

def solit(c1,a1,c2,a2):
    def sol_(x):
        return soliton1(x,c1,a1) + soliton2(x,c2,a2)
    return sol_


if __name__ == "__main__":
    time_ev(solit(20,15,40,5), 6, 0.0005, 40, 0.5)
