# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi ,exp
from scipy.stats import hypsecant
from matplotlib.animation import FuncAnimation
from matplotlib.pyplot import cm
from Ex1_2 import ZK_KdV, snaps_KdV, contour_KdV

delta = 0.022

def solit(u1,x1,u2,x2):
    Delta1 = delta/sqrt(u1/12)
    Delta2 = delta/sqrt(u2/12)
    def sol_(x):
        return u1*(pi*hypsecant.pdf((x-x1)/Delta1))**2 + u2*(pi*hypsecant.pdf((x-x2)/Delta2))**2
    return sol_

#Plot semi-animé d'Instantanés
def inst_semi_anim(u, param):
    t_span = np.arange(0,8,0.1)
    n = np.zeros(np.shape(t_span)[0])
    for t in t_span:
        n = [int(t/k)]
        if n[0] < len(U):
            plt.plot(x_range, U[n[0]], label = "$t={:2.2f}\; s$".format(t), marker ='.')
            plt.title('Instantanés dépassement de solitons sur base du schéma ZK.\n $L = {}, h = {}, k = {}, T = {}, \delta = {}$'.format(x_f, h, k, t_f, delta))
            plt.xlabel('$x$')
            plt.ylabel('Amplitude')
            plt.ylim([0,20])
            plt.legend()
            plt.show(block=False)
            plt.pause(0.05)
            plt.clf()
    plt.close()
    #Plot d'Instantanés
    t_span = [0, 0.1, 3, 5]
    n = np.zeros(np.shape(t_span)[0])
    n = [int(t_span[i]/k) for i in range(np.shape(t_span)[0])]
    for i in range(np.shape(t_span)[0]):
        plt.plot(x_range, U[n[i]], label="$t = {:2.2f}s$".format(t_span[i]), marker=',')
    plt.title('Instantanés dépassement de solitons sur base du schéma ZK.\n $L = {},T = {}, h = {}, k = {}, \delta = {}$, alpha = {}, beta = {}'.format(*param))
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
    #Plot du Meshgrid
    [xx,tt]=np.meshgrid(x_range,t_range)
    plt.contourf(xx,tt, np.array(U), levels = 10)
    plt.title('Dépassement de solitons.\n $L = {}, h = {}, k = {}, T = {}, \delta = {}$'.format(x_f, h, k, t_f, delta))
    plt.xlabel("x")
    plt.ylabel("t")
    plt.colorbar()
    plt.show()
    plt.clf()

if __name__ == "__main__":
    #Initialisation ZK dépassement solitons
    ZK = ZK_KdV(solit(0.4,0.8,0.8,0.3), 0, 2, 0.005, 0.00005, 4.51)
    param = ZK[3]
    t_span = [0, 1.5, 3, 4.5]
    U = ZK[0]
    snaps_KdV(ZK, t_span, "Zabusky-Kruskal", "\operatorname{sech}(x-a)^2 + \operatorname{sech}(x-b)^2", param)
    inst_semi_anim(U)
    contour_KdV(ZK, "Zabusky-Kruskal", "\operatorname{sech}(x-a)^2 + \operatorname{sech}(x-b)^2", param)
