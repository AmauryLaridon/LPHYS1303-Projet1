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
    """Condition initiale de deux solitons."""
    Delta1 = delta/sqrt(u1/12)
    Delta2 = delta/sqrt(u2/12)
    def sol_(x):
        return u1*(pi*hypsecant.pdf((x-x1)/Delta1))**2 + u2*(pi*hypsecant.pdf((x-x2)/Delta2))**2
    return sol_

def inst_semi_anim(u, x_range, param):
    """Plot semi-animé d'instantanées."""
    t_span = np.arange(0,param[1],0.1)
    n = np.zeros(np.shape(t_span)[0])
    for t in t_span:
        n = [int(t/param[3])]
        if n[0] < len(U):
            plt.plot(x_range, U[n[0]], label = "$t={:2.2f}\; s$".format(t), marker ='.')
            plt.title('Instantanés dépassement de solitons sur base du schéma ZK.\n $L = {},T = {}, h = {}, k = {}, \delta = {:2.6f}, alpha = {:2.6f}, beta = {:2.6f}$'.format(*param))
            plt.xlabel('$x$')
            plt.ylabel('Amplitude')
            plt.ylim([0,np.array(u).max()+0.2])
            plt.legend()
            plt.show(block=False)
            plt.pause(0.0001)
            plt.clf()
    plt.close()

if __name__ == "__main__":
    #Initialisation ZK dépassement solitons
    ZK = ZK_KdV(solit(0.4,0.8,0.8,0.3), 0, 2, 0.005, 0.00005, 4.51)
    #ZK = ZK_KdV(solit(0.4,0.8,0.8,0.3), 0, 2, 0.01, 0.0003, 2.51) essai résolution plus basse pour le plot semi anim
    x_span = ZK[1]
    param = ZK[3]
    t_span = [0, 1.5, 3, 4.5]
    U = ZK[0]
    snaps_KdV(ZK, t_span, "Zabusky-Kruskal", "\operatorname{sech}(x-a)^2 + \operatorname{sech}(x-b)^2", param)
    inst_semi_anim(U,x_span, param)
    contour_KdV(ZK, "Zabusky-Kruskal", "\operatorname{sech}(x-a)^2 + \operatorname{sech}(x-b)^2", param)
