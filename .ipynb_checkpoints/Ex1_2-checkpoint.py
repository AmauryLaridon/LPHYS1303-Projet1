# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi ,exp, cos, sin, sqrt
from scipy.stats import hypsecant

################################### Norme |kappa|^2###############################
def k(alpha, beta):
    """Calcul du modèle carré du facteur d'amplification intervenant dans la condition de stabilité du schéma Upwind."""
    rh = np.arange(0,4*pi, 0.005)
    kap = 1-alpha*(1-np.exp(-1j*rh))-beta*(np.exp(2*rh*1j)-2*np.exp(rh*1j)+2*np.exp(-rh*1j)-np.exp(-2*rh*1j))
    mod_kap_squared = np.abs(kap)**2
    return max(mod_kap_squared)
al = np.arange(-0.25,1.1, 0.005)
be = np.arange(-0.25,0.5, 0.005)

if __name__ == "__main__":
    kappa = [[k(a,b) for a in al] for b in be]

    [aa,bb]=np.meshgrid(al,be)
    levels = np.arange(0,6,0.5)
    plt.contourf(aa,bb, np.array(kappa), levels)
    plt.xlabel("$α$")
    plt.ylabel("$β$")
    plt.title("Module carré du facteur d'amplification $\kappa$ en fonction des paramètres numériques")
    plt.colorbar()
    plt.show()

    for i,ka in enumerate(kappa) :
        for j,k in enumerate(ka):
            if k > 1:
                kappa[i][j] = 1000
            else:
                kappa[i][j] = 0

    levels = [0,1,1001]
    plt.contourf(aa,bb, np.array(kappa), levels)
    plt.xlabel("$α$")
    plt.ylabel("$β$")
    plt.title("Domaine sur lequel le module carré de $\kappa$ est inférieur à 1")
    plt.colorbar()
    plt.show()

################################### Exercice 1 Schéma Upwind ###############################

delta = 0.022

def f_cos(x):
    """Première condition initiale de cosinus."""
    return np.cos(np.pi*x)

def f_sech(x):
    """Deuxième condition initiale, solution analytique des solitons."""
    u_0 = 0.5
    u_inf = 0.1
    x_0 = 0
    Delta = delta/sqrt((u_0-u_inf)/12)
    return u_inf + (u_0-u_inf)*(pi*hypsecant.pdf((x-x_0)/Delta))**2

def Upwind_KdV(u_0, x_L, x_R, h, k, T):
    """Résolution de l'équation de KdV par le schéma Upwind"""
    delta = 0.022
    x_grid = np.arange(x_L, x_R, h)
    print("\n\nRésolution numérique avec une grille spatiale de {} points".format(len(x_grid)))
    t_grid = np.arange(0, T, k)
    print("Résolution numérique avec une grille temporelle de {} points".format(len(t_grid)))
    if len(x_grid)*len(t_grid) > 3e7 :
        print("\n Attention : long temps de calcul \n")

    L = x_R - x_L
    U = []
    U.append([u_0(x) for x in x_grid])

    w0 = max(U[0])
    alpha = (k/h)*w0
    beta = (delta**2)*k/(2*(h**3))

    print("Paramètres numériques : L = {}, T = {}s, h = {:2.6f}, k = {:2.6f}, delta = {:2.6f}, alpha = {:2.6f}, beta= {:2.6f}".format(L, T, h, k, delta, alpha, beta))


    for t in t_grid[1:]:
        u1 = [U[-1][-2], U[-1][-1], *U[-1], U[-1][0], U[-1][1]] # Conditions aux bords périodiques
        nex = []
        for i in range(len(U[-1])):
            if u1[i+2] >= 0 : # Test pour rester dans un schéma Upwind
                nex.append(u1[i+2] - (k/h) *(u1[i+2] - u1[i+1])*u1[i+2] - (delta**2) * (k/(2*(h**3))) * (u1[i+4] - 2*u1[i+3] + 2*u1[i+1] - u1[i]))
            else :
                nex.append(u1[i+2] - (k/h) *(u1[i+3] - u1[i+2])*u1[i+2] - (delta**2) * (k/(2*(h**3))) * (u1[i+4] - 2*u1[i+3] + 2*u1[i+1] - u1[i]))
        U.append(nex)
    return U, x_grid, t_grid, [L, T, h, k, delta, alpha, beta]



def ZK_KdV(u_0, x_L, x_R, h, k, T):
    """Résolution de l'équation de KdV par le schéma ZK"""

    delta = 0.022
    x_grid = np.arange(x_L, x_R, h)
    print("\n\nRésolution numérique avec une grille spatiale de {} points".format(len(x_grid)))
    t_grid = np.arange(0, T, k)
    print("Résolution numérique avec une grille temporelle de {} points".format(len(t_grid)))
    print("Temps de calcul approximatif : {} s".format(len(x_grid)*len(t_grid)*2.7161e-6))
    if len(x_grid)*len(t_grid) > 3e7 :
        print("\n    Attention : long temps de calcul \n")

    L = x_R - x_L
    U = []
    U.append([u_0(x) for x in x_grid])

    w0 = max(U[0])
    alpha = (k/h)*w0
    beta = (delta**2)*k/(2*(h**3))
    print("Paramètres numériques : L = {}, T = {}s, h = {:2.6f}, k = {:2.6f}, delta = {:2.6f}, alpha = {:2.6f}, beta= {:2.6f}".format(L, T, h, k, delta, alpha, beta))


    for t in t_grid[1:]:
        u1 = [U[-1][-2], U[-1][-1], *U[-1], U[-1][0], U[-1][1]]

        if t == k: #Initialisation pour le premier temps de calcul.
            u2 = u1
        else:
            u2 = U[-2]

        U.append([u2[i] - (k/(3*h)) *(u1[i+3] + u1[i+2] + u1[i+1]) * (u1[i+3] - u1[i+1]) - (delta**2) * (k/(h**3)) * (u1[i+4] - 2*u1[i+3] + 2*u1[i+1] - u1[i])  for i in range(len(U[-1]))])
    return U, x_grid, t_grid, [L, T, h, k, delta, alpha, beta]


def snaps_KdV(U, t_range, schema, CI, parametres):
    """Plot d'instantanées pour l'équation de KdV"""
    ma = []
    mi = []
    for t in t_range:
        plt.plot(U[1], U[0][int(t/parametres[3])], label="$t = {:2.2f}s$".format(t), marker = ".")
        ma.append(max(U[0][int(t/parametres[3])]))
        mi.append(min(U[0][int(t/parametres[3])]))
    MA = max(ma)
    MI = min(mi)

    plt.ylim([MI-0.2, MA+0.2])
    plt.xlabel("$x$")
    plt.ylabel("$u(x,t)$")
    plt.title('Instantanés de la résolution de KdV par le schéma {}, CI = ${}$ ,\n L = {}, T = {}s, h = {}, k = {}, $\delta$ = {},  alpha = {:2.6f}, beta = {:2.6f}'.format(schema, CI, *parametres))
    plt.legend()
    plt.show()
    plt.close()

def contour_KdV(U, schema, CI, parametres):
    """Plot de contour plein pour l'équation de KdV"""
    [xx,tt]=np.meshgrid(U[1],U[2])
    plt.contourf(xx,tt, np.array(U[0]), levels = 30)
    plt.title('Graphes de contour de la résolution de KdV par le schéma {}, CI = ${}$ ,\n L = {}, T = {}s, h = {}, k = {}, $\delta$ = {},  alpha = {:2.6f}, beta = {:2.6f}'.format(schema, CI, *parametres))
    plt.xlabel("x")
    plt.ylabel("t")
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    #Initialisation Upwind cos
    Upwind = Upwind_KdV(f_cos, 0, 2, 0.01, 0.00001, 1.3)
    param = Upwind[3]
    t_span = [0, 1/np.pi, 3.6/np.pi]
    #Affichage
    snaps_KdV(Upwind, t_span, "Upwind", "\cos(\pi x)", param)
    contour_KdV(Upwind, "Upwind", "\cos(\pi x)", param)

    #Initialisation Upwind soliton
    Upwind = Upwind_KdV(f_sech, -0.4, 0.6, 0.008, 0.000008, 1.01)
    param = Upwind[3]
    t_span = [0, 0.25, 0.5, 1]
    snaps_KdV(Upwind, t_span, "Upwind", "\operatorname{sech}(x)^2", param)
    contour_KdV(Upwind, "Upwind", "\operatorname{sech}(x)^2", param) #fait crash mon ordi


    #Initialisation ZK cos
    ZK = ZK_KdV(f_cos, 0, 2, 0.008, 0.000008, 1.15)
    param = ZK[3]
    t_span = [0, 1/np.pi, 3.6/np.pi]
    snaps_KdV(ZK, t_span, "Zabusky-Kruskal", "\cos(\pi x)", param)
    contour_KdV(ZK, "Zabusky-Kruskal", "\cos(\pi x)", param)

    #Initialisation ZK soliton
    ZK = ZK_KdV(f_sech, -0.4, 0.6, 0.003, 0.00001, 1.01)
    param = ZK[3]
    t_span = [0, 0.25, 0.5, 1]
    snaps_KdV(ZK, t_span, "Zabusky-Kruskal", "\operatorname{sech}(x)^2", param)
    contour_KdV(ZK, "Zabusky-Kruskal", "\operatorname{sech}(x)^2", param)
