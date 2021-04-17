# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi ,exp, cos, sin
from scipy.stats import hypsecant

################################### Norme |kappa|^2###############################
def k(alpha, beta):
    rh = np.arange(0,4*pi, 0.05)
    kap = [(1+2*alpha*(sin(x/2)**2))**2 + (sin(x)**2)*(alpha - 4*beta*(sin(x/2)**2))**2 for x in rh]
    return max(kap)

al = np.arange(-1.1,0.1, 0.05)
be = np.arange(-0.7,0.3, 0.05)

kappa = [[k(a,b) for a in al] for b in be]

[aa,bb]=np.meshgrid(al,be)
levels = np.arange(0,6,0.5)
plt.contourf(aa,bb, np.array(kappa), levels)
plt.xlabel("alpha")
plt.ylabel("beta")
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
plt.xlabel("alpha")
plt.ylabel("beta")
plt.colorbar()
plt.show()






################################### Exercice 1 Schéma Upwind ###############################
def f_cos(x):
    #fonction initiale de cos
    return np.cos(np.pi*x)

def f_sech(x):
    return 0.5*(pi*hypsecant.pdf(x/2))**2



def Upwind_KdV(u_0, delta, x_L, x_R, h, k, T):
    x_grid = np.arange(x_L, x_R - h, h)
    print("Résolution numérique avec une grille spatiale de {} points".format(len(x_grid)))
    t_grid = np.arange(0, T-k, k)
    print("Résolution numérique avec une grille temporelle de {} points".format(len(t_grid)))

    L = x_R - x_L
    U = []
    U.append([u_0(x) for x in x_grid])

    w0 = max(U[0])
    alpha = (k/h)*w0
    beta = (delta**2)*k/(2*(h**3))
    print("Paramètres numérique : L = {}, T = {}s, h = {:2.4f}, k = {:2.4f}, delta = {:2.4f}, alpha = {:2.4f}, beta= {:2.4f}".format(L, T, h, k, delta, alpha, beta))


    for t in t_grid[1:]:
        u1 = [*U[-1], U[-1][0], U[-1][1]] # Conditions aux bord périodiques
        nex = []
        for i in range(len(U[-1])):
            if u1[i] >= 0 :
                nex.append(u1[i] - (k/h) *(u1[i] - u1[i-1])*u1[i] - (delta**2) * (k/(2*(h**3))) * (u1[i+2] - 2*u1[i+1] + 2*u1[i-1] - u1[i-2]))
            else :
                nex.append(u1[i] - (k/h) *(u1[i+1] - u1[i])*u1[i] - (delta**2) * (k/(2*(h**3))) * (u1[i+2] - 2*u1[i+1] + 2*u1[i-1] - u1[i-2]))
        U.append(nex)
    return U, x_grid, t_grid, [L, T, h, k, delta, alpha, beta]



def ZK_KdV(u_0, delta, x_L, x_R, h, k, T):
    x_grid = np.arange(x_L, x_R - h, h)
    print("Résolution numérique avec une grille spatiale de {} points".format(len(x_grid)))
    t_grid = np.arange(0, T-k, k)
    print("Résolution numérique avec une grille temporelle de {} points".format(len(t_grid)))

    L = x_R - x_L
    U = []
    U.append([u_0(x) for x in x_grid])

    w0 = max(U[0])
    alpha = (k/h)*w0
    beta = (delta**2)*k/(2*(h**3))
    print("Paramètres numérique : L = {}, T = {}s, h = {:2.4f}, k = {:2.4f}, delta = {:2.4f}, alpha = {:2.4f}, beta= {:2.4f}".format(L, T, h, k, delta, alpha, beta))


    for t in t_grid[1:]:
        u1 = [*U[-1], U[-1][0], U[-1][1]]

        if t == k:
            u2 = u1
        else:
            u2 = U[-2]

        U.append([u2[i] - (k/(3*h)) *(u1[i+1] + u1[i] + u1[i-1]) * (u1[i+1] - u1[i-1]) - (delta**2) * (k/(h**3)) * (u1[i+2] - 2*u1[i+1] + 2*u1[i-1] - u1[i-2])  for i in range(len(U[-1]))])
    return U, x_grid, t_grid, [L, T, h, k, delta, alpha, beta]




#Plot d'instantanné avec la CI de cos()
def snaps_KdV(U, t_range, schema, CI, parametres):
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
    plt.title('Instantanés de la résolution de KdV par le schéma {}, CI = ${}$ ,\n L = {}, T = {}s, h = {}, k = {}, $\delta$ = {},  alpha = {:2.4f}, beta = {:2.4f}'.format(schema, CI, *parametres))
    plt.legend()
    plt.show()
    plt.close()

def contour_KdV(U, schema, CI, parametres):
    #Plot 2D
    [xx,tt]=np.meshgrid(U[1],U[2])
    plt.contourf(xx,tt, np.array(U[0]))
    plt.title('Meshgrid de la résolution de KdV par le schéma {}, CI = ${}$ ,\n L = {}, T = {}s, h = {}, k = {}, $\delta$ = {},  alpha = {:2.4f}, beta = {:2.4f}'.format(schema, CI, *parametres))
    plt.xlabel("x")
    plt.ylabel("t")
    plt.colorbar()
    plt.show()


#Initialisation Upwind cos
Upwind = Upwind_KdV(f_cos, 0.022, 0, 2, 0.01, 0.00001, 1.3)
param = Upwind[3]
t_span = [0, 1/np.pi, 3.6/np.pi]
snaps_KdV(Upwind, t_span, "Upwind", "cos(\pi x)", param)
contour_KdV(Upwind, "Upwind", "cos(\pi x)", param)



#Initialisation Upwind soliton
Upwind = Upwind_KdV(f_sech, 0.022, -10, 10, 0.1739, 0.0001, 2.01)
param = Upwind[3]
t_span = [0, 0.5, 1, 2]
snaps_KdV(Upwind, t_span, "Upwind", "0.5\ \sech(x/2)^2", param)
contour_KdV(Upwind, "Upwind", "0.5\  \sech(x/2)^2", param)


#Initialisation ZK cos
ZK = ZK_KdV(f_cos, 0.022, 0, 2, 0.01, 0.00001, 1.3)
param = ZK[3]
t_span = [0, 1/np.pi, 3.6/np.pi]
snaps_KdV(ZK, t_span, "Zabusky-Kruskal", "cos(\pi x)", param)
contour_KdV(ZK, "Zabusky-Kruskal", "cos(\pi x)", param)



#Initialisation ZK soliton
ZK = ZK_KdV(f_sech, 0.022, -10, 10, 0.1739, 0.0001, 2.01)
param = ZK[3]
t_span = [0, 0.5, 1, 2]
snaps_KdV(ZK, t_span, "Zabusky-Kruskal", "0.5\ \sech(x/2)^2", param)
contour_KdV(ZK, "Zabusky-Kruskal", "0.5\  \sech(x/2)^2", param)
