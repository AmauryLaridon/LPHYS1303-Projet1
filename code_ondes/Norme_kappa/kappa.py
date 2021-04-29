# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi ,exp, cos, sin

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