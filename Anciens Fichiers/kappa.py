import matplotlib.pyplot as plt
import numpy as np
from numpy import sqrt, pi ,exp, cos, sin

def k(alpha, beta):
    rh = np.arange(0,4*pi, 0.05)
    kap = [(1+2*alpha*(sin(x/2)**2))**2 + (sin(x)**2)*(alpha - 4*beta*(sin(x/2)**2))**2 for x in rh]
    return max(kap), min(kap)


#al = np.arange(-1.1,0.1, 0.05)
al = np.arange(-0.025,0.1, 0.01)
be = np.arange(-0.3,0.3, 0.05)
#al = np.arange(-2.1,1.1, 0.05)
#be = np.arange(-1.7,1.3, 0.05)

kappa_max = [[k(a,b)[0] for a in al] for b in be]
kappa_min = [[k(a,b)[1] for a in al] for b in be]

[aa,bb]=np.meshgrid(al,be)
levels = np.arange(0.8,1.1,0.1)
for i,ka in enumerate(kappa_max) :
    for j,k in enumerate(ka):
        if k > 1:
            kappa_max[i][j] = 20
        #else:
        #    kappa_max[i][j] = 0

plt.contourf(aa,bb, np.array(kappa_max), levels)
plt.xlabel("alpha")
plt.ylabel("beta")
plt.colorbar()
plt.show()
[aa,bb]=np.meshgrid(al,be)
plt.contourf(aa,bb, np.array(kappa_min), levels = 100)
plt.xlabel("alpha")
plt.ylabel("beta")
plt.colorbar()
plt.show()


for i,ka in enumerate(kappa_max) :
    for j,k in enumerate(ka):
        if k > 1:
            kappa_max[i][j] = 1000
        else:
            kappa_max[i][j] = 0

plt.contourf(aa,bb, np.array(kappa_max))
plt.xlabel("alpha")
plt.ylabel("beta")
plt.colorbar()
plt.show()
