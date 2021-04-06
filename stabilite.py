import numpy as np
import matplotlib.pyplot as plt


def ampl_fact_sq(alpha, beta, h, r):
    kap_sq = ((1+alpha*(np.sin(r*h/2))**2)**2)*((alpha*np.cos(r*h/2)*np.sin(r*h/2)+2*beta*(np.sin(2*r*h)-2*np.sin(r*h)))**2)
    return kap_sq


alpha = np.linspace(0.1,0.4, 6)
beta = np.linspace(0.0, 0.002, 6)
r = np.arange(0,20*np.pi, 0.01)
h = 0.05
point_param = [[],[]]
for i in range(np.shape(alpha)[0]):
    point_param.append([alpha[i], beta[i]])
print(point_param)
for i in range(np.shape(alpha)[0]):
    plt.plot(r, ampl_fact_sq(alpha[i], beta[i], h, r), label = 'w{}'.format(i))
plt.legend()
plt.show()
