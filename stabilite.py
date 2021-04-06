import numpy as np
import matplotlib.pyplot as plt



def beta_1(alpha):
    lim = (alpha/4)*(1-(1/np.sqrt(alpha)))
    return lim

def beta_2(alpha):
    lim = (alpha/4)*(1+(1/np.sqrt(alpha)))
    return lim


def ampl_fact_sq1(alpha, beta,rh):
    kap_sq = ((1-2*alpha*(np.sin(rh/2))**2)**2)+np.sin(rh)**2*(alpha-4*beta*(np.sin(rh/2)**2))
    return kap_sq
def ampl_fact_sq2(alpha, beta,r, h):
    kap_sq = ((1-2*alpha*(np.sin(r*h/2))**2)**2)+np.sin(r*h)**2*(alpha-4*beta*(np.sin(r*h/2)**2))
    return kap_sq

alpha = np.arange(0., 1, 0.2)
beta_safe = np.zeros(np.shape(alpha)[0])
beta_lim_inf = np.zeros(np.shape(alpha)[0])
beta_lim_sup = np.zeros(np.shape(alpha)[0])
for i in range(np.shape(alpha)[0]):
    lim_inf = beta_1(alpha[i])
    lim_sup = beta_2(alpha[i])
    w0 = 1
    delta = 0.022
    h = 0.05
    beta = (alpha[i]/w0)*((delta**2)/2*(h**2))
    print('Pour alpha = {}, beta doit être compris entre {} et {}'.format(alpha[i], lim_inf , lim_sup))
    print('Beta par définition vaut : {}'.format(beta))
    beta_safe[i] = (lim_inf+lim_sup)/2
    beta_lim_inf[i] = lim_inf
    beta_lim_sup[i] = lim_sup


rh = np.arange(0,4*np.pi, 0.01)
#h = 0.05
#r = np.arange(0,4*np.pi, 0.01)
#rh = r*h
point_param_safe = [[],[]]
C = [np.cos(x) for x in rh]
for i in range(np.shape(alpha)[0]):
    point_param_safe.append([alpha[i], beta_safe[i]])
print(point_param_safe)
for i in range(np.shape(alpha)[0]):
    plt.plot(C, ampl_fact_sq1(alpha[i], beta_safe[i], rh), label = 'w{}'.format(i))
plt.legend()
plt.show()
