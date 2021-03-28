import numpy as np
import matplotlib.pyplot as plt

###Fonctions initiales###
def Create_Init(x,x1,x2,yy):
    y=np.zeros(np.size(x))
    for i in range(np.size(x)):
        if x[i]>x1 and x[i]<x2:
            y[i]=yy
    return(y)
def f(y,a,b):
    k=y.size
    l=b-a
    h=np.zeros(k)
    for p in range(k):
        if y[p]>=(b+a)/2-(b-a)/4 and y[p]<=(b+a)/2+(b-a)/4:
            h[p]=0.5*(1.-np.cos(2.*np.pi*((y[p]-(b+a)/2)/l)))
            h[p]=(1+np.cos(2.*np.pi*((y[p]-(b+a)/2)/(l/2))))*0.5
    return(h)
def f_cos(y,a,b):
    """fonction initiale de cos"""
    k=y.size
    l=b-a
    h=np.zeros(k)
    for p in range(k):
        h[p] = np.cos(np.pi*y[p])
    return(h)

###Méthods de différences finies###
def Leapfrog_adv(U0,a,x,tmax,lmb):
    """Schéma Leapfrog pour l'advection classique"""
    h=np.abs(x[1]-x[0])
    k=lmb*h/(abs(a))
    M=int(np.round(tmax/k))+1
    t=np.linspace(0.,tmax,M)
    k=t[1]-t[0]
    lmb=np.abs(a)*k/h
    U=np.zeros((np.shape(x)[0],np.shape(t)[0]))
    U[:,0]=U0
    for i in range(M-1):
        if i==0:
            U[1:-1,i+1]=(1-lmb)*U[1:-1,i]+lmb*U[0:-2,i]
            U[0,i+1]=0
            U[-1,i+1]=0
        else:
            U[1:-1,i+1]=U[1:-1,i-1]-lmb*(U[2:,i]-U[0:-2,i])
            U[0,i+1]=0
            U[-1,i+1]=0
    return(U,x,t)



def integration_adv(U0,a,x,tmax,lmb,int_mode):
    h=np.abs(x[1]-x[0])
    k=lmb*h/(abs(a))
    M=int(np.round(tmax/k))+1
    t=np.linspace(0.,tmax,M)
    k=t[1]-t[0]
    lmb=np.abs(a)*k/h
    U=np.zeros((np.shape(x)[0],np.shape(t)[0]))
    U[:,0]=U0
    #Upwind
    if int_mode=='Upwind':
        for i in range(M-1):
            U[1:-1,i+1]=(1-lmb)*U[1:-1,i]+lmb*U[0:-2,i]
            U[0,i+1]=0
            U[-1,i+1]=0
    if int_mode=='UpwindKDV':
        for i in range(M-1):
            
    elif int_mode=='Downwind':
        for i in range(M-1):
            U[1:-1,i+1]=(1+lmb)*U[1:-1,i]-lmb*U[2:,i]
            U[0,i+1]=0
            U[-1,i+1]=0
    elif int_mode=='Lax-Wendroff':
        A=-(1/2)*lmb*(1-lmb)
        B=(1-lmb**2)
        C=0.5*lmb*(1+lmb)
        for i in range(M-1):
            U[1:-1,i+1]=A*U[2:,i]+B*U[1:-1,i]+C*U[0:-2,i]
            U[0,i+1]=0
            U[-1,i+1]=0
    if int_mode=='LeapFrog':
        for i in range(M-1):
            if i==0:
                U[1:-1,i+1]=(1-lmb)*U[1:-1,i]+lmb*U[0:-2,i]
                U[0,i+1]=0
                U[-1,i+1]=0
            else:
                U[1:-1,i+1]=U[1:-1,i-1]-lmb*(U[2:,i]-U[0:-2,i])
                U[0,i+1]=0
                U[-1,i+1]=0
    else:
        print('Unknown integration scheme')
        U=0
    return(U,x,t)

a=1

lmb=1.
L=10
N=1001
x=np.linspace(0.,L,N+2)
tmax=100

U0_carre=Create_Init(x,4,5,1)
a_cos=0
b_cos=20
U0_cos=f(x,a_cos,b_cos)
#U_cos = f_cos(x,0, 2)



"""U_up_carre=integration_adv(U0_carre,a,x,tmax,lmb,'Upwind')
U_down_carre=integration_adv(U0_carre,a,x,tmax,lmb,'Downwind')
U_lw_carre=integration_adv(U0_carre,a,x,tmax,lmb,'Lax-Wendroff')
U_lp_carre=Leapfrog_adv(U0_carre,a,x,tmax,lmb)

U_up_cos=integration_adv(U0_cos,a,x,tmax,lmb,'Upwind')
U_down_cos=integration_adv(U0_cos,a,x,tmax,lmb,'Downwind')
U_lw_cos=integration_adv(U0_cos,a,x,tmax,lmb,'Lax-Wendroff')
U_lp_cos=Leapfrog_adv(U0_cos,a,x,tmax,lmb)

t=U_up_carre[2]
M=t.shape[0]
Mspan=[0,0.25,0.5,0.75,0.99]

fig,ax=plt.subplots(2,5,figsize=(15,7),sharey='row')
for i in range(len(Mspan)):
    n=int(np.floor(M*Mspan[i]))
    center_carre=4.5+a*t[n]
    center_cos=(b_cos+a_cos)/2+a*t[n]
    ax[0,i].plot(x,U_up_carre[0][:,n],label='Upwind',linewidth=2.)
    ax[0,i].plot(x,U_down_carre[0][:,n],label='Downwind',linewidth=2.)
    ax[0,i].plot(x,U_lw_carre[0][:,n],label='Lax-Wendroff',linewidth=2.)
    ax[0,i].plot(x,U_lp_carre[0][:,n],label='Leapfrog',linewidth=2.)
    ax[0,i].plot(x,Create_Init(x-a*t[n],4,5,1),label='Analytique')
    #ax[0,i].set_xlabel('x')
    ax[0,i].set_ylabel('y')
    ax[0,i].set_ylim((-0.5,1.5))
    ax[0,i].set_xlim((center_carre-4,center_carre+4))
    ax[0,i].set_title('Time=%.2f, lambda=%.3f'%(t[n],lmb))
    ax[0,i].legend(loc='best',fontsize='small')
    ax[1,i].plot(x,U_up_cos[0][:,n],label='Upwind',linewidth=2.)
    ax[1,i].plot(x,U_down_cos[0][:,n],label='Downwind',linewidth=2.)
    ax[1,i].plot(x,U_lw_cos[0][:,n],label='Lax-Wendroff',linewidth=2.)
    ax[1,i].plot(x,U_lp_cos[0][:,n],label='Leapfrog',linewidth=2.)
    ax[1,i].plot(x,f(x-a*t[n],a_cos,b_cos),label='Analytique')
    ax[1,i].set_xlabel('x')
    ax[1,i].set_ylabel('y')
    ax[1,i].set_ylim((-0.5,1.5))
    ax[1,i].set_xlim((center_cos-(b_cos-a_cos)/4-4,center_cos+(b_cos-a_cos)/4+4))
    ax[1,i].set_title('Time=%.2f, lambda=%.3f'%(t[n],lmb))
    ax[1,i].legend(loc='best',fontsize='small')
plt.show()

print("Nombre de points d'espace ={}, nombre de points de temps = {}".format(N+2, t.shape[0]))
fig,ax=plt.subplots(2,5,figsize=(15,7),sharey='row')
for i in range(len(Mspan)):
    n=int(np.floor(M*Mspan[i]))
    center_carre=4.5+a*t[n]
    center_cos=(b_cos+a_cos)/2+a*t[n]
    ax[0,i].plot(x,U_up_carre[0][:,n],label='Upwind',linewidth=2.)
    ax[0,i].plot(x,U_down_carre[0][:,n],label='Downwind',linewidth=2.)
    ax[0,i].plot(x,U_lw_carre[0][:,n],label='Lax-Wendroff',linewidth=2.)
    ax[0,i].plot(x,Create_Init(x-a*t[n],4,5,1),label='Analytique')
    ax[0,i].set_ylabel('y')
    ax[0,i].set_ylim((-0.5,1.5))
    ax[0,i].set_xlim((center_carre-4,center_carre+4))
    ax[0,i].set_title('Time=%.2f, lambda=%.1f'%(t[n],lmb))
    ax[0,i].legend(loc='best',fontsize='small')
    ax[1,i].plot(x,U_up_cos[0][:,n],label='Upwind',linewidth=2.)
    ax[1,i].plot(x,U_down_cos[0][:,n],label='Downwind',linewidth=2.)
    ax[1,i].plot(x,U_lw_cos[0][:,n],label='Lax-Wendroff',linewidth=2.)
    ax[1,i].plot(x,f(x-a*t[n],a_cos,b_cos),label='Analytique')
    ax[1,i].set_xlabel('x')
    ax[1,i].set_ylabel('y')
    ax[1,i].set_ylim((-0.5,1.5))
    ax[1,i].set_xlim((center_cos-4,center_cos+4))
    ax[1,i].set_title('Time=%.2f, lambda=%.1f'%(t[n],lmb))
    ax[1,i].legend(loc='best',fontsize='small')
plt.show()

###Analyse stabilité des schémas###
#Upwind
ko=np.linspace(0.01,np.pi-0.001,100)

lmb=[0.25,0.5,0.75,1.,1.0]
fig,ax=plt.subplots(1,2,figsize=(15,5),subplot_kw=dict(polar=True))
for i in range(5):
    ampl_R=1-lmb[i]+lmb[i]*np.cos(ko)
    ampl_I=-lmb[i]*np.sin(ko)
    ampl=1-lmb[i]+lmb[i]*np.cos(ko)+(1j)*(-lmb[i]*np.sin(ko))
    ratio=np.arctan(ampl_I/ampl_R)/(-lmb[i]*ko)
    ratio=np.arctan2(np.imag(ampl),np.real(ampl))/(-lmb[i]*ko)
    ax[0].plot(ko,np.sqrt(ampl_R**2+ampl_I**2),label='lambda=%.2f'%(lmb[i]))
    ax[1].plot(ko,np.abs(ratio),label='lambda=%.2f'%(lmb[i]))

fig.suptitle('Schéma Upwind')
ax[0].set_rmax(2)
ax[0].set_rgrids([0.0,0.5, 1, 1.5],labels=['0.0','0.5', '1', '1.5'],fontsize=20)
ax[0].set_frame_on(False)
ax[0].set_thetamax(180)
ax[0].set_thetagrids(angles=[0,180],labels=['rh=0','rh=$\pi$'],fontsize=15)
ax[0].legend()
ax[0].set_title(r'$|\kappa |$')
ax[1].set_rmax(2)
ax[1].set_rgrids([0.0,0.5, 1, 1.5],labels=['0.0','0.5', '1', '1.5'],fontsize=20)
ax[1].set_frame_on(False)
ax[1].set_thetamax(180)
ax[1].set_thetagrids(angles=[0,180],labels=['rh=0','rh=$\pi$'],fontsize=15)
ax[1].legend()
ax[1].set_title(r'$\frac{\theta_{d}}{\theta_{a}}$')
plt.show()

#Lax-Wendroff
ko=np.linspace(0.01,np.pi-0.001,100)

lmb=[0.25,0.5,0.75,1.]
fig,ax=plt.subplots(1,2,figsize=(15,5),subplot_kw=dict(polar=True))
for i in range(4):
    ampl_R=1-lmb[i]**2+lmb[i]**2*np.cos(ko)
    ampl_I=-lmb[i]*np.sin(ko)
    ampl=1-lmb[i]**2+lmb[i]**2*np.cos(ko)+(1j)*(-lmb[i]*np.sin(ko))
    ratio=np.arctan(ampl_I/ampl_R)/(-lmb[i]*ko)
    ratio=np.arctan2(np.imag(ampl),np.real(ampl))/(-lmb[i]*ko)
    ax[0].plot(ko,np.sqrt(ampl_R**2+ampl_I**2),label='lambda=%.2f'%(lmb[i]))
    ax[1].plot(ko,np.abs(ratio),label='lambda=%.2f'%(lmb[i]))

fig.suptitle('Schéma Lax-Wendroff')
ax[0].set_rmax(2)
ax[0].set_rgrids([0.0,0.5, 1, 1.5],labels=['0.0','0.5', '1', '1.5'],fontsize=20)
ax[0].set_frame_on(False)
ax[0].set_thetamax(180)
ax[0].set_thetagrids(angles=[0,180],labels=['rh=0','rh=$\pi$'],fontsize=15)
ax[0].legend()
ax[0].set_title(r'$|\kappa |$')
ax[1].set_rmax(2)
ax[1].set_rgrids([0.0,0.5, 1, 1.5],labels=['0.0','0.5', '1', '1.5'],fontsize=20)
ax[1].set_frame_on(False)
ax[1].set_thetamax(180)
ax[1].set_thetagrids(angles=[0,180],labels=['rh=0','rh=$\pi$'],fontsize=15)
ax[1].legend()
ax[1].set_title(r'$\frac{\theta_{d}}{\theta_{a}}$')
plt.show()


#Leapfrog
ko=np.linspace(0.01,np.pi-0.001,100)

lmb=[0.25,0.5,0.75,0.9]
fig,ax=plt.subplots(1,2,figsize=(15,5),subplot_kw=dict(polar=True))
for i in range(4):
    ampl_R=1-lmb[i]+lmb[i]*np.cos(ko)
    ampl_I=-lmb[i]*np.sin(ko)
    ampl=(-1j*lmb[i]*np.sin(ko)+np.sqrt(1-(lmb[i]*np.sin(ko))**2))
    ratio=np.arctan2(np.imag(ampl),np.real(ampl))/(-lmb[i]*ko)
    ax[0].plot(ko,np.abs(ampl),label='lambda=%.2f'%(lmb[i]))
    ax[1].plot(ko,ratio,label='lambda=%.2f'%(lmb[i]))
plt.show()
fig.suptitle('Schéma Leapfrog')
ax[0].set_rmax(2)
ax[0].set_rgrids([0.0,0.5, 1, 1.5],labels=['0.0','0.5', '1', '1.5'],fontsize=20)
ax[0].set_frame_on(False)
ax[0].set_thetamax(180)
ax[0].set_thetagrids(angles=[0,180],labels=['rh=0','rh=$\pi$'],fontsize=15)
ax[0].legend()
ax[0].set_title(r'$|\kappa |$')
ax[1].set_rmax(2)
ax[1].set_rgrids([0.0,0.5, 1, 1.5],labels=['0.0','0.5', '1', '1.5'],fontsize=20)
ax[1].set_frame_on(False)
ax[1].set_thetamax(180)
ax[1].set_thetagrids(angles=[0,180],labels=['rh=0','rh=$\pi$'],fontsize=15)
ax[1].legend()
ax[1].set_title(r'$\frac{\theta_{d}}{\theta_{a}}$')
plt.show()
"""
