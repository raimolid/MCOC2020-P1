from matplotlib.pylab import *
from scipy.integrate import odeint
from numpy import cos,sin ,pi

#Notar que eulerint es inestable
#Elegir paso estable
#Odeint aumenta subdivisiones automatico

def eulerint(zp,z0,t,Nsubdivisiones=1):
    Nt=len(t)
    Ndim=len(array([z0]))
    
    z=zeros((Nt,Ndim))
    z[0,:]=z0
    
    for i in range(1,Nt):
        t_anterior=t[i-1]
        
        dt=(t[i]-t[i-1])/Nsubdivisiones
        
        z_temp=z[i-1,:].copy()
        for k in range(Nsubdivisiones):
            z_temp+=dt*zp(z_temp,t_anterior + k*dt)
        z[i,:]=z_temp 
        
    return z

def zp(z,t):
    x=(0.29**t)*cos(352.73*t)+(0.37*0.29**t)*sin(352.73*t)
    xp=(1.02*0.29**t)*(cos(352.73*t)-6.48*sin(352.73*t))
    k=4*pi**2
    c=2.51
    m=1
    zp=(-k*x-c*xp)/m
    return zp

z0=1.
t=linspace(0,4.,100)

# Solución de Odeint
z_odeint=odeint(zp,z0,t)

# Solución Analítica 
z_real=(-41.95*0.29**t)*(cos(352.73*t)-0.05*sin(352.73*t))

# Solución de Eulerint 
z_euler1=eulerint(zp,z0,t,1)

z_euler2=eulerint(zp,z0,t,10)

z_euler3=eulerint(zp,z0,t,100)

plot(t,z_odeint,label='odeint',color='b')
plot(t,z_real,label='real',color='k',linewidth=2)
plot(t,z_euler1,label='euler N=1',color='g',linestyle='--')
plot(t,z_euler2,label='euler N=10',color='r',linestyle='--')
plot(t,z_euler3,label='euler N=100',color='orange',linestyle='--')

legend()
show()

            