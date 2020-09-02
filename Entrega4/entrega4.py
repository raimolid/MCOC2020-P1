from matplotlib.pylab import *
from scipy.integrate import odeint
from numpy import array,zeros,cos,sin,pi,exp 

# Notar que eulerint es inestable
# Elegir paso estable
# Odeint aumenta subdivisiones automatico

# Datos
M = 1.
f = 1.
e = 0.2
w = 2*pi*f 
k = M*w**2
c = 2*e*w*M

def eulerint(zp,z0,t,Nsubdivisiones = 1):
    Nt = len(t)
    Ndim = len(array(z0))
    z = zeros((Nt,Ndim))
    z[0,:] = z0[0]
    for i in range(1,Nt):
        t_anterior = t[i-1]
        dt = (t[i] - t[i-1])/Nsubdivisiones
        z_temp = z[i-1,:].copy()
        for k in range(Nsubdivisiones):
            z_temp += dt * zp(z_temp,t_anterior + k*dt)
        z[i,:] = z_temp
    return z


def zp(z,t):
    zp = zeros(2)
    zp[0] = z[1]
    z1 = z[0]
    z2 = z[1]
    zp[1]=-(c*z2+k*z1)/M
    return zp

plt.style.use('default')

# Condición Inicial
z0=[1.,1.]
t=linspace(0,4.,100)

# Solución Analítica 
z_real = (exp((-c/2)*t))*cos(w*t) 
plot(t,z_real,label='real',color='k',linewidth=2)

# Solución de Odeint
z_odeint=odeint(zp,z0,t)
x=z_odeint[:,0]
plot(t,x,label='odeint',color='b')

# Solución de Eulerint 
z_euler1=eulerint(zp,z0,t,1)
x=z_euler1[:,0]
plot(t,x,label='euler N=1',color='g',linestyle='--')

z_euler2=eulerint(zp,z0,t,10)
x=z_euler2[:,0]
plot(t,x,label='euler N=10',color='r',linestyle='--')

z_euler3=eulerint(zp,z0,t,100)
x=z_euler3[:,0]
plot(t,x,label='euler N=100',color='orange',linestyle='--')

title('Oscilador Armónico')
xlabel('Tiempo [s]')
ylabel('X(t)  [m]')
legend()
show()

            
