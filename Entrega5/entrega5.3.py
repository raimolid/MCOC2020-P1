from numpy import zeros,cos,sin,array,dot,sqrt,linspace
from scipy.integrate import odeint
from matplotlib.pylab import *
from leer_eof import leer_eof

t,x,y,z,vx,vy,vz=leer_eof('S1B_OPER_AUX_POEORB_OPOD_20200817T110752_V20200727T225942_20200729T005942.EOF')

#Condición Inicial
z0=array([x[0],y[0],z[0],vx[0],vy[0],vz[0]])

#Condición Final
zf=array([x[-1],y[-1],z[-1],vx[-1],vy[-1],vz[-1]])


Mt= 5.972e24 #kg
G= 6.67408e-11 #m3 kg-1 s-2
omega= 7.2921150e-5


J2=1.75553e25 #m5 s-2
J3=-2.61913e29 #m6 s-2

def zp(z,t):
    zp=zeros(6)
    c=cos(omega*t)
    s=sin(omega*t)
    R=array([[c,s,0],
            [-s,c,0],
            [0,0,1]])
    Rp=omega*array([[-s,c,0],
                    [-c,-s,0],
                    [0,0,0]])
    Rpp=(omega**2)*array([[-c,s,0],
                          [s,-c,0],
                          [0,0,0]])    
    z1=z[0:3]
    z2=z[3:6]
    r2=dot(z1,z1)
    r=sqrt(r2)
    
    Fg= (-G*Mt/r**2) * (R@(z1/r))
    
    zp[0:3]=z2
    zp[3:6]=R.T@(Fg-(2*(Rp@z2)+(Rpp@z1)))
    
    # Agregando corrección J2
    Fx=J2*(z1[0]/r**7)*((6*z1[2]**2)-1.5*(z1[0]**2+z1[1]**2))
    Fy=J2*(z1[1]/r**7)*((6*z1[2]**2)-1.5*(z1[0]**2+z1[1]**2))
    Fz=J2*(z1[2]/r**7)*((3*z1[2]**2)-4.5*(z1[0]**2+z1[1]**2))
    
    zp[3]+=Fx
    zp[4]+=Fy
    zp[5]+=Fz
    
    # Agregando corrección J3 
    Fx=J3*(z1[0]*z1[2]/r**9)*((10*z1[2]**2)-7.5*(z1[0]**2+z1[1]**2))
    Fy=J2*(z1[1]*z1[2]/r**9)*((10*z1[2]**2)-7.5*(z1[0]**2+z1[1]**2))
    Fz=J2*(1/r**9)*((4*z1[2]**2)*(z1[2]**2-3*(z1[0]**2+z1[1]**2))+1.5*(z1[0]**2+z1[1]**2))
    
    zp[3]+=Fx
    zp[4]+=Fy
    zp[5]+=Fz
    
    return zp


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

sol1=odeint(zp,z0,t) 
sol_ode=sol1[:,:]
x1=sol_ode[:,0]
y1=sol_ode[:,1]
z1=sol_ode[:,2]


sol2=eulerint(zp,z0,t,1)
sol_euler=sol2[:,:]
x2=sol_euler[:,0]
y2=sol_euler[:,1]
z2=sol_euler[:,2]


deriva = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
max_value=max(deriva)

figure()
plot(t/3600,deriva/1000)
grid(True)
title(f'Distancia entre prediccines odeint y eulerint \n d_max= {max_value} (m)')
ylabel('Deriva (m)')
xlabel('Tiempo (s)')
show()
    





