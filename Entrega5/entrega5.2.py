from numpy import zeros,cos,sin,array,dot,sqrt,linspace
from scipy.integrate import odeint
from matplotlib.pylab import *
from leer_eof import leer_eof
from time import perf_counter

sat_t,sat_x,sat_y,sat_z,sat_vx,sat_vy,sat_vz=leer_eof('S1B_OPER_AUX_POEORB_OPOD_20200817T110752_V20200727T225942_20200729T005942.EOF')


# Condición Inicial
z0=array([sat_x[0],sat_y[0],sat_z[0],sat_vx[0],sat_vy[0],sat_vz[0]])

# Condición Final
zf=array([sat_x[-1],sat_y[-1],sat_z[-1],sat_vx[-1],sat_vy[-1],sat_vz[-1]])


# Datos
Mt=5.972e24 #kg 
G=6.67408e-11 #m3 kg-1 s-2
omega=7.2921150e-5 #rad s-1


def zp(z,t):
    
    # Variables Internas 
    zp=zeros(6)
    z1=z[0:3]
    z2=z[3:6]
    r2=dot(z1,z1)
    r=sqrt(r2)

    # Matrices de Rotación 
    R=array([[cos(omega*t),-sin(omega*t),0],
             [sin(omega*t),cos(omega*t) ,0],
             [0           ,  0          ,1]])
    
    Rp=omega*array([[-sin(omega*t),-cos(omega*t),0],
                    [cos(omega*t) ,-sin(omega*t),0],
                    [0            ,  0          ,0]]) 
    
    Rpp=(omega**2)*array([[-cos(omega*t),sin(omega*t) ,0],
                          [-sin(omega*t),-cos(omega*t),0],
                          [0            ,  0          ,0]])
    
    
    # Vector de salida 
    zp[0:3] = z2 
    zp[3:6]= (-G*Mt/r**3)*z1  - R.T@(Rpp@z1 + 2.*Rp@z2)
    
    return zp


def eulerint(zp,z0,t,Nsubdivisiones = 1):
    Nt = len(t)
    Ndim = len(array(z0))
    z = zeros((Nt,Ndim))
    z[0,:] = z0
    for i in range(1,Nt):
        t_anterior = t[i-1]
        dt = (t[i] - t[i-1])/Nsubdivisiones
        z_temp = z[i-1,:].copy()
        for k in range(Nsubdivisiones):
            z_temp += dt * zp(z_temp,t_anterior + k*dt)
        z[i,:] = z_temp
    return z


t1=perf_counter()
sol1=odeint(zp,z0,sat_t)
t2=perf_counter()
dt=t2-t1
print(f'Tiempo Odeint: {dt}\n')
sol_ode=sol1[:,:]
x_ode=sol_ode[:,0]
y_ode=sol_ode[:,1]
z_ode=sol_ode[:,2]

t3=perf_counter()
sol2=eulerint(zp,z0,sat_t,1)
t4=perf_counter()
dt=t4-t3
print (f'Tiempo Eulerint: {dt}\n')
sol_euler=sol2[:,:]
x_euler=sol_euler[:,0]
y_euler=sol_euler[:,1]
z_euler=sol_euler[:,2]



def graficar_deriva_ode_euler():
    
    figure()
    deriva_ode= sqrt( (sat_x - x_ode)**2 + (sat_y - y_ode)**2 + (sat_z - z_ode)**2 )
    plot(sat_t/3600,deriva_ode/1000, label='Odeint')
    
    deriva_euler = sqrt( (sat_x - x_euler)**2 + (sat_y - y_euler)**2 + (sat_z - z_euler)**2 )
    plot(sat_t/3600,deriva_euler/1000, label='Eulerint')
    
    title(f'Deriva Odeint v/s Eulerint')
    ylabel('Deriva,d [Km]')
    xlabel('Tiempo, t [horas]')
    grid(True)
    tight_layout()
    legend() 
    show()
    
graficar_deriva_ode_euler()
 

# # Deriva real-ode con J2J3
# pos_final=zf-sol1[-1]

# norma_distancia_error = sqrt(pos_final[0]**2 + pos_final[1]**2 + pos_final[2]**2)
# print (norma_distancia_error)

# # Deriva real-euler con J2J3
# pos_final=zf-sol2[-1]

# norma_distancia_error = sqrt(pos_final[0]**2 + pos_final[1]**2 + pos_final[2]**2)
# print (norma_distancia_error)

