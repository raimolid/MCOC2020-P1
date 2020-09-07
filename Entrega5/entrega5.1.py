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
    
    return zp


sol_ode=odeint(zp,z0,t) 
sol=sol_ode[:,:]
  
def graficar_posicion_pred():
    # 3 Gráficos tiempo vs posición x,y,z
    figure()
    for i in range(3):
        
        subplot(3,1,1+i)
        grid(True)
        plot(t,sol[:,i])
        xlabel("Tiempo (s)")
        if i==0:
            ylabel("X (m)")
            title('Posición')
        if i==1:
            ylabel("Y (m)")
        if i==2:
            ylabel("Z (m)")
    
    show()

graficar_posicion_pred()

def graficar_posicion_real():
    # 3 Gráficos tiempo vs posición x,y,z
    figure()
    for i in range(3):
        
        subplot(3,1,1+i)
        grid(True)
        xlabel("Tiempo (s)")
        if i==0:
            plot(t,x)
            ylabel("X (m)")
            title('Posición')
        if i==1:
            plot(t,y)
            ylabel("Y (m)")
        if i==2:
            plot(t,z)
            ylabel("Z (m)")
    
    show()
    
# graficar_posicion_real()

def graficar_posicion_2():
    # 3 Gráficos tiempo vs posición x,y,z
    figure()
    subplot(3,1,1)
    title('Posición')
    grid(True)
    plot(t,x)
    plot(t,sol[:,1])
    ylabel("X (m)")
    subplot(3,1,2) 
    grid(True)
    plot(t,y)
    plot(t,sol[:,2])
    ylabel("Y (m)")    
    subplot(3,1,3) 
    grid(True)
    plot(t,z)
    plot(t,sol[:,3])
    ylabel("Z (m)")        
    xlabel("Tiempo (s)")
    show()
    
# graficar_posicion_2()




    
