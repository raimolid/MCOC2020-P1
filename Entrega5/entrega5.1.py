from numpy import zeros,cos,sin,array,dot,sqrt,linspace
from scipy.integrate import odeint
from matplotlib.pylab import *
from leer_eof import leer_eof

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

sol1=odeint(zp,z0,sat_t) 
sol_ode=sol1[:,:]
x_ode=sol_ode[:,0]
y_ode=sol_ode[:,1]
z_ode=sol_ode[:,2]


def graficar_posicion_real_pred():
    # 3 Gráficos tiempo vs posición x,y,z
    figure()
    for i in range(3):
        
        subplot(3,1,1+i)
        grid(True)
        plot(sat_t,sol_ode[:,i], label='Odeint')
        xlabel('Tiempo, t [horas]')
        if i==0:
            plot(sat_t,sat_x)
            ylabel('X [Km]')
            title('Posición real v predicha')
            xtick = [0, 18000, 36000, 54000, 72000, 90000]
            xtick1 = ["", "", "", "", "", ""]
            ytick = [-5e6, 0, 5e6]
            ytick1 = ["-5000", "0", "5000"]
            xticks(xtick,xtick1)
            yticks(ytick,ytick1)
        if i==1:
            plot(sat_t,sat_y)
            ylabel('Y [Km]')
            xtick = [0, 18000, 36000, 54000, 72000, 90000]
            xtick1 = ["", "", "", "", "", ""]
            ytick = [-5e6, 0, 5e6]
            ytick1 = ["-5000", "0", "5000"]
            xticks(xtick,xtick1)
            yticks(ytick,ytick1)
        if i==2:
            plot(sat_t,sat_z,label='Real')
            ylabel('Z [Km]')
            xtick = [0, 18000, 36000, 54000, 72000, 90000]
            xtick1 = ["0", "5", "10", "15", "20", "25"]
            ytick = [-5e6, 0, 5e6]
            ytick1 = ["-5000", "0", "5000"]
            xticks(xtick,xtick1)
            yticks(ytick,ytick1)
            legend()
    
    show()

graficar_posicion_real_pred()   

def graficar_deriva_ode_real():
    
    figure()
    deriva = sqrt( (sat_x - x_ode)**2 + (sat_y - y_ode)**2 + (sat_z - z_ode)**2 )
    plot(sat_t/3600,deriva/1000)
    max_deriva=(deriva[-1])/1000
    title(f'Deriva Oderint v/s Real\ndmax= {max_deriva} Km')
    ylabel('Deriva,d [Km]')
    xlabel('Tiempo, t [horas]')
    grid(True)
    tight_layout()
    show() 

    
graficar_deriva_ode_real() 


# # Deriva real-ode 
# pos_final=zf-sol1[-1]

# norma_distancia_error = sqrt(pos_final[0]**2 + pos_final[1]**2 + pos_final[2]**2)
# print (norma_distancia_error)








    
