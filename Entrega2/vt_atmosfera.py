from numpy import array,zeros,cos,sin,linspace
from scipy.integrate import odeint
import matplotlib.pylab as plt
import matplotlib.patches 

#Constante Gravitacional
G=6.67*10**-11  #Nm2/kg2

#Masa de la Tierra
M=5.972*10**24  #kg

#Velocidad angular de la Tierra
o=7.27*10**-5   #rad/s

#Radio Satelite
r_sat=7071000 #m

#Radio Atmósfera
r_atm=6451000 #m


def satelite(z,t):

    #Creando vector zp y guardo las velocidades
    zp=zeros(6)
    zp[0:3]=z[3:6]
    
    #Guardo los valores para ocuparlos en el calculo
    x=z[0:3] 
    xp=z[3:6]
    
    #Multiplicador que uso en z2p
    C=-(G*M)/(r_sat**3)
     
    # ot=10 #rad 
    
    #Matriz de Rotacion                                 
    R=array([[cos(o*t),-sin(o*t),0],
             [sin(o*t),cos(o*t),0],
             [0,0,1]])
    
    #Matriz de Rotacion derivada
    Rp=o*array([[-sin(o*t),-cos(o*t),0],
                [cos(o*t),-sin(o*t),0],
                [0,0,0]])
    
    #Matriz de Rotacion doble derivada
    Rpp=(o**2)*array([[-cos(o*t),sin(o*t),0],
                      [-sin(o*t),-cos(o*t),0],
                      [0,0,0]])
    
    #Calculo de z2p segun formula grande
    z2p=C*x-(R.T@(Rpp@x + 2*(Rp@xp)))
   
    #Guardo los resultados del calculo en los valores q faltaban del vector zp
    zp[3:6]=z2p[0:3]
    
    return zp


#Velocidad en y
vy=7000   #m/s

#Vector z(posicion,velocidad)
z0=array([7071000,0,0,0,vy,0])  #x,y,x,vx,vy,vz 

#Tiempo a Analizar
t=2*1.765*3600    #1.765 horas en s
T = linspace(0,t,1001)

#Posición
sol=odeint(satelite,z0,T)

x=sol[:,0]
y=sol[:,1]

plt.style.use('dark_background')

plt.plot(x,y,color='lime',label=f'V= {vy} m/s ')

ax = plt.gca()
ellipse = matplotlib.patches.Ellipse((0,0), 2*r_atm, 2*r_atm,fill=False,
                                      color='r',label='Atmósfera')
ax.add_patch(ellipse)

yTicks = [-7000*1000,-5000*1000,-3000*1000,0,3000*1000,
          5000*1000,7000*1000]
yTicks_Text = ['7000','5000','3000','0','3000','5000','7000']
    
xTicks = [-7000*1000,-5000*1000,-3000*1000,0,3000*1000,
          5000*1000,7000*1000]
xTicks_Text = ['7000','5000','3000','0','3000','5000','7000']

plt.yticks(yTicks, yTicks_Text,rotation=45)
plt.xticks(xTicks, xTicks_Text,rotation=45) 

plt.title('Trayectoria Satélite')
plt.ylabel('Y(Km)')
plt.xlabel('X(Km)')
plt.legend()
plt.grid()
plt.show()

       

