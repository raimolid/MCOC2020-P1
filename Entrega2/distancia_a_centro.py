from numpy import array,zeros,cos,sin,linspace
from scipy.integrate import odeint
import matplotlib.pylab as plt

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
z=sol[:,2]

plt.style.use('dark_background')

th=linspace(0,2*3.14,100)

X=r_atm*cos(th)
Y=r_atm*sin(th)

plt.hlines(y = 6371000, xmin = 0, xmax=t,color='cyan')
plt.hlines(y = 7071000, xmin = 0, xmax=t, color='yellow')
plt.hlines(y = -6371000, xmin = 0, xmax=t, color='cyan')
plt.hlines(y = -7071000, xmin = 0, xmax=t, color='yellow')

plt.plot(T,x,color='lime')
plt.plot(T,y,color='red')

yTicks = [-7000*1000,-5000*1000,-3000*1000,0,3000*1000,
          5000*1000,7000*1000]
yTicks_Text = ['7000','5000','3000','0','3000','5000','7000']
    
xTicks = [0,2000,4000,6000,8000,10000,12000]
xTicks_Text = [0,2000,4000,6000,8000,10000,12000]

plt.yticks(yTicks, yTicks_Text,rotation=45)
plt.xticks(xTicks, xTicks_Text,rotation=45) 

plt.title("Distancia a Tierra")
plt.ylabel("Posicion (Km)")
plt.xlabel("Tiempo (s)")
plt.legend()  
plt.show()


