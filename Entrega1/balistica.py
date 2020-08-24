import scipy as sp
from scipy.integrate import odeint
import matplotlib.pylab as plt

#Unidades base:
    
cm=0.01   #m
inch=2.54*cm
g=9.81   # m/s**2

#Coeficiente de arrastre:

p=1.225   # kg/m**3
cd=0.47
D= 8.5*inch
r=D/2
A=sp.pi*r**2
CD=0.5*p*cd*A

#Masa
m=15  # Kg

#Viento
V=0   # m/s solo en x 

#Funcion a integrar:

# z es el vector de estado

# z= [x,y,vx,vy]

# dz/dt=bala(z,t)   dz1/dt=z2

#       [z2      ]
# dz/dt=[        ]
#       [FD/m  -g]

#vector de estado
#z[0]->x
#z[1]->y
#z[2]->vx
#z[3]->vy

def bala(z,t):
    zp=sp.zeros(4)
    zp[0]=z[2]
    zp[1]=z[3]
    v=z[2:4]
    v[0]=v[0]-V
    v2=sp.dot(v,v)
    vnorm=sp.sqrt(v2)
    FD=-CD*v2*(v/vnorm)
    zp[2]=FD[0]/m
    zp[3]=FD[1]/m-g
    return zp 


plt.figure()

V=0
t=sp.linspace(0,30,1001)
vi=100*1000/3600
z0=sp.array([0,0,vi,vi])
sol=odeint(bala,z0,t)
x=sol[:,0]
y=sol[:,1]
plt.plot(x,y,label='V = 0 m/s')
V=10
t=sp.linspace(0,30,1001)
vi=100*1000/3600
z0=sp.array([0,0,vi,vi])
sol=odeint(bala,z0,t)
x=sol[:,0]
y=sol[:,1]
plt.plot(x,y,label='V = 10.0 m/s')
V=20
t=sp.linspace(0,30,1001)
vi=100*1000/3600
z0=sp.array([0,0,vi,vi])
sol=odeint(bala,z0,t)
x=sol[:,0]
y=sol[:,1]
plt.plot(x,y,label='V = 20.0 m/s')

plt.ylim(0,50)
plt.xlim(0,150)
plt.title('Trayectoria para distintos vientos')
plt.ylabel('Y(m)')
plt.xlabel('X(m)')
plt.legend()
plt.grid()
plt.savefig('balistica.png')
plt.show()
   
    
    
    