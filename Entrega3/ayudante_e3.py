# <TAI>TAI=2020-07-27T23:00:19.000000</TAI>
# <UTC>UTC=2020-07-27T22:59:42.000000</UTC>
# <UT1>UT1=2020-07-27T22:59:41.788119</UT1>
# <Absolute_Orbit>+22663</Absolute_Orbit>
# <X unit="m">-1977662.964695</X>
# <Y unit="m">4488480.565803</Y>
# <Z unit="m">5090705.333414</Z>
# <VX unit="m/s">-399.708484</VX>
# <VY unit="m/s">5614.890000</VY>
# <VZ unit="m/s">-5093.037904</VZ>
# <Quality>NOMINAL</Quality> 
      
# <TAI>TAI=2020-07-29T01:00:19.000000</TAI>
# <UTC>UTC=2020-07-29T00:59:42.000000</UTC>
# <UT1>UT1=2020-07-29T00:59:41.788931</UT1>
# <Absolute_Orbit>+22679</Absolute_Orbit>
# <X unit="m">-1732978.438394</X>
# <Y unit="m">-3123443.580625</Y>
# <Z unit="m">6098055.102419</Z>
# <VX unit="m/s">941.273635</VX>
# <VY unit="m/s">6592.256334</VY>
# <VZ unit="m/s">3635.825287</VZ>
# <Quality>NOMINAL</Quality>
from numpy import zeros,cos,sin,array,dot,sqrt,linspace,pi
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import datetime as dt

utc_EOF_format= "%Y-%m-%dT%H:%M:%S.%f"
t1=dt.datetime.strptime('2020-07-27T22:59:42.000000',utc_EOF_format)
t2=dt.datetime.strptime('2020-07-29T00:59:42.000000',utc_EOF_format)

intervalo=t2-t1
intervalo_en_segundos= intervalo.total_seconds()
# print(f'intervalo = {intervalo} s')
# print(f'intervalo_en_segundos = {intervalo_en_segundos} s') 

x_i=-1977662.964695
y_i=4488480.565803
z_i=5090705.333414

vx_i=-399.708484
vy_i=5614.890000
vz_i=-5093.037904

x_f=-1732978.438394
y_f=-3123443.580625
z_f=6098055.102419

vx_f=941.273635
vy_f=6592.256334
vz_f=3635.825287

hr=3600. #s
km= 10**3 #m
Radio= 6371.*km #km
Mt= 5.972e24 #kg
G= 6.67408e-11 #m3 kg-1 s-2
omega= 7.2921150e-5
H0= 700.*km

FgMax= G*Mt/Radio**2

def zpunto(z,t):
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

t=linspace(0,intervalo_en_segundos,9361)

#x0=Radio+H0

z0=array([x_i,y_i,z_i,vx_i,vy_i,vz_i])

sol=odeint(zpunto,z0,t)

x=sol[:,:]

pos_final=array([x_f,y_f,z_f,vx_f,vy_f,vz_f])-sol[-1]

# for el in pos_final:
#     print(el)
    
H=sqrt(x[:,0]**2+x[:,1]**2+x[:,2]**2)-Radio

def graficar():
    plt.figure()
    for i in range(3):
        plt.subplot(3,1,1+i)
        plt.grid(True)
        plt.plot(t/hr,x[:,i]) 
    
    plt.figure()
    plt.grid(True)
    plt.plot(t/hr,H/km)
    plt.axhline(80.,linestyle='--')
    plt.axhline(0.,linestyle='--')
    
    plt.figure()
    plt.grid(True)
    plt.plot(x[:,0],x[:,1])
    
    th=linspace(0,2*pi,400)
    
    plt.plot(Radio*cos(th),Radio*sin(th))
    plt.axis('equal')
    
    plt.figure()
    ax=plt.axes(projection='3d')
    ax.plot3D(x[:,0],x[:,1],x[:,2])
    plt.show()

# graficar()
    
norma_vector=sqrt(pos_final.dot(pos_final))
print (norma_vector)

