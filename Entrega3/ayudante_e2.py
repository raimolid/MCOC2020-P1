from numpy import zeros,cos,sin,array,dot,sqrt,linspace,pi
import matplotlib.pyplot as plt
from scipy.integrate import odeint

hr=3600. #s
km= 10**3 #m
Radio= 6371.*km #km
Mt= 5.972e24 #kg
G= 6.67408e-11 #m3 kg-1 s-2
omega= 7.2921150e-5
H0= 700.*km

FgMax= G*Mt/Radio**2

zp=zeros(6)
def zpunto(z,t):
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

t=linspace(0,5.*hr,1001)

x0= Radio + H0

vt=8000
z0=array([x0,0,0,0,vt,0])

sol=odeint(zpunto,z0,t)

x=sol[:,0:3]

H=sqrt(x[:,0]**2+x[:,1]**2+x[:,2]**2)-Radio

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
plt.show()

    
    