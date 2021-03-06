from numpy import zeros,cos,sin,array,dot,sqrt 
from scipy.integrate import odeint
from matplotlib.pylab import *
from leer_eof import leer_eof
from time import perf_counter

# Leer de línea de comando
from sys import argv
eofname=argv[1]

# eofname='S1B_OPER_AUX_POEORB_OPOD_20200817T110752_V20200727T225942_20200729T005942.EOF'
sat_t,sat_x,sat_y,sat_z,sat_vx,sat_vy,sat_vz=leer_eof(eofname)

eof_out=eofname.replace('.EOF','.PRED')
# print(f'Archivo de entrada: {eofname}')
# print(f'Archivo de salida: {eof_out}')


# Condición Inicial
z0=array([sat_x[0],sat_y[0],sat_z[0],sat_vx[0],sat_vy[0],sat_vz[0]])

# Condición Final
zf=array([sat_x[-1],sat_y[-1],sat_z[-1],sat_vx[-1],sat_vy[-1],sat_vz[-1]])


# Datos
Mt=5.972e24 #kg 
G=6.67408e-11 #m3 kg-1 s-2
omega=7.2921150e-5 #rad s-1

J2=1.75553e10*(1000)**5 #m5 s-2
J3=-2.61913e11*(1000)**5 #m6 s-2

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
    
    # Agregando corrección J2
    Fx=J2*(z1[0]/r**7)*((6.*z1[2]**2)-1.5*(z1[0]**2+z1[1]**2))
    Fy=J2*(z1[1]/r**7)*((6.*z1[2]**2)-1.5*(z1[0]**2+z1[1]**2))
    Fz=J2*(z1[2]/r**7)*((3.*z1[2]**2)-4.5*(z1[0]**2+z1[1]**2))
    F2 = array([Fx,Fy,Fz])
    
    # Agregando corrección J3 
    Fx=J3*(z1[0]*z1[2]/r**9)*((10.*z1[2]**2)-7.5*(z1[0]**2+z1[1]**2))
    Fy=J2*(z1[1]*z1[2]/r**9)*((10.*z1[2]**2)-7.5*(z1[0]**2+z1[1]**2))
    Fz=J2*(1/r**9)*((4.*z1[2]**2)*(z1[2]**2-3.*(z1[0]**2+z1[1]**2))+1.5*(z1[0]**2+z1[1]**2))
    F3 = array([Fx,Fy,Fz])
    
    
    # Vector de salida 
    zp[0:3] = z2 
    zp[3:6]= (-G*Mt/r**3)*z1  - R.T@(Rpp@z1 + 2.*Rp@z2) +F2+F3
    
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
# print(f'Tiempo Odeint: {dt}\n')

sol_ode=sol1[:,:]
x_ode=sol_ode[:,0]
y_ode=sol_ode[:,1]
z_ode=sol_ode[:,2]
vx_ode=sol_ode[:,3]
vy_ode=sol_ode[:,4]
vz_ode=sol_ode[:,5] 

t3=perf_counter()
sol2=eulerint(zp,z0,sat_t,1)
t4=perf_counter()
dt=t4-t3
# print (f'Tiempo Eulerint: {dt}\n')

sol_euler=sol2[:,:]
x_euler=sol_euler[:,0]
y_euler=sol_euler[:,1]
z_euler=sol_euler[:,2]

def graficar_posicion_real_pred():
    # 3 Gráficos tiempo vs posición x,y,z
    figure()
    for i in range(3):
        
        subplot(3,1,1+i)
        grid(True)
        plot(sat_t,sol_ode[:,i],label='Odeint')
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

# graficar_posicion_real_pred()

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

    
# graficar_deriva_ode_real() 
 

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
    
# graficar_deriva_ode_euler()

def graficar_deriva_euler_real():
    
    figure()

    deriva = sqrt( (sat_x - x_euler)**2 + (sat_y - y_euler)**2 + (sat_z - z_euler)**2 )
    plot(sat_t/3600,deriva/1000)
    max_deriva=(deriva[-1])/1000
    title(f'Deriva Eulerint v/s Real\ndmax= {max_deriva} Km')
    ylabel('Deriva,d [Km]')
    xlabel('Tiempo, t [horas]')
    grid(True)
    tight_layout()
    show() 

    
# graficar_deriva_euler_real()


# # Deriva real-ode con J2J3
# pos_final=zf-sol1[-1]

# norma_distancia_error = sqrt(pos_final[0]**2 + pos_final[1]**2 + pos_final[2]**2)
# print (norma_distancia_error)

# # Deriva real-euler con J2J3
# pos_final=zf-sol2[-1]

# norma_distancia_error = sqrt(pos_final[0]**2 + pos_final[1]**2 + pos_final[2]**2)
# print (norma_distancia_error)


# Creando Archivo .PRED

with open(eof_out,'w') as fout:
    fout.write('<?xml version="1.0" ?>\n'
'<Earth_Explorer_File>\n'
'  <Earth_Explorer_Header>\n'
'    <Fixed_Header>\n'
'      <File_Name>S1B_OPER_AUX_POEORB_OPOD_20200817T110752_V20200727T225942_20200729T005942</File_Name>\n'
'      <File_Description>Precise Orbit Ephemerides (POE) Orbit File</File_Description>\n'
'      <Notes></Notes>\n'
'      <Mission>Sentinel-1B</Mission>\n'
'      <File_Class>OPER</File_Class>\n'
'      <File_Type>AUX_POEORB</File_Type>\n'
'      <Validity_Period>\n'
'        <Validity_Start>UTC=2020-07-27T22:59:42</Validity_Start>\n'
'        <Validity_Stop>UTC=2020-07-29T00:59:42</Validity_Stop>\n'
'      </Validity_Period>\n'
'      <File_Version>0001</File_Version>\n'
'      <Source>\n'
'        <System>OPOD</System>\n'
'        <Creator>OPOD</Creator>\n'
'        <Creator_Version>0.0</Creator_Version>\n'
'        <Creation_Date>UTC=2020-08-17T11:07:52</Creation_Date>\n'
'      </Source>\n'
'    </Fixed_Header>\n'
'    <Variable_Header>\n'
'      <Ref_Frame>EARTH_FIXED</Ref_Frame>\n'
'      <Time_Reference>UTC</Time_Reference>\n'
'    </Variable_Header>\n'
'  </Earth_Explorer_Header>\n'
'<Data_Block type="xml">\n'
'  <List_of_OSVs count="9361">\n')
    Nt=len(sat_t)  
    for i in range(Nt):
        fout.write('    <OSV>\n'
f'      <TAI>TAI=2020-07-27T23:00:19.000000</TAI>\n'
f'      <UTC>UTC=2020-07-27T22:59:42.000000</UTC>\n'
f'      <UT1>UT1=2020-07-27T22:59:41.788119</UT1>\n'
f'      <Absolute_Orbit>+22663</Absolute_Orbit>\n'
f'      <X unit="m">{x_ode[i]}</X>\n'
f'      <Y unit="m">{y_ode[i]}</Y>\n'
f'      <Z unit="m">{z_ode[i]}</Z>\n'
f'      <VX unit="m/s">{vx_ode[i]}</VX>\n'
f'      <VY unit="m/s">{vy_ode[i]}</VY>\n'
f'      <VZ unit="m/s">{vz_ode[i]}</VZ>\n'
'      <Quality>NOMINAL</Quality>\n'
'    </OSV>\n')
    fout.write('  </List_of_OSVs>\n'
'</Data_Block>\n'
'</Earth_Explorer_File>\n')
