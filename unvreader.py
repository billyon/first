import matplotlib.pyplot as plt
import numpy as np 
import matplotlib

font = {'family' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)


name = ['F_x','F_y','F_z']

for i in range(0,3):
    f = open('/home/billykon/Desktop/OMADA 12/'+name[i]+'.lvm','r')
    l = f.readlines()
    l = l[22:]
    l = [i.split() for i in l]
    if i==0:
        df = 0.56
        
    if i==2:
        df = 0.26
    if i==1:
        df = 0
    
    l = [0.1*float(i[1])+df for i in l]
    n = len(l)
    t = np.linspace(0,n/10,num=n)

    plt.figure(i+1)
    plt.plot(t,l,linewidth=0.5)
    plt.grid()
    plt.xlabel('t [s]')
    plt.ylabel('$ '+name[i]+' [ daN ] $')
    f.close()


f = [0.1 ,0.2, 0.31]
Fc1 = [126.58*0.1, 217.66*0.1, 297.33*0.1]
Fc2 = [69.53*0.1,127.03*0.1,188.58*0.1]
Fc3 = [27.81*0.1,49.3*0.1,76.01*0.1]

plt.figure(4)
plt.loglog(f,Fc1,'-x',label='$t = 0.5 [mm]$',linewidth=0.5)
plt.loglog(f,Fc2,'-x',label='$t = 0.3 [mm]$',linewidth=0.5)
plt.loglog(f,Fc3,'-x',label='$t = 0.1 [mm]$',linewidth=0.5)
plt.xlabel('$f [mm/rev]$')
plt.ylabel('$F_c [daN]$')
plt.legend()
plt.grid(True, which="both", ls="-")
plt.axes().set_aspect('equal')

plt.figure(5)
t = [0.5, 0.3, 0.1]
fc1 = [126.58*0.1,69.53*0.1,27.81*0.1]
fc2 = [217.66*0.1,127.03*0.1,49.3*0.1]
fc3 = [297.33*0.1,188.58*0.1,76.01*0.1]
plt.loglog(t,fc1,'-x',label='$f = 0.1 [mm/rev]$',linewidth=0.5)
plt.loglog(t,fc2,'-x',label='$f = 0.2 [mm/rev]$',linewidth=0.5)
plt.loglog(t,fc3,'-x',label='$f = 0.5 [mm/rev]$',linewidth=0.5)
plt.xlabel('$t [mm]$')
plt.ylabel('$F_c [daN]$')
plt.xlim([0.1,0.6])
plt.legend()
plt.grid(True, which="both", ls="-")
plt.axes().set_aspect('equal')


plt.figure(6)
t = [0.5, 0.3, 0.1]
fx1 = [ 7.55, 3.01, 0.71]
fx2 = [ 10.21,4.26, 0.72]
fx3 = [10.21, 4.92, 0.51]
plt.loglog(t,fx1,'-x',label='$f = 0.1 [mm/rev]$',linewidth=0.5)
plt.loglog(t,fx2,'-x',label='$f = 0.2 [mm/rev]$',linewidth=0.5)
plt.loglog(t,fx3,'-x',label='$f = 0.5 [mm/rev]$',linewidth=0.5)
plt.xlabel('$t [mm]$')
plt.ylabel('$F_x [daN]$')
#plt.xlim([0.1,0.6])
plt.legend()
plt.grid(True, which="both", ls="-")
plt.axes().set_aspect('equal')


f = [0.1 ,0.2, 0.31]
Fx1 = [7.55,10.21,10.21]
Fx2 = [3.01,4.26,4.92]
Fx3 = [0.71,0.72,0.51]

plt.figure(7)
plt.loglog(f,Fx1,'-x',label='$t = 0.5 [mm]$',linewidth=0.5)
plt.loglog(f,Fx2,'-x',label='$t = 0.3 [mm]$',linewidth=0.5)
plt.loglog(f,Fx3,'-x',label='$t = 0.1 [mm]$',linewidth=0.5)
plt.xlabel('$f [mm/rev]$')
plt.ylabel('$F_x [daN]$')
plt.legend()
plt.grid(True, which="both", ls="-")
plt.axes().set_aspect('equal')


plt.figure(8)
t = [0.5, 0.3, 0.1]
fy1 = [5.91,4.73,2.61]
fy2 = [11.03,8.82,4.20]
fy3 = [15.09,13.63,6.4]
plt.loglog(t,fy1,'-x',label='$f = 0.1 [mm/rev]$',linewidth=0.5)
plt.loglog(t,fy2,'-x',label='$f = 0.2 [mm/rev]$',linewidth=0.5)
plt.loglog(t,fy3,'-x',label='$f = 0.5 [mm/rev]$',linewidth=0.5)
plt.xlabel('$t [mm]$')
plt.ylabel('$F_y [daN]$')
#plt.xlim([0.1,0.6])
plt.legend()
plt.grid(True, which="both", ls="-")
plt.axes().set_aspect('equal')



f = [0.1 ,0.2, 0.31]
Fy1 = [5.91,11.03,15.09]
Fy2 = [4.73,8.82,13.63]
Fy3 = [2.61,4.2,6.4]

plt.figure(9)
plt.loglog(f,Fy1,'-x',label='$t = 0.5 [mm]$',linewidth=0.5)
plt.loglog(f,Fy2,'-x',label='$t = 0.3 [mm]$',linewidth=0.5)
plt.loglog(f,Fy3,'-x',label='$t = 0.1 [mm]$',linewidth=0.5)
plt.xlabel('$f [mm/rev]$')
plt.ylabel('$F_y [daN]$')
plt.legend()
plt.grid(True, which="both", ls="-")
plt.axes().set_aspect('equal')
plt.show()