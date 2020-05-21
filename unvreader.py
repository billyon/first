import matplotlib.pyplot as plt
import numpy as np 

name = ['F_x','F_y','F_z']

for i in range(0,3):
    f = open('/home/billykon/Desktop/OMADA 12/'+name[i]+'.lvm','r')
    l = f.readlines()
    l = l[22:]
    l = [i.split() for i in l]
    l = [float(i[1]) for i in l]
    n = len(l)
    t = np.linspace(0,n/10,num=n)

    plt.figure(i+1)
    plt.plot(t,l,linewidth=0.5)
    plt.xlabel('t [s]')
    plt.ylabel('$ '+name[i]+' [ N ] $')
    f.close()


f = [0.1 ,0.2, 0.31]
F1 = [126.58, 217.66, 297.33]
F2 = [69.53,127.03,188.58]
F3 = [27.81,49.3,76.01]

plt.figure(4)
plt.loglog(f,F1,'-x',label='$t = 0.5 [mm]$',linewidth=0.5)
plt.loglog(f,F2,'-x',label='$t = 0.3 [mm]$',linewidth=0.5)
plt.loglog(f,F3,'-x',label='$t = 0.1 [mm]$',linewidth=0.5)
plt.xlabel('$f [mm/rev]$')
plt.ylabel('$F_c [N]$')
plt.legend()
plt.axes().set_aspect('equal')

plt.figure(5)
t = [0.5, 0.3, 0.1]
f1 = [126.58,69.53,27.81]
f2 = [217.66,127.03,49.3]
f3 = [297.33,188.58,76.01]
plt.loglog(t,f1,'-x',label='$f = 0.1 [mm/rev]$',linewidth=0.5)
plt.loglog(t,f2,'-x',label='$f = 0.2 [mm/rev]$',linewidth=0.5)
plt.loglog(t,f3,'-x',label='$f = 0.5 [mm/rev]$',linewidth=0.5)
plt.xlabel('$t [mm]$')
plt.ylabel('$F_c [N]$')
plt.legend()
plt.axes().set_aspect('equal')
plt.show()