import matplotlib.pyplot as plt
import numpy as np 
import matplotlib
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": [],           
    "font.sans-serif": [],  
})

f = ['FX_','FY_','FZ_']

plt.figure(1)
a = open('OMADA 12/'+f[0]+'1.lvm','r')
l = readlines(1)
l = l[-2:]
plt.plot(l)