#
#   PC 3  Ex 4   Eccentric Anomaly
#

#import numpy as np
from math  import sin, pi, cos
import matplotlib.pyplot as plt

period = 365.25635
omega  = 2.0 * pi / period
k      = 1000.0
a      = 1.496e6*k
e      = 0.0167
tp     = 1.0
t      = [91.0*tp , 182.0*tp , 273.0*tp]

def f(E,time,ep):
    y = E - omega * time - ep*sin(E)    
    return y

def fp(E,time,ep):
    yp = 1.0 - ep*cos(E)    
    return yp

print ('\n Part a')
print ('---------------------')

for j in range(len(t)):
    i = 0
    E0 = 10.0
    E1 = 1.0
    while ( (abs((E1 - E0)/(E1+10e-30)) > 10e-10) and (i < 100) ):
        E0 = E1
        i = i + 1
        E1 = E0 - f(E0,t[j],e)/fp(E0,t[j],e)
        plt.figure(j)
        plt.semilogy(i,abs((E1 - E0)/(E1+10e-30)),'ob')
        plt.xlabel('n')
        plt.ylabel('Relative Error')
    print ('time = ', t[j], ' E = ', E1, ' i = ', i)
    plt.savefig('plotI%i.pdf' %(j))

#part b 

print ('\n Part b')
print ('---------------------')

e = 0.9

for j in range(len(t)):
    i = 0
    E0 = 10.0
    E1 = 1.0
    while ( (abs((E1 - E0)/(E1+10e-30)) > 10e-10) and (i < 1000) ):
        E0 = E1
        i = i + 1
        E1 = E0 - f(E0,t[j],e)/fp(E0,t[j],e)
        plt.figure(j+3)
        plt.semilogy(i,abs((E1 - E0)/(E1+10e-30)),'ob')
        plt.xlabel('n')
        plt.ylabel('Relative Error')
    print ('time = ', t[j], ' E = ', E1, ' i = ', i)
    plt.savefig('plotII%i.pdf' %(j))

