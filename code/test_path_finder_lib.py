import numpy as np
import matplotlib.pyplot as plt

def analytical(x,x0,y0,b,radius):
    if b:
        return  -np.sqrt(radius**2-(x-x0)**2)+y0
    else:
        return   np.sqrt(radius**2-(x-x0)**2)+y0

def test_sine(z, N=30000, delta=0.00002, Lambda=0.1, k=1):
    radius=Lambda*k/2/np.pi
    z_test = np.zeros(N)    
    x_test = np.zeros(N) 
    z_test[0]=0
    x_test[0]=0

    for i in range(1, N):
        z_test[i] = z[i-1] 
        if np.mod((z_test[i]-radius),4*radius) < 2*radius:
            b=0
        else:
            b=1
        x_test[i]=analytical(z_test[i],int((z_test[i]+radius)/(2*radius))*2*radius,radius,b,radius)

    return z_test, x_test