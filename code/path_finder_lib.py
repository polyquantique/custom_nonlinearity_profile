import numpy as np

def find_path(f,N,delta=0.00002,d=1.0,d0=0.0,k=2,Pini=0.0): 
    """
    Function to find the path given a specified nonlinearity profile
    
    Args: 
        f(function): nonlinearity profile
        N(int): number of grid points
        delta(float): discretization in space
        d(float): susceptibility
        d0(float): angle independent susceptibility
        k(int): order angle dependency
        Pini(float): initial z position
    
    Returns: 
        tuple(array, array, array, array):(x path, z path, effective susceptibility, propagation distance)  
    """
    x = np.zeros(N)      # \rho in the manuscript
    s = np.zeros(N)      # arc length
    z = np.zeros(N)
    deff = np.zeros(N)
    theta = np.zeros(N)
    
    # initial values
    z[0] = Pini
    z[1] = Pini
    s[0] = Pini
    s[1] = Pini
    
    for i in range(2, N):   
        if -np.pi/4 < theta[i-1] < np.pi/4:
            z[i] = z[i-1]+delta;
            s[i] = s[i-1] + np.sqrt((x[i-1]-x[i-2])**2+(z[i]-z[i-1])**2);
            deff[i]=f(s[i])
            theta[i] = (1.0/k)*np.arcsin(d*f(s[i])-d0);
            dx = np.tan(theta[i]);
            x[i] = x[i-1] + delta*dx      
        else:             
            x[i] = x[i-1] + delta;
            s[i] = s[i-1] + np.sqrt((x[i]-x[i-1])**2+(z[i-1]-z[i-2])**2);
            deff[i]=f(s[i])
            theta[i] = (1.0/k)*np.arcsin(d*f(s[i])-d0);
            dz = np.tan(theta[i]);
            z[i] = z[i-1] + np.sign(theta[i])*delta/dz  
            x[i] = x[i-1] + np.sign(theta[i])*delta;
    return x,z,deff,s
