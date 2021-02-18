'''
This is a library for coordinate frame transformation functions 
it contains transformations from cartesian to spherical coordinate system and vice versa
'''
import numpy as np
import math as m



def cartesian_to_spherical(A):  #transforms given points from Cartesian to spherical coordinate system
    
    x = A[:,0]
    y = A[:,1]
    z = A[:,2]
    i = x.shape[0]
    r = []
    psi = []
    theta = []
    for a in range(i):
        r.append(m.sqrt(x[a,0]**2 + y[a,0]**2 + z[a,0]**2))
        psi.append(m.atan(y[a,0]/x[a,0]))
        theta.append(m.atan((m.sqrt(x[a,0]**2 + y[a,0]**2)/z[a,0])))
     
    r = np.matrix(r).T
    psi = np.matrix(psi).T
    theta = np.matrix(theta).T           
    return np.concatenate((r,psi,theta), axis = 1)


def spherical_to_cartesian(A):  #transforms given points from Spherical to cartesian coordinate system
     
    r = A[:,0]
    psi = A[:,1]
    theta = A[:,2]
    i = r.shape[0]
    x = []
    y = []
    z = []
    for a in range(i):
        x.append(r[a,0]*m.sin(theta[a,0])*m.cos(psi[a,0]))
        y.append(r[a,0]*m.sin(theta[a,0])*m.sin(psi[a,0]))
        z.append(r[a,0]*m.cos(theta[a,0]))
    
    x = np.matrix(x).T
    y = np.matrix(y).T
    z = np.matrix(z).T
    return np.concatenate((x,y,z), axis = 1)
