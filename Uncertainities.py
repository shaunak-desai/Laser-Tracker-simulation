'''
This is a library to calculate the uncertainies of Laser Tracker based on the distance and 
the unceretainities arising from the window imperfections
'''

import numpy as np
import math as m



#Defining Laser Tracker origin
LT_origin = np.matrix([[0], [0], [0]])


#Calculating distance of each target point from the Laser tracker origin
def dist(A):
    i = A.shape[0]
    distance = []
    
    for a in range(i):
        distance.append(m.sqrt(A[a,0]**2 + A[a,1]**2 + A[a,2]**2))
        
    return distance

#Calculating Laser trcker uncertainities
def LT_uncertainities(A, B):
    i = A.shape[0]
    
    distance = []
    
    for a in range(i):
        distance.append(m.sqrt(A[a,0]**2 + A[a,1]**2 + A[a,2]**2))

    distance = np.matrix(distance).T
    
    r_err = []
    ang_err = []
    
    if B == 0:
        for a in range(i):
            r_err.append((15e-6)/6 + ((6e-6)/6 * distance[a,0]))
            ang_err.append(m.atan(((15e-6)/6 + ((6e-6)/6 * distance[a,0])) / distance[a,0]))

    elif B == 1:
         for a in range(i):
            r_err.append(np.random.normal(0, ((15e-6)/6 + ((6e-6)/6 * distance[a,0]))))
            ang_err.append(np.random.normal(0, (m.atan(((15e-6)/6 + ((6e-6)/6 * distance[a,0])) / distance[a,0]))))
            
    r_err = np.matrix(r_err)
    ang_err = np.matrix(ang_err)        
    
    return np.concatenate((r_err, ang_err, ang_err), axis = 0)

    

def Window_uncertainities(A,B,C,D):
    i = A.shape[0]
    r_err = np.zeros((i,1))
    ang_err_theta = np.zeros((i,1))
    ang_err_phi = np.zeros((i,1))
    
    n_window = 1.51509
    n_air = 1.000273
    n_vacuum = 1
    
    C = C / 206264.806  #converts input incident arcseconds into radians
    D = D / 206264.806  #converts input incident arcseconds into radians
    
    transmit_1 = m.asin((n_window/ n_vacuum)* m.sin(C))    #calculates transmit angel based on Snell's law
    transmit_2 = m.asin((n_window/ n_vacuum)* m.sin(D))    #calculates transmit angel based on Snell's law
    
    for a in range(i):     
        r_err[a,0] = np.random.uniform(-B,B)
        ang_err_theta[a,0] = np.random.uniform(-transmit_1,transmit_1)
        ang_err_phi[a,0] = np.random.uniform(-transmit_2,transmit_2)
        
    
    return np.concatenate((r_err, ang_err_theta, ang_err_phi), axis = 1)

def Additional_ang_uncertatinities(A, B, C):
    i = A.shape[0]
    add_err = np.zeros((i,3))
    for a in range(i):
        add_err[a,1] = np.random.uniform(-B,B)
        add_err[a,2] = np.random.uniform(-C,C)
    return add_err
    

def Mechanical_uncertainities(A, B, C, D):
    i = A.shape[0]
    
    x_err = np.zeros((i,1))
    y_err = np.zeros((i,1))
    z_err = np.zeros((i,1))
    
    for a in range(i):
        x_err[a,0] = np.random.normal(0,B)
        y_err[a,0] = np.random.normal(0,C)
        z_err[a,0] = np.random.normal(0,D)
    
    return np.concatenate((x_err, y_err, z_err), axis = 1)


def vaccum_pressure(A, B, C, D):
    i = A.shape[0]
    
    slope = 2.724025e-7 #Parameteres for Pressure VS Index line
    intercept = 0.999999573
    
    n2 = (slope * D) + intercept
    
    theta = m.asin(m.sin(B)* (C/n2)) #Angular shift from change in index
    
    range_err = np.zeros((i,1))
    azimuth_err = np.zeros((i,1))
    elevation_err = np.zeros((i,1))
    
    for a in range(i):
        azimuth_err[a,0] = theta
        elevation_err[a,0] = theta
    
    return np.concatenate((range_err, azimuth_err, elevation_err), axis = 1)
    
    
    


    
    
    
    