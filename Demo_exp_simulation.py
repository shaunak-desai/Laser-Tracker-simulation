'''
Simulation for demo experiment with 4 targtes placed at an approrpirate distance to simulate the 4 corners of the 
ATHENA WFI CCD detector
'''

import numpy as np
import transformations as t
import transform 
import matplotlib.pyplot as plt
import Uncertainities as un
import math as m
from Kabsch_best_fit import best_fit

'''
###############################################################################
'''

def demo_err(window_thickness, window_centring_theta, window_centring_phi, distance, mechanical_un, pressure, glass_index,iterations):
    #Initilizing the 6-DOF arrays
    X = []
    Y = []
    Z = []
    Roll = []
    Pitch = []
    Yaw = []
    
    scaling = False
    
    for i in range(iterations):
        
        #Target coordinates in the demo_CCD cordinate frame
        CCD_target_1 = np.matrix([[-0.070],[0.070],[0],[1]])
        CCD_target_2 = np.matrix([[0.070],[0.070],[0],[1]])
        CCD_target_3 = np.matrix([[-0.070],[-0.070],[0],[1]])
        CCD_target_4 = np.matrix([[0.070],[-0.070],[0],[1]])
        
        targets_demo_CCD_frame = np.concatenate((CCD_target_1, CCD_target_2, CCD_target_3, CCD_target_4), axis = 1)
        
        #transformation matrix to convert target coordinates from demo_CCD frame to LT frame
        demo_CCD_to_LT = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,-distance],[0,0,0,1]])
        
        #Conversin of target coordinats to LT frame
        Homogeneous_modelled_points = (demo_CCD_to_LT * targets_demo_CCD_frame).T
        
        #removing the fourth dimension for ease of calculations
        modelled_points = np.delete(Homogeneous_modelled_points, 3, 1)
        
        #Various sources of errors
        LT_err = un.LT_uncertainities(modelled_points, 1).T
        Window_err = un.Window_uncertainities(modelled_points, window_thickness, window_centring_theta, window_centring_phi)
        #print (Window_err)
        
        #Vaccum_index_err = un.vaccum_pressure(modelled_points, window_centring, glass_index, pressure)
        #print(Vaccum_index_err)
        Mechanical_err = un.Mechanical_uncertainities(modelled_points, mechanical_un, mechanical_un, mechanical_un)
        
        
        modelled_points_1 = modelled_points 
        
        #Adding mechanical uncertainities of mounting the targets
        modelled_points_1 = modelled_points + Mechanical_err

        #Converting the modelled points into Laser tracker's spherical coordinate system
        spherical_points = transform.cartesian_to_spherical(modelled_points_1)

        #Adding the uncertainities from the Laser Tracker
        spherical_points = spherical_points + LT_err

        #Adding uncertainities from the window
        spherical_points = spherical_points + Window_err 
        
        #Adding angular shifts due to vaccum from Snell's law
        #spherical_points = spherical_points + Vaccum_index_err
        
        #Converting back to the cartesian coordinate system
        cartesian_points = transform.spherical_to_cartesian(spherical_points)
        
        Homogeneous_transform = t.superimposition_matrix(modelled_points, cartesian_points, usesvd=True)
        Rotation = np.matrix(Homogeneous_transform[0:3, 0:3])
        Translation = np.matrix(Homogeneous_transform[0:3,3]).T
        
        #[Rotation, Translation] = best_fit(cartesian_points, modelled_points)

#calculating homogeneous transformation function 
        #Homogeneous_transform = np.zeros((4,4))
        #Homogeneous_transform[:3,:3] = Rotation
        #Homogeneous_transform[0,3] = Translation[0]
        #Homogeneous_transform[1,3] = Translation[1]
        #Homogeneous_transform[2,3] = Translation[2]
        #Homogeneous_transform[3,3] = 1
        
        
        #Finding the euler angles from the homogeneous transformation matrix
        euler_angles = np.matrix(t.euler_from_matrix(Homogeneous_transform)).T
        
        for i in range(3):
            euler_angles[i,0] = m.degrees(euler_angles[i,0]) * 3600

        #Appending strings for 6-DOF values in each iteration
        X.append((Translation[0,0])*1e6)
        Y.append((Translation[1,0])*1e6)
        Z.append((Translation[2,0])*1e6)
        Roll.append(euler_angles[0,0])
        Pitch.append(euler_angles[1,0])
        Yaw.append(euler_angles[2,0])
        
    #calculating the standard deviation for the 6-DOF values
    X = np.std(X)
    Y = np.std(Y)
    Z = np.std(Z)
    Roll = np.std(Roll)
    Pitch = np.std(Pitch)
    Yaw = np.std(Yaw)
    
    return X, Y, Z, Roll, Pitch, Yaw


'''
###############################################################################
'''
X = []
Y = []
Z = []
Roll = []
Pitch = []
Yaw = []

# changing LT distance
       
for LT_distance in np.linspace(1.7, 2.9, num = 5):
    [x, y, z, roll, pitch, yaw] = demo_err(4e-5, 135, 135, LT_distance, 0, 0.01, 1.51289,10000)
    
    X.append(x)
    Y.append(y)
    Z.append(z)
    Roll.append(roll)
    Pitch.append(pitch)
    Yaw.append(yaw)

'''
#changing pressure in vaccum
for pressure in np.linspace(1e-6, 0.01, num = 100):
    [x, y, z, roll, pitch, yaw] = demo_err(0, 0.00016968, 2, 0, pressure, 1.51289, 100)
  
  
    X.append(x)
    Y.append(y)
    Z.append(z)
    Roll.append(roll)
    Pitch.append(pitch)
    Yaw.append(yaw)

'''  
    
'''
###############################################################################
'''

distance = np.linspace(1.7, 2.9, num =5)
    
plt.figure()
plt.plot(distance, X, 'r')
plt.plot(distance, Y, 'b')
plt.plot(distance, Z, 'g')
plt.legend(['X-err', 'Y-err', 'Z-err'], loc = 2)
plt.xlabel('Distance of LT from the Target plane (meters)')
plt.ylabel('Errors ($\mu$m)') 
plt.title('Translational error propagation with increasing distance')
plt.grid() 
plt.show  
    

plt.figure()
plt.plot(distance, Roll, 'r')
plt.plot(distance, Pitch, 'b')
plt.plot(distance, Yaw, 'g')
plt.legend(['Roll', 'Pitch', 'Yaw'], loc =2)
plt.xlabel('Distance of LT from the Target plane (meters)')
plt.ylabel('Errors (arcecond)') 
plt.title('Rotational error propagation with increasing distance')
plt.grid() 
plt.show        


'''
###############################################################################
'''
'''
pressure = np.linspace(1e-6, 0.01, num = 100)
    
plt.figure()
plt.plot(pressure, X, 'r')
plt.plot(pressure, Y, 'b')
plt.plot(pressure, Z, 'g')
plt.legend(['X-err', 'Y-err', 'Z-err'], loc = 2)
plt.xlabel('Prressure in vaccum chamber (mBar)')
plt.ylabel('Errors ($\mu$m)') 
plt.title('Translational error')
plt.grid() 
plt.show  
    

plt.figure()
plt.plot(pressure, Roll, 'r')
plt.plot(pressure, Pitch, 'b')
plt.plot(pressure, Yaw, 'g')
plt.legend(['Roll', 'Pitch', 'Yaw'], loc =2)
plt.xlabel('Prressure in vaccum chamber (mBar)')
plt.ylabel('Errors (arcecond)') 
plt.title('Rotational error')
plt.grid() 
plt.show
'''



        