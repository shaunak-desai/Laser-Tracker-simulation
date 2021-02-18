#Alternate approach/simpler code


import numpy as np
import math as m
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from transform import cartesian_to_spherical
from transform import spherical_to_cartesian
from Kabsch_best_fit import best_fit
import Uncertainities as un
import transformations as t

##def calculate(parameter):
scaling = False

####
#test = 0.00002424 * 2
####

'''
#########################################################################################################################################
'''

def MA_error(window_angular_error, window_range_error, angular_error, mechanical_error, distance, iterations):

#Defining the radius of the mirror
    r = 1.3

#defining the LT origin
    LT_origin = np.matrix([[0],[0],[0]])

#Initializing the arrays for the 6-DOF values    
    X_MA = []
    Y_MA = []
    Z_MA = []
    Roll_MA = []
    Pitch_MA = []
    Yaw_MA = []
    MA_centroid = np.zeros((3,1)) 
    
    for i in range(iterations):

#definfing target points in MACS
        target_1 = np.matrix([[r*m.sin(m.radians(60))],[r*m.cos(m.radians(60))],[0],[1]])
        target_2 = np.matrix([[r*m.sin(m.radians(120))],[r*m.cos(m.radians(120))],[0],[1]])
        target_3 = np.matrix([[r*m.sin(m.radians(180))],[r*m.cos(m.radians(180))],[0],[1]])
        target_4 = np.matrix([[r*m.sin(m.radians(240))],[r*m.cos(m.radians(240))],[0],[1]])
        target_5 = np.matrix([[r*m.sin(m.radians(300))],[r*m.cos(m.radians(300))],[0],[1]])
        target_6 = np.matrix([[r*m.sin(m.radians(0))],[r*m.cos(m.radians(0))],[0],[1]])

#Transformation to convert targets from MA frame to Laser Tracker frame
        MA_to_LT = np.matrix([[1, 0, 0, 1.3],[0,1,0,1.2],[0,0,1,distance],[0,0,0,1]])
        
        Targets_MA_frame = np.concatenate((target_1, target_2, target_3, target_4, target_5, target_6), axis = 1)

#Converting from MA to Laser tracker coordinate frame
        Homogeneous_modelled_points = (MA_to_LT * Targets_MA_frame).T

#Removing fourth dimension for ease of calculation
        modelled_points = np.delete(Homogeneous_modelled_points, 3, 1)

#Introducing Laser Tracker's uncertainities
        LT_err_MA = un.LT_uncertainities(modelled_points, 1).T
        Window_err_MA = un.Window_uncertainities(modelled_points, window_range_error, window_angular_error, window_angular_error)
        Additional_err = un.Additional_ang_uncertatinities(modelled_points, angular_error, angular_error)
        Mechanical_err = un.Mechanical_uncertainities(modelled_points, mechanical_error, mechanical_error, mechanical_error)

#Adding mechanical uncertainities of mounting the targets
        modelled_points_1 = modelled_points + Mechanical_err

#Converting the modelled points into Laser tracker's spherical coordinate system
        spherical_points = cartesian_to_spherical(modelled_points_1)

#Adding the uncertainities from the Laser Tracker
        #spherical_points = spherical_points + LT_err_MA

#Adding uncertainities from the window
        spherical_points = spherical_points + Window_err_MA 

#Adding additional angular uncertainities
        spherical_points = spherical_points + Additional_err

#Converting back to the cartesian coordinate system
        cartesian_points = -spherical_to_cartesian(spherical_points)
        
        Homogeneous_transform = t.superimposition_matrix(modelled_points, cartesian_points, usesvd=True)
        Rotation = np.matrix(Homogeneous_transform[0:3, 0:3])
        Translation = np.matrix(Homogeneous_transform[0:3,3]).T

        #[Rotation, Translation] = best_fit(cartesian_points, modelled_points, scaling)
        
        n = modelled_points.shape[0]
    
        X = modelled_points.T     
        fit_points = (Rotation * modelled_points.T) + np.tile(Translation, (1, n))
        
#Finding centroid of the fitted point cloud
        fit_point_centroid = np.mean(fit_points, axis = 1)
        
        MA_centroid = np.append(MA_centroid, fit_point_centroid, 1)
        

#Finding centroid of the modelled point cloud
        modelled_points_centroid = np.mean(modelled_points, axis = 0).T

#calculating homogeneous transformation function 
        #Homogeneous_transform = np.zeros((4,4))
        #Homogeneous_transform[:3,:3] = Rotation
        #Homogeneous_transform[0,3] = Translation[0]
        #Homogeneous_transform[1,3] = Translation[1]
        #Homogeneous_transform[2,3] = Translation[2]
        #Homogeneous_transform[3,3] = 1
        
        #Finding the euler angles from the homogeneous transformation matrix
        euler_angles_MA = np.matrix(t.euler_from_matrix(Homogeneous_transform)).T
        for i in range(3):
            euler_angles_MA[i,0] = m.degrees(euler_angles_MA[i,0]) * 3600

#Appending strings for 6-DOF values in each iteration
        X_MA.append((Translation[0,0])*1e6)
        Y_MA.append((Translation[1,0])*1e6)
        Z_MA.append((Translation[2,0])*1e6)
        Roll_MA.append(euler_angles_MA[0,0])
        Pitch_MA.append(euler_angles_MA[1,0])
        Yaw_MA.append(euler_angles_MA[2,0])
        
#calculating STD of each 6-DOF strings to find the average over the number of iterations
    X_MA = np.std(X_MA)
    Y_MA = np.std(Y_MA)
    Z_MA = np.std(Z_MA)
    Roll_MA = np.std(Roll_MA)
    Pitch_MA = np.std(Pitch_MA)
    Yaw_MA = np.std(Yaw_MA)
    
    MA_centroid = np.delete(MA_centroid, 0, 1)
    MA_centroid_STD = np.std(MA_centroid, axis = 1)
    MA_centroid = np.mean(MA_centroid, axis = 1)

    return X_MA, Y_MA, Z_MA, Roll_MA, Pitch_MA, Yaw_MA, MA_centroid, MA_centroid_STD








'''
#########################################################################################################################################
'''

def MA_LT_error(distance, iterations):
    
    #Defining the radius of the mirror
    r = 1.3 

#defining the LT origin
    LT_origin = np.matrix([[0],[0],[0]])
    
#Initializing the arrays for the 6-DOF values    
    LT_X_MA = []
    LT_Y_MA = []
    LT_Z_MA = []
    LT_Roll_MA = []
    LT_Pitch_MA = []
    LT_Yaw_MA = []
    LT_MA_centroid = np.zeros((3,1)) 
    
    for i in range(iterations):
        
        
#definfing target points in MACS
        target_1 = np.matrix([[r*m.sin(m.radians(60))],[r*m.cos(m.radians(60))],[0],[1]])
        target_2 = np.matrix([[r*m.sin(m.radians(120))],[r*m.cos(m.radians(120))],[0],[1]])
        target_3 = np.matrix([[r*m.sin(m.radians(180))],[r*m.cos(m.radians(180))],[0],[1]])
        target_4 = np.matrix([[r*m.sin(m.radians(240))],[r*m.cos(m.radians(240))],[0],[1]])
        target_5 = np.matrix([[r*m.sin(m.radians(300))],[r*m.cos(m.radians(300))],[0],[1]])
        target_6 = np.matrix([[r*m.sin(m.radians(0))],[r*m.cos(m.radians(0))],[0],[1]])

#Transformation to convert targets from MA frame to Laser Tracker frame
        MA_to_LT = np.matrix([[1, 0, 0, 1.3],[0,1,0,1.2],[0,0,1,distance],[0,0,0,1]])
        
        Targets_MA_frame = np.concatenate((target_1, target_2, target_3, target_4, target_5, target_6), axis = 1)

#Converting from MA to Laser tracker coordinate frame
        Homogeneous_modelled_points = (MA_to_LT * Targets_MA_frame).T

#removing the fourth dimension for ease of calculations
        modelled_points = np.delete(Homogeneous_modelled_points, 3, 1)
    
        LT_err_MA = un.LT_uncertainities(modelled_points, 1).T

#Converting the modelled points into Laser tracker's spherical coordinate system
        spherical_points = cartesian_to_spherical(modelled_points)

#Adding the uncertainities from the Laser Tracker
        spherical_points = spherical_points + LT_err_MA

#Converting back to the cartesian coordinate system
        cartesian_points = -spherical_to_cartesian(spherical_points)


        [Rotation, Translation] = best_fit(cartesian_points, modelled_points, scaling)

        n = modelled_points.shape[0]

#modelled_points = modelled_points.T     
        fit_points = (Rotation * modelled_points.T) + np.tile(Translation, (1, n))

#Finding centroid of the fitted point cloud
        fit_point_centroid = np.mean(fit_points, axis = 1)
        
        LT_MA_centroid = np.append(LT_MA_centroid, fit_point_centroid, 1)
        
#Finding centroid of the modelled point cloud
        modelled_points_centroid = np.mean(modelled_points, axis = 0).T

#calculating homogeneous transformation function 
        Homogeneous_transform = np.zeros((4,4))
        Homogeneous_transform[:3,:3] = Rotation
        Homogeneous_transform[0,3] = Translation[0]
        Homogeneous_transform[1,3] = Translation[1]
        Homogeneous_transform[2,3] = Translation[2]
        Homogeneous_transform[3,3] = 1
        
#Finding the euler angles from the homogeneous transformation matrix
        LT_euler_angles_MA = np.matrix(t.euler_from_matrix(Homogeneous_transform)).T
        for i in range(3):
            LT_euler_angles_MA[i,0] = m.degrees(LT_euler_angles_MA[i,0]) * 3600

#Appending strings for 6-DOF values in each iteration
        LT_X_MA.append((Translation[0,0])*1e6)
        LT_Y_MA.append((Translation[1,0])*1e6)
        LT_Z_MA.append((Translation[2,0])*1e6)
        LT_Roll_MA.append(LT_euler_angles_MA[0,0])
        LT_Pitch_MA.append(LT_euler_angles_MA[1,0])
        LT_Yaw_MA.append(LT_euler_angles_MA[2,0])
        
#calculating mean of each 6-DOF strings to find the average over the number of iterations
    LT_X_MA = np.std(LT_X_MA)
    LT_Y_MA = np.std(LT_Y_MA)
    LT_Z_MA = np.std(LT_Z_MA)
    LT_Roll_MA = np.std(LT_Roll_MA)
    LT_Pitch_MA = np.std(LT_Pitch_MA)
    LT_Yaw_MA = np.std(LT_Yaw_MA)
    
    LT_MA_centroid = np.delete(LT_MA_centroid, 0, 1)
    LT_MA_centroid_STD = np.std(LT_MA_centroid, axis = 1)
    LT_MA_centroid = np.mean(LT_MA_centroid, axis = 1)

    return LT_X_MA, LT_Y_MA, LT_Z_MA, LT_Roll_MA, LT_Pitch_MA, LT_Yaw_MA, LT_MA_centroid, LT_MA_centroid_STD
    
'''
#########################################################################################################################################
'''

def MA_mechanical_error(mechanical_error, iterations):
    
    #Defining the radius of the mirror
    r = 1.3 

#defining the LT origin
    LT_origin = np.matrix([[0],[0],[0]])
    
#Initializing the arrays for the 6-DOF values    
    ME_X_MA = []
    ME_Y_MA = []
    ME_Z_MA = []
    ME_Roll_MA = []
    ME_Pitch_MA = []
    ME_Yaw_MA = []
    ME_MA_centroid = np.zeros((3,1)) 
    
    for i in range(iterations):
        
        
#definfing target points in MACS
        target_1 = np.matrix([[r*m.sin(m.radians(60))],[r*m.cos(m.radians(60))],[0],[1]])
        target_2 = np.matrix([[r*m.sin(m.radians(120))],[r*m.cos(m.radians(120))],[0],[1]])
        target_3 = np.matrix([[r*m.sin(m.radians(180))],[r*m.cos(m.radians(180))],[0],[1]])
        target_4 = np.matrix([[r*m.sin(m.radians(240))],[r*m.cos(m.radians(240))],[0],[1]])
        target_5 = np.matrix([[r*m.sin(m.radians(300))],[r*m.cos(m.radians(300))],[0],[1]])
        target_6 = np.matrix([[r*m.sin(m.radians(0))],[r*m.cos(m.radians(0))],[0],[1]])

#Transformation to convert targets from MA frame to Laser Tracker frame
        MA_to_LT = np.matrix([[1, 0, 0, 1.3],[0,1,0,1.2],[0,0,1,-6.3],[0,0,0,1]])
        
        Targets_MA_frame = np.concatenate((target_1, target_2, target_3, target_4, target_5, target_6), axis = 1)

#Converting from MA to Laser tracker coordinate frame
        Homogeneous_modelled_points = (MA_to_LT * Targets_MA_frame).T

#removing the fourth dimension for ease of calculations
        modelled_points = np.delete(Homogeneous_modelled_points, 3, 1)
    
        LT_err_MA = un.LT_uncertainities(modelled_points, 1).T
        Mechanical_err = un.Mechanical_uncertainities(modelled_points, -mechanical_error, -mechanical_error, -mechanical_error)

#Adding mechanical uncertainities of mounting the targets
        modelled_points_1 = modelled_points + Mechanical_err

#Converting the modelled points into Laser tracker's spherical coordinate system
        spherical_points = cartesian_to_spherical(modelled_points_1)

#Adding the uncertainities from the Laser Tracker
        spherical_points = spherical_points + LT_err_MA

#Converting back to the cartesian coordinate system
        cartesian_points = -spherical_to_cartesian(spherical_points)


        [Rotation, Translation] = best_fit(cartesian_points, modelled_points, scaling)

        n = modelled_points.shape[0]

#modelled_points = modelled_points.T     
        fit_points = (Rotation * modelled_points.T) + np.tile(Translation, (1, n))

#Finding centroid of the fitted point cloud
        fit_point_centroid = np.mean(fit_points, axis = 1)
        
        ME_MA_centroid = np.append(ME_MA_centroid, fit_point_centroid, 1)
        
#Finding centroid of the modelled point cloud
        modelled_points_centroid = np.mean(modelled_points, axis = 0).T

#calculating homogeneous transformation function 
        Homogeneous_transform = np.zeros((4,4))
        Homogeneous_transform[:3,:3] = Rotation
        Homogeneous_transform[0,3] = Translation[0]
        Homogeneous_transform[1,3] = Translation[1]
        Homogeneous_transform[2,3] = Translation[2]
        Homogeneous_transform[3,3] = 1
        
#Finding the euler angles from the homogeneous transformation matrix
        ME_euler_angles_MA = np.matrix(t.euler_from_matrix(Homogeneous_transform)).T
        for i in range(3):
            ME_euler_angles_MA[i,0] = m.degrees(ME_euler_angles_MA[i,0]) * 3600

#Appending strings for 6-DOF values in each iteration
        ME_X_MA.append((Translation[0,0])*1e6)
        ME_Y_MA.append((Translation[1,0])*1e6)
        ME_Z_MA.append((Translation[2,0])*1e6)
        ME_Roll_MA.append(ME_euler_angles_MA[0,0])
        ME_Pitch_MA.append(ME_euler_angles_MA[1,0])
        ME_Yaw_MA.append(ME_euler_angles_MA[2,0])
        
#calculating mean of each 6-DOF strings to find the average over the number of iterations
    ME_X_MA = np.mean(ME_X_MA)
    ME_Y_MA = np.mean(ME_Y_MA)
    ME_Z_MA = np.mean(ME_Z_MA)
    ME_Roll_MA = np.mean(ME_Roll_MA)
    ME_Pitch_MA = np.mean(ME_Pitch_MA)
    ME_Yaw_MA = np.mean(ME_Yaw_MA)
    
    ME_MA_centroid = np.delete(ME_MA_centroid, 0, 1)
    ME_MA_centroid = np.mean(ME_MA_centroid, axis = 1)

    return ME_X_MA, ME_Y_MA, ME_Z_MA, ME_Roll_MA, ME_Pitch_MA, ME_Yaw_MA, ME_MA_centroid



