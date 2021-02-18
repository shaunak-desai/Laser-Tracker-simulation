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

scaling = False
#####
test = 0.00002424 * 2
#####

'''
#########################################################################################################################################
'''

def CCD_error(window_angular_error, window_range_error, angular_error, mechanical_error, distance, iterations):
    
#Initializing the arrays for the 6-DOF values    
    X_CCD = []
    Y_CCD = []
    Z_CCD = []
    Roll_CCD = []
    Pitch_CCD = []
    Yaw_CCD = []
    CCD_centroid = np.zeros((3,1)) 
    
   
    
    for i in range(iterations):
        
        
#the cartesian cordinates of targets in the CCD frame is known
        CCD_target_1 = np.matrix([[0.034125],[0.039125],[0],[1]])
        CCD_target_2 = np.matrix([[-0.034125],[0.039125],[0],[1]])
        CCD_target_3 = np.matrix([[-0.034125],[-0.039125],[0],[1]])
        CCD_target_4 = np.matrix([[0.034125],[-0.039125],[0],[1]])

        Targets_CCD_frame = np.concatenate((CCD_target_1, CCD_target_2, CCD_target_3, CCD_target_4), axis = 1)
        
        

#transformation matrix for converting the CCD frame to the LT frame
        CCD_to_LT = np.matrix([[1,0,0,1.3],[0,1,0,1.2],[0,0,1,distance],[0,0,0,1]])

        Homogeneous_modelled_points = (CCD_to_LT * Targets_CCD_frame).T

#removing the fourth dimension for ease of calculations
        modelled_points = np.delete(Homogeneous_modelled_points, 3, 1)
    
        LT_err_CCD = un.LT_uncertainities(modelled_points, 1).T
        Window_err_CCD = un.Window_uncertainities(modelled_points, window_range_error, window_angular_error, window_angular_error)
        Additional_err = un.Additional_ang_uncertatinities(modelled_points, angular_error, angular_error)
        Mechanical_err = un.Mechanical_uncertainities(modelled_points, mechanical_error, mechanical_error, mechanical_error)

#Adding mechanical uncertainities of mounting the targets
        modelled_points_1 = modelled_points + Mechanical_err

#Converting the modelled points into Laser tracker's spherical coordinate system
        spherical_points = cartesian_to_spherical(modelled_points_1)

#Adding the uncertainities from the Laser Tracker
        #spherical_points = spherical_points + LT_err_CCD

#Adding uncertainities from the window
        spherical_points = spherical_points + Window_err_CCD 

#Adding additional angular uncertainities
        spherical_points = spherical_points + Additional_err

#Converting back to the cartesian coordinate system
        cartesian_points = spherical_to_cartesian(spherical_points)
        
        Homogeneous_transform = t.superimposition_matrix(modelled_points, cartesian_points, usesvd=True)
        Rotation = np.matrix(Homogeneous_transform[0:3, 0:3])
        Translation = np.matrix(Homogeneous_transform[0:3,3]).T


        #[Rotation, Translation] = best_fit(cartesian_points, modelled_points, scaling)

        n = modelled_points.shape[0]

#modelled_points = modelled_points.T     
        fit_points = (Rotation * modelled_points.T) + np.tile(Translation, (1, n))

#Finding centroid of the fitted point cloud
        fit_point_centroid = np.mean(fit_points, axis = 1)
        
        CCD_centroid = np.append(CCD_centroid, fit_point_centroid, 1)
       
        
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
        euler_angles_CCD = np.matrix(t.euler_from_matrix(Homogeneous_transform)).T
        
        for i in range(3):
            euler_angles_CCD[i,0] = m.degrees(euler_angles_CCD[i,0]) * 3600

#Appending strings for 6-DOF values in each iteration
        X_CCD.append((Translation[0,0])*1e6)
        Y_CCD.append((Translation[1,0])*1e6)
        Z_CCD.append((Translation[2,0])*1e6)
        Roll_CCD.append(euler_angles_CCD[0,0])
        Pitch_CCD.append(euler_angles_CCD[1,0])
        Yaw_CCD.append(euler_angles_CCD[2,0])
        
#calculating STD of each 6-DOF strings to find the average over the number of iterations
    X_CCD = np.std(X_CCD)
    Y_CCD = np.std(Y_CCD)
    Z_CCD = np.std(Z_CCD)
    Roll_CCD = np.std(Roll_CCD)
    Pitch_CCD = np.std(Pitch_CCD)
    Yaw_CCD = np.std(Yaw_CCD)
    
    CCD_centroid = np.delete(CCD_centroid, 0, 1)
    CCD_centroid_STD = np.std(CCD_centroid, axis = 1)
    CCD_centroid = np.mean(CCD_centroid, axis = 1)

    return X_CCD, Y_CCD, Z_CCD, Roll_CCD, Pitch_CCD, Yaw_CCD, CCD_centroid, CCD_centroid_STD










'''
#########################################################################################################################################
'''

def CCD_LT_error(distance, iterations):
    
#Initializing the arrays for the 6-DOF values    
    LT_X_CCD = []
    LT_Y_CCD = []
    LT_Z_CCD = []
    LT_Roll_CCD = []
    LT_Pitch_CCD = []
    LT_Yaw_CCD = []
    LT_CCD_centroid = np.zeros((3,1)) 
    
    for i in range(iterations):
        
        
#the cartesian cordinates of targets in the CCD frame is known
        CCD_target_1 = np.matrix([[0.034125],[0.039125],[0],[1]])
        CCD_target_2 = np.matrix([[-0.034125],[0.039125],[0],[1]])
        CCD_target_3 = np.matrix([[-0.034125],[-0.039125],[0],[1]])
        CCD_target_4 = np.matrix([[0.034125],[-0.039125],[0],[1]])

        Targets_CCD_frame = np.concatenate((CCD_target_1, CCD_target_2, CCD_target_3, CCD_target_4), axis = 1)

#transformation matrix for converting the CCD frame to the LT frame
        CCD_to_LT = np.matrix([[1,0,0,1.3],[0,1,0,1.2],[0,0,1,distance],[0,0,0,1]])

        Homogeneous_modelled_points = (CCD_to_LT * Targets_CCD_frame).T

#removing the fourth dimension for ease of calculations
        modelled_points = np.delete(Homogeneous_modelled_points, 3, 1)
    
        LT_err_CCD = un.LT_uncertainities(modelled_points, 1).T

#Converting the modelled points into Laser tracker's spherical coordinate system
        spherical_points = cartesian_to_spherical(modelled_points)

#Adding the uncertainities from the Laser Tracker
        spherical_points = spherical_points + LT_err_CCD

#Converting back to the cartesian coordinate system
        cartesian_points = spherical_to_cartesian(spherical_points)


        [Rotation, Translation] = best_fit(cartesian_points, modelled_points, scaling)

        n = modelled_points.shape[0]

#modelled_points = modelled_points.T     
        fit_points = (Rotation * modelled_points.T) + np.tile(Translation, (1, n))

#Finding centroid of the fitted point cloud
        fit_point_centroid = np.mean(fit_points, axis = 1)
        
        LT_CCD_centroid = np.append(LT_CCD_centroid, fit_point_centroid, 1)
        
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
        LT_euler_angles_CCD = np.matrix(t.euler_from_matrix(Homogeneous_transform)).T
        for i in range(3):
            LT_euler_angles_CCD[i,0] = m.degrees(LT_euler_angles_CCD[i,0]) * 3600

#Appending strings for 6-DOF values in each iteration
        LT_X_CCD.append((Translation[0,0])*1e6)
        LT_Y_CCD.append((Translation[1,0])*1e6)
        LT_Z_CCD.append((Translation[2,0])*1e6)
        LT_Roll_CCD.append(LT_euler_angles_CCD[0,0])
        LT_Pitch_CCD.append(LT_euler_angles_CCD[1,0])
        LT_Yaw_CCD.append(LT_euler_angles_CCD[2,0])
        
#calculating mean of each 6-DOF strings to find the average over the number of iterations
    LT_X_CCD = np.std(LT_X_CCD)
    LT_Y_CCD = np.std(LT_Y_CCD)
    LT_Z_CCD = np.std(LT_Z_CCD)
    LT_Roll_CCD = np.std(LT_Roll_CCD)
    LT_Pitch_CCD = np.std(LT_Pitch_CCD)
    LT_Yaw_CCD = np.std(LT_Yaw_CCD)
    
    LT_CCD_centroid = np.delete(LT_CCD_centroid, 0, 1)
    LT_CCD_centroid_STD = np.std(LT_CCD_centroid, axis = 1)
    LT_CCD_centroid = np.mean(LT_CCD_centroid, axis = 1)

    return LT_X_CCD, LT_Y_CCD, LT_Z_CCD, LT_Roll_CCD, LT_Pitch_CCD, LT_Yaw_CCD, LT_CCD_centroid, LT_CCD_centroid_STD

'''
#########################################################################################################################################
'''

def CCD_mechanical_error(mechanical_error, iterations):
    
#Initializing the arrays for the 6-DOF values    
    ME_X_CCD = []
    ME_Y_CCD = []
    ME_Z_CCD = []
    ME_Roll_CCD = []
    ME_Pitch_CCD = []
    ME_Yaw_CCD = []
    ME_CCD_centroid = np.zeros((3,1)) 
    
    for i in range(iterations):
        
        
#the cartesian cordinates of targets in the CCD frame is known
        CCD_target_1 = np.matrix([[0.034125],[0.039125],[0],[1]])
        CCD_target_2 = np.matrix([[-0.034125],[0.039125],[0],[1]])
        CCD_target_3 = np.matrix([[-0.034125],[-0.039125],[0],[1]])
        CCD_target_4 = np.matrix([[0.034125],[-0.039125],[0],[1]])

        Targets_CCD_frame = np.concatenate((CCD_target_1, CCD_target_2, CCD_target_3, CCD_target_4), axis = 1)

#transformation matrix for converting the CCD frame to the LT frame
        CCD_to_LT = np.matrix([[1,0,0,1.3],[0,1,0,1.2],[0,0,1,5.7],[0,0,0,1]])

        Homogeneous_modelled_points = (CCD_to_LT * Targets_CCD_frame).T

#removing the fourth dimension for ease of calculations
        modelled_points = np.delete(Homogeneous_modelled_points, 3, 1)
    
        LT_err_CCD = un.LT_uncertainities(modelled_points, 1).T
        Mechanical_err = un.Mechanical_uncertainities(modelled_points, mechanical_error, mechanical_error, mechanical_error)

#Adding mechanical uncertainities of mounting the targets
        modelled_points_1 = modelled_points + Mechanical_err

#Converting the modelled points into Laser tracker's spherical coordinate system
        spherical_points = cartesian_to_spherical(modelled_points_1)

#Adding the uncertainities from the Laser Tracker
        spherical_points = spherical_points + LT_err_CCD

#Converting back to the cartesian coordinate system
        cartesian_points = spherical_to_cartesian(spherical_points)


        [Rotation, Translation] = best_fit(cartesian_points, modelled_points, scaling)

        n = modelled_points.shape[0]

#modelled_points = modelled_points.T     
        fit_points = (Rotation * modelled_points.T) + np.tile(Translation, (1, n))

#Finding centroid of the fitted point cloud
        fit_point_centroid = np.mean(fit_points, axis = 1)
        
        ME_CCD_centroid = np.append(ME_CCD_centroid, fit_point_centroid, 1)
        
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
        ME_euler_angles_CCD = np.matrix(t.euler_from_matrix(Homogeneous_transform)).T
        for i in range(3):
            ME_euler_angles_CCD[i,0] = m.degrees(ME_euler_angles_CCD[i,0]) * 3600

#Appending strings for 6-DOF values in each iteration
        ME_X_CCD.append((Translation[0,0])*1e6)
        ME_Y_CCD.append((Translation[1,0])*1e6)
        ME_Z_CCD.append((Translation[2,0])*1e6)
        ME_Roll_CCD.append(ME_euler_angles_CCD[0,0])
        ME_Pitch_CCD.append(ME_euler_angles_CCD[1,0])
        ME_Yaw_CCD.append(ME_euler_angles_CCD[2,0])
        
#calculating mean of each 6-DOF strings to find the average over the number of iterations
    ME_X_CCD = np.mean(ME_X_CCD)
    ME_Y_CCD = np.mean(ME_Y_CCD)
    ME_Z_CCD = np.mean(ME_Z_CCD)
    ME_Roll_CCD = np.mean(ME_Roll_CCD)
    ME_Pitch_CCD = np.mean(ME_Pitch_CCD)
    ME_Yaw_CCD = np.mean(ME_Yaw_CCD)
    
    ME_CCD_centroid = np.delete(ME_CCD_centroid, 0, 1)
    ME_CCD_centroid = np.mean(ME_CCD_centroid, axis = 1)

    return ME_X_CCD, ME_Y_CCD, ME_Z_CCD, ME_Roll_CCD, ME_Pitch_CCD, ME_Yaw_CCD, ME_CCD_centroid
