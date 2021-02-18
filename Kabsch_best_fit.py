#Kabsch algorithm for best fit
import numpy as np
#import CCD_best_fit


#A = CCD_best_fit.modelled_points
#B = CCD_best_fit.measured_points

def best_fit(A,B):
    
    assert len(A) == len(B)
    
    N = A.shape[0]; #total number of points
    
    A_centroid = np.mean(A, axis = 0)
    B_centroid = np.mean(B, axis = 0)

    #bringing the data sets to the origin (alllows to calculate just the rotation matrix)
    A_origin = A - np.tile(A_centroid, (N,1))
    B_origin = B - np.tile(B_centroid, (N,1))
    
    #calculating the covariance matrix
    H =  np.transpose(A_origin) * B_origin
    
    #Using SVD to obtain the desired rotation matrix
    [U,S,V] = np.linalg.svd(H)
    
    #rotation matrix 
    R = V * np.transpose(U)
    
     # special reflection case
    if np.linalg.det(R) < 0:
        V[2, :] *= -1
        R = V.T * U.T
    
    #finding the translation matrix
    t = A_centroid.T - R * B_centroid.T
    
    return R, t

#A = CCD_best_fit.modelled_points
#B = CCD_best_fit.measured_points

#[Rotation, Translate] =  best_fit(A, B)
