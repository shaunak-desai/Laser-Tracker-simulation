'''
Each prefix before any vairable determines the case fo the simulation i.e. :
    LT -> Only errors from LT considered (case 1: with varying distance of LT from MA)
    ME -> Mechanical mounting errors considered on top of LT errors (case 2: with varying mechanical uncertainites)
    AE -> Angular errors and window imperfections added to the previous case. (case 3: with varying angular errors)
'''

import numpy as np
import CCD_alternate_approach as CCD
import MA_alternate_approach as MA
import matplotlib.pyplot as plt
import math as m

angular_errors = np.linspace(0, 30 ,num=10)
distance_MA = np.linspace(2.3, 9.7, num=10)
mech_errors = np.linspace(0, 20, num=10)

'''
#########################################################################################################################################
'''

def focal_length_error(A, B):
    
    i = A.shape[0]
    
    focal_length = []
    error_in_focal_length = []
    
    for j in range(10):
        x = []
        y = []
        for a in range(i):
            x.append((A[a,j] - B[a,j])**2)
        
        y = np.matrix(x)
        y = np.sum(x)
        y = m.sqrt(y)
        focal_length.append(y)
        
   
        
    for i in range(10):
        error_in_focal_length.append((12 - focal_length[i])*1e6)
    
    return focal_length, error_in_focal_length, y
    
def deviation(A,B):
    standard_deviation = []
    for i in range(10):
        X = (A[2,i]**2 + B[2,i]**2)
        Y = np.matrix(X)
        Y = np.sum(Y)
        Y = np.sqrt(Y)
        standard_deviation.append(Y*1e6) 
    return standard_deviation


'''
#########################################################################################################################################
'''
'''
This block of code determines the focal length errors and the 6-DOF errors with taking into considertaion 
the uncertainities from the LT, Mechanical mountings of the targets and the window imperfections for the targets mounted on the CCD

Further on, additional angular uncertainities have also been added to simulate the effects of using a window with 
lower centring tolerance

'CCD_angular_error' function takes in the arguements as follows:
    1. The angular errors from the window imperfections (centring/parallelism tolerances)
    2. The thickness uncertainities of the window
    3. The added angular errors for depicting different cases of window with increasing centring tolerances
    4. The mechanical errors during mounting of the LT targets
    5. The number of iterations you would like the code to run
'''



AE_X_CCD = []
AE_Y_CCD = []
AE_Z_CCD = []
AE_roll_CCD = []
AE_pitch_CCD = []
AE_yaw_CCD = []

AE_CCD_centroid = np.zeros((3,1))
AE_CCD_centroid_STD = np.zeros((3,1))


for angular_error in np.linspace(0,0.00002424*6,num=10):
    
    [X, Y, Z, Roll, Pitch, Yaw, centroid, std] = CCD.CCD_error(0, 3e-5, angular_error, 0, 5.7, 10000)
    
    AE_X_CCD.append(X)
    AE_Y_CCD.append(Y)
    AE_Z_CCD.append(Z)
    AE_roll_CCD.append(Roll)
    AE_pitch_CCD.append(Pitch)
    AE_yaw_CCD.append(Yaw)
    AE_CCD_centroid = np.append(AE_CCD_centroid, centroid, 1)
    AE_CCD_centroid_STD = np.append(AE_CCD_centroid_STD, std, 1)
    
AE_CCD_centroid = np.delete(AE_CCD_centroid, 0, 1)
AE_CCD_centroid_STD = np.delete(AE_CCD_centroid_STD, 0, 1)

'''
#########################################################################################################################################
'''

'''
This block of code determines the focal length errors and the 6-DOF errors with taking into considertaion 
the uncertainities from the LT, Mechanical mountings of the targets and the window imperfections for the targets mounted on the MA

Further on, additional angular uncertainities have also been added to simulate the effects of using a window with 
lower centring tolerance

'MA_angular_error' function takes in the arguements as follows:
    1. The angular errors from the window imperfections (centring/parallelism tolerances)
    2. The thickness uncertainities of the window
    3. The added angular errors for depicting different cases of window with increasing centring tolerances
    4. The mechanical errors during mounting of the LT targets
    5. The number of iterations you would like the code to run
'''



AE_X_MA = []
AE_Y_MA = []
AE_Z_MA = []
AE_roll_MA = []
AE_pitch_MA = []
AE_yaw_MA = []

AE_MA_centroid = np.zeros((3,1))
AE_MA_centroid_STD = np.zeros((3,1))


for angular_error in np.linspace(0,0.00002424*6,num=10):
    
    [X, Y, Z, Roll, Pitch, Yaw, centroid, std] = MA.MA_error(0, 3e-5, angular_error, 0, -6.3, 10000)
    
    AE_X_MA.append(X)
    AE_Y_MA.append(Y)
    AE_Z_MA.append(Z)
    AE_roll_MA.append(Roll)
    AE_pitch_MA.append(Pitch)
    AE_yaw_MA.append(Yaw)
    AE_MA_centroid = np.append(AE_MA_centroid, centroid, 1)
    AE_MA_centroid_STD = np.append(AE_MA_centroid_STD, std, 1)
    
AE_MA_centroid = np.delete(AE_MA_centroid, 0, 1)
AE_MA_centroid_STD = np.delete(AE_MA_centroid_STD, 0, 1)

'''
#########################################################################################################################################
'''
'''
This block takes in the matrix with the centroids of MA and CCD targets of each individual cases and calculates the focal length
and errors in focal length for each case and plots the errors in focal length as well as the 6-DOF values agaainst the various cases
of the simulation
'''


#[AE_focal_length, AE_error_in_focal_length, y] = focal_length_error(AE_MA_centroid, AE_CCD_centroid)

AE_focal_length_error = deviation(AE_CCD_centroid_STD, AE_MA_centroid_STD)


plt.figure()
plt.plot(angular_errors, AE_roll_MA, 'b')
plt.plot(angular_errors, AE_pitch_MA, 'g')
plt.plot(angular_errors, AE_yaw_MA, 'y')
plt.plot(angular_errors, AE_roll_CCD, '--b')
plt.plot(angular_errors, AE_pitch_CCD, '--g')
plt.plot(angular_errors, AE_yaw_CCD, '--y')
plt.ylabel('Rotational uncertainties (arcseconds) (1-sigma)')
plt.xlabel('Incresing angular errors introduced by the window (arcseconds) (1-sigma)')
plt.legend(['AE_roll_MA (1-sigma)','AE_pitch_MA (1-sigma)', 'AE_yaw_MA (1-sigma)', 'AE_roll_CCD (1-sigma)', 'AE_pitch_CCD (1-sigma)', 'AE_yaw_CCD (1-sigma)' ], loc=2)
plt.title('Rotational errors (1-sigma) induced by adding the window imperfections')
plt.rcParams['figure.figsize'] = [16,9]
plt.grid()
plt.show   

plt.figure()
fig, ax = plt.subplots()
ax.plot(angular_errors, AE_X_MA,'b',)
ax.plot(angular_errors, AE_Y_MA,'g')
ax.plot(angular_errors, AE_Z_MA,'y')
ax.plot(angular_errors, AE_X_CCD,'--b')
ax.plot(angular_errors, AE_Y_CCD,'--g')
ax.plot(angular_errors, AE_Z_CCD,'--y')
ax.plot(angular_errors, AE_focal_length_error, color = 'r', linewidth = 4)
#ax2 = ax.twinx()
#ax2.plot(angular_errors, AE_focal_length_error, color = 'r')
ax.set_ylabel('Translational uncertainties (micron) & Focal length error (micron) (1-sigma)')
ax.set_xlabel('Increasing angular errors introduced by the window (arcseconds) (1-sigma)')
#ax2.set_ylabel('Focal length error (microns) (1-sigma)')
ax.legend(['AE_X_MA (1-sigma)', 'AE_Y_MA (1-sigma)', 'AE_Z_MA (1-sigma)','AE_X_CCD (1-sigma)', 'AE_Y_CCD (1-sigma)', 'AE_Z_CCD (1-sigma)', 'Focal length error (1-sigma)' ], loc=2)
#ax2.legend(['Focal length error (1-sigma)'], loc=4)
plt.title('Translational errors (1-sigma) and focal length errors (1-sigma) introduced by the window imperfections ')
plt.rcParams['figure.figsize'] = [16,9]
ax.grid()
plt.show

'''
#########################################################################################################################################
'''

'''
This block of code determines the focal length errors and the 6-DOF errors with taking into considertaion 
the uncertainities from the LT for the targets mounted on the CCD 

Further on, various siulation results for different distances of the LT from the CCD modules is depicted

'CCD_LT_error' function takes in the arguements as follows:
    1. The  distance of the LT from the CCD module
    2. The number of iterations you would like the code to run
'''
'''
LT_X_CCD = []
LT_Y_CCD = []
LT_Z_CCD = []
LT_roll_CCD = []
LT_pitch_CCD = []
LT_yaw_CCD = []

LT_CCD_centroid = np.zeros((3,1))
LT_CCD_centroid_STD = np.zeros((3,1))

for distance in np.linspace(9.7, 2.3, num=10):
    
    [X, Y, Z, Roll, Pitch, Yaw, centroid, std] = CCD.CCD_error(0, 0, 0, 0, distance, 10000)

    LT_X_CCD.append(X)
    LT_Y_CCD.append(Y)
    LT_Z_CCD.append(Z)
    LT_roll_CCD.append(Roll)
    LT_pitch_CCD.append(Pitch)
    LT_yaw_CCD.append(Yaw)
    LT_CCD_centroid = np.append(LT_CCD_centroid, centroid, 1)
    LT_CCD_centroid_STD = np.append(LT_CCD_centroid_STD, std, 1)
    
LT_CCD_centroid = np.delete(LT_CCD_centroid, 0, 1)
LT_CCD_centroid_STD = np.delete(LT_CCD_centroid_STD, 0, 1)
'''
'''
#########################################################################################################################################
'''

'''
This block of code determines the focal length errors and the 6-DOF errors with taking into considertaion 
the uncertainities from the LT for the targets mounted on the MA 

Further on, various siulation results for different distances of the LT from the MA modules is depicted

'MA_LT_error' function takes in the arguements as follows:
    1. The  distance of the LT from the MA module
    2. The number of iterations you would like the code to run
'''
'''
LT_X_MA = []
LT_Y_MA = []
LT_Z_MA = []
LT_roll_MA = []
LT_pitch_MA = []
LT_yaw_MA = []
LT_MA_centroid_STD = np.zeros((3,1))

LT_MA_centroid = np.zeros((3,1))

for distance in np.linspace(-2.3, -9.7, num=10):
    
    [X, Y, Z, Roll, Pitch, Yaw, centroid, std] = MA.MA_error(0, 0, 0, 0, distance, 10000)
    
    LT_X_MA.append(X)
    LT_Y_MA.append(Y)
    LT_Z_MA.append(Z)
    LT_roll_MA.append(Roll)
    LT_pitch_MA.append(Pitch)
    LT_yaw_MA.append(Yaw)
    LT_MA_centroid = np.append(LT_MA_centroid, centroid, 1)
    LT_MA_centroid_STD = np.append(LT_MA_centroid_STD, std, 1)
    
LT_MA_centroid = np.delete(LT_MA_centroid, 0, 1) 
LT_MA_centroid_STD = np.delete(LT_MA_centroid_STD, 0, 1)  
'''
    
'''
#########################################################################################################################################
'''  
'''
This block takes in the matrix with the centroids of MA and CCD targets of each individual cases and calculates the focal length
and errors in focal length for each case and plots the errors in focal length as well as the 6-DOF values agaainst the various cases
of the simulation
'''

'''
#[LT_focal_length, LT_error_in_focal_length, y] = focal_length_error(LT_MA_centroid, LT_CCD_centroid)
LT_focal_length_error = deviation(LT_CCD_centroid_STD, LT_MA_centroid_STD)

plt.figure()
plt.rcParams['figure.figsize'] = [16,9]
plt.plot(distance_MA, LT_roll_MA, 'b')
plt.plot(distance_MA, LT_pitch_MA, 'g')
plt.plot(distance_MA, LT_yaw_MA, 'y')
plt.plot(distance_MA, LT_roll_CCD, '--b')
plt.plot(distance_MA, LT_pitch_CCD, '--g')
plt.plot(distance_MA, LT_yaw_CCD, '--y')
plt.ylabel('Rotational uncertainties (arcseconds) (1-sigma)')
plt.xlabel('Incresing distance of LT from MA (metres)')
plt.legend(['LT_roll_MA (1-sigma)','LT_pitch_MA (1-sigma)', 'LT_yaw_MA (1-sigma)', 'LT_roll_CCD (1-sigma)', 'LT_pitch_CCD (1-sigma)', 'LT_yaw_CCD (1-sigma)' ], loc=2)
plt.title('Rotational errors introduced by the LT (1-sigma)')
plt.grid()
plt.show   

plt.figure()
fig, ax = plt.subplots()
ax.plot(distance_MA, LT_X_MA,'b',)
ax.plot(distance_MA, LT_Y_MA,'g')
ax.plot(distance_MA, LT_Z_MA,'y')
ax.plot(distance_MA, LT_X_CCD,'--b')
ax.plot(distance_MA, LT_Y_CCD,'--g')
ax.plot(distance_MA, LT_Z_CCD,'--y')
ax.plot(distance_MA, LT_focal_length_error, color = 'r', linewidth = 4)
#ax2 = ax.twinx()
#ax2.plot(distance_MA, LT_focal_length_error, color = 'r')
ax.set_ylabel('Translational uncertainties (micron) & Focal length error (micron) (1-sigma)')
ax.set_xlabel('Incresing distance of LT from MA (metres)')
#ax2.set_ylabel('Focal length error (microns) (1-sigma)')
ax.legend(['LT_X_MA (1-sigma)', 'LT_Y_MA (1-sigma)', 'LT_Z_MA (1-sigma)','LT_X_CCD (1-sigma)', 'LT_Y_CCD (1-sigma)', 'LT_Z_CCD (1-sigma)', 'Focal length error (1-sigma)' ], loc=2)
#ax2.legend(['Focal length error (1-sigma)'], loc=4)
plt.title('Transalational errors (1-sigma) and errors in focal length (1-sigma) introduced by the LT')
plt.rcParams['figure.figsize'] = [16,9]
ax.grid()
plt.show
'''


'''
#########################################################################################################################################
''' 


'''
This block of code determines the focal length errors and the 6-DOF errors with taking into considertaion 
the uncertainities from the LT, Mechanical mountings of the targets mounted on the CCD 

Further on, the simualtions depicts the impacts of increasing the mechanical mounting errors on the focal length errors
as well as the 6-DOF values for the CCD targets

'CCD_mechanical_errors' takes in the arguements as follows:
    1. The mechanical mounting error for the targets
    2. The number of iterations you would like the code to run
'''

'''
ME_X_CCD = []
ME_Y_CCD = []
ME_Z_CCD = []
ME_roll_CCD = []
ME_pitch_CCD = []
ME_yaw_CCD = []

ME_CCD_centroid = np.zeros((3,1))
ME_CCD_centroid_STD = np.zeros((3,1))

for mechanical_error in np.linspace(0, 20e-6, num=10):
    
     [X, Y, Z, Roll, Pitch, Yaw, centroid, std] = CCD.CCD_error(0, 0, 0, mechanical_error, 5.7,  10000)
     
     ME_X_CCD.append(X)
     ME_Y_CCD.append(Y)
     ME_Z_CCD.append(Z)
     ME_roll_CCD.append(Roll)
     ME_pitch_CCD.append(Pitch)
     ME_yaw_CCD.append(Yaw)
     ME_CCD_centroid = np.append(ME_CCD_centroid, centroid, 1)
     ME_CCD_centroid_STD = np.append(ME_CCD_centroid_STD, std, 1)
    
ME_CCD_centroid = np.delete(ME_CCD_centroid, 0, 1)
ME_CCD_centroid_STD = np.delete(ME_CCD_centroid_STD, 0, 1)
'''
   
'''
#########################################################################################################################################
'''  

'''
This block of code determines the focal length errors and the 6-DOF errors with taking into considertaion 
the uncertainities from the LT, Mechanical mountings of the targets mounted on the MA 

Further on, the simualtions depicts the impacts of increasing the mechanical mounting errors on the focal length errors
as well as the 6-DOF values for the MA targets

'MA_mechanical_errors' takes in the arguements as follows:
    1. The mechanical mounting error for the targets
    2. The number of iterations you would like the code to run
'''

'''
ME_X_MA = []
ME_Y_MA = []
ME_Z_MA = []
ME_roll_MA = []
ME_pitch_MA = []
ME_yaw_MA = []

ME_MA_centroid = np.zeros((3,1)) 
ME_MA_centroid_STD = np.zeros((3,1)) 
    
for mechanical_error in np.linspace(0, 20e-6, num=10):
    
     [X, Y, Z, Roll, Pitch, Yaw, centroid, std] = MA.MA_error(0, 0, 0, mechanical_error, -6.3,  10000) 
     ME_X_MA.append(X)
     ME_Y_MA.append(Y)
     ME_Z_MA.append(Z)
     ME_roll_MA.append(Roll)
     ME_pitch_MA.append(Pitch)
     ME_yaw_MA.append(Yaw)
     ME_MA_centroid = np.append(ME_MA_centroid, centroid, 1)
     ME_MA_centroid_STD = np.append(ME_MA_centroid_STD, std, 1)
    
ME_MA_centroid = np.delete(ME_MA_centroid, 0, 1)
ME_MA_centroid_STD = np.delete(ME_MA_centroid_STD, 0, 1)
'''

'''
#########################################################################################################################################
'''

'''
This block takes in the matrix with the centroids of MA and CCD targets of each individual cases and calculates the focal length
and errors in focal length for each case and plots the errors in focal length as well as the 6-DOF values agaainst the various cases
of the simulation
'''
'''
#[ME_focal_length, ME_error_in_focal_length, y] = focal_length_error(ME_MA_centroid, ME_CCD_centroid)
ME_focal_length_error = deviation(ME_CCD_centroid_STD, ME_MA_centroid_STD)


plt.figure(3)
plt.rcParams['figure.figsize'] = [16,9]
plt.plot(mech_errors, ME_roll_MA, 'b', label = 'ME_roll_MA')
plt.plot(mech_errors, ME_pitch_MA, 'g', label = 'ME_pitch_MA')
plt.plot(mech_errors, ME_yaw_MA, 'y', label = 'ME_yaw_MA')
plt.plot(mech_errors, ME_roll_CCD, '--b', label = 'ME_roll_CCD')
plt.plot(mech_errors, ME_pitch_CCD, '--g', label = 'ME_pitch_CCD')
plt.plot(mech_errors, ME_yaw_CCD, '--y', label = 'ME_ya_CCD')
plt.ylabel('Rotational uncertainties (arcseconds) (1-sigma)')
plt.xlabel('Incresing mechanical uncertainities (microns)(1-sigma)')
plt.legend(['ME_roll_MA (1-sigma)','ME_pitch_MA (1-sigma)', 'ME_yaw_MA (1-sigma)', 'ME_roll_CCD (1-sigma)', 'ME_pitch_CCD (1-sigma)', 'ME_yaw_CCD (1-sigma)' ], loc=2)
plt.title('Rotational errors (1-sigma) from increasing the mechanical mounting uncertainities')

plt.grid()
plt.show    


plt.figure(4)
fig, ax = plt.subplots()
ax.plot(mech_errors, ME_X_MA,'b',)
ax.plot(mech_errors, ME_Y_MA,'g')
ax.plot(mech_errors, ME_Z_MA,'y')
ax.plot(mech_errors, ME_X_CCD,'--b')
ax.plot(mech_errors, ME_Y_CCD,'--g')
ax.plot(mech_errors, ME_Z_CCD,'--y')
ax.plot(mech_errors, ME_focal_length_error, color = 'r', linewidth = 4)
#ax2 = ax.twinx()
#ax2.plot(mech_errors, ME_focal_length_error, color = 'r')
ax.set_ylabel('Translational uncertainties & Focal length error (microns) (1-sigma)')
ax.set_xlabel('Incresing mechanical uncertainities (microns) (1-sigma)')
#ax2.set_ylabel('Errors in focal length (microns) (1-sigma)')
ax.legend(['ME_X_MA (1-sigma)', 'ME_Y_MA (1-sigma)', 'ME_Z_MA (1-sigma)','ME_X_CCD (1-sigma)', 'ME_Y_CCD (1-sigma)', 'ME_Z_CCD (1-sigma)', 'Focal length error (1-sigma)' ], loc=2)
#ax2.legend(['Focal length error (1-sigma)'], loc=4)
plt.title('MA transalational errors (1-sigma) anf focal length errors (1-sigma) by increasing the mechanical mounting uncertainities ')
plt.rcParams['figure.figsize'] = [16,9]
ax.grid()
plt.show
'''
   

'''
#########################################################################################################################################
''' 



