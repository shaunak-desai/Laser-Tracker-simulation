# Laser-Tracker-simulation
Simulations of the errors sourced in the Laser Tacker metrology


The file error_main is the main file which calls the rest of the files upon execution. It calculates the overall error in measurement in each cartesian coordinate axes in the Mirror Assembly Module and the CCD module. It takes the distance of the Laser Tracker from the Mirror Assembly Module, mechanical uncertainty and the tolerences on thickness and the centring of the window as the input parameters. From these, it cacluates the errors from by calling functions from 'Uncertainities'. These errors are added to the respective targets on both the CCD and Mirror Assembly Module by calling the MA_alternate_approach and CCD_alternate_approach. These files also perform the best fit of these points with the said uncertainties to the original defined points, without the errors. After this, the errors in individual coordinates are appended for each test case by the main file and plotted. 

