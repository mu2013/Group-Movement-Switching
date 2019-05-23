


This code provide the mcmc samples for all movement diffusion parameters and switching parameters. 
User can make inference based on the real dataset by uncomment "source('Data_real.r')" in the Setup.r file. Or make inference for the simulation data by uncomment "source('Data_sim.r')" in the Setup.r file. 
In order to run the code, first run Main_Inference.r. The initial looping is 50 iterations. 

setup.r		        load group movement reindeer data or simulation data and initialise model parameter and prior
Main_Inference.r 	the main mcmc loop
RunKF.r	            run Kalman Filter to calculate the log likelihood
SdeBM.r		        calculate the diffusion parameter matrix(A in paper) and the initial P for Kalman filter
SolveSdeBM.r 	    calculate the variance matrix (\Xi in paper)

This model can deal with missing values. 

