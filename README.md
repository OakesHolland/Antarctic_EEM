# Antarctic_EEM
 Esemble ecosystem modelling of impacts of marine invasive species in Antarctic ecosystems
 
Author: Oakes Holland (o.holland@hdr.qut.edu.au)

File information:

	invasion_functions.r : 	The functions required to run the introduction of an invasive species.
							This code is produced for a situation where species can both eat, and 
							be eaten by each other. Time period for simulation to run (in years) 
							can be changed on line 22. 

	run_invasion.r  :		This is the code to run the functions. rep = the number of simulations
							to run. dist_mean can be changed, however Baker et al. 2016 found no
							difference in the outcome when this value was changed. Time functions 
							are optional to keep time of how long the simulations take to complete.
							Note: with the data provided the simulations will take a long time to 
							complete (2 - 8 hours per rep). Less species = less time to run.
	
	growth_rate.csv  :		The average growth rates for the species within each node. 
	
	interactions.csv  :		The matrix of interactions between species. 0 is no interaction, -1 is
							a negative impact (eg. prey), 1 is a positive impact (eg. predator).
	
	interactions2.csv  :	The matrix of interactions where species eat each other. Note: these
							interactions are only coded as 1, there is no negative numbers in this
							matrix.		



This code has been adapted from:
 Baker et al. 2016 
 Ensemble ecosystem modelling for predicting ecosystem repsonse to predator reintroduction
 DOI 10.1111/cobi.12798
 Original code and data: https://doi.org/10.6084/m9.figshare.c.3320556.v1
 
 

 


