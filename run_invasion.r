rep = 10 #Number of plausible ensembles created
dist_mean = 1 #The mean of the uniform distribution that the interactions are drawn from

source("invasion_functions.r") #File where functions are stored
old <- Sys.time() #Used for timing the run, optional
invasion_functions(dist_mean,rep) #Run the reintro functions

new <- Sys.time() #Used for timing the run, optional
time <- new - old #Used for timing the run, optional
print(time) #Print the run time, optional






