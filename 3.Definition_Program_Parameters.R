##################################################
### --- 1. DEFINITION OF PROGRAM PARAMETERS      #
##################################################

# Number of clusters for mixture of distributions
k <- 10
# Numero de variables a imputar
#r <- 2 # Se genera en PrepareData
# Number of samples to keep
n.keep <- 10000
# Number of jumps between samples to keep
n.gap <- 1
# Number of samples burned
n.burn <- 1000
# Required effective sample size
ta.ef.mu <- 1000
# Variable to host seed for pseudo-random numbers
seed <- 1970 
# Variable in case you want to print final values of the parameters
verbose <- FALSE


##-- Options for ordering parameter estimates for each cluster. 
##-- Label swiching problem.
sort.alpha.all <- TRUE  # With TRUE sort all estimates by the alpha parameter
sort.alpha <- FALSE # With TRUE order only the alpha parameter
sort.mu.all <- FALSE # With TRUE sort all estimates by the MU parameter

#######################################
## --  OPTIONS TO MAKE IMPUTATION OR  #
## --           ESTIMATION            #
#######################################

# Variable for imputation or estimation
imputation <- TRUE # If TRUE make imputation, otherwise estimation
# Variable to include auxiliary information (Si FALSE imputa con la media)
inf.aux <- TRUE #If imputation = TRUE and inf.aux = TRUE includes auxiliary variables, otherwise it imputes with the mean


############################################
## --END:   OPTIONS TO MAKE IMPUTATION     #
## --           OR  ESTIMATION             #
############################################


#######################################################
### --- END: 1. DEFINITION OF PROGRAM PARAMETERS      #
#######################################################
