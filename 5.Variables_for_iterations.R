#################################################
### ---5.1 CREATING VARIABLES TO SAVE ITERATIONS #
#################################################

###############################################
## Variables to save iterations in estimation #
###############################################

## --- Creation of a set of matrices and vectors to store the information

#each line of the matrix contains the mixing probabilities of the k clusters
alpha.matrix  <- matrix(nrow = n.keep, ncol = k)
#each line of the matrix has the k vectors of cluster means each of dimension p
mu.matrix     <- matrix(nrow = n.keep, ncol = k*p)
#each line of the matrix has the entries of the k cluster covariate matrices
#each size pxp
sigma2.matrix <- matrix(nrow = n.keep, ncol = k*p*p)
#Indicator matrix of the cluster to which the observations belong
z.matrix      <- matrix(nrow = n, ncol = k)
#Matrix to store the diagonal components of the Phi matrix
Phi.matrix <- matrix(nrow = n.keep, ncol = p) #It's not necesary
#Matrix to save the sizes of the components in each iteration
n.mix.matrix <- matrix(nrow = n.keep, ncol = k)


#Matrix to save the sizes of the components in each iteration in case of imputation
if(imputation){
  n.mix.imp.matrix <- matrix(nrow = n.keep, ncol = k)
}


################################################
## Variables to save log-likelihood iterations #
################################################
#Vector for posterior log-likelihood values
loglik_vector <- NA
#Vector for effective sample sizes
efectsize <- NA


#######################################################
### ---5.1 END: CREATING VARIABLES TO SAVE ITERATIONS  #
#######################################################


###############################################################
### ---5.2 CREATION OF VARIABLES TO SAVE TEMPORARY ITERATIONS #
##############################################################

###################################
### --For estimation              #
###################################

###########################################################################
## Variables to save samples of the covariances of the clusters: Wishart  #
###########################################################################

#variable to save the inverses of the covariance matrices
sig_inv <- array(dim = c(p,p,k))


#####################################################################
## Variables to store iterations for mean 'mu' and covariance 'sig' #
#####################################################################

mu  <- matrix(nrow = k,ncol = p)
sig <- array(dim = c(p,p,k))


###################################
### --For imputation              #
###################################

if(imputation){
  
  #array to save imputed databases
  bases_imp <- array(dim=c(n,p,imp.keep))
  #array to save indicator of component occupied by observation
  z_bases_imp <- array(dim=c(n,k,imp.keep))
  
  
  ## Latent indicator variable initialization for missing data
  ## Binary matrix with n-m rows and k columns
  ## Each column contains a single one and k-1 zeros
  ## Cada columna tiene 1 uno y k-1 ceros. 
  ## (m: number of obs. with missing data)
  ## Indicator matrix of component assigned for each imputed observation
  zx <- cbind(rep(1,n-m), matrix(0, nrow=n-m, ncol=k-1))
}


##########################################
### --For MAP-iteration                  #
##########################################
## Proportion of MAP mixtures
alpha_MAP <- rep(0, times = k)
## mean and covariance MAP component
mu_MAP  <- matrix(nrow = k, ncol = p)
sig_MAP <- array(dim = c(p,p,k)) 
## variables for MAP likelihood computation
loglik_MAP <- -Inf
loglik.1=loglik.2=loglik.3=loglik.4=loglik.5=0


####################################################################
### ---5.2 FIN: CREATION OF VARIABLES TO SAVE TEMPORARY ITERATIONS #
