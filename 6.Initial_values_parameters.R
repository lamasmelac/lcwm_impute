### ---6. INITIAL VALUES FOR PARAMETERS AND HYPERPARAMETERS #
#############################################################

#####################################
## --- Values for hyperparameters #
#####################################

h <- 1
a_phi <- 0.25
b_phi <- 0.25
a_eta <- 0.25
b_eta <- 0.25

mu.0 <- apply(x,2,mean)
#mu.0 <- rep(0,p)
########################################
# Initial value: matrix parameter phi  #
########################################
phi <- NA
for(i in 1:p){ #for each dimension p
  phi[i] <- rgamma(1,shape=a_phi,rate=b_phi)
}
phi_matrix <- diag(phi,nrow = p)

###############################
# Initial value: eta parameter#
###############################
eta <- rgamma(1,shape=a_eta,rate=b_eta)

############################################################
# Initial values: alpha_0, mu_0, z_0, n.mix_0 parameters   #
############################################################

## Latent indicator variable initialization
## Binary matrix with n rows, k columns
## Each column has a one (1) and k-1 zeros (0)
z <- cbind(rep(1,n), matrix(0, nrow=n, ncol=k-1))
for (i in 1:n) { z[i,] <- sample(z[i,]) }

## Clusters size initialization
n.mix <- apply(z,2,sum) 

## Proportion of mixtures: alpha_0
alpha <- apply(z, 2, mean) #istarting with stick-breaking remove

if(p>1){
  for(j in 1:k){
    sig[,,j]=cov(x[z[,j]==1,])
    mu[j,]=apply(x[z[,j]==1,],MARGIN = 2,FUN = mean)
  }
}else{
  for(j in 1:k){
    sig[,,j]=var(x[z[,j]==1,])
    mu[j,]=mean(x[z[,j]==1,])
  }
}


## --- Sort all parameters based on: alpha    #
###############################################
if(sort.alpha.all){
  Order  <- order(alpha,decreasing = T)
  alpha  <- alpha[Order]
  mu     <- matrix(mu[Order,], nrow = k, ncol = p)
  sig    <- array(sig[,,Order],dim=c(p,p,k))
  z      <- matrix(z[,Order], nrow = n, ncol = k)
  n.mix  <- n.mix[Order] #new
}


## --- Sort all parameters based on: mu    #
###############################################
if(sort.mu.all){
  
  d_mu = rep(0,k)
  for(i in 1:k){
    d_mu[i]=sqrt(sum(mu[i,]^2))
  }
  
  Order  <- order(d_mu,decreasing=T)
  alpha  <- alpha[Order]
  mu     <- matrix(mu[Order,], nrow = k, ncol = p)
  sig    <- array(sig[,,Order],dim=c(p,p,k))
  z      <- matrix(z[,Order], nrow = n, ncol = k)
  n.mix  <- n.mix[Order] #new
}

##################################################################
### ---6. END: INITIAL VALUES FOR PARAMETERS AND HYPERPARAMETERS #