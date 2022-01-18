##################
#RESULTS VECTORS#
#################

## Three lists are constructed and two more ("setup",""dat) are appended 
##    to the "mcmc" result list:
###  par: list to save iterations of the estimation of each one of the parameters and other results
###  MAP: list to save estimates in the MAP iteration and other results in this iteration
### summary: list to save posterior means of the parameters




## Convert z.matrix to ratios
z.postmean <- z.matrix / length(efectsize)


## --- List "par" to store arrays of iterations
par <- NULL
par$alpha  <- alpha.matrix[1:length(efectsize),]
par$mu     <- mu.matrix[1:length(efectsize),]
par$sigma2 <- sigma2.matrix[1:length(efectsize),]
par$z.postmean <- z.postmean
if(imputation){
  par$bases_imp <- bases_imp
  par$z_bases_imp <- z_bases_imp  
}
par$loglik_vector <- loglik_vector


## --- "MAP" list to store MAP iteration results

MAP <- NULL
MAP$alpha_MAP <- alpha_MAP
MAP$mu_MAP <- mu_MAP
MAP$sig_MAP <- sig_MAP
MAP$loglik_MAP <- loglik_MAP
MAP$z_MAP <- z_MAP
if(imputation){
  MAP$zx_MAP <- zx_MAP
  MAP$x_MAP <- x_MAP 
}



################################
## --- SUMMARY OF RESULTS   ---#
################################

## Names for parameter vectors

parnames.alpha <- paste(rep("alpha.",each=k),1:k,sep="")
parnames.mu <- paste(rep("mu.",each=k*p), paste(rep(1:k,each=p), rep(1:p,times=k), sep=""), sep="")
parnames.sigma <- paste(rep("sigma.", each=k*p*p), paste(rep(1:k,each=p*p), paste(rep(1:p, each=p), rep(1:p, times=p), sep=""), sep=""), sep="")
parnames <- c(parnames.alpha, parnames.mu, parnames.sigma)


##############################################
## Posterior means for alpha, mu and sigma   #
##############################################
postmean <- data.frame(t(c(apply(par$alpha,  2, mean),
                           apply(par$mu,     2, mean),
                           apply(par$sigma2, 2, mean))))
names(postmean)     <- parnames
row.names(postmean) <- "posterior.mean"



##########################################
## 95% credibility intervals and median #
#########################################
cred <- data.frame(cbind(apply(par$alpha,  2, quantile, prob=c(0.025, 0.5, 0.975)),
                         apply(par$mu,     2, quantile, prob=c(0.025, 0.5, 0.975)),
                         apply(par$sigma2, 2, quantile, prob=c(0.025, 0.5, 0.975))))
colnames(cred) <- parnames


## Creation of the list "summary" for storage of summary objects
summary <- NULL
summary$posterior.means <- postmean
summary$credible.intervals <- cred


## Creation of the list "mcmc" to store lists of results
mcmc <- NULL
mcmc$setup   <- setup
mcmc$dat     <- dat
mcmc$par     <- par
mcmc$summary <- summary
mcmc$MAP <- MAP


##-- Number of occupied clusters 
ind.pi = apply(n.mix.matrix[1:length(efectsize),], MAR=2, FUN=median)>0
nz.K = sum(ind.pi)
