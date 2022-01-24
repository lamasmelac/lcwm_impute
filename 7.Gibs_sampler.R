#########################################
#----------MCMC GIBBS SAMPLER start-----#
#----------MCMC GIBBS SAMPLER start-----#
#----------MCMC GIBBS SAMPLER start-----#
#########################################

iter <- 1
efsi <- 1

t <- proc.time()
while ((i.count[iter]<ta.ef.mu || efsi<ta.ef.mu) & iter<n.sim) {
  
  cat("\n Simulation ",iter,"of ", n.sim, "with effectisize", efsi) 
  
  ## --- 5. Cluster sampling from the multinomial distribution    #
  #################################################################
  
  temp=sample.zi(n,k,x,alpha,mu,sig)
  
  z <- matrix(0,nrow = n,ncol = k)
  for(i in 1:n){
    z[i,temp[i]] <- 1 
  }
  ## Variable with number of observations in each cluster
  n.mix <- as.vector(table(factor(temp,levels = 1:k)))
  
  ######################################################################
  ## --- 5. End: Cluster sampling from the multinomial distribution    #
  
  
  ## --- 4. eta parameter sampling.         #
  ##########################################
  log.alphaG <- ifelse(alpha[k]==0,-10^5,log(alpha[k]))
  eta <- rgamma(1,shape=a_eta+k-1,rate=b_eta-log.alphaG)
  ###############################################
  ## --- 4. End: eta parameter sampling.        #
  
  
  ## --- 2. Sampling of Mixture Probabilities (alpha) from Stick Breaking    #
  ###########################################################################
  alpha <- stick(n.mix,eta)
  ## --- Sort parameters based on alpha without ordering the other parameters
  if(sort.alpha){
    alpha <- sort(alpha, decreasing=T)
  }
  
  ###############################################################################
  ## --- 2. End: Sampling of Mixture Probabilities (alpha) from Stick Breaking #
  
  
  ## --- 3. Element sampling of the Phi matrix.     #
  #######################################################
  
  # start create inverses of covariance matrices
  # expr. 2.13 Paiva(2014)
  for(j in 1:k){
    sig_inv[,,j] <- solve(sig[,,j])
  }
  #end create inverses of covariance matrices
  
  #start updating phi_j components for phi matrix
  # expr. 2.13 Paiva(2014)
  for(i in 1:p){ #for each dimension p
    phi[i] <- rgamma(1,shape=a_phi+k*(p+1)/2,rate=b_phi+0.5*sum(sig_inv[i,i,]))
  }
  #}
  phi_matrix <- diag(phi,nrow = p)
  #end updating phi_j components for phi matrix
  
  ############################################################
  ## --- 3. End: Element sampling of the Phi matrix.     #
  
  
  ## --- 1. Sampling means and covariances of mixtures  #
  #######################################################
  for(j in 1:k){
    mu_n_k <- mu.0
    S <- matrix(0,nrow = p, ncol = p)
    cov.mu <- matrix(0,nrow = p, ncol = p) #new
    
    if(n.mix[j]>0){
      #Construction of the vector of means of the data
      #Paragraph before expr.2.10 Paiva(2014)
      Y_bar <- apply(matrix(x[z[,j] == 1,],nrow = n.mix[j]),2,mean)
      #Construction of matrix of sums of squares of the data
      #Paragraph before expr.2.10 Paiva(2014)
      S <- matrix(apply(matrix(x[z[,j] == 1,],ncol=p),MARGIN = 1,FUN = "-",Y_bar),nrow = p)%*%t(matrix(apply(matrix(x[z[,j] == 1,],ncol=p),MARGIN = 1,FUN = "-",Y_bar),nrow = p))
      
      #Inverse-Wishart distribution matrix update
      #Paragraph after expr.2.11 Paiva(2014)
      cov.mu <- (Y_bar - mu.0)%*%t(Y_bar - mu.0)/(1/h+1/n.mix[j])
      #Paragraph after expr.2.11 Paiva(2014)
      mu_n_k <- (h*mu.0+n.mix[j]*Y_bar)/(h+n.mix[j])
    }
    
    #Inverse-Wishart distribution matrix update
    #Paragraph after expr.2.11 Paiva(2014)
    phi_matrix_k <- phi_matrix + S + cov.mu
    #Inverse-Wishart variance and covariance matrix update
    # expr.2.10 Paiva(2014)
    sig[,,j] <- solve(rwish(p+1+n.mix[j],solve(phi_matrix_k)))
    
    #PParagraph after expr.2.11 Paiva(2014)
    Sigma_n_k <- (1/(h+n.mix[j]))*sig[,,j]
    #Normal mean vector update
    # expr.2.11 Paiva(2014)
    mu[j,] <- mvrnorm(1, mu_n_k, Sigma_n_k)
  }
  
  ############################################################
  ## --- 1. End: Sampling means and covariances of mixtures  #
  
  
  
  ## --- 6. Imputation process  #
  ###############################
  
  ### -If imputation = TRUE, the imputation process is carried out
  
  if(imputation){
    
    ## IMPUTATION START #
    #####################
    if(inf.aux){
      # - If inf.aux = TRUE, the imputation is done using the 
      #information from the auxiliary variables, otherwise the 
      #imputation is done through the "mean method"
      x.aux = matrix(data=x[(m+1):n,(r+1):ncol(x)],nrow = n-m,ncol = p-r) #es fijo puede ir en inicio del codigo
      mu.aux = matrix(data=mu[,(r+1):ncol(x)],nrow = k,ncol = p-r)
      sig.aux = array(data=sig[(r+1):ncol(x),(r+1):ncol(x),],dim = c(p-r,p-r,k))
      temp.imp=sample.zi(n-m,k,x.aux,alpha,mu.aux,sig.aux)
      
      zx=matrix(0,nrow = n-m,ncol = k)
      for(i in 1:(n-m)){
        zx[i,temp.imp[i]] <- 1 
      }
      #Matrix to save the sizes of the components in each iteration
      n.mix.imp <- as.vector(table(factor(temp.imp,levels = 1:k)))
      
      for(i in 1:(n-m)){
        j <- temp.imp[i]
        me <- mu[j,1:r] + (sig[(1:r),(r+1):p,j]%*%solve(sig[(r+1):p,(r+1):p,j]))%*%(x[m+i,(r+1):p] - mu[j,(r+1):p])
        si <- sig[(1:r),(1:r),j] - (sig[(1:r),(r+1):p,j]%*%solve(sig[(r+1):p,(r+1):p,j]))%*%sig[(r+1):p,(1:r),j]
        x[(m+i),1:r] <- mvrnorm(1, mu = me, Sigma = si)
      }
    }else{
      temp.imp=sample(1:k,size = n-m, replace = TRUE, prob = alpha)
      zx=matrix(0,nrow = n-m,ncol = k)
      #Matrix to save the sizes of the components in each iteration
      n.mix.imp <- as.vector(table(factor(temp.imp,levels = 1:k)))
      
      for(i in 1:(n-m)){
        zx[i,temp.imp[i]] <- 1 
      }
      for(i in 1:(n-m)){
        j <- temp.imp[i]
        me <- mu[j,1:r]
        si <- sig[1:r,1:r,j]
        x[(m+i),1:r] <- mvrnorm(1, mu = me, Sigma = si) 
      }
    }
    
    ##################  
    ## END IMPUTATION#
    
    
    ## -Save the last imputed databases   #
    ######################################
    if(any(iter==i.imp.keep)){
      # last "imp.keep" imputed databases to keep
      bases_imp[,,imp.count[iter]] <- x
      z_bases_imp[,,imp.count[iter]] <- z
    }
    ###########################################
    ## -End: Save the last imputed databases #
    
  }
  
  #######################################
  ## --- 6. End: Imputation process  #
  
  ## --- Sort all parameters based on: alpha    #
  ###############################################

  if(sort.alpha.all){
    Order  <- order(alpha,decreasing = T)
    alpha  <- alpha[Order]
    mu     <- matrix(mu[Order,], nrow = k, ncol = p)
    sig    <- array(sig[,,Order],dim=c(p,p,k))
    z      <- matrix(z[,Order], nrow = n, ncol = k)
    n.mix  <- n.mix[Order] #new
    if(imputation){
      n.mix.imp <- n.mix.imp[Order]
    }
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
    n.mix  <- n.mix[Order] #Nuevo
  }
  
  ###############################################
  ##-- End: Sort all parameters based on: alpha #
  
  
  
  ## --- 7. Save estimation process results #
  ###########################################
  
  if (any(iter == i.keep)){
    ## Increment of the z matrix
    z.matrix  <- z.matrix + z
    
    ## Storing the mixture parameters
    alpha.matrix[i.count[iter],]  <- alpha
    mu.matrix[i.count[iter],]     <- c(t(mu))
    sigma2.matrix[i.count[iter],] <- c(sig)
    Phi.matrix[i.count[iter],] <- phi  #new
    n.mix.matrix[i.count[iter],]  <- n.mix
    if (imputation) {
      n.mix.imp.matrix[i.count[iter],]  <- n.mix.imp
    }
    ## Print parameters if verbose = TRUE
    if (verbose) { print (c(i.count[iter], alpha, mu, sigma2)) }
  }
  
  #################################################
  ## --- 7. End: Save estimation process results #
  
  
  
  
  
  ## --- 8. Log likelihood computation and MAP iteration #
  #######################################################
  
  if (any(iter == i.keep)) {
    
    temp.1 <- matrix(data = 0, nrow = n, ncol = k)
    temp.2 <- rep(0, k)
    
    for(j in 1:k){
      # part 1: log likelihood
      temp.1[,j] = alpha[j]*dmvnorm(x,mean = mu[j,],sigma = as.matrix(sig[,,j]),log=FALSE)
      # part 2: priors for mu and InvSig
      temp.2[j] = dmvnorm(matrix(mu[j,],1,p),mu.0,(1/h)*as.matrix(sig[,,j]),log=TRUE) +
        log(dwish(solve(sig[,,j]),p+1,solve(phi_matrix)))
    }
    loglik.1 = sum(log(apply(temp.1,MAR=1,FUN=sum)))
    loglik.2 = sum(temp.2)
    
    # part 3: prior for Phi
    loglik.3 = sum(log(dgamma(phi,shape=a_phi,rate=b_phi)))
    
    # part 4: prior for Alpha
    
    part.sum <- numeric(k-1)
    for(j in 1:(k-1)){
      part.sum[j] = sum(alpha[j:k]) }
    
    log.part.sum <- numeric(k-1)
    for(j in 1:(k-1)){
      log.part.sum[j]<- ifelse(part.sum[j]==0,-10^5,log(part.sum[j]))
    }
    
    loglik.4 = (k-1)*log(eta) + (k-1)*log.alphaG - sum(log.part.sum)  
    
    # part 5: prior for eta
    loglik.5 = dgamma(eta, shape=a_eta, rate=b_eta, log=TRUE)
    
    loglik=loglik.1+loglik.2+loglik.3+loglik.4+loglik.5
    
    loglik_vector[i.count[iter]] <- loglik
    
    if(i.count[iter]>1){
      efsi <- effectiveSize(loglik_vector)
      efectsize[i.count[iter]] <- efsi
    }
    
    if(loglik>loglik_MAP){
      alpha_MAP  <- alpha
      mu_MAP     <- mu
      sig_MAP    <- sig
      loglik_MAP <- loglik
      z_MAP      <- z
      if(imputation){
        zx_MAP     <- zx  
      }
      x_MAP      <- x
    }
  }
  
  ##############################################################
  ## --- 8. End: Log likelihood computation and MAP iteration #
  
  iter <- iter + 1
} #end of MCMC
proc.time() - t



###########################
### ---- END GIBS SAMPLER #
### ---- END GIBS SAMPLER #
### ---- END GIBS SAMPLER #
###########################