#########################################################
##----AUXILIARY FUNCTIONS REQUIRED TO RUN THE PROGRAM  #
########################################################


##########################################################################
## --- 1. Sampling Mixing Probabilities (alpha) from Stick Breaking     #
#########################################################################


stick = function(counts,alpha){
  
  # counts:  vector with the current number of observations on each component
  # alpha:   current sample of the dispersion parameter
  
  K = length(counts)  # number of components
  n = sum(counts)     # total number of observations
  
  ncur = csumlog = 0           # to store some cumulative quantities
  logprob = prob = numeric(K)  # to return the probabilities by component
  
  for(h in 1:(K-1)){
    ncur = ncur + counts[h]              # cumulative sum of counts
    a = 1 + counts[h]                    # 1st parameter of the beta distribution for the stick-breaking weights
    b = alpha + n - ncur                 # 2nd parameter of the beta distribution for the stick-breaking weights
    lga = log(rgamma(1,shape=a,rate=1))  # log(X) where X~Gamma(a,1)
    lgb = log(rgamma(1,shape=b,rate=1))  # log(Y) where Y~Gamma(b,1)
    lgsum = log(exp(lga)+exp(lgb))       # log(X+Y)
    
    logprob[h] = lga - lgsum + csumlog   # log(pi_h) = log(v_h) + sum_{g<h} log(1-v_g)   where log(v_h)=log(X)-log(X+Y)
    csumlog = csumlog + lgb - lgsum      # sum_{g<=h} log(1-v_g)   where log(1-v_g)=log(Y)-log(X+Y)
    prob[h] = exp(logprob[h])            # exp(log(pi_h))
  }
  
  logprob[K] = csumlog                   # pi_K = v_K * prod_{g<K}(1-v_h) and since v_K=1, log(pi_K)=sum_{g<K} log(1-v_g)
  prob[K] = exp(logprob[K])              # exp(log(pi_K))
  
  return(prob)        # return the vector of probabilities
  
}

#################################################################
## --- 2. Sampling the cluster from the multinomial distribution#
#################################################################

sample.zi = function(n,K,Y,cur.Pi,cur.Mu,cur.S){
  
  # n:      current sample size
  # K:      number of components
  # Y:      current set of observations
  # cur.Pi: current sample of Pi (Pi[] size: K )
  # cur.Mu: current sample of Mu (Mu[,] size: K x p )
  # cur.S: current sample of Sig (Sig[,,]) (size: p x p x K)
  
  # compute the updated multinomial probabilities for each obs and each group
  log.pi = matrix(0, nrow = n, ncol = K) # matrix (n x K) of the log of the numerator
  for(k in 1:K){
    # log of the probability by component
    #     log.pi[,k] = log(cur.Pi[k]) + dmvn_prec(Y,cur.Mu[k,],cur.IS[,,k],iflog=T)
    log.pi[,k] = log(cur.Pi[k]) + dmvnorm(Y,cur.Mu[k,],as.matrix(cur.S[,,k]),log=T)
  }
  # subtract the max by row to take care of extreme cases
  log.pi = sweep(log.pi, MAR=1, apply(log.pi,MAR=1,FUN=max), FUN="-")
  # normalizing by row and taking exponential
  pi.norm = sweep(exp(log.pi), MAR=1, rowSums(exp(log.pi)), FUN="/")
  
  # sample the zi's (avoid a loop over i=1,...,n)
  cum.pi = t(apply(pi.norm,MAR=1,FUN=cumsum))      # cumulative probabilities
  Z.out = max.col( (cum.pi - runif(n))>0, "first") # select which category
  # a random number falls into
  return(as.vector(Z.out))
  
  
  ###################################################
  ##--Function to generate colors with transparency #
  ###################################################
  
  
  t_col <- function(color, percent = 50, name = NULL) {
    #      color = color name
    #    percent = % transparency
    #       name = an optional name for the color
    
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    
    ## Save the color
    invisible(t.col)
  }
  
  rojo_t = t_col("red",50)
}

