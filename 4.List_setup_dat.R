################################################################################
### ---4.1. CREATION OF LIST 'setup' WITH VARIABLES FOR NUMBER OF ITERATIONS   #
################################################################################

setup <- NULL

## Seed setting for pseudo-random numbers
if ( is.null(seed) ) { seed <- ceiling(runif(1, 0, 10000)) }
set.seed(seed)
setup$seed <- seed

# Number of simulations
n.sim   <- n.burn + n.gap * n.keep

# Save program variables to list: "setup"
setup$n.sim  <- n.sim
setup$n.keep <- n.keep
setup$n.gap  <- n.gap
setup$n.burn <- n.burn

## Determination of iterations keep
i.keep  <- seq(from = (n.burn+n.gap), to = n.sim, by = n.gap)
i.count <- rep(0, times = n.sim)
i.count[i.keep] = 1:n.keep

#################################################################
## --IN CASE OF IMPUTATION: VARIABLES FOR NUMBER OF ITERATIONS #
################################################################

if(imputation){
  ## Determination of which imputed data matrices to keep
  ## (last "imp.keep" in jumps of "imp.gap")
  
  # last imputed matrices to keep
  imp.keep  <- 1
  # jump between those last matrices to keep
  imp.gap   <- 1
  i.imp.keep <- seq(from = (n.sim-(imp.keep-1)*imp.gap), to = n.sim, by = imp.gap)
  imp.count <- rep(0, times = n.sim)
  imp.count[i.imp.keep] = 1:imp.keep
  
  
  # Save program variables to list: "setup"
  setup$imp.keep <- imp.keep
  setup$imp.gap <- imp.gap
}

######################################################################
## --END: IN CASE OF IMPUTATION: VARIABLES FOR NUMBER OF ITERATIONS  #
######################################################################

#################################################################################
### ---4.1. END: CREATION OF LIST 'setup' WITH VARIABLES FOR NUMBER OF ITERATIONS #
#################################################################################


###################################################################
### ---4.2. CREATION OF 'dat' LIST WITH DATABASE CHARACTERISTICS #
##################################################################

dat <- NULL
dat$x <- x
dat$n <- n
dat$k <- k
dat$p <- p

dat$r <- r

#######################################################################
### ---4.2. END: CREATION OF 'dat' LIST WITH DATABASE CHARACTERISTICS #
######################################################################
