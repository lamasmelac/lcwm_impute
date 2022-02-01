########################
##--PLOTS PLOTS PLOTS  #
########################

###--- Plot 1
###---TRACE PLOTS DE CLUSTERS OCUPADOS PARA ALPHA Y MU #
########################################################

##--Trace plots for alpha

mis.colores <- colorRampPalette(c("blue","brown","forestgreen","chocolate1","mediumvioletred","deepskyblue"))

par(mfrow=c(1,2))
j=1
plot(mcmc$par$alpha[,j], main=TeX('Trace-plots for $\\alpha$'), type="l", col=mis.colores(nz.K*2)[j],ylab="",ylim = c(0,1))
r=nz.K
for(j in 2:r){
  lines (mcmc$par$alpha[,j], type="l", ylab="", col=mis.colores(nz.K)[j], ylim = c(0,1))
}

##--Trace plots for mu
aux = range(mcmc$par$mu[,1:p*nz.K])

j=1
plot(mcmc$par$mu[,j], main=TeX('Trace-plots for $\\mu$'), type="l", col=mis.colores(2*nz.K)[j], ylab="",ylim = aux)
for(j in 2:(p*nz.K)){
  lines(mcmc$par$mu[,j], main=TeX('Trace-plots for $\\mu$'), type="l", col=mis.colores(p*nz.K)[j], ylab="")
}

dev.off()


###--- Plot 2
###---PAIR PLOT FOR OBSERVED AND IMPUTED DATA
#################################################


color_pairs = c("black",rojo_t)
pairs(x, pch = 19, upper.panel = NULL, col=color_pairs[R_ind.mis+1])
par(xpd=TRUE)
legend(x=0.7, y=1.10, legend ="observed", pch=19, col = "black", cex=1.0, bty = "n", horiz = FALSE)
legend(x=0.7, y=1.05, legend ="imputed", pch=19, col = rojo_t, cex=1.0, bty = "n", horiz = FALSE)

dev.off()

########################
##--    ESTIMATION     #
########################

##-- MAP iteration in the estimation case
##-- to calculate the KL divergence
#########################################

mcmc$MAP$alpha_MAP
mcmc$MAP$mu_MAP
mcmc$MAP$sig_MAP
