################################################################################
## 1) Prepare Data:           
##
##      First make sure your data are arranged as a matrix X (say)
##      with rows representing "cases" and columns variables, 
##      and represent missing values as NA.
##        
##        The code delivers a database called "x" with the first
##        variables with missing data and the last with fully observed 
##        information. It also establishes values for the variables 
##        n number of cases and r number of variables. The variable 
##        d number of variables with missing data and the variable p 
##        number of variables with complete information.
################################################################################

## Load required packages
require(norm)
require(dplyr)


require(mvtnorm)
require(MASS)
require(MCMCpack)
require(car)
require(latex2exp)
require(scatterplot3d)
require(Hmisc)



################################################################################

## Ask for input data
data.file = readline("Enter file name of input data: ")


## Load data
datos = read.csv(data.file, sep=";")
datos_matrix = as.matrix(datos)

## Data dimensions
n = dim(datos_matrix)[1]  # number of observations
p = dim(datos_matrix)[2]  # number of variables

## Permutations are made between rows, also between columns in such 
## a way that the database that is delivered has the first columns 
## with missing data and the last ones with complete information. 
## The missing data are in the last rows.
## The matrix "x" is returned with the data thus organized
s = prelim.norm(datos_matrix)
r0 <- s$ro
r1 <- 1:n
perm <- cbind(r0,r1)
perm <- perm[order(perm[,"r0"]),]
x <- as.matrix(datos[perm[,2],sort(s$nmis, index.return=TRUE, decreasing = TRUE)$ix])
## matrix of original data sorted according to NA
x_NA <- x


s = prelim.norm(x)

## Number of input variables
d = sum(s$nmis==0) 
## Number of output variables
r = p-d 
## Number of individuals with missing data
m = sum(is.na(x[,1])) 
## Indicator vector of individuals with missing information 
## (0 if information is complete, 1 if information is missing)
R_ind.mis = c(rep(0,n-m),rep(1,m))

## Imputation of NAs by the average of values observed in each variable
## Initial imputations
for(i in 1:p){
  x[is.na(x[,i]),i] <- mean(x[,i],na.rm=TRUE)
}

## matrix of original data imputing NAs with the mean of observed data
x_imp.ini <- x
