##########################################################################################################
## Bayesain Markov Random Field (BMRF)
## function: BMRF()
## BMRF is a function to generate probabilistic estimation of network structure and partial correlation
##########################################################################################################
## Input for BMRF function
## 1. Dataset; 2. Est.Index; 3. Prior.Edge
##########################################################################################################
## 1. Dataset: a n by p dataset (n represents # of samples, p represents # of variables)
##########################################################################################################
## 2. Est.Index: a (p*(p-1)/2)*2 matrix which indicates the upper off-diagonal elements of precision matrix
## list the upper off-diagonal elements of precision matrix by "row"
## In BMRF, we use Est.Index to assign "Beta" to the corresponding element in precision matrix
## take p = 3 for example, we need to estimate 3 edges
## and the corresponding parameters are Beta1, Beta2, Beta3 in BMRF
## Est.Index is a 3*2 matrix with 1st row: (1,2); 2nd row:(1,3); 3rd row: (2,3)
## we will assign the (1,2) element to Beta1, (1,3) element to Beta2, and (2,3) element to Beta3
##########################################################################################################
## 3. Prior.Edge: a p*(p-1)/2 vector which contains prior knowledge for each possible edge in network
## each element in Prior.Edge is a binary value (0/1)
## Pior.Edge[i] = 1: edge i might have a higher probability to exist in network
## set prior prob for that edge exisiting to 0.8 (on avaerage)
## Pior.Edge[i] = 0: set prior probabilisity of P(edge = 1) = 0.5 on average
##########################################################################################################
##########################################################################################################
## The output of BMRF function
##########################################################################################################
## Because we use the "coda" package to store the posterior samples, 
## the output of the BMRF function is a "MCMC list" object.
## It contains "T" posterior samples for each parameter in the generating list object.
## "T" can be set by users.
## Users can also define which parameters they are interested in making an inference.
## The default setting is generating 5000 posterior samples for both "Beta" and "Gamma" in the BMRF model
##########################################################################################################


BMRF <- function(Dataset, Est.Index, Prior.Edge){

## standerdize Dataset by col (gene)
S.Data <- apply(Dataset, 2, function(input) (input-mean(input))/sd(input))

## check the input data is a matrix object
S.Data <- as.matrix(S.Data)

## N: # of samples
N <- nrow(S.Data)

## P: # of genes
P <- ncol(S.Data)

## Num: number of all possible edges in the network
Num <- choose(P, 2)

## Start the OpenBUGS code
Model <- function(){

#########################################################
## spike-and slab lasso prior
## slab distribution: DE(0, 2)
## spike distribution: DE(0, 20)
## Beta-Bernoulli conjugate prior for Gamma
## if Prior.Edge[i] == 1,
## set mean of Gamma = 1 to 0.8 (on avaerage)
## if Prior.Edge[i] == 0,
## set mean of Gamma = 1 to 0.5 (on avaerage)
#########################################################

for (k in 1:Num){

Gamma[k] ~ dbern(p[k])

ALPHA[k] <- (30*(equals(Prior.Edge[k], 1))) + 10

p[k] ~ dbeta(ALPHA[k], 10)

tauprior[k] <- (18*(equals(Gamma[k], 0))) + 2

Beta[k] ~ ddexp(0,tauprior[k])

}

#########################################################
## construct the Beta matrix for estimation
## Beta.Matrix is a p*p matrix
## each element is the corresponding Beta value
## Beta: strength of edge
## take p = 3 for example
## Beta.Matrix[1,2] = Beta[1]
## Beta.Matrix[1,3] = Beta[2]
## Beta.Matrix[2,3] = Beta[3]
## add symmetric constraint
## Beta.Matrix[i,j] = Beta.Matrix[j,i]
#########################################################

for (S in 1:Num){  
	Beta.Matrix[Est.Index[S,1], Est.Index[S,2]] <- Beta[S]
	Beta.Matrix[Est.Index[S,2], Est.Index[S,1]] <- Beta.Matrix[Est.Index[S,1], Est.Index[S,2]]
}

#########################################################
## set diagonal elements of Beta.Matrix to zero
## no self-loop in the network
#########################################################

for (Diag in 1:P){
	Beta.Matrix[Diag, Diag] <- 0
}

#########################################################
## CAR model
#########################################################

for (col in 1:P){

for (row in 1:N){

S.Data[row, col] ~ dnorm(mu[row, col], tau)

mu[row, col] <- inprod(S.Data[row, ], Beta.Matrix[,col])

}

}

#########################################################
## prior dist. for varaince in CAR model
## inverse Gamma conjugate prior for variance
## with mean approximate 1
#########################################################

tau ~ dgamma(8, 8)

} 

#########################################################
## set the working directory
## set the model.file used in OpenBUGS
#########################################################

Direction <- getwd()

OpenBUGS.file <- gsub(" ", "", paste(Direction,"/Model.odc"))

write.model(Model, con = OpenBUGS.file)

my.data <- list("S.Data", "N", "P", "Num", "Est.Index", "Prior.Edge")


inits <- function() {
list(Beta = rep(0, Num), Gamma = rep(0, Num))
}

params <- c("Beta", "Gamma")

#########################################################
## start sampler
## for each parameter, generated 5000 posterior samples
## thin = 10
## burn in = 5000
## iteration = 10000
## users can change the settings here!
#########################################################

out <- bugs(data = my.data, inits = inits, parameters.to.save = params, 
model.file = OpenBUGS.file, codaPkg = TRUE,
n.iter = 10000, n.chains = 1, n.burnin = 5000, n.thin = 10, debug = F)

out.coda <- read.bugs(out)

return(out.coda)

}





