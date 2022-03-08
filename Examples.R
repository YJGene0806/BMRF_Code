#############################################
## install packages and load into R first
#############################################
## install.packages("R2OpenBUGS")
## install.packages("MASS")
## install.packages("coda")
## install.packages("gdata")
## install.packages("ograph")
library(R2OpenBUGS)
library(MASS)
library(coda)
library(gdata)
library(igraph)
#############################################

#############################################
## Preload BMRF_function.R
#############################################
source("BMRF_Function.R")
#############################################
## Preload Performance_Criteria.R
#############################################
source("Performance_Criteria.R")
#############################################

#############################################
## Example
## a star network structure
## p = 5
## n = 50
#############################################
#############################################
## Generate data from MVN
#############################################

Sample.Size <- 50

N.Genes <- 5

Total.Edges <- choose(N.Genes, 2)

## sett true precision matrix
PREC.M <- matrix(
c(1, -0.25, -0.25, -0.25, -0.25,  
-0.25, 1, 0, 0, 0, 
-0.25, 0, 1, 0, 0, 
-0.25, 0, 0, 1, 0,
-0.25, 0, 0, 0, 1), 5, 5, byrow=T)

COR.M <- cov2cor(solve(PREC.M))

Data.1 <- mvrnorm(Sample.Size, rep(0, N.Genes), COR.M)

True.Network <- apply(PREC.M, 2, function(input) ifelse(input != 0, 1, 0))

diag(True.Network) <- 0

True.Edge <- upperTriangle(True.Network, byrow = T)

G <- graph_from_adjacency_matrix(True.Network, mode = "undirected", diag = F)

## visualize true network
plot(G)



#############################################
## setting the input "Est.Edge" 
#############################################
#############################################
## XX() is a function
## which can transfer the number of Beta
## to the upper off-diagonal elements in the precision matrix
#############################################

XX <- function(P, cell) {

N.Cell.each.row <- rev(1:(P-1))

Cum.Sum <- c(0, cumsum(N.Cell.each.row))

CUM.Index <- min(which(Cum.Sum >= cell))

Row.Index <- (CUM.Index - 1)

Col.Index <- Row.Index + (cell - Cum.Sum[CUM.Index-1])

return(c(Row.Index, Col.Index))

}


#############################################
## Edge.Number is a 10*2 matrix
## indicates the upper off-diagonal elements in precision matrix
#############################################

Edge.Number <- matrix(0, Total.Edges, 2)

for (i in 1:Total.Edges){

Edge.Number[i, ] <- XX(N.Genes, i)

}


#Edge.Number


#############################################
## set non-informative prior
#############################################

Prior.Edge <- rep(0, Total.Edges)

#############################################
## run BMRF
#############################################

star <- Sys.time()

Output <- BMRF(Data.1, Edge.Number, Prior.Edge)

as.numeric(Sys.time()-star)


#############################################
## outputs
#############################################

## posterior samples
Post.Samp <- Output[[1]]

## generate 5000 posterior samples
## 10 Beta + 10 Gamma + 1 deviance
dim(Post.Samp)

colnames(Post.Samp)

## posterior mean of each variable
summary(Output)

## estimated prob(edge)
Prob.edge <- summary(Output)[[1]][(Total.Edges+1):(2*Total.Edges), 1]

## posterior mean of Beta
Est.ParCor <- summary(Output)[[1]][1:Total.Edges, 1]

## select edges 
## prob(edg) >= 0.5 or not
Est.Edge <- ifelse(Prob.edge >= 0.5, 1, 0)

Est.Network <- matrix(0, N.Genes, N.Genes)

upperTriangle(Est.Network, byrow = T) <- Est.Edge

Est.Network <- Est.Network + t(Est.Network)

Est.G <- graph_from_adjacency_matrix(Est.Network, mode = "undirected", diag = F)

## visualize the estimated network
plot(Est.G)


#############################################
## performance
#############################################

ACC.Criterion(rbind(Est.Edge, True.Edge))


####################################
## set informative prior
####################################

Prior.Edge2 <- True.Edge

####################################
## run BMRF
####################################

star <- Sys.time()

Output <- BMRF(Data.1, Edge.Number, Prior.Edge2)

as.numeric(Sys.time()-star)


####################################
## outputs
####################################

## summary from OpenBUGS
summary(Output)

## estimated prob(edge)
Prob.edge <- summary(Output)[[1]][(Total.Edges+1):(2*Total.Edges), 1]

## posterior mean of Beta
Est.ParCor <- summary(Output)[[1]][1:Total.Edges, 1]

## select edges 
## prob(edg) >= 0.5
Est.Edge <- ifelse(Prob.edge >= 0.5, 1, 0)


Est.Network <- matrix(0, N.Genes, N.Genes)

upperTriangle(Est.Network, byrow = T) <- Est.Edge

Est.Network <- Est.Network + t(Est.Network)

Est.G <- graph_from_adjacency_matrix(Est.Network, mode = "undirected", diag = F)

## visualize the estimated network
plot(Est.G)


####################################
## performance
####################################

ACC.Criterion(rbind(Est.Edge, True.Edge))






