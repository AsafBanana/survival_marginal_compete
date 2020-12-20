source("competing_risks_model.R")

library(mvtnorm)
library(abind)
set.seed(323214)
N_QUESTION = 12
N_EVENTS = 3
#N_EVENTS = 2

N=2000
wb_k = 1
wb_g = 1
MAX_ITR = 150
P = 0.5
C = 3
beta = t(matrix(c(rep(c(1,1.5,0.5),4),rep(c(1.5,1,0.6),4),rep(c(0.75,1,0.5),4)),ncol=N_QUESTION))

r12 = 0.2
r13 = 0.3
r23 = 0.5

sigma = 1 
omega_varcov = matrix(c(sigma,r12,r13,r12,sigma,r23,r13,r23,sigma),ncol=3)


X = rbinom(N*N_QUESTION, 1, P)
dim(X) <- c(N,N_QUESTION,1)


base_rate = 1

simulate_survial_times <- function(beta,X,n=N,rv,base_rate=1) {
	proph = exp(t(matrix(rep(beta[,1],n),ncol=n)) * X) * exp(matrix(rep(t(rv),N_QUESTION),ncol=N_QUESTION))
	T = matrix(rep(0,N*N_QUESTION),nrow=N)
	for (i in seq(1,N)) {
		for (j in seq(1,N_QUESTION)) {
			T[i,j] = rexp(base_rate, rate = base_rate*proph[i,j])
		}
	}
	return(T)
}


sim_RV <- function(n=N,Sigma=omega_varcov) {
	return(mvtnorm::rmvnorm(n, rep(0,N_EVENTS), Sigma))
}


omega = sim_RV()
T1 = simulate_survial_times(matrix(beta[,1]),X[,,1],rv=omega[,1],base_rate=1)
T2 = simulate_survial_times(matrix(beta[,2]),X[,,1],rv=omega[,2],base_rate=1)
T3 = simulate_survial_times(matrix(beta[,3]),X[,,1],rv=omega[,3],base_rate=1)

T = pmin(T1,T2,T3,C)
delta1 = (T == T1)
delta2 = (T == T2)
delta3 = (T == T3)
delta= abind(delta1,delta2,delta3,along=3)


source("competing_risks_model.R")
comprsk = CompetingRisksModel()
comprsk$fit(T,delta,X)

