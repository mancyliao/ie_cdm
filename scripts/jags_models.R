library(R2jags)
library(sirt)

# set dir -----------------------------------------------------------------
work_dir <- getwd()
data_dir <- paste(work_dir,"/data",sep = "")
output_dir <- paste(work_dir,"/output",sep = "")
dir.create(output_dir, recursive = T, showWarnings = F)
# model_out_dir <- paste(work_dir,"/Model_Output",sep = "")
# result_out_dir <- paste(work_dir,"/Result_Output",sep = "")
# code_dir <- paste(work_dir,"/JAGS_code",sep = "")



# load data--------------
response_data <- read.csv(paste(data_dir,"/data_released_to_use_medium_more_country.csv",sep=""))
response_data <- as.matrix(response_data)
qmatrix <- read.csv(paste(data_dir,"/q_matrix_medium.csv",sep=""))
item_feature <- read.csv(paste(data_dir,"/item_data_to_use_unstd.csv",sep=""))
item_feature_use <- as.matrix(item_feature[,setdiff(names(item_feature),"Item_ID")])


# metadata ----------------------------------------------------------------

I=ncol(response_data) # number of items
N=nrow(response_data) # number of persons
K=ncol(qmatrix) # number of attributes
M=ncol(item_feature_use) # number of features

#MCMC chain prms
n.iter = 20000
n.burnin = 10000
n.thin=1
use_parallel = T ## wheter or not use parallel to run JAGS

# model to run ------------------------------------------------------------
## values of `condition` could be one of: "HO_DINA_g", "HO_DINA_s", "HO_DINA_g_resid","HO_DINA_s_resid"
condition <- "HO_DINA_g"


# JAGS model (IE_HO_DINAs)-------------------------------------------

if(condition=="HO_DINA_g"){
  bayes.mod <- function(){
    #Standardize x's
    for (m in 1:M){
      for (i in 1:I){
        Z[i,m]<-(X[i,m]-mean(X[,m]))/sd(X[,m])}}
    
    for (n in 1:N) {
      for (i in 1:I) {
        for (k in 1:K) {w[n, i, k] <- pow(alpha[n, k], Q[i, k])}
        eta[n, i] <- prod(w[n, i, 1:K])
        p[n, i] <- pow((1 - s[i]), eta[n, i]) * pow(g[i], (1 - eta[n, i]))
        Y[n, i] ~ dbern(p[n, i])}}
    for(n in 1:N){
      for(k in 1:K){
        logit(prob.a[n, k]) <- xi[k] * theta[n] - beta[k]
        alpha[n, k] ~ dbern(prob.a[n, k])}
      theta[n] ~ dnorm(0,1)}
    for(k in 1:K){
      beta[k] ~ dnorm(0, 0.5)
      xi[k] ~ dnorm(0, 0.5) %_% T(0,)}
    
    for (i in 1:I) {
      logit(g[i]) <- gamma0+gamma[1]*Z[i,1]+gamma[2]*Z[i,2]+gamma[3]*Z[i,3]+gamma[4]*Z[i,4]+
        gamma[5]*Z[i,5]+gamma[6]*Z[i,6]+gamma[7]*Z[i,7]+gamma[8]*Z[i,8]
      #s[i] ~ dbeta(1, 1)
      s[i] ~ dbeta(1, 1) %_% T(, 1 - g[i])
    }
    #Priors
    gamma0~dnorm(0,1.0E-6)
    for (m in 1:M){
      gamma[m]~dnorm(0,1.0E-6) #independent coefficients
    } 
  }
}


if(condition=="HO_DINA_g_resid"){
  bayes.mod <- function(){
    #Standardize x's
    for (m in 1:M){
      for (i in 1:I){
        Z[i,m]<-(X[i,m]-mean(X[,m]))/sd(X[,m])}}
    
    for (n in 1:N) {
      for (i in 1:I) {
        for (k in 1:K) {w[n, i, k] <- pow(alpha[n, k], Q[i, k])}
        eta[n, i] <- prod(w[n, i, 1:K])
        p[n, i] <- pow((1 - s[i]), eta[n, i]) * pow(g[i], (1 - eta[n, i]))
        Y[n, i] ~ dbern(p[n, i])}}
    for(n in 1:N){
      for(k in 1:K){
        logit(prob.a[n, k]) <- xi[k] * theta[n] - beta[k]
        alpha[n, k] ~ dbern(prob.a[n, k])}
      theta[n] ~ dnorm(0,1)}
    for(k in 1:K){
      beta[k] ~ dnorm(0, 0.5)
      xi[k] ~ dnorm(0, 0.5) %_% T(0,)}
    
    for (i in 1:I) {
      logit(g[i]) <- gamma0+gamma[1]*Z[i,1]+gamma[2]*Z[i,2]+gamma[3]*Z[i,3]+gamma[4]*Z[i,4]+
        gamma[5]*Z[i,5]+gamma[6]*Z[i,6]+gamma[7]*Z[i,7]+gamma[8]*Z[i,8]+resid[i]
      resid[i]~dnorm(0,tau)
      #s[i] ~ dbeta(1, 1)
      s[i] ~ dbeta(1, 1) %_% T(, 1 - g[i])
    }
    #Priors
    tau~dgamma(1,1)
    resid_sigma <- 1/tau
    gamma0~dnorm(0,1.0E-6)
    for (m in 1:M){
      gamma[m]~dnorm(0,1.0E-6) #independent coefficients
    } 
  }
}


if(condition=="HO_DINA_s"){
  bayes.mod <- function(){
    #Standardize x's
    for (m in 1:M){
      for (i in 1:I){
        Z[i,m]<-(X[i,m]-mean(X[,m]))/sd(X[,m])}}
    
    for (n in 1:N) {
      for (i in 1:I) {
        for (k in 1:K) {w[n, i, k] <- pow(alpha[n, k], Q[i, k])}
        eta[n, i] <- prod(w[n, i, 1:K])
        p[n, i] <- pow((1 - s[i]), eta[n, i]) * pow(g[i], (1 - eta[n, i]))
        Y[n, i] ~ dbern(p[n, i])}}
    for(n in 1:N){
      for(k in 1:K){
        logit(prob.a[n, k]) <- xi[k] * theta[n] - beta[k]
        alpha[n, k] ~ dbern(prob.a[n, k])}
      theta[n] ~ dnorm(0,1)}
    for(k in 1:K){
      beta[k] ~ dnorm(0, 0.5)
      xi[k] ~ dnorm(0, 0.5) %_% T(0,)}
    
    for (i in 1:I) {
      logit(s[i]) <- gamma0+gamma[1]*Z[i,1]+gamma[2]*Z[i,2]+gamma[3]*Z[i,3]+gamma[4]*Z[i,4]+
        gamma[5]*Z[i,5]+gamma[6]*Z[i,6]+gamma[7]*Z[i,7]+gamma[8]*Z[i,8]
      g[i] ~ dbeta(1, 1) %_% T(, 1 - s[i])
    }
    #Priors
    gamma0~dnorm(0,1.0E-6)
    for (m in 1:M){
      gamma[m]~dnorm(0,1.0E-6) #independent coefficients
    } 
  }
}


if(condition=="HO_DINA_s_resid"){
  bayes.mod <- function(){
    #Standardize x's
    for (m in 1:M){
      for (i in 1:I){
        Z[i,m]<-(X[i,m]-mean(X[,m]))/sd(X[,m])}}
    
    for (n in 1:N) {
      for (i in 1:I) {
        for (k in 1:K) {w[n, i, k] <- pow(alpha[n, k], Q[i, k])}
        eta[n, i] <- prod(w[n, i, 1:K])
        p[n, i] <- pow((1 - s[i]), eta[n, i]) * pow(g[i], (1 - eta[n, i]))
        Y[n, i] ~ dbern(p[n, i])}}
    for(n in 1:N){
      for(k in 1:K){
        logit(prob.a[n, k]) <- xi[k] * theta[n] - beta[k]
        alpha[n, k] ~ dbern(prob.a[n, k])}
      theta[n] ~ dnorm(0,1)}
    for(k in 1:K){
      beta[k] ~ dnorm(0, 0.5)
      xi[k] ~ dnorm(0, 0.5) %_% T(0,)}
    
    for (i in 1:I) {
      logit(s[i]) <- gamma0+gamma[1]*Z[i,1]+gamma[2]*Z[i,2]+gamma[3]*Z[i,3]+gamma[4]*Z[i,4]+
        gamma[5]*Z[i,5]+gamma[6]*Z[i,6]+gamma[7]*Z[i,7]+gamma[8]*Z[i,8]+resid[i]
      resid[i]~dnorm(0,tau)
      g[i] ~ dbeta(1, 1) %_% T(, 1 - s[i])
    }
    #Priors
    tau~dgamma(1,1)
    resid_sigma <- 1/tau
    gamma0~dnorm(0,1.0E-6)
    for (m in 1:M){
      gamma[m]~dnorm(0,1.0E-6) #independent coefficients
    } 
  }
}



# jags model prms ---------------------------------------------------------

if(condition=="HO_DINA_g"|condition=="HO_DINA_s"){
  bayes.mod.params <- c("s","g","alpha","gamma0","gamma")
}

if(condition=="HO_DINA_g_resid"|condition=="HO_DINA_s_resid"){
  bayes.mod.params <- c("s","g","alpha","gamma0","gamma","resid_sigma")
}


dat.jags <- list(N=N, I=I, K=K,M=M,X=item_feature_use,
                 Y=response_data,Q=qmatrix)




# run model ---------------------------------------------------------------

if(use_parallel==F){
  bayes.mod.fit <- jags(data = dat.jags, 
                        parameters.to.save = bayes.mod.params, 
                        n.chains = 2, 
                        n.iter = n.iter,
                        n.burnin = n.burnin,n.thin = n.thin,
                        model.file = bayes.mod,
                        DIC=T,
                        refresh=n.iter/10)
}else if(use_parallel==T){
  bayes.mod.fit <- jags.parallel(data = dat.jags,
                                 parameters.to.save = bayes.mod.params,
                                 n.chains = 2, n.iter = 20000,
                                 n.burnin = 10000,n.thin = 1,
                                 model.file = bayes.mod,DIC=T)
}
