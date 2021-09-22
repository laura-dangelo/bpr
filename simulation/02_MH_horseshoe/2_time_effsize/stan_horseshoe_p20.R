rm(list=ls())
setwd("~/Documents/Dottorato/3.12 Poisson Regression Sim")
library('rstan')
library("coda")
model <- stan_model("poisreg_horseshoe.stan")

p = 20
set.seed(9)
sim_mix <- function(n, p, m, v)
{
  k <- sample(1:length(m), n, replace = TRUE, prob = p)
  return(rnorm(n, m[k], v[k]))
}
beta_p20 <- sim_mix(p, c(0.3,0.5,0.4,0.4), c(0.6,-0.1,2.5,-2.4), c(0.7,0.7,0.7,0.7))


est_tau = function(n, beta)
{
  p_n = sum(abs(beta) > 0.05)
  return( p_n/n * sqrt(log(n/p_n)) )
}


### n = 25 
stan_horseshoe_p20_n25 <- list()
tmp = list()
seeds = 1:50

tau = est_tau(25, beta_p20)

for(seed in seeds)
{
  filenameload = paste0("Data_gen/Data/X_p20_n25_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p20_n25_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p20_n25
  y = y_p20_n25
  rm(X_p20_n25)
  rm(y_p20_n25)
  p = ncol(X)
  n = nrow(X)
  
  nsim = 10
  nrep = 10000
  burnin = 1:5000
  
  set.seed(seed)
  start <- Sys.time()
  run <-  sampling(model,
                   data = list(N = n, p = p, X = X, y = y, tau = tau),
                   chains = 1, iter = nrep, warmup = 5000, verbose = FALSE)
  end <- Sys.time()
  time <- end - start
  tmp$run <- cbind(run@sim$samples[[1]]$`beta[1]`, run@sim$samples[[1]]$`beta[2]`,
                   run@sim$samples[[1]]$`beta[3]`, run@sim$samples[[1]]$`beta[4]`,
                   run@sim$samples[[1]]$`beta[5]`, run@sim$samples[[1]]$`beta[6]`,
                   run@sim$samples[[1]]$`beta[7]`, run@sim$samples[[1]]$`beta[8]`,
                   run@sim$samples[[1]]$`beta[9]`, run@sim$samples[[1]]$`beta[10]`,
                   run@sim$samples[[1]]$`beta[11]`, run@sim$samples[[1]]$`beta[12]`,
                   run@sim$samples[[1]]$`beta[13]`, run@sim$samples[[1]]$`beta[14]`,
                   run@sim$samples[[1]]$`beta[15]`, run@sim$samples[[1]]$`beta[16]`,
                   run@sim$samples[[1]]$`beta[17]`, run@sim$samples[[1]]$`beta[18]`,
                   run@sim$samples[[1]]$`beta[19]`, run@sim$samples[[1]]$`beta[20]`
                   ) ### modificare al variare di p
  effSize <- apply(tmp$run[-burnin,], 2, effectiveSize)
  tmp$time <- time
  tmp$effSize <- effSize
  stan_horseshoe_p20_n25[[seed]] <- tmp
  
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("02_MH_horseshoe/2_time_effsize/stan_horseshoe_toess_p", p, "_n", n, ".Rdata")
save(stan_horseshoe_p20_n25, file = filenamesave)


### n = 50
stan_horseshoe_p20_n50 <- list()
tmp = list()
seeds = 1:50

tau = est_tau(50, beta_p20)
tau

for(seed in seeds)
{
  filenameload = paste0("Data_gen/Data/X_p20_n50_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p20_n50_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p20_n50
  y = y_p20_n50
  rm(X_p20_n50)
  rm(y_p20_n50)
  p = ncol(X)
  n = nrow(X)
  
  nsim = 10
  nrep = 10000
  burnin = 1:5000
  
  set.seed(seed)
  start <- Sys.time()
  run <-  sampling(model,
                   data = list(N = n, p = p, X = X, y = y, tau = tau),
                   chains = 1, iter = nrep, warmup = 5000, verbose = FALSE)
  end <- Sys.time()
  time <- end - start
  tmp$run <- cbind(run@sim$samples[[1]]$`beta[1]`, run@sim$samples[[1]]$`beta[2]`,
                   run@sim$samples[[1]]$`beta[3]`, run@sim$samples[[1]]$`beta[4]`,
                   run@sim$samples[[1]]$`beta[5]`, run@sim$samples[[1]]$`beta[6]`,
                   run@sim$samples[[1]]$`beta[7]`, run@sim$samples[[1]]$`beta[8]`,
                   run@sim$samples[[1]]$`beta[9]`, run@sim$samples[[1]]$`beta[10]`,
                   run@sim$samples[[1]]$`beta[11]`, run@sim$samples[[1]]$`beta[12]`,
                   run@sim$samples[[1]]$`beta[13]`, run@sim$samples[[1]]$`beta[14]`,
                   run@sim$samples[[1]]$`beta[15]`, run@sim$samples[[1]]$`beta[16]`,
                   run@sim$samples[[1]]$`beta[17]`, run@sim$samples[[1]]$`beta[18]`,
                   run@sim$samples[[1]]$`beta[19]`, run@sim$samples[[1]]$`beta[20]`
  ) ### modificare al variare di p
  effSize <- apply(tmp$run[-burnin,], 2, effectiveSize)
  tmp$time <- time
  tmp$effSize <- effSize
  stan_horseshoe_p20_n50[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("02_MH_horseshoe/2_time_effsize/stan_horseshoe_toess_p", p, "_n", n, ".Rdata")
save(stan_horseshoe_p20_n50, file = filenamesave)


### n = 100
stan_horseshoe_p20_n100 <- list()
tmp = list()
seeds = 1:50

tau = est_tau(100, beta_p20)
tau

for(seed in seeds)
{
  filenameload = paste0("Data_gen/Data/X_p20_n100_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p20_n100_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p20_n100
  y = y_p20_n100
  rm(X_p20_n100)
  rm(y_p20_n100)
  p = ncol(X)
  n = nrow(X)
  
  nsim = 10
  nrep = 10000
  burnin = 1:5000
  
  set.seed(seed)
  start <- Sys.time()
  run <-  sampling(model,
                   data = list(N = n, p = p, X = X, y = y, tau = tau),
                   chains = 1, iter = nrep, warmup = 5000, verbose = FALSE)
  end <- Sys.time()
  time <- end - start
  tmp$run <- cbind(run@sim$samples[[1]]$`beta[1]`, run@sim$samples[[1]]$`beta[2]`,
                   run@sim$samples[[1]]$`beta[3]`, run@sim$samples[[1]]$`beta[4]`,
                   run@sim$samples[[1]]$`beta[5]`, run@sim$samples[[1]]$`beta[6]`,
                   run@sim$samples[[1]]$`beta[7]`, run@sim$samples[[1]]$`beta[8]`,
                   run@sim$samples[[1]]$`beta[9]`, run@sim$samples[[1]]$`beta[10]`,
                   run@sim$samples[[1]]$`beta[11]`, run@sim$samples[[1]]$`beta[12]`,
                   run@sim$samples[[1]]$`beta[13]`, run@sim$samples[[1]]$`beta[14]`,
                   run@sim$samples[[1]]$`beta[15]`, run@sim$samples[[1]]$`beta[16]`,
                   run@sim$samples[[1]]$`beta[17]`, run@sim$samples[[1]]$`beta[18]`,
                   run@sim$samples[[1]]$`beta[19]`, run@sim$samples[[1]]$`beta[20]`
  ) ### modificare al variare di p
  effSize <- apply(tmp$run[-burnin,], 2, effectiveSize)
  tmp$time <- time
  tmp$effSize <- effSize
  stan_horseshoe_p20_n100[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("02_MH_horseshoe/2_time_effsize/stan_horseshoe_toess_p", p, "_n", n, ".Rdata")
save(stan_horseshoe_p20_n100, file = filenamesave)


### n = 200
stan_horseshoe_p20_n200 <- list()
tmp = list()
seeds = 1:50

tau = est_tau(200, beta_p20)
tau


for(seed in seeds)
{
  filenameload = paste0("Data_gen/Data/X_p20_n200_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p20_n200_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p20_n200
  y = y_p20_n200
  rm(X_p20_n200)
  rm(y_p20_n200)
  p = ncol(X)
  n = nrow(X)
  
  nsim = 10
  nrep = 10000
  burnin = 1:5000
  
  set.seed(seed)
  start <- Sys.time()
  run <-  sampling(model,
                   data = list(N = n, p = p, X = X, y = y, tau = tau),
                   chains = 1, iter = nrep, warmup = 5000, verbose = FALSE)
  end <- Sys.time()
  time <- end - start
  tmp$run <- cbind(run@sim$samples[[1]]$`beta[1]`, run@sim$samples[[1]]$`beta[2]`,
                   run@sim$samples[[1]]$`beta[3]`, run@sim$samples[[1]]$`beta[4]`,
                   run@sim$samples[[1]]$`beta[5]`, run@sim$samples[[1]]$`beta[6]`,
                   run@sim$samples[[1]]$`beta[7]`, run@sim$samples[[1]]$`beta[8]`,
                   run@sim$samples[[1]]$`beta[9]`, run@sim$samples[[1]]$`beta[10]`,
                   run@sim$samples[[1]]$`beta[11]`, run@sim$samples[[1]]$`beta[12]`,
                   run@sim$samples[[1]]$`beta[13]`, run@sim$samples[[1]]$`beta[14]`,
                   run@sim$samples[[1]]$`beta[15]`, run@sim$samples[[1]]$`beta[16]`,
                   run@sim$samples[[1]]$`beta[17]`, run@sim$samples[[1]]$`beta[18]`,
                   run@sim$samples[[1]]$`beta[19]`, run@sim$samples[[1]]$`beta[20]`
  ) ### modificare al variare di p
  effSize <- apply(tmp$run[-burnin,], 2, effectiveSize)
  tmp$time <- time
  tmp$effSize <- effSize
  stan_horseshoe_p20_n200[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("02_MH_horseshoe/2_time_effsize/stan_horseshoe_toess_p", p, "_n", n, ".Rdata")
save(stan_horseshoe_p20_n200, file = filenamesave)
