rm(list=ls())
setwd("~/Documents/Dottorato/3.12 Poisson Regression Sim")
library('rstan')
model <- stan_model("poisreg.stan")
library("coda")


### n = 25 
stan_p50_n25 <- list()
tmp = list()
seeds = 1:50

for(seed in seeds)
{
  filenameload = paste0("Data_gen/Data/X_p50_n25_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p50_n25_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p50_n25
  y = y_p50_n25
  rm(X_p50_n25)
  rm(y_p50_n25)
  p = ncol(X)
  n = nrow(X)
  
  nsim = 10
  nrep = 10000
  burnin = 1:5000
  
  set.seed(seed)
  start <- Sys.time()
  run <-  sampling(model,
                   data = list(N = n, p = p, X = X, y = y),
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
                   run@sim$samples[[1]]$`beta[19]`, run@sim$samples[[1]]$`beta[20]`,
                   run@sim$samples[[1]]$`beta[21]`, run@sim$samples[[1]]$`beta[22]`,
                   run@sim$samples[[1]]$`beta[23]`, run@sim$samples[[1]]$`beta[24]`,
                   run@sim$samples[[1]]$`beta[25]`, run@sim$samples[[1]]$`beta[26]`,
                   run@sim$samples[[1]]$`beta[27]`, run@sim$samples[[1]]$`beta[28]`,
                   run@sim$samples[[1]]$`beta[29]`, run@sim$samples[[1]]$`beta[30]`,
                   run@sim$samples[[1]]$`beta[31]`, run@sim$samples[[1]]$`beta[32]`,
                   run@sim$samples[[1]]$`beta[33]`, run@sim$samples[[1]]$`beta[34]`,
                   run@sim$samples[[1]]$`beta[35]`, run@sim$samples[[1]]$`beta[36]`,
                   run@sim$samples[[1]]$`beta[37]`, run@sim$samples[[1]]$`beta[38]`,
                   run@sim$samples[[1]]$`beta[39]`, run@sim$samples[[1]]$`beta[40]`,
                   run@sim$samples[[1]]$`beta[41]`, run@sim$samples[[1]]$`beta[42]`,
                   run@sim$samples[[1]]$`beta[43]`, run@sim$samples[[1]]$`beta[44]`,
                   run@sim$samples[[1]]$`beta[45]`, run@sim$samples[[1]]$`beta[46]`,
                   run@sim$samples[[1]]$`beta[47]`, run@sim$samples[[1]]$`beta[48]`,
                   run@sim$samples[[1]]$`beta[49]`, run@sim$samples[[1]]$`beta[50]`
                   ) ### modificare al variare di p
  effSize <- apply(tmp$run[-burnin,], 2, effectiveSize)
  tmp$time <- time
  tmp$effSize <- effSize
  stan_p50_n25[[seed]] <- tmp
  
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("01_MH_Gaussian/2_time_effsize/stan_toess_p", p, "_n", n, ".Rdata")
save(stan_p50_n25, file = filenamesave)


### n = 50
stan_p50_n50 <- list()
tmp = list()
seeds = 1:50

for(seed in seeds)
{
  filenameload = paste0("Data_gen/Data/X_p50_n50_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p50_n50_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p50_n50
  y = y_p50_n50
  rm(X_p50_n50)
  rm(y_p50_n50)
  p = ncol(X)
  n = nrow(X)
  
  nsim = 10
  nrep = 10000
  burnin = 1:5000
  
  set.seed(seed)
  start <- Sys.time()
  run <-  sampling(model,
                   data = list(N = n, p = p, X = X, y = y),
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
                   run@sim$samples[[1]]$`beta[19]`, run@sim$samples[[1]]$`beta[20]`,
                   run@sim$samples[[1]]$`beta[21]`, run@sim$samples[[1]]$`beta[22]`,
                   run@sim$samples[[1]]$`beta[23]`, run@sim$samples[[1]]$`beta[24]`,
                   run@sim$samples[[1]]$`beta[25]`, run@sim$samples[[1]]$`beta[26]`,
                   run@sim$samples[[1]]$`beta[27]`, run@sim$samples[[1]]$`beta[28]`,
                   run@sim$samples[[1]]$`beta[29]`, run@sim$samples[[1]]$`beta[30]`,
                   run@sim$samples[[1]]$`beta[31]`, run@sim$samples[[1]]$`beta[32]`,
                   run@sim$samples[[1]]$`beta[33]`, run@sim$samples[[1]]$`beta[34]`,
                   run@sim$samples[[1]]$`beta[35]`, run@sim$samples[[1]]$`beta[36]`,
                   run@sim$samples[[1]]$`beta[37]`, run@sim$samples[[1]]$`beta[38]`,
                   run@sim$samples[[1]]$`beta[39]`, run@sim$samples[[1]]$`beta[40]`,
                   run@sim$samples[[1]]$`beta[41]`, run@sim$samples[[1]]$`beta[42]`,
                   run@sim$samples[[1]]$`beta[43]`, run@sim$samples[[1]]$`beta[44]`,
                   run@sim$samples[[1]]$`beta[45]`, run@sim$samples[[1]]$`beta[46]`,
                   run@sim$samples[[1]]$`beta[47]`, run@sim$samples[[1]]$`beta[48]`,
                   run@sim$samples[[1]]$`beta[49]`, run@sim$samples[[1]]$`beta[50]`
  ) ### modificare al variare di p
  effSize <- apply(tmp$run[-burnin,], 2, effectiveSize)
  tmp$time <- time
  tmp$effSize <- effSize
  stan_p50_n50[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("01_MH_Gaussian/2_time_effsize/stan_toess_p", p, "_n", n, ".Rdata")
save(stan_p50_n50, file = filenamesave)


### n = 100
stan_p50_n100 <- list()
tmp = list()
seeds = 1:50

for(seed in seeds)
{
  filenameload = paste0("Data_gen/Data/X_p50_n100_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p50_n100_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p50_n100
  y = y_p50_n100
  rm(X_p50_n100)
  rm(y_p50_n100)
  p = ncol(X)
  n = nrow(X)
  
  nsim = 10
  nrep = 10000
  burnin = 1:5000
  
  set.seed(seed)
  start <- Sys.time()
  run <-  sampling(model,
                   data = list(N = n, p = p, X = X, y = y),
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
                   run@sim$samples[[1]]$`beta[19]`, run@sim$samples[[1]]$`beta[20]`,
                   run@sim$samples[[1]]$`beta[21]`, run@sim$samples[[1]]$`beta[22]`,
                   run@sim$samples[[1]]$`beta[23]`, run@sim$samples[[1]]$`beta[24]`,
                   run@sim$samples[[1]]$`beta[25]`, run@sim$samples[[1]]$`beta[26]`,
                   run@sim$samples[[1]]$`beta[27]`, run@sim$samples[[1]]$`beta[28]`,
                   run@sim$samples[[1]]$`beta[29]`, run@sim$samples[[1]]$`beta[30]`,
                   run@sim$samples[[1]]$`beta[31]`, run@sim$samples[[1]]$`beta[32]`,
                   run@sim$samples[[1]]$`beta[33]`, run@sim$samples[[1]]$`beta[34]`,
                   run@sim$samples[[1]]$`beta[35]`, run@sim$samples[[1]]$`beta[36]`,
                   run@sim$samples[[1]]$`beta[37]`, run@sim$samples[[1]]$`beta[38]`,
                   run@sim$samples[[1]]$`beta[39]`, run@sim$samples[[1]]$`beta[40]`,
                   run@sim$samples[[1]]$`beta[41]`, run@sim$samples[[1]]$`beta[42]`,
                   run@sim$samples[[1]]$`beta[43]`, run@sim$samples[[1]]$`beta[44]`,
                   run@sim$samples[[1]]$`beta[45]`, run@sim$samples[[1]]$`beta[46]`,
                   run@sim$samples[[1]]$`beta[47]`, run@sim$samples[[1]]$`beta[48]`,
                   run@sim$samples[[1]]$`beta[49]`, run@sim$samples[[1]]$`beta[50]`
  ) ### modificare al variare di p
  effSize <- apply(tmp$run[-burnin,], 2, effectiveSize)
  tmp$time <- time
  tmp$effSize <- effSize
  stan_p50_n100[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("01_MH_Gaussian/2_time_effsize/stan_toess_p", p, "_n", n, ".Rdata")
save(stan_p50_n100, file = filenamesave)


### n = 200
stan_p50_n200 <- list()
tmp = list()
seeds = 1:50

for(seed in seeds)
{
  filenameload = paste0("Data_gen/Data/X_p50_n200_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p50_n200_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p50_n200
  y = y_p50_n200
  rm(X_p50_n200)
  rm(y_p50_n200)
  p = ncol(X)
  n = nrow(X)
  
  nsim = 10
  nrep = 10000
  burnin = 1:5000
  
  set.seed(seed)
  start <- Sys.time()
  run <-  sampling(model,
                   data = list(N = n, p = p, X = X, y = y),
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
                   run@sim$samples[[1]]$`beta[19]`, run@sim$samples[[1]]$`beta[20]`,
                   run@sim$samples[[1]]$`beta[21]`, run@sim$samples[[1]]$`beta[22]`,
                   run@sim$samples[[1]]$`beta[23]`, run@sim$samples[[1]]$`beta[24]`,
                   run@sim$samples[[1]]$`beta[25]`, run@sim$samples[[1]]$`beta[26]`,
                   run@sim$samples[[1]]$`beta[27]`, run@sim$samples[[1]]$`beta[28]`,
                   run@sim$samples[[1]]$`beta[29]`, run@sim$samples[[1]]$`beta[30]`,
                   run@sim$samples[[1]]$`beta[31]`, run@sim$samples[[1]]$`beta[32]`,
                   run@sim$samples[[1]]$`beta[33]`, run@sim$samples[[1]]$`beta[34]`,
                   run@sim$samples[[1]]$`beta[35]`, run@sim$samples[[1]]$`beta[36]`,
                   run@sim$samples[[1]]$`beta[37]`, run@sim$samples[[1]]$`beta[38]`,
                   run@sim$samples[[1]]$`beta[39]`, run@sim$samples[[1]]$`beta[40]`,
                   run@sim$samples[[1]]$`beta[41]`, run@sim$samples[[1]]$`beta[42]`,
                   run@sim$samples[[1]]$`beta[43]`, run@sim$samples[[1]]$`beta[44]`,
                   run@sim$samples[[1]]$`beta[45]`, run@sim$samples[[1]]$`beta[46]`,
                   run@sim$samples[[1]]$`beta[47]`, run@sim$samples[[1]]$`beta[48]`,
                   run@sim$samples[[1]]$`beta[49]`, run@sim$samples[[1]]$`beta[50]`
  ) ### modificare al variare di p
  effSize <- apply(tmp$run[-burnin,], 2, effectiveSize)
  tmp$time <- time
  tmp$effSize <- effSize
  stan_p50_n200[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("01_MH_Gaussian/2_time_effsize/stan_toess_p", p, "_n", n, ".Rdata")
save(stan_p50_n200, file = filenamesave)
