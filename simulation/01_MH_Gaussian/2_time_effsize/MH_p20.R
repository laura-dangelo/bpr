rm(list=ls())
setwd("~/Documents/Dottorato/3.12 Poisson Regression Sim")
# install.packages("/home/laura/Documents/Dottorato/3.10 Pacchetto Poisson/bpr_0.1.tar.gz", repos = NULL, type = "source")
library("bpr")
library("coda")


### n = 25 
MH_p20_n25 <- list()
tmp = list()
load("01_MH_Gaussian/1_tuning/dist_best_MH_p20_n25.Rdata")
seeds = 1:50

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
  run <- sample_bpr( y ~ X-1, data = data.frame(y,X),
                     iter = nrep, perc_burnin = 0.5,
                     pars = list(max_dist = best_dist[seed]),
                     seed = 1,
                     verbose = FALSE)
  end <- Sys.time()
  time <- end - start
  effSize <- apply(run$sim$beta[-burnin,], 2, effectiveSize)
  tmp$run <- run$sim
  tmp$time <- time
  tmp$effSize <- effSize
  MH_p20_n25[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("01_MH_Gaussian/2_time_effsize/MH_toess_p", p, "_n", n, ".Rdata")
save(MH_p20_n25, file = filenamesave)


### n = 50
MH_p20_n50 <- list()
tmp = list()
load("01_MH_Gaussian/1_tuning/dist_best_MH_p20_n50.Rdata")
seeds = 1:50

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
  run <- sample_bpr( y ~ X-1, data = data.frame(y,X),
                     iter = nrep, perc_burnin = 0.5,
                     pars = list(max_dist = best_dist[seed]),
                     seed = 1,
                     verbose = FALSE)
  end <- Sys.time()
  time <- end - start
  effSize <- apply(run$sim$beta[-burnin,], 2, effectiveSize)
  tmp$run <- run$sim
  tmp$time <- time
  tmp$effSize <- effSize
  MH_p20_n50[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("01_MH_Gaussian/2_time_effsize/MH_toess_p", p, "_n", n, ".Rdata")
save(MH_p20_n50, file = filenamesave)


### n = 100
MH_p20_n100 <- list()
tmp = list()
load("01_MH_Gaussian/1_tuning/dist_best_MH_p20_n100.Rdata")
seeds = 1:50

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
  run <- sample_bpr( y ~ X-1, data = data.frame(y,X),
                     iter = nrep, perc_burnin = 0.5,
                     pars = list(max_dist = best_dist[seed]),
                     seed = 1,
                     verbose = FALSE)
  end <- Sys.time()
  time <- end - start
  effSize <- apply(run$sim$beta[-burnin,], 2, effectiveSize)
  tmp$run <- run$sim
  tmp$time <- time
  tmp$effSize <- effSize
  MH_p20_n100[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("01_MH_Gaussian/2_time_effsize/MH_toess_p", p, "_n", n, ".Rdata")
save(MH_p20_n100, file = filenamesave)


### n = 200
rm(list=ls())
MH_p20_n200 <- list()
tmp = list()
load("01_MH_Gaussian/1_tuning/dist_best_MH_p20_n200.Rdata")
seeds = 1:50

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
  run <- sample_bpr( y ~ X-1, data = data.frame(y,X),
                     iter = nrep, perc_burnin = 0.5,
                     pars = list(max_dist = best_dist[seed]),
                     seed = 1,
                     verbose = FALSE)
  end <- Sys.time()
  time <- end - start
  effSize <- apply(run$sim$beta[-burnin,], 2, effectiveSize)
  tmp$run <- run$sim
  tmp$time <- time
  tmp$effSize <- effSize
  MH_p20_n200[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("01_MH_Gaussian/2_time_effsize/MH_toess_p", p, "_n", n, ".Rdata")
save(MH_p20_n200, file = filenamesave)
