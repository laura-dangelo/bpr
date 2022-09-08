rm(list=ls())
setwd("~/Documents/Dottorato/3.12 Poisson Regression Sim")
# install.packages("/home/laura/Documents/Dottorato/3.10 Pacchetto Poisson/bpr_0.1.tar.gz", repos = NULL, type = "source")
library("bpr")
effectiveSize = function(x) sum(x)^2 / sum(x^2)

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
  p_n = sum(beta != 0)
  return( p_n/n * sqrt(log(n/p_n)) )
}


### n = 25 
IS_horseshoe_p20_n25 <- list()
tmp = list()
seeds = 1:50
load("04_IS_horseshoe/1_tuning/dist_best_IS_horseshoe_p20_n25.Rdata")
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
  run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                     iter = nrep,
                     prior = list(type = "horseshoe", tau = tau),
                     pars = list( method = "IS",
                                  max_dist = best_dist[seed] ),
                     verbose = FALSE,
                     seed = 1)
  end <- Sys.time()
  time <- end - start
  meanw = mean(run$sim$logw[-burnin])
  w = exp( run$sim$logw - meanw )
  w = w / mean(w)
  effSize <- effectiveSize(w)
  tmp$run <- run$sim
  tmp$time <- time
  tmp$effSize <- effSize
  
  IS_horseshoe_p20_n25[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("04_IS_horseshoe/2_time_effsize/IS_horseshoe_toess_p", p, "_n", n, ".Rdata")
save(IS_horseshoe_p20_n25, file = filenamesave)
rm(IS_horseshoe_p20_n25)


### n = 50
IS_horseshoe_p20_n50 <- list()
tmp = list()
seeds = 1:50
load("04_IS_horseshoe/1_tuning/dist_best_IS_horseshoe_p20_n50.Rdata")
tau = est_tau(50, beta_p20)

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
  run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                     iter = nrep,
                     prior = list(type = "horseshoe", tau = tau),
                     pars = list( method = "IS",
                                  max_dist = best_dist[seed] ),
                     verbose = FALSE,
                     seed = 1)
  end <- Sys.time()
  time <- end - start
  meanw = mean(run$sim$logw[-burnin])
  w = exp( run$sim$logw - meanw )
  w = w / mean(w)
  effSize <- effectiveSize(w)
  tmp$run <- run$sim
  tmp$time <- time
  tmp$effSize <- effSize
  
  IS_horseshoe_p20_n50[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("04_IS_horseshoe/2_time_effsize/IS_horseshoe_toess_p", p, "_n", n, ".Rdata")
save(IS_horseshoe_p20_n50, file = filenamesave)
rm(IS_horseshoe_p20_n50)

### n = 100
IS_horseshoe_p20_n100 <- list()
tmp = list()
seeds = 1:50
load("04_IS_horseshoe/1_tuning/dist_best_IS_horseshoe_p20_n100.Rdata")
tau = est_tau(100, beta_p20)

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
  run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                     iter = nrep,
                     prior = list(type = "horseshoe", tau = tau),
                     pars = list( method = "IS",
                                  max_dist = best_dist[seed] ),
                     verbose = FALSE,
                     seed = 1)
  end <- Sys.time()
  time <- end - start
  meanw = mean(run$sim$logw[-burnin])
  w = exp( run$sim$logw - meanw )
  w = w / mean(w)
  effSize <- effectiveSize(w)
  tmp$run <- run$sim
  tmp$time <- time
  tmp$effSize <- effSize
  
  IS_horseshoe_p20_n100[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("04_IS_horseshoe/2_time_effsize/IS_horseshoe_toess_p", p, "_n", n, ".Rdata")
save(IS_horseshoe_p20_n100, file = filenamesave)
rm(IS_horseshoe_p20_n100)

### n = 200

IS_horseshoe_p20_n200 <- list()
tmp = list()
seeds = 1:50
load("04_IS_horseshoe/1_tuning/dist_best_IS_horseshoe_p20_n200.Rdata")
tau = est_tau(200, beta_p20)

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
  run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                     iter = nrep,
                     prior = list(type = "horseshoe", tau = tau),
                     pars = list( method = "IS",
                                  max_dist = best_dist[seed] ),
                     verbose = FALSE,
                     seed = 1)
  end <- Sys.time()
  time <- end - start
  meanw = mean(run$sim$logw[-burnin])
  w = exp( run$sim$logw - meanw )
  w = w / mean(w)
  effSize <- effectiveSize(w)
  tmp$run <- run$sim
  tmp$time <- time
  tmp$effSize <- effSize
  
  IS_horseshoe_p20_n200[[seed]] <- tmp
  
  if(seed%%10==0) print(seed)
}

filenamesave = paste0("04_IS_horseshoe/2_time_effsize/IS_horseshoe_toess_p", p, "_n", n, ".Rdata")
save(IS_horseshoe_p20_n200, file = filenamesave)
