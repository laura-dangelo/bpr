rm(list=ls())
setwd("~/Documents/Dottorato/3.12 Poisson Regression Sim")
# install.packages("/home/laura/Documents/Dottorato/3.10 Pacchetto Poisson/bpr_0.1.tar.gz", repos = NULL, type = "source")
library("bpr")
library("coda")


#-------# #-------# #-------# #-------# #-------# #-------# #-------# 
#-------#          Sample with MH polya adaptive r          #-------# 
#-------# #-------# #-------# #-------# #-------# #-------# #-------# 

p = 5
set.seed(9)
sim_mix <- function(n, p, m, v)
{
  k <- sample(1:length(m), n, replace = TRUE, prob = p)
  return(rnorm(n, m[k], v[k]))
}
beta_p5 <- sim_mix(p, c(0.3,0.5,0.4,0.4), c(0.6,-0.1,2.5,-2.4), c(0.7,0.7,0.7,0.7))


est_tau = function(n, beta)
{
  p_n = sum(beta != 0)
  return( p_n/n * sqrt(log(n/p_n)) )
}




#-------# #-------# n = 25

nrep = 1000
burnin <- c(1:500)
dist_seq = sort(rep(seq(1, 300, length.out = 5),3))
best_dist = numeric(50)

tau = est_tau(25, beta_p5)

for(seed in 1:50)
{
  filenameload = paste0("Data_gen/Data/X_p5_n25_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p5_n25_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p5_n25
  y = y_p5_n25
  rm(X_p5_n25)
  rm(y_p5_n25)
  p = ncol(X)
  n = nrow(X)
  
  effic_sim_mh_p5_n25 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep, perc_burnin = 0.5,
                       prior = list(type = "horseshoe", tau = tau),
                       pars = list(max_dist = dist_seq[d]),
                       verbose = FALSE)
    
    effic_sim_mh_p5_n25[[d]] <- run$sim$beta[-burnin,]
    
  }

  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_mh_p5_n25, function(x) mean(apply(x, 2, effectiveSize ) )) ) )
  
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )

  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="02_MH_horseshoe/1_tuning/dist_best_MH_p5_n25.Rdata")





#-------# #-------# n = 50

nrep = 1000
burnin <- c(1:500)
dist_seq = sort(rep(seq(1, 300, length.out = 5),3))
best_dist = numeric(50)

tau = est_tau(50, beta_p5)

for(seed in 1:50)
{
  filenameload = paste0("Data_gen/Data/X_p5_n50_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p5_n50_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p5_n50
  y = y_p5_n50
  rm(X_p5_n50)
  rm(y_p5_n50)
  p = ncol(X)
  n = nrow(X)
  
  effic_sim_mh_p5_n50 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep, perc_burnin = 0.5,
                       prior = list(type = "horseshoe", tau = tau),
                       pars = list(max_dist = dist_seq[d]),
                       verbose = FALSE)
    
    effic_sim_mh_p5_n50[[d]] <- run$sim$beta[-burnin,]
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_mh_p5_n50, function(x) mean(apply(x, 2, effectiveSize ) )) ) )
  
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="02_MH_horseshoe/1_tuning/dist_best_MH_p5_n50.Rdata")




#-------# #-------# n = 100
nrep = 1000
burnin <- c(1:500)
dist_seq = sort(rep(seq(1, 300, length.out = 5),3))
best_dist = numeric(50)

tau = est_tau(100, beta_p5)

for(seed in 1:50)
{
  filenameload = paste0("Data_gen/Data/X_p5_n100_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p5_n100_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p5_n100
  y = y_p5_n100
  rm(X_p5_n100)
  rm(y_p5_n100)
  p = ncol(X)
  n = nrow(X)
  
  effic_sim_mh_p5_n100 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep, perc_burnin = 0.5,
                       prior = list(type = "horseshoe", tau = tau),
                       pars = list(max_dist = dist_seq[d]),
                       verbose = FALSE)
    
    effic_sim_mh_p5_n100[[d]] <- run$sim$beta[-burnin,]
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_mh_p5_n100, function(x) mean(apply(x, 2, effectiveSize ) )) ) )
  
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="02_MH_horseshoe/1_tuning/dist_best_MH_p5_n100.Rdata")





#-------# #-------# n = 200
nrep = 1000
burnin <- c(1:500)
dist_seq = sort(rep(seq(1, 300, length.out = 5),3))
best_dist = numeric(50)

tau = est_tau(200, beta_p5)

for(seed in 1:50)
{
  filenameload = paste0("Data_gen/Data/X_p5_n200_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p5_n200_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p5_n200
  y = y_p5_n200
  rm(X_p5_n200)
  rm(y_p5_n200)
  p = ncol(X)
  n = nrow(X)
  
  effic_sim_mh_p5_n200 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep, perc_burnin = 0.5,
                       prior = list(type = "horseshoe", tau = tau),
                       pars = list(max_dist = dist_seq[d]),
                       verbose = FALSE)
    
    effic_sim_mh_p5_n200[[d]] <- run$sim$beta[-burnin,]
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_mh_p5_n200, function(x) mean(apply(x, 2, effectiveSize ) )) ) )
  
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="02_MH_horseshoe/1_tuning/dist_best_MH_p5_n200.Rdata")

