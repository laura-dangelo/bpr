rm(list=ls())
setwd("~/Documents/Dottorato/3.12 Poisson Regression Sim")
# install.packages("/home/laura/Documents/Dottorato/3.10 Pacchetto Poisson/bpr_0.1.tar.gz", repos = NULL, type = "source")
library("bpr")

effectiveSize = function(x) sum(x)^2 / sum(x^2)

#-------# #-------# #-------# #-------# #-------# #-------# #-------# 
#-------#          Sample with MH polya adaptive r          #-------# 
#-------# #-------# #-------# #-------# #-------# #-------# #-------# 

#-------# #-------# n = 25

nrep = 2000
burnin <- c(1:1000)
dist_seq = sort(rep(seq(50, 500, length.out = 10),4))
best_dist = numeric(50)


for(seed in 1:50)
{
  filenameload = paste0("Data_gen/Data/X_p10_n25_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p10_n25_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p10_n25
  y = y_p10_n25
  rm(X_p10_n25)
  rm(y_p10_n25)
  p = ncol(X)
  n = nrow(X)
  
  effic_sim_is_p10_n25 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep,
                       pars = list( method = "IS",
                                    max_dist = dist_seq[d] ),
                       verbose = FALSE)
    
    w = exp( run$sim$logw[-burnin] - mean(run$sim$logw[-burnin]) )
    w = w / mean(w)
    effic_sim_is_p10_n25[[d]] <- w
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_is_p10_n25, effectiveSize  )) ) 
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="03_IS_Gaussian/1_tuning/dist_best_IS_p10_n25.Rdata")





#-------# #-------# n = 50

best_dist = numeric(50)

for(seed in 1:50)
{
  filenameload = paste0("Data_gen/Data/X_p10_n50_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p10_n50_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p10_n50
  y = y_p10_n50
  rm(X_p10_n50)
  rm(y_p10_n50)
  p = ncol(X)
  n = nrow(X)
  
  effic_sim_is_p10_n50 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep,
                       pars = list( method = "IS",
                                    max_dist = dist_seq[d] ),
                       verbose = FALSE)
    
    w = exp( run$sim$logw[-burnin] - mean(run$sim$logw[-burnin]) )
    w = w / mean(w)
    effic_sim_is_p10_n50[[d]] <- w
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_is_p10_n50, effectiveSize  )) ) 
  if(var(df$value) == 0) print(seed)
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="03_IS_Gaussian/1_tuning/dist_best_IS_p10_n50.Rdata")




#-------# #-------# n = 100

best_dist = numeric(50)


for(seed in 1:50)
{
  filenameload = paste0("Data_gen/Data/X_p10_n100_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p10_n100_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p10_n100
  y = y_p10_n100
  rm(X_p10_n100)
  rm(y_p10_n100)
  p = ncol(X)
  n = nrow(X)
  
  effic_sim_is_p10_n100 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep,
                       pars = list( method = "IS",
                                    max_dist = dist_seq[d] ),
                       verbose = FALSE)
    
    w = exp( run$sim$logw[-burnin] - mean(run$sim$logw[-burnin]) )
    w = w / mean(w)
    effic_sim_is_p10_n100[[d]] <- w
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_is_p10_n100, effectiveSize  )) ) 
  if(var(df$value) == 0) print(seed)
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="03_IS_Gaussian/1_tuning/dist_best_IS_p10_n100.Rdata")



#-------# #-------# n = 200

best_dist = numeric(50)

for(seed in 1:50)
{
  filenameload = paste0("Data_gen/Data/X_p10_n200_seed", seed, ".Rdata")
  load(filenameload)
  filenameload = paste0("Data_gen/Data/y_p10_n200_seed", seed, ".Rdata")
  load(filenameload)
  X = X_p10_n200
  y = y_p10_n200
  rm(X_p10_n200)
  rm(y_p10_n200)
  p = ncol(X)
  n = nrow(X)
  
  effic_sim_is_p10_n200 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep,
                       pars = list( method = "IS",
                                    max_dist = dist_seq[d] ),
                       verbose = FALSE)
    
    w = exp( run$sim$logw[-burnin] - mean(run$sim$logw[-burnin]) )
    w = w / mean(w)
    effic_sim_is_p10_n200[[d]] <- w
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_is_p10_n200, effectiveSize  )) ) 
  if(var(df$value) == 0) print(seed)
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="03_IS_Gaussian/1_tuning/dist_best_IS_p10_n200.Rdata")

