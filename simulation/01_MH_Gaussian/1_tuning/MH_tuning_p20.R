rm(list=ls())
setwd("~/Documents/Dottorato/3.12 Poisson Regression Sim")
# install.packages("/home/laura/Documents/Dottorato/3.10 Pacchetto Poisson/bpr_0.1.tar.gz", repos = NULL, type = "source")
library("bpr")
library("coda")


#-------# #-------# #-------# #-------# #-------# #-------# #-------# 
#-------#          Sample with MH polya adaptive r          #-------# 
#-------# #-------# #-------# #-------# #-------# #-------# #-------# 

#-------# #-------# n = 25

nrep = 2000
burnin <- c(1:1000)
dist_seq = sort(rep(seq(10, 300, length.out = 5),4))
best_dist = numeric(50)


for(seed in 1:50)
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
  
  effic_sim_mh_p20_n25 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep,
                       pars = list(max_dist = dist_seq[d]),
                       verbose = FALSE)
    if(run$sim$acceptance_rate < 0.05) print(seed)
    effic_sim_mh_p20_n25[[d]] <- run$sim$beta[-burnin,]
    
  }

  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_mh_p20_n25, function(x) mean(apply(x, 2, effectiveSize ) )) ) )
  
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )

  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="01_MH_Gaussian/1_tuning/dist_best_MH_p20_n25.Rdata")



#-------# #-------# n = 50
rm(list=ls())
nrep = 2000
burnin <- c(1:1000)
dist_seq = sort(rep(seq(10, 300, length.out = 5),4))
best_dist = numeric(50)


for(seed in 1:50)
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
  
  effic_sim_mh_p20_n50 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep,
                       pars = list(max_dist = dist_seq[d]),
                       verbose = FALSE)
    
    effic_sim_mh_p20_n50[[d]] <- run$sim$beta[-burnin,]
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_mh_p20_n50, function(x) mean(apply(x, 2, effectiveSize ) )) ) )
  
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="01_MH_Gaussian/1_tuning/dist_best_MH_p20_n50.Rdata")



#-------# #-------# n = 100
rm(list=ls())
nrep = 2000
burnin <- c(1:1000)
dist_seq = sort(rep(seq(10, 300, length.out = 5),4))
best_dist = numeric(50)


for(seed in 1:50)
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
  
  effic_sim_mh_p20_n100 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep,
                       pars = list(max_dist = dist_seq[d]),
                       verbose = FALSE)
    
    effic_sim_mh_p20_n100[[d]] <- run$sim$beta[-burnin,]
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_mh_p20_n100, function(x) mean(apply(x, 2, effectiveSize ) )) ) )
  
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="01_MH_Gaussian/1_tuning/dist_best_MH_p20_n100.Rdata")



#-------# #-------# n = 200
rm(list=ls())
nrep = 2000
burnin <- c(1:1000)
dist_seq = sort(rep(seq(50, 500, length.out = 5),4))
best_dist = numeric(50)


for(seed in 1:50)
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
  
  effic_sim_mh_p20_n200 <- list()
  for(d in 1:length(dist_seq))
  {
    run <- sample_bpr( y ~ X - 1, data = data.frame("y" = y, X), 
                       iter = nrep,
                       pars = list(max_dist = dist_seq[d]),
                       verbose = FALSE)
    
    effic_sim_mh_p20_n200[[d]] <- run$sim$beta[-burnin,]
    
  }
  
  df = data.frame(dist = dist_seq, 
                  value = unlist(lapply(effic_sim_mh_p20_n200, function(x) mean(apply(x, 2, effectiveSize ) )) ) )
  
  best_dist[seed] = as.numeric( attr(which.max(tapply(df$value, df$dist, mean)), "names") )
  
  if(seed%%5 == 0) { print(seed) }
}

best_dist   = round(best_dist)

save(best_dist, 
     file="01_MH_Gaussian/1_tuning/dist_best_MH_p20_n200.Rdata")

