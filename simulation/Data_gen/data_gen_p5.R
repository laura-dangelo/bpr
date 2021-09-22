p = 5
set.seed(9)
sim_mix <- function(n, p, m, v)
{
  k <- sample(1:length(m), n, replace = TRUE, prob = p)
  return(rnorm(n, m[k], v[k]))
}
beta_p5 <- sim_mix(p, c(0.3,0.5,0.4,0.4), c(0.6,-0.1,2.5,-2.4), c(0.7,0.7,0.7,0.7))



seeds = 1:50

for(seed in seeds)
{
  n = 1000
  set.seed(seed)
  X <- matrix(rep(NA, n*p), n, p)
  X[,1] <- 1
  X[,2] <- runif(n, -2, 1)
  X[,3] <- rgamma(n, 2,2)
  X[,4] <- rnorm(n, 2,2)
  X[,5] = rexp(n, 2) + rpois(n, 1)
  
  for(j in c(2:5))
  {
    X[,j] = (X[,j] - mean(X[,j]))/sd(X[,j])
  }
  
  ix = ( exp(X%*%beta_p5) > 1 ) & ( exp(X%*%beta_p5) < 200 )
  if(sum(ix) < 200) stop()
  X = X[ix,]
  
  #-------------#   n=25   #-------------#
  n = 25
  ixn = sample.int(sum(ix),n)
  X_p5_n25 <- X[ixn,]
  
  set.seed(seed)
  y_p5_n25 <- rpois(n, exp(X_p5_n25%*%beta_p5))
  
  filenameX = paste0("X_p5_n", n, "_seed", seed, ".Rdata")
  filenamey = paste0("y_p5_n", n, "_seed", seed, ".Rdata")
  save(X_p5_n25, file=filenameX)
  save(y_p5_n25, file=filenamey)
  
  #-------------#   n=50   #-------------#
  n = 50
  ixn = sample.int(sum(ix),n)
  X_p5_n50 <- X[ixn,]
  
  set.seed(seed)
  y_p5_n50 <- rpois(n, exp(X_p5_n50%*%beta_p5))
  filenameX = paste0("X_p5_n", n, "_seed", seed, ".Rdata")
  filenamey = paste0("y_p5_n", n, "_seed", seed, ".Rdata")
  save(X_p5_n50, file=filenameX)
  save(y_p5_n50, file=filenamey)
  
  
  #-------------#   n=100   #-------------#
  n = 100
  ixn = sample.int(sum(ix),n)
  X_p5_n100 <- X[ixn,]
  
  set.seed(seed)
  y_p5_n100 <- rpois(n, exp(X_p5_n100%*%beta_p5))
  filenameX = paste0("X_p5_n", n, "_seed", seed, ".Rdata")
  filenamey = paste0("y_p5_n", n, "_seed", seed, ".Rdata")
  save(X_p5_n100, file=filenameX)
  save(y_p5_n100, file=filenamey)
  
  
  #-------------#   n=200   #-------------#
  n = 200
  ixn = sample.int(sum(ix),n)
  X_p5_n200 <- X[ixn,]
  
  set.seed(seed)
  y_p5_n200 <- rpois(n, exp(X_p5_n200%*%beta_p5))
  filenameX = paste0("X_p5_n", n, "_seed", seed, ".Rdata")
  filenamey = paste0("y_p5_n", n, "_seed", seed, ".Rdata")
  save(X_p5_n200, file=filenameX)
  save(y_p5_n200, file=filenamey)
  
  
}



