rm(list=ls())
library("coda")
library(Rcpp)
library(RcppDist)
load("data_calcium_imaging_for_poisson.RData")
str(data)


#### quadratic on depth
data$depth2 = data$depth^2
X = model.matrix(~ ., data = data[,c(2:4,6,9)])
str(X)
p = ncol(X)
n = nrow(X)
y = data$n_spikes
str(y)

mle.start <- c(summary(glm(y ~ X - 1, family = "poisson"(link = "log")))$coef[,1])
str(mle.start)

library(bpr)
nrep = 10000
burnin = 1:5000
data = data.frame(y=y, X)


run = bpr::sample_bpr(y ~ . - 1, data = data,
                      iter = nrep, burnin = max(burnin),
                      prior = list(type="gaussian", b = rep(0,p), B = diag(p)*2), 
                      pars = list(max_dist = 1e+6),
                      state = mle.start)
run$sim$acceptance_rate
effSize_mh <- apply(run$sim$beta[-burnin,], 2, effectiveSize)
effSize_mh
run$sim$time


#----------------------------------------------------------------------------------#



library(TeachingDemos)
post_mean_mh = colMeans(run$sim$beta[-burnin,])
CI_mh = apply(run$sim$beta[-burnin,], 2, function(x) emp.hpd(x, conf=0.95))

var_area = 2:6
var_cre = 7:18
var_depth = c(19,23)
var_experiment = 20:22

### CRE line ###
order_cre = sort(post_mean_mh[var_cre], index.return = T, decreasing = T)
labels_cre = sapply(colnames(X)[var_cre], function(x) strsplit(x, "CRE.line")[[1]][2] ) 
labels_cre[1] = "Emx1"
labels_cre[2] = "Fezf2"
labels_cre[3] = "Nr5a1"
labels_cre[4] = "Ntsr1"
labels_cre[5] = "Pvalb"
labels_cre[6] = "Rbp4"
labels_cre[7] = "Rorb"
labels_cre[8] = "Scnn1a"
labels_cre[9] = "Slc17a7"
labels_cre[10] = "Sst"
labels_cre[11] = "Tlx3"
labels_cre[12] = "Vip"

### Area ###
order_area = sort(post_mean_mh[var_area], index.return = T, decreasing = T)

### Experiment ###
order_experiment = sort(post_mean_mh[var_experiment], index.return = T, decreasing = T)



library(ggplot2)
library(ggridges)
library(grid)
library(reshape2)
library(ggstance)
library(dplyr)


#----------------------# CRE #----------------------#
# Generating the dataset
x = data.frame(run$sim$beta[-burnin,var_cre[order_cre$ix]])
colnames(x) = labels_cre[order_cre$ix]
str(x)

library(tidyr)
data <- gather(x)
str(data)

# Summary statistics
mean_vec <- colMeans(x)
n        <- nrow(x)
left  <- CI_mh[1,var_cre[order_cre$ix]]
right <- CI_mh[2,var_cre[order_cre$ix]]
sum_stat <- data.frame("label" = labels_cre[order_cre$ix], left, mean_vec, right)
sum_stat$label = factor(sum_stat$label, levels = labels_cre[order_cre$ix])
data$key = factor(data$key, levels = labels_cre[order_cre$ix])


plot_cre =  ggplot() +
  geom_density_ridges(data=data, aes(x = value, y = key), alpha=0.2, scale=3.9, size = 0.2, bandwidth = 0.005) +
  geom_point(data = sum_stat, aes(mean_vec, label, color = mean_vec)) +
  geom_segment(data = sum_stat, aes(x = left, y = label, xend = right, yend = label, color = mean_vec),
               size = 0.8, 
               linetype = 1) +
  scale_colour_gradient( low = "blue", high = "red") +
  theme_light() +
  xlab("Coefficient") +
  ylab("Cre-line") +
  theme(legend.position = "none")
plot_cre



#----------------------# AREA #----------------------#
# Generating the dataset
labels_area = sapply(colnames(X)[var_area], function(x) strsplit(x, "area")[[1]][2] ) 
x = data.frame(run$sim$beta[-burnin,var_area[order_area$ix]])
colnames(x) = labels_area[order_area$ix]
str(x)

data <- gather(x)
str(data)
# Summary statistics
mean_vec <- colMeans(x)
n        <- nrow(x)
left  <- CI_mh[1,var_area[order_area$ix]]
right <- CI_mh[2,var_area[order_area$ix]]
sum_stat <- data.frame("label" = labels_area[order_area$ix], left, mean_vec, right)
sum_stat$label = factor(sum_stat$label, levels = labels_area[order_area$ix])
data$key = factor(data$key, levels = labels_area[order_area$ix])

plot_area =  ggplot() +
  geom_density_ridges(data=data, aes(x = value, y = key), alpha=0.2, scale=3.9, size = 0.2, bandwidth = 0.005) +
  geom_point(data = sum_stat, aes(mean_vec, label, color = mean_vec)) +
  geom_segment(data = sum_stat, aes(x = left, y = label, xend = right, yend = label, color = mean_vec),
               size = 0.8, 
               linetype = 1) +
  theme_light() +
  xlab("") +
  ylab("Area") +
  theme(legend.position = "none") +
  scale_colour_gradient2( low = "blue", mid = "#A1015D", high = "red", midpoint = 0)
plot_area




#----------------------# EXPERIMENT #----------------------#
# Generating the dataset
labels_experiment = sapply(colnames(X)[var_experiment], function(x) strsplit(x, "experiment")[[1]][2] ) 
x = data.frame(run$sim$beta[-burnin,var_experiment[order_experiment$ix]])
colnames(x) = labels_experiment[order_experiment$ix]
str(x)

data <- gather(x)
str(data)

# Summary statistics
mean_vec <- colMeans(x)
n        <- nrow(x)
left  <- CI_mh[1,var_experiment[order_experiment$ix]]
right <- CI_mh[2,var_experiment[order_experiment$ix]]
sum_stat <- data.frame("label" = labels_experiment[order_experiment$ix], left, mean_vec, right)
sum_stat$label = factor(sum_stat$label, levels = labels_experiment[order_experiment$ix])
data$key = factor(data$key, levels = labels_experiment[order_experiment$ix])

plot_experiment =  ggplot() +
  geom_density_ridges(data=data, aes(x = value, y = key), alpha=0.2, scale=3.9, size = 0.2, bandwidth = 0.005) +
  geom_point(data = sum_stat, aes(mean_vec, label, color = mean_vec)) +
  geom_segment(data = sum_stat, aes(x = left, y = label, xend = right, yend = label, color = mean_vec),
               size = 0.8, 
               linetype = 1) +
  theme_light() +
  xlab("") +
  ylab("Experiment") +
  theme(legend.position = "none") +
  scale_colour_gradient2( low = "blue", mid = "#A1015D", high = "red", midpoint = 0)
plot_experiment


