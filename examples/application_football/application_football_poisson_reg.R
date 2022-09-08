rm(list=ls())
library(Rcpp)
library(RcppDist)
library(coda)
library(TeachingDemos)
library(bpr)

data = read.csv("20-21.csv", header = TRUE)

y_1 = data$FTHG
y_2 = data$FTAG
y = c(y_1,y_2)

X_team = factor(c(data$HomeTeam, data$AwayTeam))
X_team_tmp = model.matrix(~ X_team)
X_tmp = data[,-c(1:23)] # solo scommesse

away = c(rep(0, nrow(X_tmp)), rep(1, nrow(X_tmp)))
X = cbind(X_team_tmp, rbind(X_tmp, X_tmp)) 

sum(is.na(y))
sum(is.na(X))
X[is.na(X)] = 0

p = ncol(X)
n = nrow(X)
n
p

table(y) 
X = as.matrix(X)

est_tau = function(n)
{
  p_n = 80 #sum(beta != 0)
  return( p_n/n * sqrt(log(n/p_n)) )
}
tau = est_tau(n)

nrep = 10000
burnin = 1:5000

data = data.frame("y" = y, X)


########### horseshoe ###########
run <- bpr::sample_bpr(y ~ . -1, data = data, iter = nrep, 
                       burnin = 5000, 
                       prior = list(type = "horseshoe", tau = tau),
                       pars = list(max_dist = 0.7, max_dist_burnin = 0.5))

run$sim$acceptance_rate
colMeans(run$sim$beta[-burnin,])
CI_mh = apply(run$sim$beta[-burnin,], 2, function(x) emp.hpd(x, conf = 0.9))

which(apply(CI_mh, 2, function(x) ((x[1]<0)&(x[2]<0))|((x[1]>0)&(x[2]>0)) ) )
nonzero_h = which(apply(CI_mh, 2, function(x) ((x[1]<0)&(x[2]<0))|((x[1]>0)&(x[2]>0)) ) )



########### gaussian ###########
run2 <- bpr::sample_bpr(y ~ . -1, data = data, iter = nrep, state = colMeans(run$sim$beta[-burnin,]),
                       prior = list(type = "gaussian", B = diag(p)*2),
                       pars = list(max_dist = 0.1))

run2$sim$acceptance_rate
colMeans(run2$sim$beta[-burnin,])

CI_mh2 = apply(run2$sim$beta[-burnin,], 2, function(x) emp.hpd(x, conf = 0.9))
which(apply(CI_mh2, 2, function(x) ((x[1]<0)&(x[2]<0))|((x[1]>0)&(x[2]>0)) ) )
nonzero_nh = which(apply(CI_mh2, 2, function(x) ((x[1]<0)&(x[2]<0))|((x[1]>0)&(x[2]>0)) ) )



## plot
order = sort(colMeans(run2$sim$beta[-burnin,]), index.return=T, decreasing = T)
run$sim$beta = run$sim$beta[-burnin,order$ix]
run2$sim$beta = run2$sim$beta[-burnin,order$ix]

df = data.frame("prior" = c(rep("Horseshoe",p),rep("Gaussian",p)) ,
                "variable" = c(1:p, (1:p)+0.2),
                "lower" = c(CI_mh[1,],CI_mh2[1,]),
                "upper" = c(CI_mh[2,],CI_mh2[2,]),
                "mean" = c(colMeans(run$sim$beta), colMeans(run2$sim$beta)),
                "iszero" = rep(0,2*p)
)
df$iszero[df$prior=="Horseshoe"][nonzero_h] = 1
df$iszero[df$prior=="Gaussian"][nonzero_nh] = 1


library(ggplot2)

selezione = (df$variable>p/2)&(df$variable<p)
p1 = ggplot(df[selezione,]) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_pointrange(aes(x = variable, y = mean, 
                      ymax = upper, ymin = lower, color = prior, shape = prior)) +
  geom_point(data = df[(selezione)&(df$iszero==1)&(df$prior=="Gaussian"),], aes(x = variable, y = mean, color = prior), 
             shape = 16, size=2 ) +
  geom_point(data = df[(selezione)&(df$iszero==1)&(df$prior=="Horseshoe"),], aes(x = variable, y = mean, color = prior), 
             shape = 17, size = 2 ) +
  scale_shape_manual(values=c(1,2)) +
  coord_flip() +
  theme_light() +
  xlab("Variable") +
  ylab("Value") +
  theme(legend.position = "none") +
  theme(
    # panel.grid.major.y = element_blank()
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  ) +
  scale_x_continuous(breaks = sort(unique(round(df[(selezione)&(df$iszero==1),]$variable))) + 0.1, 
                     labels = colnames(X)[order$ix][sort(unique(round(df[(selezione)&(df$iszero==1),]$variable)))]) 


selezione = df$variable<p/2
p3 = ggplot(df[selezione,]) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_pointrange(aes(x = variable, y = mean, 
                      ymax = upper, ymin = lower, color = prior, shape = prior)) +
  geom_point(data = df[(selezione)&(df$iszero==1)&(df$prior=="Gaussian"),], aes(x = variable, y = mean, color = prior), 
             shape = 16, size = 2 ) +
  geom_point(data = df[(selezione)&(df$iszero==1)&(df$prior=="Horseshoe"),], aes(x = variable, y = mean, color = prior), 
             shape = 17, size =2 ) +
  scale_shape_manual(values=c(1,2)) +
  coord_flip() +
  theme_light() +
  xlab("") +
  ylab("Value") +
  theme(legend.position = "none") +
  theme(
    # panel.grid.major.y = element_blank()
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")
  ) +
  scale_x_continuous(breaks = sort(unique(round(df[(selezione)&(df$iszero==1),]$variable))) + 0.1, 
                     labels = colnames(X)[order$ix][sort(unique(round(df[(selezione)&(df$iszero==1),]$variable)))]) 


library(gridExtra)
grid.arrange(p1,p3, ncol = 2)








