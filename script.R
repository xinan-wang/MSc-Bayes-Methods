library(MCMCpack)
library(lattice)

obs.dat=dget("CompRank21.txt")
n <- length(obs.dat) #134
num_player <- 24

# plist contains the data for each player
plist <- vector(mode = "list", length = num_player)
for (i in 1:num_player){
  plist[[i]] <-  list('G'=c(),'S'=c(), 'RS'=c(), 'R'=c())
}
par(mfrow=c(1,1))
num_per_game <- c()
for (i in 1:134){
  game_i <- obs.dat[[i]]
  num_player_i <- length(game_i[[1]])
  for (j in 1:num_player_i){
    # the jth player in game i
    player_j <- game_i[[1]][j]
    # Record game index, seniority, relative seniority and result rank
    plist[[player_j]][['G']] <- c(plist[[player_j]][['G']],i)
    plist[[player_j]][['S']] <- c(plist[[player_j]][['S']],game_i[[2]][j])
    plist[[player_j]][['RS']] <- c(plist[[player_j]][['RS']],
                                   rank(game_i[[2]],ties.method = "min")[j])
    plist[[player_j]][['R']] <- c(plist[[player_j]][['R']],j)
  }
  num_per_game <- c(num_per_game,num_player_i)
}

# Obtain number of games each player participated
num_per_player <- c()
for (i in 1:num_player){
  num_per_player[i] = length(plist[[i]][[1]])
}

# Diagnostics
barplot(num_per_player,names.arg = 1:24, xlab = 'Player',
        main = 'Number of Games Played by Each Competitor')
hist(num_per_game, main = 'Distributions of Number of Players',
     xlab='Number of Players', ylab='Number of Games Count')
hist(num_per_player,breaks=seq(0,80,10), main = 'Distributions of Number of Games',
     xlab='Number of Games', ylab='Number of Player Count')


## Line chart of games player by player i
plot_player <- function(i){
  num_game <- length(plist[[i]][[1]])
  rank_range <- range(c(plist[[i]][['R']],plist[[i]][['RS']]))
  plot(plist[[i]][['R']], ylim = rev(c(rank_range[1], rank_range[2])),
       col='orange',lwd=2, type="b" , xlab ='Game Index', 
       ylab="Ranks", axes=FALSE, 
       main=paste('Ranks of Player',i, 'at', num_game, 'Games'))
  lines(plist[[i]][['RS']], col='blue',lwd=2, type="b")
  axis(1, at=1:num_game, labels=plist[[i]][['G']]) 
  axis(2, at=rank_range[1]:rank_range[2])
  legend('bottomleft', c("Game Result", "Relative Seniority"), box.lwd = 0.2, 
         col = c("blue", "orange"),bty = 'n', lwd=c(2,2))
  
}

plot_player(5)
plot_player(7)
plot_player(10)
plot_player(15)


## Scatter plot
a = 0.3

plot(jitter(plist[[1]][['RS']], amount=a),jitter(plist[[1]][['R']], amount=a), 
     pch=19, ylab = 'Rank of Game Result', xlab = 'Relative Seniority Rank',
     ylim = c(1,10), xlim = c(1,10), axes = FALSE,
     main = 'Game Result v.s. Relative Seniority Rank')
for (k in 2:24){
  points(jitter(plist[[k]][['RS']],amount=a),
         jitter(plist[[k]][['R']],amount=a), pch=19)
}
axis(1,at=1:10)
axis(2,at=1:10)

plot(jitter(plist[[1]][['S']], amount=a),jitter(plist[[1]][['R']], amount=a), 
     pch=19, ylab = 'Rank of Game Result', xlab = 'Absolute Seniority Rank',
     xlim = c(1,16), ylim = c(1,10), axes = FALSE,
     main = 'Game Result v.s. Absolute Seniority Rank')
for (k in 2:24){
  points(jitter(plist[[k]][['S']],amount=a),
         jitter(plist[[k]][['R']],amount=a), pch=19)
}
axis(1,at=1:16)
axis(2,at=1:10)

## Grouped barplot
rank_table <- matrix(rep(0,9),nrow=3,ncol=3)
for (i in 1:length(plist)){
  playeri <- plist[[i]]
  for (j in 1:2){
    for (k in 1:2){
      rank_table[j,k] = rank_table[j,k] +
        sum((playeri$RS==j)&(playeri$R==k))
    }
  }
  
  for (l in 1:2){
    rank_table[l,3] = rank_table[l,3] +
      sum((playeri$RS==l)&(playeri$R>=3))
    rank_table[3,l] = rank_table[3,l] +
      sum((playeri$RS>=3)&(playeri$R==l))
  }
  rank_table[3,3] = rank_table[3,3] + 
    sum((playeri$RS>=3)&(playeri$R>=3))
}
rank_table

x<-barplot(rank_table, beside=T, ylim = c(0,70), xlab="Relative Seniority Rank",
        ylab = 'Count', names.arg = c('1st','2nd','>=3rd'), legend=c('1st','2nd','>=3rd'), 
        col = c('azure4','antiquewhite3','cornsilk2'),
        args.legend = list(title="Game Result", cex=0.8),
        main = 'Number of Players Grouped by Seniority Rank')
text(x,rank_table+3,labels=as.character(rank_table))


## Observational model
obs_prob <- function(d, lambda, beta = rep(0,24)) {
 out <- 1
 for (i in 1:length(d)) {
   players <- d[[i]]$o
   srank <- d[[i]]$e
   mi <- length(players)
   for (j in 1:mi){
     # Multiply by likelihood
     out <- out * exp(lambda[players[j]]+beta[srank[j]])/
       sum(exp(lambda[players[j:mi]]+beta[srank[j:mi]]))
   }
 }
 return(out)
}

## Construct bridge estimator
b_est <- function(d, l1, l2, b1 = matrix(rep(0,24*1000),ncol=24), 
                  b2 = matrix(rep(0,24*1000),ncol=24)){
  n1 = 0; n2 = 0
  for (i in 1:1000){
    n1 = n1 + sqrt(obs_prob(d, l1[i,], b1[i,])) # numerator
    n2 = n2 + (sqrt(obs_prob(d, l2[i,], b2[i,])))^(-1) # Denominator
  }
  out <- n1/n2
  return(out)
}


#initialize state for MCMC with lambda and beta
lambda <- rnorm(24)
beta <- rbeta(24, shape1 = 1 +1/(1:24), shape2 = 2 - 1/(1:24))
oldp=obs_prob(obs.dat, lambda, beta)

K=1e5; SS=100; 
out.mcmc=matrix(NA,K/SS,24*2)
for (k in 1:K) {
  # Choose state to update
  state = sample(1:48,1)
  if (state<=24){
    lambdap = lambda
    lambdap[state] = rnorm(1)
    betap = beta
  }
  else {
    betap = beta
    betap[state-24] = rbeta(1, shape1 = 1 +1/(1:24), shape2 = 2 - 1/(1:24))
    lambdap = lambda
  }

  newp=obs_prob(obs.dat, lambdap, betap)
   
  # Update probability  
  alpha = newp/oldp
  if (runif(1)<alpha) {  
    lambda=lambdap; beta=betap; 
    oldp=newp
  }
  # Subsampling
  if (k%%SS==0) {
    out.mcmc[k/SS,]=c(lambda,beta)
  }
}

# sample prior, obtain posterior
lambda1 = matrix(rnorm(24*1000),ncol=24)
lambda2 = out.mcmc[,1:24]
beta1 = matrix(nrow=1000,ncol=24)
for (j in 1:1000){
  beta1[j,] = rbeta(24, shape1 = 1 +1/(1:24), shape2 = 2 - 1/(1:24))
}
beta2 = out.mcmc[,25:48]

#Compute bridge estimates
p1 <- b_est(obs.dat, lambda1, lambda2, beta1, beta2)
p1
# lambda  normal
# beta  Beta(5-rank/6, 1+rank/6)  Beta((5-rank/6)/2, (1+rank/6)/2) 
#         1.672759e-84             8.034321e-87
#        Beta(1+1/rank, 2-1/rank)   Beta(2+2/rank, 4-2/rank)
#        5.941653e-84               8.116373e-89


out.mcmc <- as.mcmc(out.mcmc)
par(mfrow=c(1,1))
#correlation can be a problem, only mild correlation here
levelplot(out.mcmc,scales = list(x = list(rot = 45)), main='Correlation Plot') 
splom(as.matrix(out.mcmc[,c(7,9,10,11)]))  #some of the more strongly correlated cases
effectiveSize(out.mcmc)  #these could be bigger closer to 1000
hist(effectiveSize(out.mcmc), breaks=seq(0,1000,100),
     xlab='Effective Size',ylab='Count',main='Histogram of Effective Size')
par(mfrow=c(2,4),oma=c(1,1,1,1)); #acf plots
for (i in 1:dim(out.mcmc)[2]) {
  par(mai=0.2*c(1,1,1,1)); 
  plot(acf(out.mcmc[,i],lag.max=200,plot=F),
       type='l',ann=F,xaxp=c(0,200,2),yaxp=c(0,1,1)); 
  text(50,0.8,colnames(out.mcmc)[i])
}

#initialize state for MCMC with lambda only
lambda0 <- rnorm(24)
oldp0=obs_prob(obs.dat, lambda0)

K=1e5; SS=100; 
out.mcmc0=matrix(NA,K/SS,24)
for (k in 1:K) {
  state = sample(1:24,1)
  lambdap0 = lambda0
  lambdap0[state] = rnorm(1)
  
  newp0=obs_prob(obs.dat, lambdap0)
  
  alpha = newp0/oldp0
  if (runif(1)<alpha) {  
    lambda0=lambdap0
    oldp0=newp0
  }
  if (k%%SS==0) {
    out.mcmc0[k/SS,]=lambda0
  }

}



lambda01 = matrix(rnorm(24*1000),ncol=24)
lambda02 = out.mcmc0
p0 <- b_est(obs.dat, lambda01, lambda02)
p0
# lambda normal         t(2)                  t(5)          
#        1.347753e-82   1.006573e-92          2.860277e-84

out.mcmc0 <- as.mcmc(out.mcmc0)
par(mfrow=c(1,1))
levelplot(out.mcmc0,scales = list(x = list(rot = 45)), main='Correlation Plot')
splom(as.matrix(out.mcmc0[,c(7,10,11)]))  #some of the more strongly correlated cases
effectiveSize(out.mcmc0) 
hist(effectiveSize(out.mcmc0), breaks=seq(0,1000,100),
     xlab='Effective Size',ylab='Count',main='Histogram of Effective Size')
par(mfrow=c(2,4),oma=c(1,1,1,1));   #acf plots
for (i in 1:dim(out.mcmc0)[2]) {
  par(mai=0.2*c(1,1,1,1)); 
  plot(acf(out.mcmc0[,i],lag.max=200,plot=F),
       type='l',ann=F,xaxp=c(0,200,2),yaxp=c(0,1,1)); 
  text(50,0.8,colnames(out.mcmc0)[i])
}

## Predictive Distribution

library(combinat)

game <- obs.dat[[68]]
players <- game[[1]] # get player list
n_player <- length(players)
players <- players[-2]
# remove player 7 first
comb <-  permn(players)
prob = 0
for (i in 1:length(comb)){
  # Add player 7 as first
  ranki = c(7, comb[[i]])
  for (j in 1:1000){
    out = 1e-3
    # Posterior sample of lambda
    lambdaj = lambda02[j,]
    for (k in 1:n_player){
      # Multiply by the kth component
      out <- out * exp(lambdaj[ranki[k]])/
        sum(exp(lambdaj[ranki[k:n_player]]))
    }
    prob = prob + out
  }
}

prob
# 0.1783173

