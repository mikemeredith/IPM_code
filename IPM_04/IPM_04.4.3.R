# Schaub & Kéry (2022) Integrated Population Models
# Chapter 4 : Components of integrated population models
# ------------------------------------------------------

# Run time approx. 3 mins

library(IPMbook) ; library(jagsUI)

# ~~~ this needs data created in section 4.4.1 ~~~
nbrood <- 1000           # Number of broods with young counted
brood.mean <- 1.5        # Average brood size
sd.brood <- 0.3          # log-linear brood random effect
set.seed(24)
expNyoung <- exp(log(brood.mean) + rnorm(nbrood, 0, sd.brood))
C <- rpois(nbrood, expNyoung)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 4.4 Models for productivity surveys
# ===================================

# 4.4.3 Zero-truncation in brood size data
# ----------------------------------------

# Create a variant of the data with partial zero truncation
set.seed(68)
C1 <- C                                   # Copy data set
table(C1)                                 # Remind ourselves of the actual data
prop.uncertain <- 0.6                     # Proportion of zeros dropped
zeros <- which(C1==0)
toss <- zeros[rbinom(n=length(zeros), size=1, prob=prop.uncertain) == 1]
C1 <- C1[-toss]                           # Toss out some of the zeros
table(C1)                                 # Frequency dist. of new brood size data
# C1
#  0   1   2   3   4   5   6   7   8
# 89 330 237 134  39  16   7   3   3

mean(C); mean(C1)                         # Compare sample means
# [1] 1.529
# [1] 1.782051

# Make another copy of C1 data and kick out remaining zeroes as well
C2 <- C1[C1 > 0]                          # Make a copy and toss out remaining zeroes

# Data bundle
jags.data <- list(C1=C1, C2=C2, nC1=length(C1), nC2=length(C2))
str(jags.data)
# List of 4
# $ C1 : int [1:858] 0 0 1 1 4 1 2 0 1 3 ...
# $ C2 : int [1:769] 1 1 4 1 2 1 3 1 1 1 ...
# $ nC1: int 858
# $ nC2: int 769

# Write JAGS model file
cat(file="model10.txt", "
model {
  # Priors and linear models
  rho1 ~ dunif(0, 3)                      # Expected brood size in model 1
  rho2 ~ dunif(0, 3)                      # Expected brood size in model 2

  # Likelihoods
  # Model 1: Standard Poisson GLM
  for (i in 1:nC1){
    C1[i] ~ dpois(rho1)
  }

  # Model 2: Poisson GLM with response truncated at one (i.e., no zeroes)
  for (i in 1:nC2){
    C2[i] ~ dpois(rho2)T(1,)
  }
}
")

# Initial values
inits <- function(){list(rho1=runif(1, 0.5, 2.5), rho2=runif(1, 0.5, 2.5))}

# Parameters monitored
parameters <- c("rho1", "rho2")

# MCMC settings
ni <- 60000; nb <- 10000; nc <- 3; nt <- 10; na <- 1000

# Call JAGS from R (ART 3 min), check convergence and summarize posteriors
out13 <- jags(jags.data, inits, parameters, "model10.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out13) # Not shown
print(out13, 3)
#              mean    sd     2.5%      50%    97.5% overlap0 f Rhat n.eff
# rho1        1.784 0.046    1.694    1.784    1.874    FALSE 1    1  9064
# rho2        1.580 0.053    1.481    1.579    1.685    FALSE 1    1  5698

# ~~~~ Compare the estimates of the two analyses (Figure 4.11) ~~~~
library(scales)
co <- viridis_pal(option='E')(20)[c(2, 18)]
cl <- c(alpha(co[1], 0.5), alpha(co[2], 0.5))
breaks <- hist(c(out13$sims.list$rho1, out13$sims.list$rho2),
  breaks=60, plot=FALSE)$breaks
op <- par(las=1, mar=c(4.5,5,1,1), cex=1.2)
hist(out13$sims.list$rho1, breaks=breaks, col=cl[1],
    xlab=expression(paste('Estimate of mean productivity (', rho,')')),
    xlim=c(1.35, 2), border=NA, ylim=c(0, 1500), main=NA)
hist(out13$sims.list$rho2, breaks=breaks, col=cl[2], border=NA, add=TRUE)
segments(1.5, 0, 1.5, 1300, lwd=2, col='red')
legend(x=1.325, y=1550, pch=rep(15,2), col=cl,
    legend=c('Ignoring truncation', 'Accounting for truncation'), bty='n', pt.cex=1.2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Number of zeroes in data set corrupted by partial zero-truncation
sum(C1==0)
# [1] 89

# Expected number of zeroes from model accommodating zero-truncation
round(nbrood * dpois(0, out13$mean$rho2))
# [1] 206

# True number of zeroes in the original data set C
sum(C==0)
# [1] 231
