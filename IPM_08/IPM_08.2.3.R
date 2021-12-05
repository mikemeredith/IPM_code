# Schaub & Kéry (2022) Integrated Population Models
# Chapter 8 : Integrated population models with density-dependence
# ----------------------------------------------------------------

library(IPMbook) ; library(jagsUI)

# ~~~ this uses the output from section 8.2.2 ~~~
load("Shrike-mod1.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 8.2 Density-dependence in red-backed shrikes
# ============================================

# 8.2.3 Assessing density-dependence at the population level
# ----------------------------------------------------------

N <- out1$sims.list$N
T <- ncol(N)
n <- nrow(N)

r <- sigma.r <- alpha.star <- beta.star <- sigma.star <- alpha <- beta <- sigma <- numeric(n)
N.star <- matrix(NA, nrow=n, ncol=T)
for (i in 1:n){                                       # loop over all MCMC draw
  if(round(i %% 500 ) == 0) {cat(paste('** Processing draw', i, '**\n'))}
  # 1. Estimate the growth parameters under the exponential model (i.e., without DD)
  z <- lm(log(N[i,-1] / N[i,-T]) ~ 1)
  r[i] <- coef(z)
  sigma.r[i] <- summary(z)$sigma
  # 2. Simulate population trajectories under the exponential growth model
  N.star[i,1] <- N[i,1]
  for (t in 1:(T-1)){
    N.star[i,t+1] <- N.star[i,t] * exp(r[i] + rnorm(1, 0, sigma.r[i]))
  } #t
  # 3. Fit the DD model to simulated population sizes from step 2
  z <- lm(log(N.star[i,-1]) ~ log(N.star[i,-T]))
  alpha.star[i] <- coef(z)[1]
  beta.star[i] <- coef(z)[2]-1
  sigma.star[i] <- summary(z)$sigma
  # 4. Fit DD model to estimated population sizes (from the IPM)
  z <- lm(log(N[i,-1]) ~ log(N[i,-T]))
  alpha [i] <- coef(z)[1]
  beta[i] <- coef(z)[2]-1
  sigma[i] <- summary(z)$sigma
} #i


# ~~~~ Code for Fig. 8.4 ~~~~
plot(density(beta), ylim=c(0, 20), xlim=c(-0.6, 0.2), main=NA,
    xlab=expression('Strength of density-dependence ('*beta*')'),
    ylab="Density", col="red", axes=FALSE)
lines(density(beta.star), col="black")
legend("topleft", lwd=c(1, 1), col=c("red", "black"),
    legend=c(expression('Observed ('*beta*')'),
    expression('Expected under exponential growth model (bias) ('*beta^'*'*')')), bty="n")
axis(1)
axis(2, las=1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mean(beta < beta.star)
# [1] 0.8665333
