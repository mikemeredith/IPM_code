# Schaub & Kéry (2022) Integrated Population Models
# Chapter 2 : Bayesian statistical modeling using JAGS
# ----------------------------------------------------

# Run time approx. 1 min

library(IPMbook) ; library(jagsUI)

# ~~~ need following code from 2.7.1 ~~~~~~~~
original.elev <- c(500, 400, 700, 500, 600, 800, 1000, 900, 1000, 900)
elev <- (original.elev - mean(original.elev)) / 1000    # Center and scale
C <- c(6, 10, 2, 7, 4, 1, 1, 2, 0, 0)                   # Counts
times.bigger <- 100                                     # Results in sample size of 1000
newC <- rep(C, times.bigger)
newElev <- rep(elev, times.bigger)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.7 Using JAGS to fit simple statistical models from R: GLMs and GLMMs
# ======================================================================

# 2.7.2 Bernoulli generalized linear models
# -----------------------------------------

# Turn counts into detection/nondetection data
y <- as.numeric(newC > 0); table(y)       # Not shown

# Bundle data
jags.data <- list(y=y, elev=newElev, n=length(y))
str(jags.data)

# List of 3
# $ y   : num [1:1000] 1 1 1 1 1 1 1 1 0 0 ...
# $ elev: num [1:1000] -0.23 -0.33 -0.03 -0.23 -0.13 0.07 0.27 0.17...
# $ n   : int 1000

# Write JAGS model file
cat(file="model2.txt", "
model {
  # Priors and linear models
  alpha <- logit(mean.theta)                      # Intercept on logit link scale
  mean.theta ~ dunif(0, 1)                        # Intercept on prob. scale
  beta ~ dnorm(0, 1.0E-06)                        # Slope on logit link scale

  # Likelihood of the Bernoulli GLM
  for (i in 1:n){
    y[i] ~ dbern(theta[i])                        # Stochastic part of response
    logit(theta[i]) <- alpha + beta * elev[i]     # Link function and lin. pred.
    # theta[i] <- ilogit(alpha + beta * elev[i])  # Same written differently
  }
}
")

# Initial values
inits <- function() list(mean.theta=runif(1), beta=rnorm(1))

# Parameters monitored
parameters <- c("alpha", "beta", "mean.theta")

# MCMC settings
ni <- 50000; nb <- 10000; nc <- 3; nt <- 10; na <- 1000

# Call JAGS from R (ART 1 min), check convergence and summarize posteriors
out3 <- jags(jags.data, inits, parameters, "model2.txt", n.iter=ni, n.burnin=nb, n.chains=nc,
    n.thin=nt, n.adapt=na, parallel=TRUE)
traceplot(out3)                 # May need jagsUI::traceplot
print(out3, 3)
#               mean    sd    2.5%     50%   97.5% overlap0 f Rhat n.eff
# alpha        2.972 0.231   2.536   2.967   3.441    FALSE 1    1 12000
# beta       -12.717 1.081 -14.890 -12.699 -10.686    FALSE 1    1  8781
# mean.theta   0.950 0.011   0.927   0.951   0.969    FALSE 1    1 12000
# deviance   632.603 1.997 630.651 631.980 638.029    FALSE 1    1  6411

# Get frequentist MLEs for comparison
summary(glm(y ~ newElev, family='binomial'))

# Coefficients:
#             Estimate Std. Error z value Pr(>|z|)
# (Intercept)   2.9922     0.2363   12.66   <2e-16 ***
# newElev     -12.8055     1.1039  -11.60   <2e-16 ***

# Inspect the posterior draws available after our run of JAGS
str(out3$sims.list)
# List of 4
# $ alpha     : num [1:12000] 2.83 2.54 2.95 2.79 2.82 ...
# $ beta      : num [1:12000] -11.5 -11.9 -12 -11.9 -11.7 ...
# $ mean.theta: num [1:12000] 0.944 0.927 0.95 0.942 0.944 ...
# $ deviance  : num [1:12000] 633 640 632 631 632 ...

# Prepare covariate for prediction
# Apply identical scaling as for the covariate used in the model fitting
mean(original.elev)                                     # Remind ourselves of the mean in the original cov.
original.elev.pred <- seq(400, 1000, length.out=1000)   # For prediction
elev.pred <- (original.elev.pred - mean(original.elev)) / 1000

# Prepare an array for holding the predictions of theta
pred1 <- array(NA, dim=c(1000, 12000))
# Fill the array: get posterior distribution of predictions
for (i in 1:12000){
  pred1[,i] <- plogis(out3$sims.list$alpha[i] + out3$sims.list$beta[i] * elev.pred)
}
# Summarize posterior distribution of theta_hat
pm <- apply(pred1, 1, mean)                             # Posterior mean
psd <- apply(pred1, 1, sd)                              # Posterior SD (like prediction SE)
CRI <- apply(pred1, 1, quantile, probs=c(0.025, 0.975)) # 95% CRI

# ~~~~ extra code for Figure 2.13 ~~~~
op <- par(las=1, mar=c(4.1, 4.5, 2, 2))
plot(original.elev.pred, pm, xlab='Site elevation (m)',
    ylab=expression(paste('Apparent occupancy (', hat(theta), ')')),
    type='n', ylim=c(0.3, 1))
polygon(c(original.elev.pred, rev(original.elev.pred)), c(CRI[1,],
    rev(CRI[2,])), col='grey', border=NA)
lines(original.elev.pred, pm, type='l', lty=1, lwd=3, col='blue')
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Prepare another covariate with different range
mean(original.elev)                                     # Remind ourselves of the mean in the original cov.
alt.elev.pred <- 601:1600                               # For prediction
elev.pred <- (alt.elev.pred - mean(original.elev)) / 1000

# Get posterior distribution of predictions of theta
pred2 <- array(NA, dim=c(length(elev.pred), 12000))
for(i in 1:12000){
  pred2[,i] <- plogis(out3$sims.list$alpha[i] + out3$sims.list$beta[i] * elev.pred)
}

# Get posterior mean of theta (theta_hat)
pm <- apply(pred2, 1, mean)                             # Posterior mean
# Compute probability that theta_hat < 0.01
prob <- apply(pred2, 1, function(x) mean((x - 0.01) < 0))
# Compute minimum elevation at which prob > 0.95
min.elev <- min(alt.elev.pred[prob > 0.95])
print(min.elev)
# [1] 1391

# ~~~~ additional code for Figure 2.14 ~~~~
plot(alt.elev.pred, pm, type='l', lwd=2, col='blue', frame=FALSE,
    xlab='Site elevation (m)', ylab='Probability')
lines(alt.elev.pred, prob, lwd=2, col='black')
abline(h=0.95, lwd=1, lty=3)
abline(v=min.elev, lwd=2, lty=1, col='red')
legend(x=600, y=0.3, lwd=rep(2, 2), col=c('blue', 'black'), bty='n',
    legend=c(expression(paste('Occupancy (', hat(theta),')')),
        expression(paste(italic(p),'(', hat(theta), ')<0.01'))))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
