# Schaub & Kéry (2022) Integrated Population Models
# Chapter 2 : Bayesian statistical modeling using JAGS
# ----------------------------------------------------

library(IPMbook)

# 2.2 Parametric statistical modeling
# ===================================

# 2.2.2 Parametric statistical models for inference about chance processes
# ------------------------------------------------------------------------

# Read the data for the tadpole experiment into R
N <- 50     # Number of tadpoles released initially
y <- 20     # Number of tadpoles detected later

# 2.3 Maximum likelihood estimation in a nutshell
# ===============================================

# Data for Fig. 2.3
all.possible <- 0:50                              # All possible values
pmf1 <- dbinom(all.possible, size=50, prob=0.1)   # pmf 1
pmf2 <- dbinom(all.possible, size=50, prob=0.5)   # pmf 2
pmf3 <- dbinom(all.possible, size=50, prob=0.9)   # pmf 3

# ~~~~ additional code for the plot ~~~~
op <- par(mfrow=c(1, 3), mar=c(5, 5, 4, 1), cex.lab=1.5, cex.axis=1.5, cex.main=2, las=1)
plot(all.possible, pmf1, type="h", lend="butt", lwd=3, frame=FALSE,
    xlab="Counts (y)", ylab="Probability of y", main=expression(paste(theta, " = 0.1")))
abline(v=20, col="blue", lwd=2)
plot(all.possible, pmf2, type="h", lend="butt", lwd=3, frame=FALSE,
    xlab="Counts (y)", ylab="", main=expression(paste(theta, " = 0.5")))
abline(v=20, col="blue", lwd=2)
plot(all.possible, pmf3, type="h", lend="butt", lwd=3, frame=FALSE,
    xlab="Counts (y)", ylab="", main=expression(paste(theta, " = 0.9")))
abline(v=20, col="blue", lwd=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Use RNG to obtain binomial density for pmf3
hist(rbinom(10^6, size=50, prob=0.9), freq=FALSE, xlim=c(0, 50)) # Not shown

# Brute-force search for MLE of theta for the tadpole data set (Fig. 2.4)
try.theta <- seq(0, 1, by=0.01)                     # Values to try out
like <- dbinom(20, 50, try.theta, log=FALSE)        # Likelihood
loglike <- dbinom(20, 50, try.theta, log=TRUE)      # Log-Likelihood
negloglike <- -dbinom(20, 50, try.theta, log=TRUE)  # Negative log-likelihood (NLL)

op <- par(mfrow=c(1, 3), mar=c(5, 5, 4, 1), cex.lab=1.5, cex.axis=1.5)
plot(x=try.theta, y=like, xlab=expression(paste("Detection probability (", theta, ")")),
    ylab="Likelihood", frame=FALSE, type="p", pch=16, col="black")
abline(v=try.theta[which(like == max(like))], col="red", lwd=3)
plot(x=try.theta, y=loglike, xlab=expression(paste("Detection probability (", theta, ")")),
    ylab="Log-Likelihood", frame=FALSE, type="p", pch=16, col="black")
abline(v=try.theta[which(loglike == max(loglike))], col="red", lwd=3)
plot(x=try.theta, y=negloglike, xlab=expression(paste("Detection probability (", theta, ")")),
    ylab="Negative log-Likelihood", frame=FALSE, type="p", pch=16, col="black")
abline(v=try.theta[which(negloglike == min(negloglike))], col="red", lwd=3)
par(op)

# 2.4 Bayesian inference
# ======================

theta.vals <- seq(0, 1, 0.001)

# Define likelihood function (same for all four analyses)
like <- dbinom(20, 50, theta.vals)

# Define four prior distributions
prior0 <- dbeta(theta.vals, 1, 1)
prior1 <- dbeta(theta.vals, 4, 6)
prior2 <- dbeta(theta.vals, 40, 60)
prior3 <- dbeta(theta.vals, 60, 40)

# Derive four posterior distributions
post0 <- dbeta(theta.vals, 20 + 1, 30 + 1)
post1 <- dbeta(theta.vals, 20 + 4, 30 + 6)
post2 <- dbeta(theta.vals, 20 + 40, 30 + 60)
post3 <- dbeta(theta.vals, 20 + 60, 30 + 40)


# ~~~~ additional code for plotting Fig 2.5 ~~~~
sc.like <- like * (50 + 1)     # Scale likelihood. Multiplied by n + 1
    # because the area under the curve for trial size 1 is 1/(n + 1).
library(scales)
co <- viridis_pal(option="E")(20)[c(18, 11, 2)]
lwd <- 3; cx <- 1.5
op <- par(mfrow=c(2, 2), mar=c(5, 5, 4, 2), cex.axis=cx, cex.lab=cx, cex.main=cx)

# Analysis 1 with vague prior
plot(theta.vals, post0, type ="l", col=co[3], xlab="",
    ylab="Scaled likelihood or density", las=1, frame=FALSE, lwd=2, ylim=c(0, 10))
mtext("Vague prior", side=3, line=0.5, font=2)
lines(theta.vals, sc.like, lty=2, lwd=2, col=co[2])
lines(theta.vals, prior0, lwd=2, col=co[1])
legend(0.5, 10, c("Prior dist.", "Likelihood function", "Posterior dist."),
    col=co, lty=1, lwd=2, bty="n")

# Analysis 2 with informative prior 1
plot(theta.vals, post1, type="l", lwd=2, col=co[3], xlab="", ylab="", las=1,
    frame=FALSE, ylim=c(0, 10))
mtext("Informative prior 1", side=3, line=0.5, font=2)
lines(theta.vals, sc.like, lwd=2, col=co[2])
lines(theta.vals, prior1, lty=1, lwd=2, col=co[1])

# Analysis 3 with informative prior 2
plot(theta.vals, post2, type="l", lwd=2, col=co[3], xlab= expression(theta),
    ylab="Scaled likelihood or density", las=1, frame=FALSE, ylim=c(0, 10))
mtext("Informative prior 2", side=3, line=0.5, font=2)
lines(theta.vals, sc.like, lwd=2, col=co[2])
lines(theta.vals, prior2, lty=1, lwd=2, col=co[1])

# Analysis 4 with informative prior 3
plot(theta.vals, post3, type="l", lwd=2, col=co[3], xlab= expression(theta),
    ylab='', las=1, frame=FALSE, ylim=c(0, 10))
mtext("Informative prior 3", side=3, line=0.5, font=2)
lines(theta.vals, sc.like, lwd=2, col=co[2])
lines(theta.vals, prior3, lty=1, lwd= 2, col=co[1])
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ additional code  for Fig 2.6 ~~~~
# How information in the data changes our state of knowledge
co <- c("grey92", "grey60", "black")
xylim <- 0:1
# plot(theta.vals, sc.prior1, type="l", lwd=5, col=co[1],  xlab=expression(theta),
plot(theta.vals, prior0, type="l", lwd=5, col=co[1],  xlab=expression(theta),
    ylab="", xlim=xylim, ylim=xylim, las=1, frame=FALSE, axes=FALSE)
axis(1)
sc.post1 <- post1/max(post1)
lines(theta.vals, sc.post1, col=co[2], lwd=5)
abline(v=0.4, col=co[3], lwd=5)
legend(0.6, 0.9, c("Complete ignorance", "Improved state of knowledge",
    "Certainty (fixed parameter)"), col=co, lty=1, lwd=5, bty="n")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 2.5 Bayesian computation
# ========================

# Choose initial value for logit(theta) and tuning parameters
ltheta1 <- 1                      # Initial value for tadpole detection prob.
sigma_prop <- 1                   # SD of Gaussian proposal distribution

# Array to hold the MCMC samples
ltheta <- numeric()

# Initial value becomes first (and 'current') value in the chain
ltheta[1] <- ltheta1
ltheta                            # Our posterior sample up to now (not shown)

# Randomly perturb the current value
set.seed(1)                       # To initalize your RNGs identically to ours
( ltheta_star <- rnorm(1, ltheta[1], sigma_prop) )
# [1] 0.3735462

# Compute likelihood times prior evaluated for the proposed new value of ltheta
( pd_star <- dbinom(20, 50, plogis(ltheta_star)) * dbeta(plogis(ltheta_star), 1, 1) )
# [1] 0.002716919

# Compute likelihood times prior evaluated for the current value of ltheta
( pd_1 <- dbinom(20, 50, plogis(ltheta[1])) * dbeta(plogis(ltheta[1]), 1, 1) )
# [1] 6.951277e-07

# Compute posterior density ratio R
( R <- pd_star / pd_1 )
# [1] 3908.518

# Add theta_star into MCMC sample
ltheta[2] <- ltheta_star
ltheta                            # Our posterior sample up to now (not shown)

( ltheta_star <- rnorm(1, ltheta[2], sigma_prop) )
# [1] 0.5571895

pd_star <- dbinom(20, 50, plogis(ltheta_star)) * dbeta(plogis(ltheta_star), 1, 1)
pd_t <- dbinom(20, 50, plogis(ltheta[2])) * dbeta(plogis(ltheta[2]), 1, 1)
( R <- pd_star / pd_t )
# [1] 0.1398872

( keep.ltheta_star <- rbinom(1, 1, R) )
# [1] 0

ltheta[3] <- ltheta[2]
ltheta                            # Our posterior sample up to now (not shown)

# Iteration 4 to T of RW-MH algorithm
T <- 60000                        # Choose chain length
for (t in 4:T){                   # Continue where we left off
  ltheta_star <- rnorm(1, ltheta[t-1], sigma_prop)
  pd_star <- dbinom(20, 50, plogis(ltheta_star)) * dbeta(plogis(ltheta_star), 1, 1)
  pd_t <- dbinom(20, 50, plogis(ltheta[t-1])) * dbeta(plogis(ltheta[t-1]), 1, 1)
  R <- min(1, pd_star / pd_t)     # Note more general solution here
  keep.ltheta_star <- rbinom(1, 1, R)
  ltheta[t] <- ifelse(keep.ltheta_star == 1, ltheta_star, ltheta[t-1])
}


# ~~~~ code for Fig. 2.7 (left) ~~~~
op <- par(mfrow=c(1, 3), mar=c(6, 7, 6, 3), cex.lab=2, cex.axis=2, cex.main=2, las=1)
plot(1:10, plogis(ltheta[1:10]), xlab='Iteration', ylab=expression(theta),
    type='l', frame=FALSE, lwd=3, main='First ten iterations')
abline(h=0.4, col='red', lwd=2)             # The maximum likelihood estimate
abline(h=mean(plogis(ltheta[1:10])), lty=2, col='blue', lwd=2)

# Update trace-plot of time series of posterior draws (Fig. 2-7 middle)
plot(1:T, plogis(ltheta), xlab='Iteration', ylab= expression(theta), type='l',
    frame=FALSE, lwd=1, main='All iterations')
abline(h=0.4, col='red', lwd=3)             # The maximum likelihood estimate
abline(h=mean(plogis(ltheta)), lty=2, col='blue', lwd=3)

# Plot histogram of posterior samples of tadpole detection probability
# Fig. 2-7 right
hist(plogis(ltheta), breaks=50, col='lightgray', xlab=expression(theta),
    main=expression(bold(paste('Posterior distribution of ', theta))), border=NA)
abline(v=0.4, col='red', lwd=3)             # The maximum likelihood estimate
abline(v=mean(plogis(ltheta)), lty=2, col='blue', lwd=3) # Posterior mean
par(op)
# ~~~~~~~~~~~~~~~~~~~~~

library(IPMbook)
out <- demoMCMC(y=20, N=50, niter=25000, mu.ltheta=0, sd.ltheta=100, prop.sd=1, init=0)
    # produces Fig 2.8

# Show convergence
tmp <- demoMCMC(y=20, N=50, niter=2500, mu.ltheta=0, sd.ltheta=100, prop.sd=0.1, init=10)

# No convergence within 2500 iterations
tmp <- demoMCMC(y=20, N=50, niter=2500, mu.ltheta=0, sd.ltheta=100, prop.sd=0.1, init=100)

# But convergence is reached after about 3k iterations
tmp <- demoMCMC(y=20, N=50, niter=25000, mu.ltheta=0, sd.ltheta=100, prop.sd=0.1, init=100)

# ... and you get convergence within 2500 iters with longer step length
tmp <- demoMCMC(y=20, N=50, niter=2500, mu.ltheta=0, sd.ltheta=100, prop.sd=1, init=100)


# ~~~~ extra code to explore step size ~~~~
# Very, very small step size: very inefficient MCMC
str(out <- demoMCMC(prop.s = 0.01))

# Very small step size: fairly inefficient
str(out <- demoMCMC(prop.s = 0.1))

# Larger than default step size: efficiency goes down
str(out <- demoMCMC(prop.s = 10))

# Much larger step size..... brrrrr!
str(out <- demoMCMC(prop.s = 100))

# Brutally large step size..... ACH!
str(out <- demoMCMC(prop.s = 1000))

# Default step size: pretty good for this case
str(out <- demoMCMC(prop.s = 1))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

str(out)

# List of 7
# $ y        : num 20
# $ N        : num 50
# $ mu.ltheta: num 0
# $ sd.ltheta: num 100
# $ prop.sd  : num 1
# $ ltheta   : num [1:25000] 0 0 -0.0519 -0.0519 -0.7637 ...
# $ acc.prob : num 0.34

# Measures of central tendency to serve as point estimates
library(MCMCglmm)                                   # For function for mode
mean(plogis(out$ltheta))                            # Posterior mean
# [1] 0.4007489                                     # Your result will differ slightly

median(plogis(out$ltheta))                          # Posterior median
# [1] 0.4001582                                     # Your result will differ

posterior.mode(mcmc(plogis(out$ltheta)))            # Posterior mode (ditto)
#      var1
# 0.4152145

# Measures of spread:
# - Bayesian 'variant' of standard error (= posterior SD)
# - two Bayesian credible intervals (CRI and HPDI)
sd(plogis(out$ltheta))                              # Posterior SD
# [1] 0.06950494                                    # Your result will differ

quantile(plogis(out$ltheta), prob=c(0.025, 0.975))  # Symmetrical Bayesian CRI (your result will differ)
#      2.5%     97.5%
# 0.2705029 0.5383629

HPDinterval(mcmc(plogis(out$ltheta)))               # Highest posterior density credible
                                                    # interval(HPDI); your result will differ
#          lower     upper
# var1 0.2672082 0.5349586
# attr(,"Probability")
# [1] 0.95


# Compute p(theta > 0.5)
mean(plogis(out$ltheta) > 0.5)
# [1] 0.07584
