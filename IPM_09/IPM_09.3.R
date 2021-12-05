# Schaub & Kéry (2022) Integrated Population Models
# Chapter 9 : Retrospective population analyses
# ---------------------------------------------

library(IPMbook) ; library(jagsUI)

# ~~~ this needs 'out1' from section 9.2 ~~~
load("ResultsHoopoe.Rdata")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 9.3 Life-table response experiments
# ===================================

# Calculation of growth rate sensitivities
delta <- 0.001                                    # (Small) size of perturbation
T <- 5                                            # Number of years with projections
n.draws <- out1$mcmc.info$n.samples               # Number of MCMC draws

# Define arrays to store output
N.ref <- N.star.phij <- N.star.phia <- N.star.om <- N.star.f1 <- N.star.f2 <- array(NA,
    dim=c(n.draws, 3, T))
N.ref[,,1] <- N.star.phij[,,1] <- N.star.phia[,,1] <- N.star.om[,,1] <- N.star.f1[,,1] <-
    N.star.f2[,,1] <- 1
r.annual.ref <- r.annual.star.phij <- r.annual.star.phia <- r.annual.star.om <- r.annual.star.f1 <-
    r.annual.star.f2 <- matrix(NA, nrow=n.draws, ncol=T)
lambda <- numeric()
sens <- matrix(NA, nrow=n.draws, ncol=5)
draws <- out1$sims.list                           # store MCMC draws in a new object to simplify calculation

# Loop over all MCMC draws and project in time
for (s in 1:n.draws){                             # Loop over all MCMC draws
  for (t in 1:(T-1)){                             # Loop over all time steps
    # Calculate asymptotic population growth rate based on stage-specific population sizes
    N.ref[s,1,t+1] <- draws$mean.phij[s] * (draws$mean.f1[s] *
        (N.ref[s,1,t] + N.ref[s,3,t]) + draws$mean.f2[s] * N.ref[s,2,t]) / 2
    N.ref[s,2,t+1] <- draws$mean.phia[s] * sum(N.ref[s,,t])
    N.ref[s,3,t+1] <- draws$mean.omega[s] * sum(N.ref[s,,t])
    r.annual.ref[s,t] <- log(sum(N.ref[s,,t+1])) - log(sum(N.ref[s,,t]))

    # Sensitivity with respect to juvenile survival
    N.star.phij[s,1,t+1] <- (draws$mean.phij[s] + delta) / 2 * (draws$mean.f1[s] *
        (N.star.phij[s,1,t] + N.star.phij[s,3,t]) + draws$mean.f2[s] * N.star.phij[s,2,t])
    N.star.phij[s,2,t+1] <- draws$mean.phia[s] * sum(N.star.phij[s,,t])
    N.star.phij[s,3,t+1] <- draws$mean.omega[s] * sum(N.star.phij[s,,t])
    r.annual.star.phij[s,t] <- log(sum(N.star.phij[s,,t+1])) - log(sum(N.star.phij[s,,t]))

    # Sensitivity with respect to adult survival
    N.star.phia[s,1,t+1] <- draws$mean.phij[s] / 2 * (draws$mean.f1[s] *
        (N.star.phia[s,1,t] + N.star.phia[s,3,t]) + draws$mean.f2[s] * N.star.phia[s,2,t])
    N.star.phia[s,2,t+1] <- (draws$mean.phia[s] + delta) * sum(N.star.phia[s,,t])
    N.star.phia[s,3,t+1] <- draws$mean.omega[s] * sum(N.star.phia[s,,t])
    r.annual.star.phia[s,t] <- log(sum(N.star.phia[s,,t+1])) - log(sum(N.star.phia[s,,t]))

    # Sensitivity with respect to immigration
    N.star.om[s,1,t+1] <- draws$mean.phij[s] / 2 * (draws$mean.f1[s] * (N.star.om[s,1,t] +
        N.star.om[s,3,t]) + draws$mean.f2[s] * N.star.om[s,2,t])
    N.star.om[s,2,t+1] <- draws$mean.phia[s] * sum(N.star.om[s,,t])
    N.star.om[s,3,t+1] <- (draws$mean.omega[s] + delta) * sum(N.star.om[s,,t])
    r.annual.star.om[s,t] <- log(sum(N.star.om[s,,t+1])) - log(sum(N.star.om[s,,t]))

    # Sensitivity with respect to productivity 1y
    N.star.f1[s,1,t+1] <- draws$mean.phij[s] / 2 * ((draws$mean.f1[s] + delta) *
        (N.star.f1[s,1,t] + N.star.f1[s,3,t]) + draws$mean.f2[s] * N.star.f1[s,2,t])
    N.star.f1[s,2,t+1] <- draws$mean.phia[s] * sum(N.star.f1[s,,t])
    N.star.f1[s,3,t+1] <- draws$mean.omega[s] * sum(N.star.f1[s,,t])
    r.annual.star.f1[s,t] <- log(sum(N.star.f1[s,,t+1])) - log(sum(N.star.f1[s,,t]))
    # Sensitivity with respect to productivity ad
    N.star.f2[s,1,t+1] <- draws$mean.phij[s] / 2 * (draws$mean.f1[s] * (N.star.f2[s,1,t] +
        N.star.f2[s,3,t]) + (draws$mean.f2[s] + delta) * N.star.f2[s,2,t])
    N.star.f2[s,2,t+1] <- draws$mean.phia[s] * sum(N.star.f2[s,,t])
    N.star.f2[s,3,t+1] <- draws$mean.omega[s] * sum(N.star.f2[s,,t])
    r.annual.star.f2[s,t] <- log(sum(N.star.f2[s,,t+1])) - log(sum(N.star.f2[s,,t]))
  } #t

  lambda[s] <- exp(r.annual.ref[s,T-1])

  # Growth rate sensitivities
  sens[s,1] <- (exp(r.annual.star.phij[s,T-1]) - lambda[s]) / delta
  sens[s,2] <- (exp(r.annual.star.phia[s,T-1]) - lambda[s]) / delta
  sens[s,3] <- (exp(r.annual.star.om[s,T-1]) - lambda[s]) / delta
  sens[s,4] <- (exp(r.annual.star.f1[s,T-1]) - lambda[s]) / delta
  sens[s,5] <- (exp(r.annual.star.f2[s,T-1]) - lambda[s]) / delta
} #s

# Calculate variance-covariance matrix of the demographic rates
vcov <- array(NA, dim=c(n.draws, 5, 5))
for (s in 1:n.draws){
  vcov[s,,] <- var(cbind(draws$phij[s,], draws$phia[s,], draws$omega[s,], draws$f1[s,1:15],
      draws$f2[s,1:15]))
}

# Calculate contribution matrix for every posterior draw
cont <- array(NA, dim=c(n.draws, 5, 5))
for (s in 1:n.draws){
  cont[s,,] <- vcov[s,,] * outer(sens[s,], sens[s,])
}

# ~~~~ Make figure 9.3 ~~~~
CRI <- apply(apply(cont, c(1,2), sum), 2, quantile, probs=c(0.025, 0.975))
op <- par(cex=1.5)
names.arg=c(expression(atop("Juvenile", "survival ("*phi[italic(j)]*")")),
            expression(atop("Adult", "survival ("*phi[italic(a)]*")")),
            expression(atop("Immigration", "("*omega*")")),
            expression(atop("Productivity", "1y ("*italic(f)[1]*")")),
            expression(atop("Productivity", "ad ("*italic(f)[2]*")")))
a <- barplot(apply(apply(cont, c(1,2), sum), 2, mean), ylim=range(CRI),
    names.arg=names.arg, ylab="Contribution", xlab="", las=1, col="grey65", border="grey65")
segments(a, CRI[1,], a, CRI[2,], lwd=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~ code for figure 9.4 ~~~~
mean.cont <- apply(cont, c(2, 3), mean)
mean.cont2 <- rbind(diag(mean.cont), rowSums(mean.cont) - diag(mean.cont))
library(RColorBrewer)
co <- brewer.pal(n=8, name='Blues')[c(7,3)]

op <- par(mar=c(4, 5, 1, 1))
names.arg=c(expression(atop("Juvenile", "survival ("*phi[italic(j)]*")")),
            expression(atop("Adult", "survival ("*phi[italic(a)]*")")),
            expression(atop("Immigration", "("*omega*")")),
            expression(atop("Productivity", "1y ("*italic(f)[1]*")")),
            expression(atop("Productivity", "ad ("*italic(f)[2]*")")))
a <- barplot(mean.cont2, names.arg=names.arg, ylab=NA, xlab=NA, col=co,
    border=co, axes=FALSE)
axis(2, las=1)
mtext("Contribution", side=2, line=4)
legend("topright", legend=c("Covariance", "Variance"), col=co[2:1],
    pch=c(15, 15), bty="n", pt.cex=2)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
