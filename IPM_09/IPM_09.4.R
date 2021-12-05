# Schaub & Kéry (2022) Integrated Population Models
# Chapter 9 : Retrospective population analyses
# ---------------------------------------------

library(IPMbook) ; library(jagsUI)

# ~~~ this needs 'out1' from section 9.2 ~~~
load("ResultsHoopoe.Rdata")
# ~~~ and this from section 9.3 ~~~
n.draws <- out1$mcmc.info$n.samples
draws <- out1$sims.list
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 9.4 Transient life-table response experiments
# =============================================

lambda <- expression(((n1 + n3) * 0.5 * phij * f1 + n2 * 0.5 * phij * f2 + (phia + om) *
    (n1 + n2 + n3)) / (n1 + n2 + n3))             # Define lambda

D(lambda, "phij") # Get derivative for juvenile apparent survival
# ((n1 + n3) * 0.5 * f1 + n2 * 0.5 * f2)/(n1 + n2 + n3)

# Calculate proportional population sizes
n.years <- length(out1$mean$phij)
n1 <- draws$N[, 1, 1:n.years] / draws$Ntot[, 1:n.years]
n2 <- draws$N[, 2, 1:n.years] / draws$Ntot[, 1:n.years]
n3 <- draws$N[, 3, 1:n.years] / draws$Ntot[, 1:n.years]

# Extract the mean demographic rates and population sizes and store them in a list
mu <- list(phij=draws$mean.phij, phia=draws$mean.phia, om=draws$mean.omega, f1=draws$mean.f1,
    f2=draws$mean.f2, n1=apply(n1, 1, mean), n2=apply(n2, 1, mean), n3=apply(n3, 1, mean))

# Calculate sensitivities
sens <- matrix(NA, n.draws, 8)
colnames(sens) <- c("phij", "phia", "om", "f1", "f2", "n1", "n2", "n3")
sens[,"phij"] <- eval(D(lambda, "phij"), envir=mu)
sens[,"phia"] <- eval(D(lambda, "phia"), envir=mu)
sens[,"om"] <- eval(D(lambda, "om"), envir=mu)
sens[,"f1"] <- eval(D(lambda, "f1"), envir=mu)
sens[,"f2"] <- eval(D(lambda, "f2"), envir=mu)
sens[,"n1"] <- eval(D(lambda, "n1"), envir=mu)
sens[,"n2"] <- eval(D(lambda, "n2"), envir=mu)
sens[,"n3"] <- eval(D(lambda, "n3"), envir=mu)

# Define matrix to store results
cont <- matrix(NA, nrow=n.draws, ncol=8)
colnames(cont) <- c("phij", "phia", "om", "f1", "f2", "n1", "n2", "n3")

# Calculate contributions for each demographic rate and stage-structured population size at each MCMC draw
for (s in 1:n.draws){
  dp_stoch <- cbind(draws$phij[s,], draws$phia[s,], draws$omega[s,], draws$f1[s,1:n.years],
      draws$f2[s,1:n.years], n1[s,], n2[s,], n3[s,])
  # Derive process variance and covariance among demographic parameters using shrunk estimates of
  #     demographic rates and proportional pop. sizes
  dp_varcov <- var(dp_stoch)
  sensvec <- sens[s, ]
  # Calculate demographic contributions
  contmatrix <- dp_varcov * outer(sensvec, sensvec)
  cont[s, ] <- rowSums(contmatrix)
}

# ~~~~ code for figure 9.6 ~~~~
CRI <- apply(cont, 2, quantile, probs=c(0.025, 0.975))
names.arg <- c(expression(atop("Juvenile", "survival ("*phi[italic(j)]*")")),
    expression(atop("Adult", "survival ("*phi[italic(a)]*")")),
    expression(atop("Immigration", "("*omega*")")),
    expression(atop("Productivity", "1y ("*italic(f)[1]*")")),
    expression(atop("Productivity", "ad ("*italic(f)[2]*")")),
    expression(italic(n)[1]), expression(italic(n)[2]), expression(italic(n)[3]))
a <- barplot(colMeans(cont), names.arg=names.arg, ylab="Contribution", las=1,
    ylim=range(CRI), col=c(rep("grey65", 5), rep("black", 3)),
    border=c(rep("grey65", 5), rep("black", 3)))
segments(a, CRI[1,], a, CRI[2,], lwd=1.5)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Compute differences and means of demographic rates and population structure between successive years
diff <- array(NA, dim=c(n.draws, n.years-1, 8))
dimnames(diff)[[3]] <- c("phij", "phia", "om", "f1", "f2", "n1", "n2", "n3")

# Function to compute differences over successive time steps
getDiff <- function(x) x[,2:n.years] - x[,1:(n.years-1)]

diff[,,"phij"] <- getDiff(draws$phij)
diff[,,"phia"] <- getDiff(draws$phia)
diff[,,"om"] <- getDiff(draws$om)
diff[,,"f1"] <- getDiff(draws$f1)
diff[,,"f2"] <- getDiff(draws$f2)
diff[,,"n1"] <- getDiff(draws$N[,1,])
diff[,,"n2"] <- getDiff(draws$N[,2,])
diff[,,"n3"] <- getDiff(draws$N[,3,])

# Function to compute means over successive time steps, store them in a list
getMn <- function(x) (x[,2:n.years] + x[,1:(n.years-1)]) / 2

means <- list(phij=getMn(draws$phij), phia=getMn(draws$phia), om=getMn(draws$om),
    f1=getMn(draws$f1), f2=getMn(draws$f2), n1=getMn(draws$N[,1,]), n2=getMn(draws$N[,2,]),
    n3=getMn(draws$N[,3,]))

# Compute sensitivities
senss <- array(NA, dim=c(n.draws, n.years-1, 8))
dimnames(senss)[[3]] <- c("phij", "phia", "om", "f1", "f2", "n1", "n2", "n3")
senss[,,"phij"] <- eval(D(lambda, "phij"), envir=means)
senss[,,"phia"] <- eval(D(lambda, "phia"), envir=means)
senss[,,"om"] <- eval(D(lambda, "om"), envir=means)
senss[,,"f1"] <- eval(D(lambda, "f1"), envir=means)
senss[,,"f2"] <- eval(D(lambda, "f2"), envir=means)
senss[,,"n1"] <- eval(D(lambda, "n1"), envir=means)
senss[,,"n2"] <- eval(D(lambda, "n2"), envir=means)
senss[,,"n3"] <- eval(D(lambda, "n3"), envir=means)

conts <- diff * senss


# ~~~~ figure 9.7 ~~~~
# Create a plot for the sequential differences in finite growth rates
n.draws <- out1$mcmc.info$n.samples
n.years <- length(out1$mean$phij)
diff.lam <- draws$lambda[, 2:n.years] - draws$lambda[, 1:(n.years-1)]
differences <- cbind(2002:2015, apply(diff.lam, 2, mean))

# Mean contributions
V <- apply(conts, 3:2, mean)
V1 <- V2 <- V
V1[V1<0] <- 0
V2[V2>0] <- 0

# Make figure
colours <- c(colorRampPalette(c("blue4", "lightblue"))(2),
    colorRampPalette(c("darkred", "yellow"))(3),
    colorRampPalette(c("lightgreen", "darkgreen"))(3))
op <- par(mfrow=c(2,1), mar=c(4, 5, 1, 1))
barplot(differences[,2], ylim=c(-0.35, 0.35),
    ylab="Difference in population growth", axes=FALSE, border=NA)
axis(2, las=1)
abline(h=0)
legend.text <- c(expression(italic(phi)[italic(j)]),
                expression(italic(phi)[italic(a)]),
                expression(omega), expression(italic(f)[1]),
                expression(italic(f)[2]), expression(italic(n)[3]),
                expression(italic(n)[2]), expression(italic(n)[1]))
barplot(V1, ylim=c(-0.35, 0.35), col=colours, ylab="Contribution",
    xlab="Beginning year", axes=FALSE, names.arg=2002:2015, legend.text=legend.text,
        args.legend=list(bty="n", x="top", ncol=3, border=NA), border=NA)
barplot(V2, add=TRUE, col=colours, axes=FALSE, border=NA)
axis(2, las=1)
abline(h=0)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
