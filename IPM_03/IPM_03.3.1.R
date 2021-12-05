# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------

library(IPMbook) ; library(jagsUI)

# 3.3 Classical analysis of a matrix population model
# ===================================================

# 3.3.1 Analysis of a matrix population model without stochasticity
#       and parameter uncertainty
# -----------------------------------------------------------------

# Define the demographic rates
sj <- 0.3                                           # Juvenile survival
sa <- 0.55                                          # Adult survival
f1 <- 1.3                                           # Number of female fledglings per 1-year old female
fa <- 1.8                                           # Number of female fledglings per adult female

# Define the transition matrix A (inspired by woodchat shrike demography)
A <- matrix(c(f1 * sj, fa * sj, sa, sa), ncol=2, byrow=TRUE)
dimnames(A) <- list(c("1y", "ad"), c("1y", "ad"))
A                                                   # print the matrix
#      1y   ad
# 1y 0.39 0.54
# ad 0.55 0.55

# Asymptotic population growth rate:
# ''''''''''''''''''''''''''''''''''

N1 <- c(10, 1)

T <- 7
N <- matrix(NA, nrow=2, ncol=T)                     # Population sizes
gr <- matrix(NA, nrow=3, ncol=T-1)                  # Growth rates
N[,1] <- N1
for (t in 2:T){
  N[,t] <- A %*% N[,t-1]
  gr[1,t-1] <- N[1,t] / N[1,t-1]
  gr[2,t-1] <- N[2,t] / N[2,t-1]
  gr[3,t-1] <- (N[1,t] + N[2,t]) / (N[1,t-1] + N[2,t-1])
}
sr <- N[1, ] / colSums(N)                           # State distribution (proportion of 1y)

# ~~~~ extra code for Figure 3.8 ~~~~
op <- par(las=1, "mfrow")
layout(matrix(1:3, 1, 3,byrow=TRUE), widths=c(1, 1, 1), heights=1, TRUE)
# Plot stage-specific population size
plot(N[1,], type="b", pch=16, ylim=range(c(min(N), colSums(N))), axes=FALSE,
    ylab="Population size", xlab="Time", col="red")
points(N[2,], type="b", pch=16, col="orange")
points(colSums(N), type="b", pch=16, col="blue")
axis(1)
axis(2)
mtext("A", at=1, line=0.5, cex=1.5)
legend("bottomright", legend=c("1y", "ad", "total"), bty="n", pch=rep(16, 3),
    lty=rep(1, 3), col=c("red", "orange", "blue"))

# Plot stage-specific growth rates
plot(x=(1:(T-1)) + 0.5, y=gr[1,], type="b", pch=16, ylim=range(gr),
    xlim=c(1, T), axes=FALSE, ylab="Population growth rate", xlab="Time", col="red")
points(x=(1:(T-1)) + 0.5, y=gr[2,], type="b", pch=16, col="orange")
points(x=(1:(T-1)) + 0.5, y=gr[3,], type="b", pch=16, col="blue")
axis(1, at=1:T, labels=1:T)
axis(2)
abline(h=max(Re(eigen(A)$values)), lty=2)
mtext("B", at=1, line=0.5, cex=1.5)

# Plot actual stage distribution
plot(sr, type="b", pch=16, ylab="Proportion of stage classes", axes=FALSE,
    col="red", xlab="Time", ylim=range(c(sr, 1-sr)))
points(1-sr, type="b", pch=16, col="orange")
axis(1)
axis(2)
u <- which.max(Re(eigen(A)$values))
revec <- Re(eigen(A)$vectors[,u])
abline(h=revec[1]/sum(revec), col="red", lty=2)
abline(h=revec[2]/sum(revec), col="orange", lty=2)
mtext("C", at=1, line=0.5, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Compute asymptotic population growth rate (lambda)
lambda <- max(Re(eigen(A)$values))
# [1] 1.020818

# ~~~~ extra code for Figure 3.9 ~~~~
# Define the demographic rates
s1 <- 0.5     # Juvenile survival
s2 <- 0.8     # 1y survival
s3 <- 0.9     # Ad survival
rho1 <- 0.5   # Productivity of 2-years old females
rho2 <- 0.9   # Productivity of adult females

# Define the transition matrix A
Abat <- matrix(c(0, rho1*s1/2, rho2*s1/2,
              s2, 0, 0,
              0, s3, s3), ncol=3, byrow=TRUE)

# Define R structures and do computations
T <- 9
N1 <- c(10, 5, 5)                   # Initial population sizes
N <- matrix(NA, nrow=3, ncol=T)     # Population size
gr <- matrix(NA, nrow=4, ncol=T-1)  # Growth rate
N[,1] <- N1
for (t in 2:T){
  N[,t] <- Abat %*% N[,t-1]  # stage-specific population sizes
  gr[1,t-1] <- N[1,t] / N[1,t-1]
  gr[2,t-1] <- N[2,t] / N[2,t-1]
  gr[3,t-1] <- N[3,t] / N[3,t-1]
  gr[4,t-1] <- (N[1,t] + N[2,t] + N[3,t]) / (N[1,t-1] + N[2,t-1] + N[3,t-1])
}
sr <- sweep(N, 2, colSums(N), "/")   # Stage distribution

# Plots
co <- c("red", "orange", "tomato2")
op <- par(las=1, cex=1.1, "mfrow")
layout(matrix(1:3, 1, 3, byrow=TRUE), widths=c(1, 1, 1), heights=1, TRUE)
plot(N[1,], type="b", pch=16, ylim=range(c(min(N), colSums(N))), axes=FALSE,
    ylab="Population size", xlab="Time", col=co[1])
points(N[2,], type="b", pch=16, col=co[2])
points(N[3,], type="b", pch=16, col=co[3])
points(colSums(N), type="b", pch=16, col="blue")
axis(1, at=1:T)
axis(2)
mtext("A", at=1, line=0.5, cex=1.5)

plot(x=(1:(T-1)) + 0.5, y=gr[1,], type="b", pch=16, ylim=range(gr), xlim=c(1, T),
    axes=FALSE, ylab="Population growth rate", xlab="Time", col=co[1])
points(x=(1:(T-1)) + 0.5, y=gr[2,], type="b", pch=16, col=co[2])
points(x=(1:(T-1)) + 0.5, y=gr[3,], type="b", pch=16, col=co[3])
points(x=(1:(T-1)) + 0.5, y=gr[4,], type="b", pch=16, col="blue")
axis(1, at=1:T)
axis(2)
abline(h=max(Re(eigen(Abat)$values)), lty=2)
mtext("B", at=1, line=0.5, cex=1.5)
legend("topright", legend=c("1y", "2y", "ad", "total"), bty="n", pch=rep(16,4),
    lty=rep(1, 4), col=c(co, "blue"))

plot(sr[1,], type="b", pch=16, ylab="Proportion of stage classes", axes=FALSE,
    col=co[1], ylim=c(0, 0.8), xlab="Time")
points(sr[2,], type="b", pch=16, col=co[2])
points(sr[3,], type="b", pch=16, col=co[3])
axis(1, at=1:T)
axis(2)
u <- which.max(Re(eigen(Abat)$values))
revec <- Re(eigen(Abat)$vectors[,u])
abline(h=revec[1]/sum(revec), col=co[1], lty=2)
abline(h=revec[2]/sum(revec), col=co[2], lty=2)
abline(h=revec[3]/sum(revec), col=co[3], lty=2)
mtext("C", at=1, line=0.5, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Stable stage distribution:
# ''''''''''''''''''''''''''

u <- which.max(Re(eigen(A)$values))
# [1] 1
revec <- Re(eigen(A)$vectors[,u])
# [1] -0.6503047 -0.7596734

revec/sum(revec)
# [1] 0.4612162 0.5387838

# Stage-specific reproductive values:
# '''''''''''''''''''''''''''''''''''

u <- which.max(Re(eigen(A)$values))
levec <- Re((eigen(A)$vectors)[u,])
levec / sum(levec)
# [1] 0.4631655 0.5368345

# Net reproductive rate:
# ''''''''''''''''''''''

i <- 1:100                                          # 100 as our approximation to infinity
R0 <- sj * f1 + sj * fa * sum(sa^i)
# [1] 1.05

# Generation time
# '''''''''''''''

Q <- 100                                            # Our approximation to infinity
i <- 2:Q
G <- sj * f1 / lambda + sj * fa * sum(i * sa^(i-1) * lambda^(-i))
# [1] 2.339834

GT <- log(R0) / log(lambda)
# [1] 2.368012


# Sensitivity and elasticity (perturbation analysis):
# '''''''''''''''''''''''''''''''''''''''''''''''''''

senmat <- levec %*% t(revec)
#           [,1]      [,2]
# [1,] 0.4228963 0.4940192
# [2,] 0.4901602 0.5725957

elasmat <- senmat * A / lambda
#           1y        ad
# 1y 0.1615661 0.2613301
# ad 0.2640904 0.3085053

derivmat <- matrix(c(f1, fa, 0, 0), ncol=2, byrow=TRUE)
ssj <- sum(senmat * derivmat)                       # Sensitivity
esj <- sj / lambda * ssj                            # Elasticity
ssj; esj
# [1] 1.439
# [1] 0.4228963

# ~~~~ extra code for Figure 3.10 ~~~~
# Define the demographic rates
sj <- 0.3      # Juvenile survival
sa <- 0.55     # Adult survival
f1 <- 1.3      # Number of female fledglings per 1-year old female and year
fa <- 1.8      # Number of female fledglings per adult female and year

# Sensitivity
ds <- seq(-0.1, 0.1, length.out=501)
lams <- matrix(NA, ncol=4, nrow=501)
for (i in 1:501){
  A <- matrix(c((f1+ds[i])*sj, fa*sj, sa, sa), ncol=2, byrow=TRUE)
  lams[i,1] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj, (fa+ds[i])*sj, sa, sa), ncol=2, byrow=TRUE)
  lams[i,2] <- max(eigen(A)$values)
  A <- matrix(c(f1*(sj+ds[i]), fa*(sj+ds[i]), sa, sa), ncol=2, byrow=TRUE)
  lams[i,3] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj, fa*sj, sa+ds[i], sa+ds[i]), ncol=2, byrow=TRUE)
  lams[i,4] <- max(eigen(A)$values)
}

# Elasticity
es <- seq(0.9, 1.1, length.out=501)
lame <- matrix(NA, ncol=4, nrow=501)
for (i in 1:501){
  A <- matrix(c(f1*es[i]*sj, fa*sj, sa, sa), ncol=2, byrow=TRUE)
  lame[i,1] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj, fa*es[i]*sj, sa, sa), ncol=2, byrow=TRUE)
  lame[i,2] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj*es[i], fa*sj*es[i], sa, sa), ncol=2, byrow=TRUE)
  lame[i,3] <- max(eigen(A)$values)
  A <- matrix(c(f1*sj, fa*sj, sa*es[i], sa*es[i]), ncol=2, byrow=TRUE)
  lame[i,4] <- max(eigen(A)$values)
}

# Plots
op <- par("mfrow", "las")
layout(matrix(1:2, 1, 2, byrow=TRUE), widths=c(1, 1), heights=1, TRUE)
co <- c("red", "orange", "darkgreen", "blue")
lw <- 2
plot(y=lams[,1], x=ds, type="l", ylim=range(cbind(lams, lame)), axes=FALSE,
    ylab="Population growth rate", xlab="Absolute change in demographic rate",
    col=co[1], lwd=lw)
axis(1)
axis(1, at=c(-0.075, -0.025, 0.025, 0.075), labels=NA, tcl=-0.25)
axis(2, las=1)
lines(y=lams[,2], x=ds, col=co[2], lwd=lw)
lines(y=lams[,3], x=ds, col=co[3], lwd=lw)
lines(y=lams[,4], x=ds, col=co[4], lwd=lw)
mtext("A", at=-0.1, line=0.5, cex=1.5)
legend("topleft", lwd=rep(lw, 4), col=co, legend=c(expression(italic(f)[1]),
    expression(italic(f)[italic(a)]), expression(italic(s)[italic(j)]),
    expression(italic(s)[italic(a)])), bty="n")

plot(y=lame[,1], x=es, type="l", ylim=range(cbind(lams, lame)), axes=FALSE,
    ylab="", xlab="Relative change in demographic rate (%)", col=co[1], lwd=lw)
axis(1, at=c(0.925, 0.975, 1.025, 1.075), labels=NA, tcl=-0.25)
axis(1, at=c(0.9, 0.95, 1, 1.05, 1.1), labels=c(-10, -5, 0, 5, 10))
axis(2, las=1)
lines(y=lame[,2], x=es, col=co[2], lwd=lw)
lines(y=lame[,3], x=es, col=co[3], lwd=lw)
lines(y=lame[,4], x=es, col=co[4], lwd=lw)
mtext("B", at=0.9, line=0.5, cex=1.5)
par(op)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
