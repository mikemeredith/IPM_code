# Schaub & Kéry (2022) Integrated Population Models
# Chapter 3 : Introduction to stage-structured population models
# ----------------------------------------------------
# Code from final MS.

# 3.2 Age- and stage-structured population models
# ===============================================

# 3.2.1 From a life-cycle graph to population equations
# -----------------------------------------------------

# ~~~~ Code for Fig. 3.3 ~~~~
sj <- 0.3      # Juvenile survival
sa <- 0.55     # Adult survival
f1 <- 1.3      # Number of female fledglings per 1-year old female and year
fa <- 1.8      # Number of female fledglings per adult female and year

N <- numeric(14)
N[1] <- 10
for (i in 2:14){
  if (i %% 2 == 0){   # even numbers
    N[i] <- N[i-1] / 2 * f1 + N[i-1] / 2 * fa
  } else {               # odd numbers
    N[i] <- N[i-2] * sa + N[i-1] * sj
  }
}
x <- c(rbind(1:7, 1:7 + 0.1))

plot(x=x, y=N, type="l", ylab="Total population size", ylim=c(7, 18),
    col="blue", lwd=2, axes=FALSE, xlab="Year")
axis(2, las=1)
axis(1)
points(x, N, pch=c(22, 1), cex=2, col=c("red", "darkgreen"), lwd=2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 3.2.2 Age-structured, pre-birth-pulse model (no code)

# 3.2.3 Stage-structured, pre-birth-pulse model (no code)

# 3.2.4 Age-structured, post-birth-pulse model (no code)

# 3.2.5 Stage-structured, post-birth-pulse model (no code)
