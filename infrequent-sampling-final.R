rm(list=ls())

library(magrittr)
library(data.table)
library(ggplot2)

# Get analytic vaccine efficacy -- truth
analytic.ve <- function(nu, lambda){
  return(1-nu*lambda)
}

# Get inferred SAR based on parameters
get.mu.term <- function(tau, rho, c, k){
  if(k <= rho - c){
    term <- rho
  } else if(k <= rho + c) {
    Sk.num <- 2*k*(rho + c) - (rho - c)^2 - k^2
    Sk.den <- 4 * c * k
    Sk <- Sk.num / Sk.den

    term.num <- 3*k*(rho + c)^2 - k^3 - 2*(rho - c)^3
    term.den <- 12 * c * k * Sk
    term <- term.num / term.den
  } else {
    term <- rho + 1/(3*rho) * c^2
  }
  return(tau * term)
}

# Get inferred VE based on the inferred SAR
inferred.ve <- function(nu, lambda, tau0, rho0, c, k){

  rho1 <- rho0 * lambda
  tau1 <- tau0 * nu

  mu1 <- get.mu.term(tau=tau1, rho=rho1, c=c, k=k)
  mu0 <- get.mu.term(tau=tau0, rho=rho0, c=c, k=k)

  mu <- mu1/mu0
  return(1-mu)
}

# Set parameters
ks <- seq(1, 25, by=1)
lambdas <- 8/14
nus <- seq(0, 1, by=0.1)
RHO0 <- 14

tau0s <- c(0.01)
C <- 7

# Generate parameter grid
params <- expand.grid(lambda=lambdas, nu=nus, tau0=tau0s,
                      rho0=RHO0, c=C, k=ks)
setnames(params, c("lambda", "nu", "tau0", "rho0", "c", "k"))

# Calculate analytic and inferred VE for each parameter set
analytic <- mapply(FUN=analytic.ve, lambda=params$lambda,
                   nu=params$nu)
inferred <- mapply(FUN=inferred.ve, lambda=params$lambda,
                   nu=params$nu, tau0=params$tau0,
                   rho0=RHO0, c=C, k=params$k)

data <- data.table(params)
data[, ve := analytic]
data[, ve_inf := inferred]
data <- data[, .(lambda, nu, tau0, rho0, c, k, analytic, inferred)]
setnames(data, "analytic", "ve")

# Reduction scale for lamda and nu
data[, red.lambda := 1-lambda]
data[, red.nu := 1-nu]

df <- copy(data)
df <- data.frame(df)

# FIGURE 1 IN THE APPENDIX
pdf("fig-infrequent-monotone.pdf", height=8, width=8)
ggplot(data=data,
       aes(x=k, y=inferred, color=ve, group=ve)) +
  geom_line() +
  labs(color="VE", x="Testing Interval", y="VE Inferred") +
  theme(legend.position="bottom")
dev.off()

# Which k's (testing interval) to plot
k.breaks <- c(1, 3, 5, 7, 9, 11, 14)

# FIGURE 2(B) IN THE MAIN TEXT
pdf("bias-plot-duration.pdf", height=6, width=6)
p1 <- ggplot(data=data[k %in% k.breaks],
             aes(x=ve, y=inferred, color=k, group=k)) +
  geom_line() +
  labs(color="Testing Interval (days)", y="VE Inferred", x="VE") +
  theme(legend.position="bottom") +
  ggtitle("(b) Infrequent Sampling") +
  geom_abline(intercept=0, slope=1, linetype='dashed') +
  scale_color_gradient(
    breaks = k.breaks,
    labels = k.breaks,
  )
p1
dev.off()
