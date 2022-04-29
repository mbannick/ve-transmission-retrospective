rm(list=ls())

library(RColorBrewer)
library(magrittr)
library(data.table)
library(ggplot2)

# Get inferred vaccine efficacy
# In this scenario, we fix mu, the true vaccine efficacy
inferred.ve <- function(rho, mu, delta, lambda){
  denominator <- (lambda*rho + delta * (1-lambda*rho))
  numerator <- mu * (rho + (1-rho)*delta)
  return(1-numerator/denominator)
}

# Set parameter grid
probs <- seq(0.0, 1.0, by=0.001)
prob.grid <- expand.grid(t0=0.5,
                         s0=0.5,
                         delta=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
                         mu=probs,
                         lambda=0.2) %>% data.table

results <- mapply(FUN=inferred.ve,
                  rho=prob.grid$s0,
                  delta=prob.grid$delta,
                  mu=prob.grid$mu,
                  lambda=prob.grid$lambda)

data <- cbind(prob.grid, results)
data[, ve := 1-mu]
data[, ve_inf := results]
data[, results := NULL]

# Restrict to scenarios where
# inferred vaccine efficacy is positive
data[, valid := (mu <= (lambda * s0 + delta*(1-lambda*s0))/(s0 + delta*(1-s0)))]

# FIGURE 2(A) IN THE MAIN TEXT
pdf("bias-plot-symptoms.pdf", height=6, width=6)
p2 <- ggplot(data=data[valid==T,],
             aes(x=ve, y=ve_inf, color=1-delta, group=delta)) +
  geom_line() +
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  labs(x="VE",
       y="VE Inferred",
       color="Reduction in transmission given asymptomatic") +
  ggtitle("(a) Symptomatic Sampling") +
  theme(legend.position="bottom")
p2
dev.off()

# IN ORDER TO CREATE THE FULL FIGURE 2 IN THE MAIN TEXT
# RUN THIS, WITH P1 FROM THE INFREQUENT-SAMPLING-FINAL.R SCRIPT

# library(gridExtra)
# pdf("fig-1.pdf", height=6, width=10)
# grid.arrange(p2, p1, nrow=1)
# dev.off()

