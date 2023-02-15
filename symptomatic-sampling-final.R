# rm(list=ls())

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
pdf("figure-1-A.pdf", height=6, width=8, useDingbats=FALSE)

data <- data.table(data)
d.breaks <- rev(data$delta %>% unique())

cols <- sequential_hcl(n=length(d.breaks), palette="viridis")
par(oma=c(0, 5, 0, 0), xpd=FALSE)
for(i in 1:length(d.breaks)){
  dat.k <- data[(delta == d.breaks[i]) & (valid == T)]
  if(i == 1){
    plot(dat.k$ve_inf ~ dat.k$ve, type='l', col=cols[i], lwd=2,
         ylim=c(0, 1),
         ylab="Vaccine Efficacy Inferred",
         xlab="Vaccine Efficacy",
         cex.lab=1.3, las=1, axes=F)
    box(bty="l")
    axis(1, cex.axis=1.3, las=1)
    axis(2, cex.axis=1.3, las=1)
    par(xpd=FALSE)
    abline(a=0, b=1, lty=2)
    mtext(text="A)", side=3, line=-2, outer=T, adj=-0.1, cex=1.3)
  } else {
    lines(dat.k$ve_inf ~ dat.k$ve, col=cols[i], lwd=2)
  }
}
pos <- legend("topleft", inset=0.08, cex=1.3, title.adj=0.5, col=cols, lwd=2, ncol=2, legend=c(
  paste0(sprintf("%2.0f", (1-d.breaks)*100))), bty="n",
  title=expression(underline("Reduction in Transmission\n     Given Asymptomatic, %"))
)
xleft <- pos$rect[["left"]]
ytop <- pos$rect[["top"]]
ybottom <- ytop - pos$rect[["h"]]
xright <- xleft + pos$rect[["w"]]
rect(xleft, ybottom+0.025, xright, ytop+0.08)

# p2 <- ggplot(data=data[valid==T,],
#              aes(x=ve, y=ve_inf, color=1-delta, group=delta)) +
#   geom_line() +
#   geom_abline(slope=1, intercept=0, linetype='dashed') +
#   labs(x="VE",
#        y="VE Inferred",
#        color="Reduction in transmission given asymptomatic") +
#   ggtitle("(a) Symptomatic Sampling") +
#   theme(legend.position="bottom")
# p2
dev.off()

# IN ORDER TO CREATE THE FULL FIGURE 2 IN THE MAIN TEXT
# RUN THIS, WITH data.inf FROM THE INFREQUENT-SAMPLING-FINAL.R SCRIPT

# library(gridExtra)
# pdf("fig-1.pdf", height=6, width=10)
# grid.arrange(p2, p1, nrow=1)
# dev.off()
#

pdf("figure-1.pdf", height=12, width=8, useDingbats=FALSE)
par(mfrow=c(2, 1), oma=c(0, 5, 0, 0), xpd=F)
data <- data.table(data)
d.breaks <- rev(data$delta %>% unique())

cols <- sequential_hcl(n=length(d.breaks), palette="viridis")

for(i in 1:length(d.breaks)){
  dat.k <- data[(delta == d.breaks[i]) & (valid == T)]
  if(i == 1){
    plot(dat.k$ve_inf ~ dat.k$ve, type='l', col=cols[i], lwd=2,
         ylim=c(0, 1),
         ylab="Vaccine Efficacy Inferred",
         xlab="Vaccine Efficacy",
         cex.lab=1.3, las=1, axes=F)
    box(bty="l")
    axis(1, cex.axis=1.3, las=1)
    axis(2, cex.axis=1.3, las=1)
    par(xpd=FALSE)
    abline(a=0, b=1, lty=2)
    mtext(text="A)", side=3, line=-2, outer=T, adj=-0.1, cex=1.3)
  } else {
    lines(dat.k$ve_inf ~ dat.k$ve, col=cols[i], lwd=2)
  }
}
pos <- legend("topleft", inset=0.08, cex=1.3, title.adj=0.5, col=cols, lwd=2, ncol=2, legend=c(
  paste0(sprintf("%2.0f", (1-d.breaks)*100))), bty="n",
  title=expression(underline("Reduction in Transmission\n     Given Asymptomatic, %"))
)
xleft <- pos$rect[["left"]]
ytop <- pos$rect[["top"]]
ybottom <- ytop - pos$rect[["h"]]
xright <- xleft + pos$rect[["w"]]
rect(xleft, ybottom+0.025, xright, ytop+0.08)

data <- data.inf
cols <- sequential_hcl(n=length(k.breaks), palette="viridis")
for(i in 1:length(k.breaks)){
  dat.k <- data[k == k.breaks[i]]
  if(i == 1){
    plot(dat.k$inferred ~ dat.k$ve, type='l', col=cols[i], lwd=2,
         ylim=c(0.3, 1.0),
         ylab="Vaccine Efficacy Inferred",
         xlab="Vaccine Efficacy",
         las=1, axes=F, cex.lab=1.3)
    box(bty="l")
    axis(1, cex.axis=1.3, las=1)
    axis(2, cex.axis=1.3, las=1)
    par(xpd=FALSE)
    abline(a=0, b=1, lty=2)
    mtext(text="B)",side=3,line=-32,outer=TRUE,adj=-0.1, cex=1.3)
  } else {
    lines(dat.k$inferred ~ dat.k$ve, col=cols[i], lwd=2)
  }
}
legend("topleft", inset=0.08, col=cols, lwd=2, ncol=2, legend=c(
  "1", paste0(k.breaks[2:length(k.breaks)])),
  border=1, title=expression(underline("Testing Interval, days")), cex=1.3)
dev.off()
