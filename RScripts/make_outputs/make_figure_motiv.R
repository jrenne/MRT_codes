# ============================================
#       PLOTS DISTRIBUTION MOTIVATION
# ============================================

path <- paste0(path_graph,"tail.risks.pdf")
pdf(path, width = 8, height = 4)

par(mfrow=c(1,2))
par(plt=c(.1,.95,.1,.85))

# plot(seq(-5,7,0.01),dnorm(seq(-5,7,0.01),mean=1.89,sd=1.8), type="l", ylim=c(0,0.25), col="red", xlab="", ylab="", lwd=2, main="Anticipating no shock", cex.main = 0.9)
# lines(seq(-5,7,0.01),dnorm(seq(-5,7,0.01),mean=1.59, sd=2.2), col="black", lwd=2)
# abline(v=moments.incompl.beta(c(1.5,5,0,6),-4.7647,14-4.7647)$Mean, lty=2, col="red")
# abline(v=moments.incompl.beta(c(6,2.5,-6,5),-10,4)$Mean, col="black", lty=2)
# legend("topleft",legend=c(expression( mu^(pi)*"=1.9"),
#                           expression( sigma^(pi)*"=1.8"),
#                           expression( mu[3]^(pi)*"=0.0")),
#        cex=1, bty = "n", y.intersp=1.2, text.col="red")
# legend("topright", legend=c(expression( mu^(Delta*y)*"=1.6"),
#                             expression( sigma^(Delta*y)*"=2.2"),
#                             expression( mu[3]^(Delta*y)*"=0.0")), 
#        cex=1, bty = "n", y.intersp=1.2)

x <- seq(-11,9,0.01)
y <- pdf.incompl.beta(x,c(5,1.5,-6,0),-7+0.775,7+0.775)
y[is.na(y)] <- 0
plot(x,y,
     type="l",
     ylim=c(0,0.25),
     col="red", xlab="", ylab="", lwd=2,  main="Large (negative) demand shock",
     cex.main = 0.9)
y <- pdf.incompl.beta(x,c(6,2.5,-6,5),-10+0.924,4+0.924)
y[is.na(y)] <- 0
lines(x,y,
      col="black", lwd=2)
abline(v=moments.incompl.beta(c(5,1.5,-6,0),-7+0.775,7+0.775)$Mean, lty=2, col="red")
abline(v=moments.incompl.beta(c(6,2.5,-6,5),-10+0.924,4+0.924)$Mean, col="black", lty=2)
legend("topleft",legend=c(expression( mu^(pi)*"=1.9"),
                          expression( sigma^(pi)*"=8.7"),
                          expression( mu[3]^(pi)*"=-0.5")),
       cex=1, bty = "n", y.intersp=1.2, text.col="red")
legend("topright", legend=c(expression( mu^(Delta*y)*"=1.6"),
                            expression( sigma^(Delta*y)*"=12.1"),
                            expression( mu[3]^(Delta*y)*"=-0.4")), 
       cex=1, bty = "n", y.intersp=1.2)

legend(x=-10,y=.15,
       c("Inflation", 'GDP growth') ,
       col=c("red", "black"),
       lty=1,bg="white",lwd=2,
       cex=0.9)

x <- seq(-7,14,0.01)
y <- pdf.incompl.beta(x,c(1.5,5,0,6),-4.7647+0.775,14-4.7647+0.775)
y[is.na(y)] <- 0
plot(x,y, type="l", ylim=c(0,0.25), col="red", xlab="", ylab="", lwd=2, main="Large (negative) supply shock", cex.main = 0.9)
y <- pdf.incompl.beta(x,c(6,2.5,-6,5),-10+0.924,4+0.924)
y[is.na(y)] <- 0
lines(x,y, col="black", lwd=2)
abline(v=moments.incompl.beta(c(1.5,5,0,6),-4.7647+0.775,14-4.7647+0.775)$Mean, lty=2, col="red")
abline(v=moments.incompl.beta(c(6,2.5,-6,5),-10+0.924,4+0.924)$Mean, col="black", lty=2)
legend("topleft",legend=c(expression( mu^(pi)*"=1.9"),
                          expression( sigma^(pi)*"=8.7"),
                          expression( mu[3]^(pi)*"=0.5")),
       cex=1, bty = "n", y.intersp=1.2, text.col="red")
legend("topright", legend=c(expression( mu^(Delta*y)*"=1.6"),
                            expression( sigma^(Delta*y)*"=12.1"),
                            expression( mu[3]^(Delta*y)*"=-0.4")), 
       cex=1, bty = "n", y.intersp=1.2)


dev.off()

moments.incompl.beta(c(5,1.5,-6,0),-7+0.775,7+0.775)
moments.incompl.beta(c(1.5,5,0,6),-4.7647+0.775,14-4.7647+0.775)
moments.incompl.beta(c(6,2.5,-6,5),-10+0.924,4+0.924)
