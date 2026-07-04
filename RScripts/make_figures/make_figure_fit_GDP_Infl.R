# ============================================
#         INFLATION AND GDP EVOLUTION
# ============================================


### INFLATION

horizon.of.interest <- "4Q"
date.of.interest <- "2020-04-15"

data.of.interest <- survey.DATA.US.with.param %>%
  dplyr::select(date, contains(horizon.of.interest)) %>%
  filter(date > "1999-01-01")

# Inflation Data
path <- paste0(path_graph,"y.o.y.inflation.pdf")
pdf(path, width = 8, height = 5)
par(plt=c(.1,.95,.1,.95))
plot(DATA.US$date[!is.na(DATA.US$US.infl)],
     DATA.US$US.infl[!is.na(DATA.US$US.infl)],type='l',
     #col="#ff8243",
     las=1,
     ylab="",xlab="", lwd=2)
make_recessions()
grid()
dev.off()

# GDP data
path <- paste0(path_graph,"annual.gdp.growth.pdf")
pdf(path, width = 8, height = 5)
par(plt=c(.1,.95,.1,.95))
plot(DATA.G.US$date[!is.na(DATA.G.US$US.gdp.growth)],
     DATA.G.US$US.gdp.growth[!is.na(DATA.G.US$US.gdp.growth)],type='l',
     #col="#56B4E9",
     las=1,
     ylab="",xlab="", lwd=2)
make_recessions()
grid()
dev.off()


# Figure with approximate linear growth rates:

GDP_raw       <- read.csv("./data/US/raw/GDPC1.csv")
GDP_raw$logGDPC1 <- log(GDP_raw$GDPC1)
GDP_raw$DATE  <- as.Date(GDP_raw$DATE)

T <- dim(GDP_raw)[1]

k <- 1
GDP_raw$dy1 <- c(rep(NaN,k),GDP_raw$logGDPC1[(k+1):T] - GDP_raw$logGDPC1[1:(T-k)])
# plot(GDP_raw$DATE,GDP_raw$dy1,type="l")
# lines(as.Date(observables.with.dates$date),
#       observables.with.dates$US.growth/100,col="red")

k <- 4
GDP_raw$dy4 <- c(rep(NaN,k),GDP_raw$logGDPC1[(k+1):T] - GDP_raw$logGDPC1[1:(T-k)])
GDP_raw$dy4_geom <- c(rep(NaN,k),GDP_raw$GDPC1[(k+1):T]/GDP_raw$GDPC1[1:(T-k)]-1)
plot(GDP_raw$DATE,GDP_raw$dy4,type="l")
lines(GDP_raw$DATE,GDP_raw$dy4_geom,col="grey40")


GDP_raw$GDP4 <- lag(GDP_raw$GDPC1,0) +
  lag(GDP_raw$GDPC1,1) +
  lag(GDP_raw$GDPC1,2) +
  lag(GDP_raw$GDPC1,3)
GDP_raw$annualGDPgrowth <- GDP_raw$GDP4/lag(GDP_raw$GDP4,4) - 1

GDP_raw$dy4_approx <- 1/4 * (
  1 * lag(GDP_raw$dy1,0) +
    2 * lag(GDP_raw$dy1,1) +
    3 * lag(GDP_raw$dy1,2) +
    4 * lag(GDP_raw$dy1,3) +
    3 * lag(GDP_raw$dy1,4) +
    2 * lag(GDP_raw$dy1,5) +
    1 * lag(GDP_raw$dy1,6))

path <- paste0(path_graph,"figure_approx.pdf")
pdf(path, width = 6, height = 4)
par(plt=c(.1,.95,.1,.95))

plot(GDP_raw$DATE,GDP_raw$annualGDPgrowth,type="l",lwd=2,
     xlab="",ylab="",las=1)
grid()
lines(GDP_raw$DATE,GDP_raw$dy4_approx,lty=3,col="grey",lwd=4)

legend("topright",
       lty=c(1,3),
       c("y-o-y growth rate of annual GDP", "Affine approximation with MA(6)") ,
       col = c("black","grey"),
       pt.cex=1.5,
       seg.len=3,
       bg="white",
       lwd=c(2,3),
       cex=1)

dev.off()
