# ==============================================================================
#    DEMAND SUPPLY DECOMPOSITION
# ==============================================================================

# res_TrendCycle  <- make_TrendCycle(Model.final,
#                                    KF.result4)
# 
# Trend.inflation <- res_TrendCycle$Trend.inflation
# Trend.gdp       <- res_TrendCycle$Trend.gdp
# Cycle.inflation <- res_TrendCycle$Cycle.inflation
# Cycle.gdp       <- res_TrendCycle$Cycle.gdp
# model.p.t       <- res_TrendCycle$model.p.t
# model.gdp.t     <- res_TrendCycle$model.gdp.t

# Compute growth rates:
res_growth_rates <- compute_growth_rates(Model.final,
                                         KF.result4)
model.base.inflation.quarter <- res_growth_rates$model.base.inflation.quarter
model.base.gdp.quarter       <- res_growth_rates$model.base.gdp.quarter
model.inflation.annual       <- res_growth_rates$model.inflation.annual
model.gdp.annual             <- res_growth_rates$model.gdp.annual

# Operate trend/cycle decomposition:
res_TrendCycle  <- compute_TrendCycle(Model.final,
                                      KF.result4)
Trend.inflation <- res_TrendCycle$Trend.inflation
Trend.gdp       <- res_TrendCycle$Trend.gdp
Cycle.inflation <- res_TrendCycle$Cycle.inflation
Cycle.gdp       <- res_TrendCycle$Cycle.gdp
model.p.t       <- res_TrendCycle$model.p.t
model.gdp.t     <- res_TrendCycle$model.gdp.t

# Operate demand/supply decompositions:
res_decompositions <- compute_decompositions(Model.final,
                                             KF.result4)

Cycle.inflation.supply <- res_decompositions$Cycle.inflation.supply
Cycle.inflation.demand <- res_decompositions$Cycle.inflation.demand
Cycle.gdp.supply <- res_decompositions$Cycle.gdp.supply
Cycle.gdp.demand <- res_decompositions$Cycle.gdp.demand
Inflation.y.o.y.supply <- res_decompositions$Inflation.y.o.y.supply
Inflation.y.o.y.demand <- res_decompositions$Inflation.y.o.y.demand
gdp.annual.demand <- res_decompositions$gdp.annual.demand
gdp.annual.supply <- res_decompositions$gdp.annual.supply
Cycle.inflation.decomposition <- res_decompositions$Cycle.inflation.decomposition
Cycle.gdp.decomposition       <- res_decompositions$Cycle.gdp.decomposition
Trend.inflation.decomposition <- res_decompositions$Trend.inflation.decomposition
Trend.gdp.decomposition       <- res_decompositions$Trend.gdp.decomposition
Inflation.y.o.y.decomposition <- res_decompositions$Inflation.y.o.y.decomposition
gdp.annual.decomposition      <- res_decompositions$gdp.annual.decomposition
Inflation.q.o.q.decomposition <- res_decompositions$Inflation.q.o.q.decomposition
gdp.q.o.q.decomposition      <- res_decompositions$gdp.q.o.q.decomposition

# Compute delta.inf and delta.gdp for IC
delta.inf  <- matrix(c(Model.final$delta[,1],rep(0,q)),1,Model.final$n+Model.final$q) %x% 
  matrix(c(Model.final$delta[,1],rep(0,q)),1,Model.final$n+Model.final$q)

delta.gdp   <- matrix(c(Model.final$delta[,2],rep(0,q)),1,Model.final$n+Model.final$q) %x% 
  matrix(c(Model.final$delta[,2],rep(0,q)),1,Model.final$n+Model.final$q)

# Compute variance of q.o.q inflation and gdp
## Get Variance of latent factors
P.tT <- KF.result4$P.tT
var.cov.inf <- delta.inf%*%t(P.tT)
var.cov.gdp <- delta.gdp%*%t(P.tT)

#### PLOTS
# Cycle HCPI

library(tis)
recession_fill <- "grey90"
recession_alpha <- 0.45
if(area=="US"){
  start_date <- as.Date(as.character(nberDates()[,1]), format = "%Y%m%d")
  start_date <- start_date[start_date >= "1981-07-15"]
  end_date <- as.Date(as.character(nberDates()[,2]), format = "%Y%m%d")
  end_date <- end_date[end_date >= "1981-07-15"]
  plot_data_rec <- data.frame(start_date=start_date,end_date=end_date)
  plot_data_rec$start_date <- as.POSIXct(plot_data_rec$start_date)
  plot_data_rec$end_date <- as.POSIXct(plot_data_rec$end_date)
}
if(area=="EA"){
  plot_data_rec <- get_recessions_date() %>%
    dplyr::filter(start_date >= "1999-01-01") %>%
    mutate(start_date = as.POSIXct(start_date, tz = "UTC"),
           end_date = as.POSIXct(end_date, tz = "UTC"))
}


path <- paste0(path_graph,"log.hcpi.cyclical.decomposition.pdf")
pdf(path, width = 7, height = 3.5)
#par(plt=c(.1,.95,.15,.95))
ggplot(data=Cycle.inflation.decomposition, aes(fill=group, y=value, x=date)) + 
  geom_rect(data = plot_data_rec, aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = recession_fill, alpha = recession_alpha) +
  geom_bar(position=position_stack(), stat="identity") +
  #scale_fill_viridis(discrete=TRUE, name="",labels=c("Demand","Supply")) +
  ylab("") +
  xlab("") +
  geom_line(aes(x=date, y=Cycle.inflation, col="Cycle"), lwd=1) +
  geom_hline(yintercept=0, lwd=0.5, col="black", linetype=2) +
  ggtitle(label = "") +
  scale_fill_manual('', values=c( "#ff8243", "#56B4E9"), labels=c("Demand","Supply")) +
  scale_color_manual('', values=c("black")) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom")
dev.off()

# Cycle GDP
path <- paste0(path_graph,"log.gdp.cyclical.decomposition.pdf")
pdf(path, width = 7, height = 3.5)
par(plt=c(.1,.95,.15,.95))
ggplot(data=Cycle.gdp.decomposition, aes(fill=group, y=value, x=date)) + 
  geom_rect(data = plot_data_rec, aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = recession_fill, alpha = recession_alpha) +
  geom_bar(position=position_stack(), stat="identity") +
  ylab("") +
  xlab("") +
  geom_line(aes(x=date, y=Cycle.gdp, color = 'Cycle'), lwd=1) +
  geom_hline(yintercept=0, lwd=0.5, col="black", linetype=2) +
  ggtitle(label = "") +
  scale_fill_manual('', values=c( "#ff8243", "#56B4E9"), labels=c("Demand","Supply")) +
  scale_color_manual('', values=c("black")) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom") 
dev.off()


#### PLOTS

# Y.o.y inflation 
path <- paste0(path_graph,"y.o.y.inflation.decomposition.pdf")
pdf(path, width = 7, height = 3.5)
par(plt=c(.1,.95,.15,.95))
ggplot(data=Inflation.y.o.y.decomposition, aes(fill=group, y=value, x=date)) + 
  geom_rect(data = plot_data_rec, aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = recession_fill, alpha = recession_alpha) +
  geom_bar(position=position_stack(), stat="identity") +
  #scale_fill_viridis(discrete=TRUE, name="",labels=c("Demand","Supply")) +
  ylab("") +
  xlab("") +
  geom_line(aes(x=date, y=Inflation.y.o.y, col="Inflation", linetype="Inflation"), lwd=1) +
  geom_hline(yintercept=0, lwd=0.5, col="black", linetype=2) +
  geom_hline(aes(yintercept=4*Model.final$pi.bar[1], color="Average inflation", linetype="Average inflation"), lwd=0.5) +
  ggtitle(label = "") +
  scale_fill_manual('', values=c( "#ff8243", "#56B4E9"), labels=c("Demand","Supply")) +
  scale_color_manual('', values=c("darkgrey","black")) +
  scale_linetype_manual('',values=c("Average inflation"="dotted",'Inflation'="solid")) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10)#,
        #legend.title = element_text(size = 10) # Adjust the size as needed
  ) 

dev.off()

#  Y.o.y GDP
path <- paste0(path_graph,"annual.gdp.growth.decomposition.pdf")
pdf(path, width = 7, height = 3.5)
par(plt=c(.1,.95,.15,.95))
ggplot(data=gdp.annual.decomposition, aes(fill=group, y=value, x=date)) + 
  geom_rect(data = plot_data_rec, aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = recession_fill, alpha = recession_alpha) +
  geom_bar(position=position_stack(), stat="identity") +
  ylab("") +
  xlab("") +
  geom_line(aes(x=date, y=gdp.annual, color ='GDP growth', linetype="GDP growth"), lwd=1) +
  geom_hline(yintercept=0, lwd=0.5, col="black", linetype=2) +
  geom_hline(aes(yintercept=4*Model.final$pi.bar[2], color="Average GDP growth", linetype="Average GDP growth"), lwd=0.5) +
  ggtitle(label = "") +
  scale_fill_manual('', values=c( "#ff8243", "#56B4E9"), labels=c("Demand","Supply")) +
  scale_color_manual('', values=c("darkgrey","black")) +
  scale_linetype_manual('',values=c("Average GDP growth"="dotted",'GDP growth'="solid")) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10)
  ) 
dev.off()




# q.o.q inflation 
path <- paste0(path_graph,"q.o.q.inflation.decomposition.pdf")
pdf(path, width = 7, height = 3.5)
par(plt=c(.1,.95,.15,.95))
ggplot(data=Inflation.q.o.q.decomposition, aes(fill=group, y=value, x=date)) + 
  geom_rect(data = plot_data_rec, aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = recession_fill, alpha = recession_alpha) +
  geom_bar(position=position_stack(), stat="identity") +
  #scale_fill_viridis(discrete=TRUE, name="",labels=c("Demand","Supply")) +
  ylab("") +
  xlab("") +
  geom_line(aes(x=date, y=Inflation.q.o.q, col="q.o.q Inflation", linetype="q.o.q Inflation"), lwd=1) +
  geom_hline(yintercept=0, lwd=0.5, col="black", linetype=2) +
  geom_hline(aes(yintercept=Model.final$pi.bar[1], color="Average inflation", linetype="Average inflation"), lwd=0.5) +
  ggtitle(label = "") +
  scale_fill_manual('', values=c( "#ff8243", "#56B4E9"), labels=c("Demand","Supply")) +
  scale_color_manual('', values=c("darkgrey","black")) +
  scale_linetype_manual('',values=c("Average q.o.q inflation"="dotted",'q.o.q Inflation'="solid")) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10)#,
        #legend.title = element_text(size = 10) # Adjust the size as needed
  ) 

dev.off()

# q.o.q inflation with confidence interval
path <- paste0(path_graph,"q.o.q.inflation.decomposition.with.ic.pdf")
pdf(path, width = 7, height = 3.5)

m <- matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE)

layout(mat = m, heights = c(1))
par(plt=c(.1, .95, .2, .85))

Inflation.q.o.q.supply <- Inflation.q.o.q.decomposition %>% filter(group=="Inflation.q.o.q.supply")
Inflation.q.o.q.supply.ic <- cbind(c(Inflation.q.o.q.supply$value-1.96*sqrt(var.cov.inf)),
                                   c(Inflation.q.o.q.supply$value+1.96*sqrt(var.cov.inf)))

Inflation.q.o.q.demand <- Inflation.q.o.q.decomposition %>% filter(group=="Inflation.q.o.q.demand")
Inflation.q.o.q.demand.ic <- cbind(c(Inflation.q.o.q.demand$value-1.96*sqrt(var.cov.inf)),
                                   c(Inflation.q.o.q.demand$value+1.96*sqrt(var.cov.inf)))

Inflation.q.o.q.ic <- cbind(Inflation.q.o.q.supply.ic,
                            Inflation.q.o.q.demand.ic)

# Define the color with reduced alpha (e.g., 35% opacity)
mycol1 <- rgb(86 / 255, 180 / 255, 233 / 255, alpha = 0.35) #56B4E9
mycol2 <- rgb(255 / 255, 130 / 255, 67 / 255, alpha = 0.35) #FF8243

#Plot supply
plot(as.Date(Inflation.q.o.q.supply$date), Inflation.q.o.q.supply$value, type="l", col ="#56B4E9", xlab = "", ylab="", lwd=2, ylim=c(-2,4),#ylim=c(min(Inflation.q.o.q.ic), max(Inflation.q.o.q.ic)),
     las=1,)
polygon(c(as.Date(Inflation.q.o.q.supply$date), rev(as.Date(Inflation.q.o.q.supply$date))), c(Inflation.q.o.q.supply$value,rev(Inflation.q.o.q.supply.ic[,1])),col=mycol1, border=NA)                                          
polygon(c(as.Date(Inflation.q.o.q.supply$date), rev(as.Date(Inflation.q.o.q.supply$date))), c(Inflation.q.o.q.supply$value,rev(Inflation.q.o.q.supply.ic[,2])),col=mycol1, border=NA)                                          
grid()
make_recessions()

#Plot supply
plot(as.Date(Inflation.q.o.q.demand$date), Inflation.q.o.q.demand$value, type="l", col ="#FF8243", xlab = "", ylab="", lwd=2, ylim=c(-2,4),#ylim=c(min(Inflation.q.o.q.ic), max(Inflation.q.o.q.ic)),
     las=1,)
polygon(c(as.Date(Inflation.q.o.q.demand$date), rev(as.Date(Inflation.q.o.q.demand$date))), c(Inflation.q.o.q.demand$value,rev(Inflation.q.o.q.demand.ic[,1])),col=mycol2, border=NA)                                          
polygon(c(as.Date(Inflation.q.o.q.demand$date), rev(as.Date(Inflation.q.o.q.demand$date))), c(Inflation.q.o.q.demand$value,rev(Inflation.q.o.q.demand.ic[,2])),col=mycol2, border=NA)                                          
grid()
make_recessions()

dev.off()


#  q.o.q GDP
path <- paste0(path_graph,"q.o.q.gdp.growth.decomposition.pdf")
pdf(path, width = 7, height = 3.5)
par(plt=c(.1,.95,.15,.95))
ggplot(data=gdp.q.o.q.decomposition, aes(fill=group, y=value, x=date)) + 
  geom_rect(data = plot_data_rec, aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = recession_fill, alpha = recession_alpha) +
  geom_bar(position=position_stack(), stat="identity") +
  ylab("") +
  xlab("") +
  geom_line(aes(x=date, y=gdp.q.o.q, color ='q.o.q GDP growth', linetype="q.o.q GDP growth"), lwd=1) +
  geom_hline(yintercept=0, lwd=0.5, col="black", linetype=2) +
  geom_hline(aes(yintercept=Model.final$pi.bar[2], color="Average q.o.q GDP growth", linetype="Average q.o.q GDP growth"), lwd=0.5) +
  ggtitle(label = "") +
  scale_fill_manual('', values=c( "#ff8243", "#56B4E9"), labels=c("Demand","Supply")) +
  scale_color_manual('', values=c("darkgrey","black")) +
  scale_linetype_manual('',values=c("Average q.o.q GDP growth"="dotted",'q.o.q GDP growth'="solid")) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 10)
  ) 
dev.off()



# q.o.q gdo with confidence interval
path <- paste0(path_graph,"q.o.q.gdp.growth.decomposition.with.ic.pdf")
pdf(path, width = 7, height = 3.5)

m <- matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE)

layout(mat = m, heights = c(1))

par(plt=c(.1, .95, .2, .85))

gdp.q.o.q.supply <- gdp.q.o.q.decomposition %>% filter(group=="gdp.q.o.q.supply")
gdp.q.o.q.supply.ic <- cbind(c(gdp.q.o.q.supply$value-1.96*sqrt(var.cov.gdp)),
                                   c(gdp.q.o.q.supply$value+1.96*sqrt(var.cov.gdp)))

gdp.q.o.q.demand <- gdp.q.o.q.decomposition %>% filter(group=="gdp.q.o.q.demand")
gdp.q.o.q.demand.ic <- cbind(c(gdp.q.o.q.demand$value-1.96*sqrt(var.cov.gdp)),
                                   c(gdp.q.o.q.demand$value+1.96*sqrt(var.cov.gdp)))

gdp.q.o.q.ic <- cbind(gdp.q.o.q.supply.ic,
                            gdp.q.o.q.demand.ic)

# Define the color with reduced alpha (e.g., 35% opacity)
mycol1 <- rgb(86 / 255, 180 / 255, 233 / 255, alpha = 0.35) #56B4E9
mycol2 <- rgb(255 / 255, 130 / 255, 67 / 255, alpha = 0.35) #FF8243

#Plot supply
plot(as.Date(gdp.q.o.q.supply$date), gdp.q.o.q.supply$value, type="l", col ="#56B4E9", xlab = "", ylab="", lwd=2, ylim=c(-8,6), #ylim=c(min(gdp.q.o.q.ic), max(gdp.q.o.q.ic)),
     las=1)
polygon(c(as.Date(gdp.q.o.q.supply$date), rev(as.Date(gdp.q.o.q.supply$date))), c(gdp.q.o.q.supply$value,rev(gdp.q.o.q.supply.ic[,1])),col=mycol1, border=NA)                                          
polygon(c(as.Date(gdp.q.o.q.supply$date), rev(as.Date(gdp.q.o.q.supply$date))), c(gdp.q.o.q.supply$value,rev(gdp.q.o.q.supply.ic[,2])),col=mycol1, border=NA)                                          
grid()
make_recessions()

#Plot supply
plot(as.Date(gdp.q.o.q.demand$date), gdp.q.o.q.demand$value, type="l", col ="#FF8243", xlab = "", ylab="", lwd=2, ylim=c(-8,6), #ylim=c(min(gdp.q.o.q.ic), max(gdp.q.o.q.ic)),
     las=1,)
polygon(c(as.Date(gdp.q.o.q.demand$date), rev(as.Date(gdp.q.o.q.demand$date))), c(gdp.q.o.q.demand$value,rev(gdp.q.o.q.demand.ic[,1])),col=mycol2, border=NA)                                          
polygon(c(as.Date(gdp.q.o.q.demand$date), rev(as.Date(gdp.q.o.q.demand$date))), c(gdp.q.o.q.demand$value,rev(gdp.q.o.q.demand.ic[,2])),col=mycol2, border=NA)                                          
grid()
make_recessions()

dev.off()
