# Extract matrices from the list "all.parameters
H <- KF.load.final$all.parameters$H
A <- KF.load.final$all.parameters$A

X <- KF.load.final$X
matrix.xi.tT <- KF.result4$xi.tT


indic.var.k3rd <- grep("var|k3rd", colnames(observables))
indic.supply <- c(1,2)
indic.demand <- c(3,4)
indic.constant <- c(5)

# Compute fitted values for observables of Skewness and Kurtosis
fitted.obs.var.k3rd <- X %*% A[,indic.var.k3rd] + matrix.xi.tT%*%H[,indic.var.k3rd] # Fitted observables

# Check if fitted obs correct
#sum(fitted.obs.var.k3rd ==KF.result4$fitted.obs.smoothed[,indic.var.k3rd])

# Perform decomposition
fitted.obs.var.k3rd.supply <- matrix.xi.tT[,((n+1):(n+q))[indic.supply]]%*%H[((n+1):(n+q))[indic.supply],indic.var.k3rd] 
fitted.obs.var.k3rd.demand <- matrix.xi.tT[,((n+1):(n+q))[indic.demand]]%*%H[((n+1):(n+q))[indic.demand],indic.var.k3rd] 
fitted.obs.var.k3rd.constant <- X %*% A[,indic.var.k3rd] + matrix.xi.tT[,((n+1):(n+q))[indic.constant], drop = FALSE]%*%H[((n+1):(n+q))[indic.constant],indic.var.k3rd, drop = FALSE]

# Check if decomposition is equal to fitted obs
#sum(round(fitted.obs.var.k3rd.supply + fitted.obs.var.k3rd.demand + fitted.obs.var.k3rd.constant,6) == round(fitted.obs.var.k3rd,6)) 

fitted.obs.var.k3rd.supply <- as.data.frame(fitted.obs.var.k3rd.supply)
colnames(fitted.obs.var.k3rd.supply) <- colnames(observables)[indic.var.k3rd]
fitted.obs.var.k3rd.supply <- fitted.obs.var.k3rd.supply %>%
  mutate(date=vec.dates$date,
         group="supply") %>%
  dplyr::select(date, group, everything())

fitted.obs.var.k3rd.demand <- as.data.frame(fitted.obs.var.k3rd.demand)
colnames(fitted.obs.var.k3rd.demand) <- colnames(observables)[indic.var.k3rd]
fitted.obs.var.k3rd.demand <- fitted.obs.var.k3rd.demand %>%
  mutate(date=vec.dates$date,
         group="demand") %>%
  dplyr::select(date, group, everything())

fitted.obs.var.k3rd.constant <- as.data.frame(fitted.obs.var.k3rd.constant)
colnames(fitted.obs.var.k3rd.constant) <- colnames(observables)[indic.var.k3rd]
fitted.obs.var.k3rd.constant <- fitted.obs.var.k3rd.constant %>%
  mutate(date=vec.dates$date,
         group="constant") %>%
  dplyr::select(date, group, everything())

fitted.obs.var.k3rd <- as.data.frame(fitted.obs.var.k3rd)
colnames(fitted.obs.var.k3rd) <- colnames(observables)[indic.var.k3rd]
fitted.obs.var.k3rd <- fitted.obs.var.k3rd %>%
  mutate(date=vec.dates$date) %>%
  dplyr::select(date, everything())

higher.order.moments.decomposition <- rbind(fitted.obs.var.k3rd.supply,
                                           fitted.obs.var.k3rd.demand,
                                           fitted.obs.var.k3rd.constant) %>%
  arrange(date)

# data.to.consider.old <- higher.order.moments.decomposition %>%
#   dplyr::select(date,group,SPF.US.7Q.beta.var) %>%
#   rename(value = SPF.US.7Q.beta.var) %>%
#   full_join(fitted.obs.var.k3rd %>%
#               dplyr::select(date,SPF.US.7Q.beta.var) %>%
#               rename(fit.value = SPF.US.7Q.beta.var),
#             by="date")


for(i in 3:(dim(higher.order.moments.decomposition)[2])){

  name.to.consider <- colnames(higher.order.moments.decomposition)[i]
  
  path <- paste0(path_graph,name.to.consider,".decomposition.pdf")
  pdf(path, width = 7, height = 3.5)
  
  data.to.consider <- build_data_to_consider(name.to.consider)
  par(plt=c(.1,.95,.15,.95))
  p <- ggplot(data=data.to.consider, aes(fill=group, y=value, x=date)) + 
    geom_rect(data = plot_data_rec, aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey90", alpha = 0.45) +
    geom_bar(position=position_stack(), stat="identity") +
    ylab("") +
    xlab("") +
    geom_line(aes(x=date, y=fit.value, color ='var', linetype="var"), lwd=1, show.legend = FALSE) +
    geom_hline(yintercept=0, lwd=0.5, col="black", linetype=2) +
    ggtitle(label = "") +
    scale_fill_manual('', values=c( "darkgrey", "#ff8243", "#56B4E9"), labels=c("Residual", "Demand", "Supply")) +
    scale_color_manual('', values=c("black")) +
    scale_linetype_manual('',values=c('var'="solid")) +
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 10)
    ) 
  print(p)
  
 dev.off() 
 
}
