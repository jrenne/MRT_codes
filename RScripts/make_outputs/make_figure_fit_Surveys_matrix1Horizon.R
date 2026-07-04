# Convert in matrix form for simplicity to plot
observables <- as.matrix(observables)
indic.to.chart <- c(1,grep(paste(c("US.4Q."), collapse = "|"), colnames(observables))) # Inflation 4Q
#indic.to.chart <- c(22,grep(paste(c("G.4Q."), collapse = "|"), colnames(observables))) # GDP 4Q
#indic.to.chart <- c(1,2,5,8) # 1y + inflation
#indic.to.chart <- c(11,12,15,18) #1 + growth
#indic.to.chart <- c(3,4,13,14) # inflation and growth survey PEs 2-5y, 
#indic.to.chart <- c(6,7,16,17) # inflation and growth survey var 2-5y,
#indic.to.chart <- c(9,10,19,20) # inflation and growth survey k3rd 2-5y,
#indic.to.chart <- c(21,22) # fit HP filters,
# Plot fit obtained with filter (all):

path <- paste0(path_graph,"Fitted.1.year.inflation" ,".pdf")
#path <- paste0(path_graph,"Fitted.1.year.gdp.growth" ,".pdf")
#path <- paste0(path_graph,"Fitted.PE.2.5" ,".pdf")
#path <- paste0(path_graph,"Fitted.Var.2.5" ,".pdf")
#path <- paste0(path_graph,"Fitted.k3rd.2.5" ,".pdf")
#path <- paste0(path_graph,"Fitted.HP.filter" ,".pdf")

pdf(path, width = 5, height = 4)
b <- length(indic.to.chart)
#m.1 <- matrix(c(1,2,3,3), nrow = (2),ncol = (2),byrow = TRUE)
#layout(mat = m.1, heights = c(rep(0.8/(b/2),b/2),0.2))
m.1 <- matrix(c(seq(1,b),rep(b+1,b/2)), nrow = (b/2+1),ncol = (b/2),byrow = TRUE)
layout(mat = m.1, heights = c(rep(0.86/(b/2),b/2),0.14))
#layout.show(5)

for(i in indic.to.chart){
  
  h <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(observables)))
  
  par(mar=c(2, 3, 1, 1))
  
  plot(vec.dates$date, KF.result4$fitted.obs.smoothed[,i],type='l',col='darkgrey', lwd=2, ylab="", xlab="", main="",
       ylim=c(min(KF.result4$fitted.obs.smoothed[,i], observables[,i],na.rm=TRUE)-0.15,
              max(KF.result4$fitted.obs.smoothed[,i], observables[,i],na.rm=TRUE)+2.15))
  points(vec.dates$date, observables[,i], col="black", lwd=0.25, pch=19, cex = 0.75)
  
  # Legend condition
  ## Condition for Inflation
  if(grepl("infl", colnames(observables)[i])){
    legend("topleft",legend=bquote(pi[t*"-1,"*t]), cex=1, bty = "n")
  } 
  
  ## Condition for PE
  if(grepl("pe", colnames(observables)[i]) & !grepl("G", colnames(observables)[i])){
    legend("topleft",legend=bquote(E[t]*"("*pi[t*","*t*"+"*.(h[i])]*")"), cex=1, bty = "n")
  } 
  
  ## Condition for Var
  if(grepl("var", colnames(observables)[i]) & !grepl("G", colnames(observables)[i])){
    legend("topleft",legend=bquote(Var[t]*"("*pi[t*","*t*"+"*.(h[i])]*")"), cex=1, bty = "n")
  } 
  
  ## Condition for 3rd Cumulant
  if(grepl("k3rd", colnames(observables)[i]) & !grepl("G", colnames(observables)[i])){
    legend("topleft",legend=bquote(mu["3,"*t]*"("*pi[t*","*t*"+"*.(h[i])]*")"), cex=1, bty = "n")
  } 
  
  ## Condition for 4th Cumulant
  if(grepl("k4th", colnames(observables)[i]) & !grepl("G", colnames(observables)[i])){
    legend("topleft",legend=bquote(mu["4,"*t]*"("*pi[t*","*t*"+"*.(h[i])]*")"), cex=1, bty = "n")
  } 
  
  ## Condition for GDP
  if(grepl("growth", colnames(observables)[i])){
    legend("topleft",legend=bquote(Delta*y[t*"-1,"*t]), cex=1, bty = "n")
  } 
  
  ## Condition for PE GDP
  if(grepl("pe", colnames(observables)[i]) & grepl("G", colnames(observables)[i])){
    legend("topleft",legend=bquote(E[t]*"("*Delta*y[t*","*t*"+"*.(h[i])]*")"), cex=1, bty = "n")
  } 
  
  ## Condition for Var GDP
  if(grepl("var", colnames(observables)[i]) & grepl("G", colnames(observables)[i])){
    legend("topleft",legend=bquote(Var[t]*"("*Delta*y[t*","*t*"+"*.(h[i])]*")"), cex=1, bty = "n")
  } 
  
  ## Condition for 3rd Cumulant GDP
  if(grepl("k3rd", colnames(observables)[i]) & grepl("G", colnames(observables)[i])){
    legend("topleft",legend=bquote(mu["3,"*t]*"("*Delta*y[t*","*t*"+"*.(h[i])]*")"), cex=1, bty = "n")
  } 
  
  ## Condition for 4th Cumulant GDP
  if(grepl("k4th", colnames(observables)[i]) & grepl("G", colnames(observables)[i])){
    legend("topleft",legend=bquote(mu["4,"*t]*"("*Delta*y[t*","*t*"+"*.(h[i])]*")"), cex=1, bty = "n")
  } 
  
}

# Plot legend below
#expression(paste(E[t](pi[t*","*t+paste(h)])))
par(mar=c(0.1, 0.1, 0.1, 0.1))
plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",inset = 0, title="Variables:", c("Model-implied","Survey observations"), 
       fill=c("darkgrey","black"), cex=0.86, horiz = FALSE)

dev.off()
