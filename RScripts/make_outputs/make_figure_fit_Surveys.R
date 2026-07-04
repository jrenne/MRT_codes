# ==============================================================================
#             PLOTS FIT SURVEY
# ==============================================================================

# Convert in matrix form for simplicity to plot
observables <- as.matrix(observables)
# Plot fit obtained with filter (all):

for(i in 1:dim(observables)[2]){
  path <- paste0(path_graph,i,".Fitted.", colnames(observables)[i],".pdf")
  pdf(path, width = 5, height = 2.5)
  
  h <- as.numeric(gsub(".*?([0-9]+).*", "\\1", colnames(observables)))
  
  ## Condition for Inflation
  if(grepl("infl", colnames(observables)[i])){
    titletopleft <- bquote(pi[t*"-1,"*t])
    
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    observed_variable      <- observables[,i]
  } 
  ## Condition for PE & Inflation
  if(grepl("pe", colnames(observables)[i]) & !grepl("G", colnames(observables)[i])){
    titletopleft <- bquote(E[t]*"("*pi[t*","*t*"+"*.(h[i])]*")")
    
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    observed_variable      <- observables[,i]
  } 
  ## Condition for Var & Inflation
  if(grepl("var", colnames(observables)[i]) & !grepl("G", colnames(observables)[i])){
    titletopleft <- bquote(Var[t]*"("*pi[t*","*t*"+"*.(h[i])]*")")
    
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    observed_variable      <- observables[,i]
  } 
  ## Condition for 3rd Cumulant & Inflation
  if(grepl("k3rd", colnames(observables)[i]) & !grepl("G", colnames(observables)[i])){
    titletopleft <- bquote(mu["3,"*t]*"("*pi[t*","*t*"+"*.(h[i])]*")")
    # Detect variance:
    indic_variance <- which(grepl("var", colnames(observables))&
                              !grepl("G", colnames(observables))&
                              (h[i]==h))
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    model_implied_variance <- KF.result4$fitted.obs.smoothed[,indic_variance]
    model_implied_variable <- model_implied_variable/model_implied_variance^(3/2)
    
    observed_variable <- observables[,i]
    observed_variance <- observables[,indic_variance]
    observed_variable <- observed_variable/observed_variance^(3/2)
  }
  ## Condition for 4th Cumulant & Inflation
  if(grepl("k4th", colnames(observables)[i]) & !grepl("G", colnames(observables)[i])){
    titletopleft <- bquote(mu["4,"*t]*"("*pi[t*","*t*"+"*.(h[i])]*")")
    
    # Detect variance:
    indic_variance <- which(grepl("var", colnames(observables))&
                              !grepl("G", colnames(observables))&
                              (h[i]==h))
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    model_implied_variance <- KF.result4$fitted.obs.smoothed[,indic_variance]
    model_implied_variable <- model_implied_variable + 3*model_implied_variance^2
    model_implied_variable <- model_implied_variable/model_implied_variance^2
      
    observed_variable <- observables[,i]
    observed_variance <- observables[,indic_variance]
    observed_variable <- observed_variable + 3*observed_variance^2
    observed_variable <- observed_variable/observed_variance^2
  } 
  ## Condition for GDP
  if(grepl("growth", colnames(observables)[i])){
    titletopleft <- bquote(Delta*y[t*"-1,"*t])
    
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    observed_variable      <- observables[,i]
  } 
  ## Condition for PE GDP
  if(grepl("pe", colnames(observables)[i]) & grepl("G", colnames(observables)[i])){
    titletopleft <- bquote(E[t]*"("*Delta*y[t*","*t*"+"*.(h[i])]*")")
    
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    observed_variable      <- observables[,i]
  } 
  ## Condition for Var GDP
  if(grepl("var", colnames(observables)[i]) & grepl("G", colnames(observables)[i])){
    titletopleft <- bquote(Var[t]*"("*Delta*y[t*","*t*"+"*.(h[i])]*")")
    
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    observed_variable      <- observables[,i]
  } 
  ## Condition for 3rd Cumulant GDP
  if(grepl("k3rd", colnames(observables)[i]) & grepl("G", colnames(observables)[i])){
    titletopleft <- bquote(mu["3,"*t]*"("*Delta*y[t*","*t*"+"*.(h[i])]*")")
    
    # Detect variance:
    indic_variance <- which(grepl("var", colnames(observables))&
                              grepl("G", colnames(observables))&
                              (h[i]==h))
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    model_implied_variance <- KF.result4$fitted.obs.smoothed[,indic_variance]
    model_implied_variable <- model_implied_variable/model_implied_variance^(3/2)
    
    observed_variable <- observables[,i]
    observed_variance <- observables[,indic_variance]
    observed_variable <- observed_variable/observed_variance^(3/2)
  }
  ## Condition for 4th Cumulant GDP
  if(grepl("k4th", colnames(observables)[i]) & grepl("G", colnames(observables)[i])){
    titletopleft <- bquote(mu["4,"*t]*"("*Delta*y[t*","*t*"+"*.(h[i])]*")")
    
    # Detect variance:
    indic_variance <- which(grepl("var", colnames(observables))&
                              grepl("G", colnames(observables))&
                              (h[i]==h))
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    model_implied_variance <- KF.result4$fitted.obs.smoothed[,indic_variance]
    model_implied_variable <- model_implied_variable + 3*model_implied_variance^2
    model_implied_variable <- model_implied_variable/model_implied_variance^2
    
    observed_variable <- observables[,i]
    observed_variance <- observables[,indic_variance]
    observed_variable <- observed_variable + 3*observed_variance^2
    observed_variable <- observed_variable/observed_variance^2
  }
  
  ## Condition for Inflation
  if(grepl("cycle", colnames(observables)[i])){
    
    if(grepl("cpi", colnames(observables)[i])){
      titletopleft <- bquote(C[t]^{(pi)})
    }
    
    if(grepl("gdp", colnames(observables)[i])){
      titletopleft <- bquote(C[t]^{(Delta*y)})
    }
    
    model_implied_variable <- KF.result4$fitted.obs.smoothed[,i]
    observed_variable      <- observables[,i]
  } 
  
  par(plt=c(.1,.95,.15,.95))
  
  y.lim <- c(min(model_implied_variable, observed_variable,na.rm=TRUE),
             max(model_implied_variable, observed_variable,na.rm=TRUE))
  
  rge <- y.lim[2] - y.lim[1]
  y.lim[2] <- y.lim[2] + .3*rge
  
  plot(vec.dates$date, model_implied_variable,type='l',col='darkgrey', lwd=3, ylab="", xlab="", main="",
       ylim=y.lim,las=1)
  grid()
  
  if(grepl("infl", colnames(observables)[i]) | grepl("growth", colnames(observables)[i])){
    points(vec.dates$date, observed_variable, col="black", lwd=1, pch=19, cex = 0.5)
  } else{
    points(vec.dates$date, observed_variable, col="black", lwd=1, pch=19)
  }
    
  # Legend condition
  legend("topleft",legend=titletopleft, cex=1, bty = "n")
  
  legend("topright",
         c("Model","Data"), 
         lty=c(1,NaN),pch=c(NaN,19),
         col=c("darkgrey","black"),lwd=2,
         cex=0.86,bg="white")
  
  dev.off()
}
