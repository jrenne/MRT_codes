# ==============================================================================
#       PLOTS DATA
# ==============================================================================

path <- paste0("graphs/US_2024/figure_data.pdf")
pdf(path, width = 11, height = 7)
par(plt=c(.1,.95,.15,.8))

Horizons <- 5:8
pch <- c(3,3,5,5)
col <- c("black","dark grey","black","dark grey")

Variables <- c("beta.pe","beta.var","beta.k3rd","beta.k4th")

par(mfrow=c(2,4))

for(macro in c("infl","gdp")){
  
  count_variable <- 0
  for(variable in Variables){
    count_variable <-   count_variable + 1
    
    count_horizon <- 0
    
    main.t = paste("(",ifelse(macro=="infl","a","b"),".",count_variable,
                   ") ",ifelse(macro=="infl","Inflation","GDP growth"),
                   " - Cumulant of order ",count_variable,sep="")
    
    for(horizon in Horizons){
      count_horizon <-   count_horizon + 1
      eval(parse(text = gsub(" ","",paste("observations <- observables.with.dates$SPF.US.",
                                          ifelse(macro=="infl","","G."),
                                          horizon,"Q.",variable, sep="")))) 
      if(count_horizon == 1){
        range <- range(observations,na.rm=TRUE)
        ylim =  c(range[1] - .1 * (range[2] - range[1]),
                  range[2] + .2 * (range[2] - range[1]))
        plot(as.Date(observables.with.dates$date),
             observations,
             pch=pch[count_horizon],
             col = col[count_horizon],
             lwd=1,xlab = "",ylab="",las=1,
             main=main.t,
             ylim=ylim)
        grid()
        
        if((variable == Variables[1])&(macro == "infl")){
          legend("topright",
                 paste(Horizons,"Q",sep=""), 
                 lty=NaN,
                 pch=c(3,3,5,5),
                 col=c("black","darkgrey"),
                 lwd=1,
                 cex=1,bg="white")
        }
      }else{
        points(as.Date(observables.with.dates$date),
               observations,
               pch=pch[count_horizon],
               col = col[count_horizon],
               lwd=1)
      }
    }
  }
}

dev.off()

