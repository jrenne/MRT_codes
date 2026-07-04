# ==============================================================================
#  PLOT IMPLIED AND ORIGINAL DITRIBUTIONS 
# ==============================================================================


# Plot the distribution 

if(area=="US"){
  set.of.dates.gdp <- c("1987-10-15",
                        "2007-10-15",
                        "2021-10-15")
  
  set.of.dates.infl <- c("1987-10-15",
                         "2007-10-15",
                         "2021-10-15")
  
  set.of.dates.gdp <- c("1985-10-15",
                        "2007-10-15",
                        "2017-10-15")
  
  set.of.dates.infl <- c("1985-10-15",
                         "2007-10-15",
                         "2017-10-15")
}


if(area=="EA"){
  set.of.dates.gdp <- c("1999-11-15",
                        "2014-11-15",
                        "2021-10-15")
  
  set.of.dates.infl <- c("1999-11-15",
                         "2014-11-15",
                         "2022-10-15")
}

# set.of.dates.gdp <- c("2020-10-15",
#                       "2021-10-15",
#                       "2022-10-15")
# set.of.dates.infl <- c("2007-10-15",
#                        "2008-10-15",
#                        "2009-10-15")
#set.of.dates.infl <- set.of.dates.gdp

# set.of.dates.gdp <- c("2015-10-15",
#                       "2016-10-15",
#                       "2017-10-15")
# set.of.dates.infl <- set.of.dates.gdp


horizon <- 5

for(type.var.implied in c("gdp","infl")){
  
  if(type.var.implied == "gdp"){
    first.year <- format(as.Date(set.of.dates.gdp[1]),"%Y")
  }else{
    first.year <- format(as.Date(set.of.dates.infl[1]),"%Y")
  }
  
  for(display.implied in c(TRUE,FALSE)){
    
    path <- paste(path_graph,"figure_distributions_",
                  ifelse(!display.implied,"noModeled_",""),
                  type.var.implied,"_",
                  first.year,
                  ".pdf",sep="")
    
    if(type.var.implied=="gdp"){
      set.of.dates <- set.of.dates.gdp
    }else{
      set.of.dates <- set.of.dates.infl
    }
    
    pdf(path, width = 9, height = 6)
    par(mfrow=c(2,2))
    par(plt=c(.15,.95,.25,.95))
    make.model.implied.plot(set.of.dates[1],type.var.implied,horizon,
                            display.legend = FALSE,
                            display.implied = display.implied)
    make.model.implied.plot(set.of.dates[2],type.var.implied,horizon,
                            display.legend = FALSE,
                            display.implied = display.implied)
    make.model.implied.plot(set.of.dates[3],type.var.implied,horizon,
                            display.legend = FALSE,
                            display.implied = display.implied)
    
    plot.new()
    
    if(display.implied){
      legend("topleft",
             lty=c(NaN,3,1),
             title="Distributions:",
             c("Survey (observed)", 'Survey ("mixture-smoothed")', 'Modeled') ,
             col = c("grey","black", "red"),
             pt.cex=1.5,
             seg.len=3,
             lwd=2,
             pch=c(15,NaN,NaN),
             cex=1.1)
    }else{
      legend("topleft",
             lty=c(NaN,3,1),
             title="Distributions:",
             c("Survey (observed)", 'Survey ("mixture-smoothed")') ,
             col = c("grey","black"),
             pt.cex=1.5,
             seg.len=3,
             lwd=2,
             pch=c(15,NaN),
             cex=1.1)
    }
    
    dev.off()
    
    # path <- paste0("graphs/US_2024/figure_distributions.pdf")
    # pdf(path, width = 7, height = 4)
    # par(mfrow=c(1,1))
    # par(plt=c(.15,.95,.2,.95))
    # make.model.implied.plot("2020-10-15","gdp",5,TRUE,TRUE,display.legend = F)
    # dev.off()
    
  }

}

