# ==============================================================================
#  PLOT IMPLIED AND ORIGINAL DITRIBUTIONS 
# ==============================================================================


# Plot the distribution 

if(area=="US"){
  distribution_figure_sets <- list(
    list(
      first.year = "1985",
      display.implied = c(TRUE, FALSE),
      gdp = c("1985-10-15", "2007-10-15", "2017-10-15"),
      infl = c("1985-10-15", "2007-10-15", "2017-10-15")
    ),
    list(
      first.year = "2015",
      display.implied = FALSE,
      gdp = c("2015-10-15", "2016-10-15", "2017-10-15"),
      infl = c("2015-10-15", "2016-10-15", "2017-10-15")
    ),
    list(
      first.year = "2020",
      display.implied = FALSE,
      gdp = c("2020-10-15", "2021-10-15", "2022-10-15"),
      infl = c("2020-10-15", "2021-10-15", "2022-10-15")
    )
  )
}


if(area=="EA"){
  distribution_figure_sets <- list(
    list(
      first.year = "1999",
      display.implied = c(TRUE, FALSE),
      gdp = c("1999-11-15", "2014-11-15", "2021-10-15"),
      infl = c("1999-11-15", "2014-11-15", "2022-10-15")
    )
  )
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

plot_distribution_set <- function(type.var.implied, set.of.dates, first.year, display.implied) {
  path <- paste(path_graph,"figure_distributions_",
                ifelse(!display.implied,"noModeled_",""),
                type.var.implied,"_",
                first.year,
                ".pdf",sep="")
  
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
           col = c("grey","black", "grey35"),
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
}

for(figure_set in distribution_figure_sets){
  for(type.var.implied in c("gdp","infl")){
    set.of.dates <- figure_set[[type.var.implied]]
    for(display.implied in figure_set$display.implied){
      plot_distribution_set(
        type.var.implied,
        set.of.dates,
        figure_set$first.year,
        display.implied
      )
    }
  }
}
