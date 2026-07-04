# ==============================================================================
#          MODEL IMPLIED DISTRIBUTIONS 
# --------------------------------------------
#                 * ALL *
# ==============================================================================

# ==============================================================================
# This program derives different analytical 
# distributions and, if requested, std dev associated with them
# ==============================================================================

indic.make.implied.distri <- "FALSE"

par(plt=c(.1,.9,.1,.9))

# HHH is a vector of considered horizons, expressed in quarters:
HHH <- seq(4,8,by=1)
# Select.if.compute.stdv is a matrix, whose rows correspond
# to the horizons in HHH, the number of columns corresponds to the 
# different countries. contains 1 if one wants to compute stdv of PDF/CDF (0 otherwise)


# ==============================================================================
# Specify model:
n <- n # number of Y factors
q <- q # number of z factors
# Select Area:
#Area <- area


Model <- Model.final
estimated.Model <- Model
observables <- observables
select.inflation.types <- select.inflation.types
H <- H
KF.res <- KF.result4
print(KF.res$log.lik)

Select.if.compute.stdv <- matrix(0,length(HHH),estimated.Model$r)

# Compute smoothed results:
KF.res.0 <- KF.result4

# Get smoothed factors:
X <- KF.res.0$xi.tT

Indic.5y.in.Xy <- 0


# Initialize the condition for the y at which the distribution will be evaluated
# FYI:   vec.y <- seq(min.bx,max.bx,by=step.distri), where 
## min.bx = center.distri - estimated.Model$pi.bar[area] - width.distri/2
## max.bx = center.distri - estimated.Model$pi.bar[area] + width.distri/2
step.distri <- .05
center.distri <- 2
width.distri <- 30 

# ==============================================================================

if(indic.make.implied.distri=="TRUE"){
  
  # Initialize the max for v and its increment. (Riemann sum)
  # the start is exp(-10)
  # vector is seq(-10, log(50), by=0.02) => take exponential of that
  step.4.integral <- .02
  max.v.4.integral <- 50
  
  # Number of considered abscissa, i.e nbr of y at which the PDF is evaluated
  nb.abscissa <- width.distri/step.distri
  
  # PDF (and CDF) results will be in array of dimension:
  #      T x (nb points where pdf evaluated) x (nb horizons) X (nb countries)
  #### WARNINGS : nb.abscissa-1 before
  PDF.all <- array(NaN,c(dim(X)[1],nb.abscissa,length(HHH),estimated.Model$r))
  CDF.all <- array(NaN,c(dim(X)[1],nb.abscissa,length(HHH),estimated.Model$r))
  
  PDF.stdv.all <- array(NaN,c(dim(X)[1],nb.abscissa,length(HHH),estimated.Model$r))
  CDF.stdv.all <- array(NaN,c(dim(X)[1],nb.abscissa,length(HHH),estimated.Model$r))
  
  PDF.average.all <- array(NaN,c(dim(X)[1],nb.abscissa,length(HHH),estimated.Model$r))
  CDF.average.all <- array(NaN,c(dim(X)[1],nb.abscissa,length(HHH),estimated.Model$r))
  
  PDF.average.stdv.all <- array(NaN,c(dim(X)[1],nb.abscissa,length(HHH),estimated.Model$r))
  CDF.average.stdv.all <- array(NaN,c(dim(X)[1],nb.abscissa,length(HHH),estimated.Model$r))
  
  PDF.x.all <- NULL
  
  # ==============================================================================
  
  ### Loop on the indicator of average... 0 by default in the AB function
  
  for(indic.average in 2:1){
    # ============================================================================
    
    for(area.var in 1:estimated.Model$r){
      b <- c(estimated.Model$delta[,area.var],rep(0,estimated.Model$q))
      # ==========================================================================
      
      # Parameterzing v:
      # v is the point at which the CDF is evaluated
      
      ## Max value for v
      max.x <- log(max.v.4.integral)
      
      ## Vectors with all the considered v => points at which the CDF is evaluated!!
      ## Here we use exponential from -10 to 4
      ## Hence v goes from 0 to 49.
      v <- matrix(exp(seq(-10,max.x,by=step.4.integral)),nrow=1)
      
      ## Find the arguments to compute A_h and B_h functions
      u_Y <- 1i * matrix(b[1:estimated.Model$n],ncol=1) %*% v  # our paper, section 7.3 Appendix i*v*gamma for Y
      u_z <- 1i * matrix(b[(estimated.Model$n+1):(estimated.Model$n+estimated.Model$q)],ncol=1) %*% v  # our paper, section 7.3 Appendix i*v*gamma for z
      AB.list <- compute_AB(estimated.Model,u_Y,u_z,HHH,indic.average)
      
      ## Compute the min and the max value of the distribution
      min.bx <- center.distri - estimated.Model$pi.bar[area.var]*4 - width.distri/2
      max.bx <- center.distri - estimated.Model$pi.bar[area.var]*4 + width.distri/2
      
      ## Begin the loop to compute the PDF
      count <- 0
      for(HH in HHH){
        
        ### Increment count in the loop
        count <- count + 1
        
        ### Condition to extract P.tT if needed. 
        if(Select.if.compute.stdv[count,area.var]==1){# want to compute stdev of distributions
          Sigma.X <- KF.res.0$P.tT
        }else{
          Sigma.X <- 0
        }
        
        ### Find the position in the vector HHH
        indic.matur <- which(HH==HHH)
        
        ### Extract the AB for the considered horizon and put it in an array of final dimension 1 => should change it maybe
        AB.list.external <- list(A=array(AB.list$A[,,indic.matur],c(dim(AB.list$A)[1],dim(AB.list$A)[2],1)),
                                 B=array(AB.list$B[,,indic.matur],c(dim(AB.list$B)[1],dim(AB.list$B)[2],1)))
        
        res.analytical.computation <-
          compute.distri.plus.stdv(estimated.Model,X,b,HH,
                                   min.bx,max.bx,step.distri,
                                   step.4.integral,
                                   max.v.4.integral,
                                   indic.average,
                                   Sigma.X,
                                   AB.list.external)
        
        if(indic.average==2){
          PDF.all[,,count,area.var] <- res.analytical.computation$PDF
          CDF.all[,,count,area.var] <- res.analytical.computation$CDF
          if(Select.if.compute.stdv[count,area.var]==1){# want to compute stdev of distributions
            PDF.stdv.all[,,count,area.var] <- res.analytical.computation$PDF.stdv
            CDF.stdv.all[,,count,area.var] <- res.analytical.computation$CDF.stdv
          }
        }
        if(indic.average==1){
          PDF.average.all[,,count,area.var] <- res.analytical.computation$PDF
          CDF.average.all[,,count,area.var] <- res.analytical.computation$CDF
          if(Select.if.compute.stdv[count,area.var]==1){# want to compute stdev of distributions
            PDF.average.stdv.all[,,count,area.var] <- res.analytical.computation$PDF.stdv
            CDF.average.stdv.all[,,count,area.var] <- res.analytical.computation$CDF.stdv
          }
        }
      }
      PDF.x.all <- rbind(PDF.x.all,res.analytical.computation$PDF.x)
    }
  }
  
  # implied.distribution.data <- list("PDF.all"=PDF.all,"PDF.average.all"=PDF.average.all,"PDF.x.all"=PDF.x.all)
  # saveRDS(implied.distribution.data, file="./results/US/US.Model.Implied.Distribution.RData")
  # ?saveRDS
  
} else{
  implied.distribution.data <- readRDS(file="results/US/US.Model.Implied.Distribution.RData", refhook = NULL)
  
  PDF.all <- implied.distribution.data$PDF.all
  PDF.average.all <- implied.distribution.data$PDF.average.all
  PDF.x.all <- implied.distribution.data$PDF.x.all
}



#NOTE:
## PDF.x,all is the matrix that contains the values for the x-axis of the integral.
## lines 1 and 2 are for area.var 1 and 2 and indic.average=0.
## PDF.average.all is the array that contains the values for the y-axis of the integral.
## 1st dimension is time, 2nd is nbr of point at which the PDF is evaluated, 3rd is horizons considered length(HHH)
## 4th is the area.var; PDF.all[time,,horizon,area.var]

# ============================================
#   B. PLOT IMPLIED AND ORIGINAL DITRIBUTIONS 
# ============================================

# FOR INFLATION
# Plot the distribution for the last 6 observations.



for(horizon in 4:8){
area.var <- 1
#horizon <- 8  
#horizon <- 8
horizon.y <- horizon-(freq-1)
#period <- "1981_1985"
#period <- "1985_1992"
period <- "1992_2014"
#period <- "post_2014"

## Parameters in order to diplay well the chart
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)

layout(mat = m, heights = c(0.87/3,0.87/3,0.87/3,0.13))

#layout.show(7)
par(mar=c(4, 4.5, 1, 1))

eval(parse(text = gsub(" ","", paste('survey.bins.distri <- US.SPF.DISTRI.',
                                     paste(horizon,sep = ""), 'Q.', paste(period,sep = ""), ' %>% left_join(survey.DATA.US.with.param %>% dplyr::select(date, contains("',
                                     paste(horizon,sep = ""), 'Q" )), by="date")',
                                     sep="")))) 

eval(parse(text = gsub(" ","", paste('nbr.bins <- dim(US.SPF.DISTRI.',
                                     paste(horizon,sep = ""), 'Q.', paste(period,sep = ""),')[2]', sep="")))) 


survey.bins.distri.without.na <- na.omit(survey.bins.distri) 
max.m <- dim(survey.bins.distri.without.na)[1]
seq.m <- seq(max.m-min(8,max.m-1),max.m, by=1)

for (m in seq.m) {
  
  all.data.new <- survey.bins.distri.without.na[m,]
  Date <- as.Date(all.data.new[,1])
  data1 <- all.data.new[,2:nbr.bins]*100
  data1 <- data1[length(data1):1] # reverse order
  time <- which(Date==as.matrix(vec.dates))
  
  param1 <- as.numeric(all.data.new[,grep("param", colnames(as.data.frame(all.data.new)))])
  #max.c1 <- as.numeric(all.data.new[,grep("max.c", colnames(as.data.frame(all.data.new)))])
  #min.d1 <- as.numeric(all.data.new[,grep("min.d", colnames(as.data.frame(all.data.new)))])
  max.bin <- as.numeric(all.data.new[,grep("max.bin", colnames(as.data.frame(all.data.new)))])
  min.bin <- as.numeric(all.data.new[,grep("min.bin", colnames(as.data.frame(all.data.new)))])
  x1 <- seq(min.bin-0.5,max.bin+0.5,by=0.05)
  
  q <- format(Date, "%m")
  t <- format(Date, "%Y")
  
  lower.bound.interval <- round(sort(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",colnames(data1))))))*2)/2
  lower.bound.interval[1] <- min.bin
  upper.bound.interval <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(data1)))))*2)/2 
  upper.bound.interval[length(upper.bound.interval)] <- max.bin
  all.intervals <- upper.bound.interval - lower.bound.interval
  
  max.bin.fit <- max(upper.bound.interval)
  min.bin.fit <- min(lower.bound.interval)
  
  cdf.values.new <- c(0,cumsum(c(data1)))
  x.new <- c(lower.bound.interval[1], upper.bound.interval)      
  min.sigma <- max(min(diff(x.new))/2,0.25)
  
  max.bin.fit <- x.new[min(which(cdf.values.new == max(cdf.values.new)))]
  min.bin.fit <- x.new[max(which(cdf.values.new == min(cdf.values.new)))]
  
  
  if(substr(colnames(all.data.new)[nbr.bins+1], 1, 3)=="PDS"){# PDS
    mean.class <- c(0.75,seq(1.25, 2.75, 0.5),3.5)
    break2 <- sort(c(c(0.5,4), seq(1, 3, 0.5), seq(0.95, 2.95, by =0.5)))
    
  }else if(substr(colnames(all.data.new)[nbr.bins+1], 1, 3)=="SPF"){# US SPF
    
    mean.class <- lower.bound.interval + all.intervals/2
    break2 <- sort(c(upper.bound.interval, lower.bound.interval, lower.bound.interval[-order(lower.bound.interval)[1]] - 0.0))
    break2 <- break2[!duplicated(break2)]
  }
  
  plot.fit.survey.distribution.mixture(data1,x1,param1,min.bin=min.bin.fit,max.bin=max.bin.fit,min.sigma=min.sigma, mean.class, break2)
  lines(PDF.x.all[area.var,], PDF.all[time,,horizon.y,area.var], col="red", lwd=1.5)
  legend("topleft",legend=paste(t,":",q, sep=""), cex=0.6)
  legend("topright",legend=bquote(E[t]*"("*pi[t*","*t*"+"*.(horizon)]*")"), cex=0.6, bty = "n")
  
}

par(mar=c(0.5, 5, 0.5, 1))
plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",inset = 0, title="Variables:", lty=c('blank','solid','solid'),
       c("Survey (observed)", 'Survey ("Mixture" smoothed)', 'Modeled') ,
       col = c("grey","black", "red"), fill=c("grey", "white", "white"), border = c("grey", "white", "white"),
       cex=0.9, horiz = TRUE)


# FOR GDP
# Plot the distribution for the last 6 observations.

area.var <- 2
horizon.y <- horizon-(freq-1)
#period <- "1981_1992"
period <- "1992_2009"
#period <- "2009_2020"
#period <- "post_2020"

## Parameters in order to diplay well the chart
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)

layout(mat = m, heights = c(0.87/3,0.87/3,0.87/3,0.13))

#layout.show(7)
par(mar=c(4, 4.5, 1, 1))

eval(parse(text = gsub(" ","", paste('survey.bins.distri <- US.G.SPF.DISTRI.',
                                     paste(horizon,sep = ""), 'Q.', paste(period,sep = ""), ' %>% left_join(survey.DATA.US.G.with.param %>% dplyr::select(date, contains(".',
                                     paste(horizon,sep = ""), 'Q" )), by="date")',
                                     sep="")))) 

eval(parse(text = gsub(" ","", paste('nbr.bins <- dim(US.G.SPF.DISTRI.',
                                     paste(horizon,sep = ""), 'Q.', paste(period,sep = ""),')[2]', sep="")))) 

survey.bins.distri.without.na <- na.omit(survey.bins.distri) 
max.m <- dim(survey.bins.distri.without.na)[1]
seq.m <- seq(max.m-min(8,max.m-1),max.m, by=1)

## Parameters in order to diplay well the chart
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)

layout(mat = m, heights = c(0.87/3,0.87/3,0.87/3,0.13))


if(max.m < 6){
  
  if(max.m < 4){
    m <- matrix(c(seq(1,max.m,1),rep(max.m+1,min(length(c(seq(1,max.m,1))),3))),
                nrow = (ceiling(max.m/3)+1),ncol = min(length(c(seq(1,max.m,1))),3),byrow = TRUE)
    
    layout(mat = m, heights = c(0.42,0.58))
    
  }else if(max.m == 4){m <- matrix(c(1,2,3,4,5,5), ncol=2,nrow = 3,byrow = TRUE)
  layout(mat = m, heights = c(0.84/2,0.84/2,0.16))
  }else{ m <- matrix(c(1,2,3,4,5,5,6,6,6), ncol=3,nrow = 3,byrow = TRUE )
  layout(mat = m, heights = c(0.84/2,0.84/2,0.16))}
  
  
}

for (m in seq.m) {
  
  
  all.data.new <- survey.bins.distri.without.na[m,]
  Date <- as.Date(all.data.new[,1])
  data1 <- all.data.new[,2:nbr.bins]*100
  data1 <- data1[length(data1):1] # reverse order
  time <- which(Date==as.matrix(vec.dates))
  
  param1 <- as.numeric(all.data.new[,grep("param", colnames(as.data.frame(all.data.new)))])
  #max.c1 <- as.numeric(all.data.new[,grep("max.c", colnames(as.data.frame(all.data.new)))])
  #min.d1 <- as.numeric(all.data.new[,grep("min.d", colnames(as.data.frame(all.data.new)))])
  max.bin <- as.numeric(all.data.new[,grep("max.bin", colnames(as.data.frame(all.data.new)))])
  min.bin <- as.numeric(all.data.new[,grep("min.bin", colnames(as.data.frame(all.data.new)))])
  x1 <- seq(min.bin-0.5,max.bin+0.5,by=0.05)
  
  q <- format(Date, "%m")
  t <- format(Date, "%Y")
  
  lower.bound.interval <- round(sort(as.numeric(sub(".*?([-+]?\\d*\\.?\\d+).*", "\\1", gsub("_",".",gsub("N","-",colnames(data1))))))*2)/2
  lower.bound.interval[1] <- min.bin
  upper.bound.interval <- round(as.numeric(gsub(".*?(-?\\d+(?:\\.\\d+)?)[^0-9]*$", "\\1", gsub("_",".",gsub("N","-",colnames(data1)))))*2)/2 
  upper.bound.interval[length(upper.bound.interval)] <- max.bin
  all.intervals <- upper.bound.interval - lower.bound.interval
  
  max.bin.fit <- max(upper.bound.interval)
  min.bin.fit <- min(lower.bound.interval)
  
  cdf.values.new <- c(0,cumsum(c(data1)))
  x.new <- c(lower.bound.interval[1], upper.bound.interval)      
  min.sigma <- max(min(diff(x.new))/2,0.25)
  
  max.bin.fit <- x.new[min(which(cdf.values.new == max(cdf.values.new)))]
  min.bin.fit <- x.new[max(which(cdf.values.new == min(cdf.values.new)))]
  
  
  if(substr(colnames(all.data.new)[nbr.bins+1], 1, 3)=="PDS"){# PDS
    mean.class <- c(0.75,seq(1.25, 2.75, 0.5),3.5)
    break2 <- sort(c(c(0.5,4), seq(1, 3, 0.5), seq(0.95, 2.95, by =0.5)))
    
  }else if(substr(colnames(all.data.new)[nbr.bins+1], 1, 3)=="SPF"){# US SPF
    
    mean.class <- lower.bound.interval + all.intervals/2
    break2 <- sort(c(upper.bound.interval, lower.bound.interval, lower.bound.interval[-order(lower.bound.interval)[1]] - 0.0))
    break2 <- break2[!duplicated(break2)]
  }
  
  plot.fit.survey.distribution.mixture(data1,x1,param1,min.bin=min.bin.fit,max.bin=max.bin.fit,min.sigma=min.sigma, mean.class, break2)
  lines(PDF.x.all[area.var,], PDF.all[time,,horizon.y,area.var], col="red", lwd=1.5)
  legend("topleft",legend=paste(t,":",q, sep=""), cex=0.6)
  legend("topright",legend=bquote(E[t]*"("*Delta*y[t*","*t*"+"*.(horizon)]*")"), cex=0.6, bty = "n")
  
}

par(mar=c(0.5, 5, 0.5, 1))
plot(1,1, type = "n", axes=FALSE, xlab="", ylab="")
legend("bottom",inset = 0, title="Variables:", lty=c('blank','solid'),
       c("Survey (observed)", 'Survey ("Mixture" smoothed)') ,
       col = c("grey","black"), fill=c("grey", "white"), border = c("grey", "white"),
       cex=0.9, horiz = TRUE)
}
