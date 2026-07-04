# ==============================================================================
#          MODEL IMPLIED DISTRIBUTIONS 
# --------------------------------------------
#                 * ALL *
# ==============================================================================

# ==============================================================================
# This program derives different analytical 
# distributions and, if requested, std dev associated with them
# ==============================================================================

#indic.make.implied.distri <- "FALSE"
if (!exists("show_implied_distribution_progress")) show_implied_distribution_progress <- TRUE

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
if (show_implied_distribution_progress) {
  message("MRT replication: Kalman log-likelihood for implied distributions: ", round(KF.res$log.lik, 3))
}

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
  if (show_implied_distribution_progress) {
    message("MRT replication: recomputing model-implied distributions.")
  }
  implied_distribution_start_time <- Sys.time()
  
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
  
  # ============================================================================
  
  ### Loop on the indicator of average... 0 by default in the AB function
  total_implied_distribution_steps <- length(2:1) * estimated.Model$r * length(HHH)
  implied_distribution_step <- 0
  
  for(indic.average in 2:1){
    # ==========================================================================
    
    for(area.var in 1:estimated.Model$r){
      b <- c(estimated.Model$delta[,area.var],rep(0,estimated.Model$q))
      # ========================================================================
      
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
        implied_distribution_step <- implied_distribution_step + 1
        if (show_implied_distribution_progress) {
          message(
            "MRT replication: implied distribution step ",
            implied_distribution_step, "/", total_implied_distribution_steps,
            " (average=", indic.average,
            ", variable=", area.var,
            ", horizon=", HH, "Q)."
          )
        }
        
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
  
  implied.distribution.data <- list("PDF.all"=PDF.all,"PDF.average.all"=PDF.average.all,"PDF.x.all"=PDF.x.all)
  if (exists("path_implied_distribution")) {
    dir.create(dirname(path_implied_distribution), recursive = TRUE, showWarnings = FALSE)
    if (show_implied_distribution_progress) {
      message("MRT replication: saving model-implied distributions to ", path_implied_distribution)
    }
    saveRDS(implied.distribution.data, file = path_implied_distribution, compress = "xz")
  }
  if (show_implied_distribution_progress) {
    message(
      "MRT replication: model-implied distributions completed in ",
      round(as.numeric(difftime(Sys.time(), implied_distribution_start_time, units = "secs")), 1),
      " seconds."
    )
  }
  
} else{
  if (!exists("path_implied_distribution")) {
    path_implied_distribution <- "results/US/US.Model.Implied.Distribution.RData"
  }
  if (show_implied_distribution_progress) {
    message("MRT replication: loading model-implied distributions from ", path_implied_distribution)
  }
  implied_distribution_start_time <- Sys.time()
  implied.distribution.data <- readRDS(file = path_implied_distribution, refhook = NULL)
  if (show_implied_distribution_progress) {
    message(
      "MRT replication: model-implied distributions loaded in ",
      round(as.numeric(difftime(Sys.time(), implied_distribution_start_time, units = "secs")), 1),
      " seconds."
    )
  }
  
  PDF.all <- implied.distribution.data$PDF.all
  PDF.average.all <- implied.distribution.data$PDF.average.all
  PDF.x.all <- implied.distribution.data$PDF.x.all
}
