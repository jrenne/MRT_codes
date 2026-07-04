# ==============================================================================
#   A. PREPARE THE MODEL
# ==============================================================================
indic.cycle     <- TRUE
indic.cycle.use <- FALSE
indic.3rd       <- TRUE
indic.4th       <- TRUE
indic.4th.use   <- FALSE 
if (!exists("indic.estimate")) indic.estimate <- FALSE
if (!exists("estimation_start")) estimation_start <- "saved"
if (!exists("save_estimation_results")) save_estimation_results <- indic.estimate
indic.no.4Q <- TRUE


estimate.var.only= FALSE

all_mat <- NULL
if(indic.4th){
  all_mat <- readRDS(file="data/processed/all_mat_sparse_4_5_new.RData", refhook = NULL)
}

# ==============================================================================
# Specify the data to be use

if(!indic.compute.mixture){
  
  observables.with.dates <- read_excel("data/processed/Output.EA.xlsx")
  
  # Read data from Sheet1
  survey.DATA.new.with.param <- read_excel("data/processed/survey.DATA.EA.with.param.xlsx", sheet = "survey.DATA.EA.with.param")
  
  # Read data from Sheet2
  survey.DATA.G.new.with.param <- read_excel("data/processed/survey.DATA.EA.with.param.xlsx", sheet = "survey.DATA.EA.G.with.param")
  
} else{
  
  # Specify the data to be use
  observables.with.dates <- observables.with.dates.EA %>%
    full_join(observables.with.dates.EA.G, by = "date") %>%
    filter(date >= as.Date("1999-01-15") & date <= as.Date("2024-01-15")) %>%
    filter(
      (date < as.Date("2015-01-01") & substr(date, 6, 7) %in% c("02", "05", "08", "11")) |
        (date >= as.Date("2015-01-01") & substr(date, 6, 7) %in% c("01", "04", "07", "10"))
    )
  
  observables.with.dates$date <- as.POSIXct(observables.with.dates$date)
  
  
  if(indic.cycle){
    observables.with.dates <- observables.with.dates  %>%
      left_join(DATA %>% dplyr::select(date,EA.log.hcpi.deseasonalized.cycle), by="date") %>%
      left_join(DATA.G %>%  
                  filter(date != as.Date("2014-12-15")) %>%
                  mutate(date = if_else(year(date) < 2015, 
                                        date + months(1),  # Subtract one month
                                        date)) %>% 
                  dplyr::select(date,EA.log.gdp.deseasonalized.cycle), by="date") 
  }
  
}

# Store initial observables
observables.with.dates.all <- observables.with.dates


#### remove when PDS data available
H.EA <- matrix(c(12, 15, 18, 21, 24),5,1) # number of columns = number of areas
H.EA.G <- matrix(c(12, 15, 18, 21, 24),5,1)

#H.EA <- H.EA[1:min(5,length(c(H.EA)))]
select.inflation.types.EA <- array(0,c(length(c(H.EA)),1,4))
select.inflation.types.EA[1:min(5,length(c(H.EA))),1,3] <- 11 

if(indic.3rd){
  select.inflation.types.EA[1:min(5,length(c(H.EA))),1,3] <- 111  
}

if(indic.4th){
  select.inflation.types.EA[1:min(5,length(c(H.EA))),1,3] <- 1111  
}

#H.EA.G <- H.EA.G[1:min(5,length(c(H.EA.G)))]
select.inflation.types.EA.G <- array(0,c(length(c(H.EA.G)),1,4))
select.inflation.types.EA.G[1:length(c(H.EA.G)),1,3] <- 11 

if(indic.3rd){
  select.inflation.types.EA.G[1:length(c(H.EA.G)),1,3] <- 111 
}

if(indic.4th){
  select.inflation.types.EA.G[1:length(c(H.EA.G)),1,3] <- 1111  
}

observables <- observables.with.dates %>% dplyr::select(-date)

# ================================
# Specify the data to be use

if(!indic.3rd){
  observables.with.dates <- observables.with.dates %>% 
    dplyr::select(-contains(c("k3rd")))
  
}

if(!indic.4th){
  observables.with.dates <- observables.with.dates %>% 
    dplyr::select(-contains(c("k4th")))
}

if(indic.4th){
  if(!indic.4th.use){
    columns_to_replace <- c(".4th")
    observables.with.dates[, grepl(paste(columns_to_replace, collapse = "|"), colnames(observables.with.dates))] <- NA
  }
}

# Extract the frequency
freq.count.all <- observables.with.dates %>%
  mutate(year =format(as.Date(date, format="%d/%m/%Y"),"%Y")) %>%
  group_by(year) %>% count
freq <- max(freq.count.all$n)

# ================================
# Specify Matrix of horizons:
max.H.nb <- max(length(H.EA),length(H.EA.G))
H <- matrix(NaN,max.H.nb,2)
H[1:length(H.EA),1] <- H.EA
H[1:length(H.EA.G),2] <- H.EA.G

# Specify the horizon of interest taking into account the frequency
H <- H/(12/freq)

# Specify select.inflation.types array
select.inflation.types <- array(NaN,c(max.H.nb,2,4))
select.inflation.types[1:length(H.EA),1,] <- select.inflation.types.EA
select.inflation.types[1:length(H.EA.G),2,] <- select.inflation.types.EA.G

## Delete 4Q data
if(indic.no.4Q){
  columns_to_replace <- c("4Q.beta.pe","4Q.beta.var","4Q.beta.k3rd","4Q.beta.k4th")
  observables[, grepl(paste(columns_to_replace, collapse = "|"), colnames(observables))] <- NA
}

# ================================
# Specify the type of inflation for h > 60
## It means that for horizons higher than 5 years, consider 5y-in-Xy inflation 
## types (for case-1 inflation rates)
Indic.5y.in.Xy <- 1 
indic.smoother <- 1
indic.charts=1
indic.suffix=""

# ================================
# Specify model:
m <- 4 # number of Y factors
n <- 2*m + 2*(freq -1)
q <- 5 # number of z factors
r <- 2 # number of key observables
nbr.horizon.max <- max.H.nb #max number of horizon observed
# Select Area:
area <- "EA"
indic.observed <- "q.o.q"
var.type <- "infl.gdp"

## Calculate delta
delta.aux <- matrix(1,m,r)
aux <- apply(delta.aux,1,function(x){sum(x^2)}) #squared element of the column

## Estimate stdv for observables (delta)
stdv.measur <- make.stdv.measure(observables.with.dates,select.inflation.types,r)


if(!indic.cycle.use){
  observables[,dim(observables)[2] - c(1,0)] <- NA
  observables.with.dates[,dim(observables)[2] - c(1,0)+1] <- NA
}

if(indic.4th){
  if(!indic.4th.use){
    #columns_to_replace <- c(".3rd",".4th")
    columns_to_replace <- c(".4th")
    observables[, grepl(paste(columns_to_replace, collapse = "|"), colnames(observables))] <- NA
  }
}

# Convert in matrix form for simplicity to plot
observables <- as.matrix(observables)

if(indic.estimate){
  
  if (estimation_start == "saved") {
    Model.initial <- readRDS(file="results/EA/EA.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.no.4Q.RData", refhook = NULL)  
    
    model <- list(pi.bar = matrix(c(mean(observables.with.dates$EA.infl,na.rm=TRUE),
                                    mean(observables.with.dates$EA.growth,na.rm=TRUE)), 1,r),
                  delta.t = Model.initial$delta.t,
                  delta.c = Model.initial$delta.c,
                  Phi.Y.r = Model.initial$Phi.Y.r,
                  Theta = Model.initial$Theta,
                  Gamma.Y0.r = Model.initial$Gamma.Y0.r, #0
                  Gamma.Y1.r =  Model.initial$Gamma.Y1.r,
                  nu = Model.initial$nu, # 1st parameter of the non centered gamma process (AGP(nu,phi,mu))
                  phi = Model.initial$phi,
                  mu = Model.initial$mu, # 3rd parameter of the AGP
                  sigma.av = stdv.measur$sigma.av,
                  sigma.var = stdv.measur$sigma.var #, sigma.k3rd = stdv.measur$sigma.k3rd, sigma.k4th = stdv.measur$sigma.k4th
                  
    )
  } else if (estimation_start == "generic") {
    model <- list(pi.bar = matrix(c(mean(observables.with.dates$EA.infl,na.rm=TRUE),
                                    mean(observables.with.dates$EA.growth,na.rm=TRUE)), 1,r),
                  delta.t = matrix(0, m, r),
                  delta.c = matrix(0, m, r),
                  Phi.Y.r = diag(0.7, m),
                  Theta = matrix(0, n, q),
                  Gamma.Y0.r = matrix(0.05, m, 1),
                  Gamma.Y1.r = matrix(0, q, m),
                  nu = matrix(0.05, q, 1),
                  phi = diag(0.7, q),
                  mu = matrix(1, q, 1),
                  sigma.av = stdv.measur$sigma.av,
                  sigma.var = stdv.measur$sigma.var)
    model$delta.c[1, 1] <- -1
    model$delta.c[min(2, m), 2] <- 1
  } else {
    stop("estimation_start must be either 'saved' or 'generic'.")
  }
  
  
  model$sigma.var <- model$sigma.av
  model$sigma.av <- model$sigma.av/2
  
  if(indic.3rd){
    model$sigma.k3rd = stdv.measur$sigma.k3rd
  }
  
  if(indic.4th){
    model$sigma.k3rd = stdv.measur$sigma.k3rd
    model$sigma.k4th = stdv.measur$sigma.k4th
    model$sigma.k4th[,,3][,1] <- model$sigma.k4th[,,3][,1]*20
    model$sigma.k4th[,,3][,2] <- model$sigma.k4th[,,3][,2]*10
    
  }
  
  
  # Replace element theta[1,2] and theta[2,4] by 0
  if(q>1){
    model$Theta[1,2] <- 0}
  if(q>3){model$Theta[2,4] <- 0}
  
  
  # Create vector with estimates to consider
  estimated.Model <- Make.thetas.indicator.trend.cycle.model(model, delta.t.s=TRUE,  delta.c.s=TRUE, Phi.Y.r.s=TRUE, Gamma.Y0.r.s=TRUE,
                                                             Gamma.Y1.r.s=TRUE, nu.s=TRUE, phi.s=TRUE, Theta.s =TRUE)
  
  # Estimate only delta.t and delta.c for the second column (GDP growth)
  #step 2 and final
  # TO UNCOMMENT
  if(q==2 & m==1){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[(m+1):(2*m)],2] <- 0
  }
  if(q==1 & m==2){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[(m+1):(2*m)],2] <- 0
  }
  if(q==4 & m==2){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[c(1,(m+2):(2*m))],2] <- 0
  }
  if(q==3 & m==3){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[c((m+1):(2*m))],2] <- 0
  }
  if(q==5 & m==4){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[c(1,(m+2):(2*m))],2] <- 0
  }
  

  # TO UNCOMMENT
  if(q==1 & m==2){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[(1:2)],3] <- c(2,3)
  }
  if(q==4 & m==2){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[c(1:(m+1))],3] <- c(3,2,2)
  }
  if(q==3 & m==3){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[c(1:(m))],3] <- c(2,2,3)
  }
  if(q==5 & m==4){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[c(1:(m+1))],3] <- c(3,2,2,3,2)
  }
  
  
  nbr.to.excl <- 0
  estimated.Model[grep("Phi.Y.r", rownames(as.data.frame(estimated.Model)))[seq(1,m*m-nbr.to.excl*(m+1),by=(m+1))],2] <- 1
  estimated.Model[grep("Phi.Y.r", rownames(as.data.frame(estimated.Model)))[seq(1,m*m-nbr.to.excl*(m+1),by=(m+1))],3] <- 4
  
  if(estimate.var.only){
    estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[c(2,5)],2] <- 0
    estimated.Model[grep("Phi.Y.r", rownames(as.data.frame(estimated.Model)))[c(1,6)],2] <- 0
    estimated.Model[grep("Gamma.Y0.r", rownames(as.data.frame(estimated.Model)))[c(1,2)],2] <- 0
    estimated.Model[grep("nu", rownames(as.data.frame(estimated.Model)))[c(1:4)],2] <- 0
    estimated.Model[grep("phi", rownames(as.data.frame(estimated.Model)))[c(1,7,13,19)],2] <- 0
  }
  
  ### Define the vector of parameters to be estimated. 
  ### Check solution when initial guess consists of true simulation parameters.
  ### Initial parameters with no transformation  
  estimated.Model.initial <- estimated.Model[,1][ which(!estimated.Model[,2] == 0)] 
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  ### Use Inverse Mapping function to insure stationarity. This will transform
  ### some parameters for which we imposed restrictions.
  thetas <- Mapping.function.inverse(estimated.Model.initial, indicator.transformation)
  
  ### Refind initial parameters
  estimated.Model.list <- Retrieve.initial.par.trend.cycle.model(estimated.Model, n, m, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by -theta[1,1] and theta[2,4] by -theta[2,3]
  if(q>1){
    estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]}
  if(q>3){
    estimated.Model.list$Theta[2,4] <- - estimated.Model.list$Theta[2,3]}
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type, all_mat)
  
  ## Add other information to "estimated.Model"
  estimated.Model.list$areas <- area
  
  # Store all
  Model <- estimated.Model.list
  
  # Run the function that allows to create all the parameters of the KF for the 
  # defined model. 
  KF.load <- prepare.KF.model(Model, observables)
  
  KF.result1 <- Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X, 
                              KF.load$xi.00, KF.load$P.00, S="default", 
                              indic.pos.z=KF.load$indic_pos)
  if(indic.smoother==1){
    KF.result2 <- Kalman.smoother(KF.load$all.parameters, KF.load$Y, KF.load$X,
                                  KF.load$xi.00, KF.load$P.00, S="default", 
                                  indic.pos.z=KF.load$indic_pos)
  }

  
  # ============================================
  #   B. CALIBRATE THE MODEL
  # ============================================
  
  # Give results of the first guess
  fit.log.lik.trend.cycle.joint.model.3.1(thetas, estimated.Model, n, m, q, r, nbr.horizon.max, var.type, all_mat)
  
  #prepare for optimization 
  thetas.inter <- thetas
  thetas.inter.all <- thetas
  
  ###  Perform optimization
  for(j in 1:4){
    
    solution1 <- optim(fn = fit.log.lik.trend.cycle.joint.model.3.1, par = thetas.inter, estimated.Model=estimated.Model,
                       n=n, m.Y=m, q=q, r=r, nbr.horizon.max=nbr.horizon.max, var.type=var.type, all_mat=all_mat, method=c("Nelder-Mead"), 
                       control=list(trace=TRUE, maxit=1500), hessian = FALSE)
    
    #saveRDS(solution1, file="results/US/last_solution_optim.RData")
    #?saveRDS
    # USE THIS: solution1 <- readRDS(file="results/US/last_solution_optim.RData", refhook = NULL)
    
    thetas.inter <- solution1$par
    thetas.inter.all <- rbind(thetas.inter.all, thetas.inter)
    
    #note: remove also kkt=FALSE if wants to estimate hessian
    solution1 <- optimx(fn = fit.log.lik.trend.cycle.joint.model.3.1, par = thetas.inter, estimated.Model=estimated.Model,
                        n=n, m.Y=m, q=q, r=r, nbr.horizon.max=nbr.horizon.max, var.type=var.type, all_mat=all_mat, method=c("nlminb"), 
                        control=list(trace=TRUE, maxit=50, kkt=FALSE), hessian = FALSE)
    
    #saveRDS(solution1, file="results/US/last_solution_optim.RData")
    #?saveRDS
    # USE THIS: solution1 <- readRDS(file="results/US/last_solution_optim.RData", refhook = NULL)
    
    thetas.inter <- as.matrix(solution1)[,1:length(thetas)]
    thetas.inter.all <- rbind(thetas.inter.all, thetas.inter)
    
    solution1 <- optim(fn = fit.log.lik.trend.cycle.joint.model.3.1, par = thetas.inter, estimated.Model=estimated.Model,
                       n=n, m.Y=m, q=q, r=r, nbr.horizon.max=nbr.horizon.max, var.type=var.type, all_mat=all_mat, method=c("Nelder-Mead"),
                       control=list(trace=TRUE, maxit=1500), hessian = FALSE)
    
    #saveRDS(solution1, file="results/US/last_solution_optim.RData")
    #?saveRDS
    # USE THIS: solution1 <- readRDS(file="results/US/last_solution_optim.RData", refhook = NULL)
    
    thetas.inter <- solution1$par
    thetas.inter.all <- rbind(thetas.inter.all, thetas.inter)
    
  }
  
  #solution_par <- solution$par
  indic <- tail(which(rowSums(!is.na(thetas.inter.all)) > 0), 1)
  solution_par1 <- thetas.inter.all[indic,]
  
  # Recover the optimal parameters of interest
  coeff_KF_mc <- Mapping.function(solution_par1, indicator.transformation)
  
  # Create the vector of parameters with thetas
  all.thetas.final <- estimated.Model
  all.thetas.final[rownames(as.data.frame(coeff_KF_mc)),1] <- coeff_KF_mc
  
  all.thetas.final.list <- Retrieve.initial.par.trend.cycle.model(all.thetas.final, n, m, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by -theta[1,1] and theta[2,4] by -theta[2,3]
  if(q>1){
    all.thetas.final.list$Theta[1,2] <- - all.thetas.final.list$Theta[1,1]}
  if(q>3){
    all.thetas.final.list$Theta[2,4] <- - all.thetas.final.list$Theta[2,3]}
  
  ## Compute the parameters of the Model
  all.thetas.final.list <- make.parameters.model(all.thetas.final.list, var.type, all_mat)
  
  Model.final <- all.thetas.final.list                
  
} else{
  
  if(!indic.model.var.only){
    if(indic.4th.use){ #all moments
      Model.final <- readRDS(file="results/EA/EA.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.no.4Q.mean.higher.errors.RData", refhook = NULL)
    } else{ #until k3rd
      Model.final <- readRDS(file="results/EA/EA.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.no.4Q.mean.higher.errors.RData", refhook = NULL)
    }
  } else{ #no higer order moments
    Model.final <- readRDS(file="results/EA/EA.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.no.4Q.mean.higher.errors.RData", refhook = NULL)
  }
}


# ============================================
#   D. PLOT THE RESULTS
# ============================================


KF.load.final <- prepare.KF.model(Model.final, observables)


KF.result3 <- Kalman.filter(KF.load.final$all.parameters, KF.load.final$Y, KF.load.final$X, 
                            KF.load.final$xi.00, KF.load.final$P.00, S="default", 
                            indic.pos.z=KF.load.final$indic_pos)
if(indic.smoother==1){
  KF.result4 <- Kalman.smoother(KF.load.final$all.parameters, KF.load.final$Y, KF.load.final$X,
                                KF.load.final$xi.00, KF.load.final$P.00, S="default",indic.pos.z=KF.load.final$indic_pos)
}

vec.dates <- data.frame(date = observables.with.dates[[1]])


KF.load <- KF.load.final

# Store initial observables
observables <- observables.with.dates.all %>% dplyr::select(-date)
observables <- as.matrix(observables)
