
# ==============================================================================
#   A. PREPARE THE MODEL
# ==============================================================================
indic.cycle     <- TRUE
indic.cycle.use <- FALSE
indic.4th       <- TRUE
#indic.4th.use   <- TRUE 
if (!exists("indic.estimate")) indic.estimate <- FALSE
if (!exists("estimation_start")) estimation_start <- "saved"
if (!exists("save_estimation_results")) save_estimation_results <- indic.estimate
indic.no.4Q <- TRUE

all_mat <- NULL
if(indic.4th){
  all_mat <- readRDS(file="data/processed/all_mat_sparse_4_5_new.RData", refhook = NULL)
}

# ==============================================================================
# Specify the data to be use

if(!indic.compute.mixture){
  
  observables.with.dates <- read_excel("data/processed/Output.US.xlsx")
  
  # Read data from Sheet1
  survey.DATA.US.with.param <- read_excel("data/processed/survey.DATA.US.with.param.xlsx", sheet = "survey.DATA.US.with.param")
  
  # Read data from Sheet2
  survey.DATA.US.G.with.param <- read_excel("data/processed/survey.DATA.US.with.param.xlsx", sheet = "survey.DATA.US.G.with.param")
  
} else{
  
  # Specify the data to be use
  observables.with.dates <- observables.with.dates.US %>%
    full_join(observables.with.dates.US.G, by="date") %>% 
    dplyr::select(-contains("PDS")) %>%
    filter(date > "1981-06-15") %>%
    filter(substr(date,6,7) %in% c("01","04","07","10")) %>%
    dplyr::select(-contains(c("9Q.", "10Q.", "11Q.", "12Q.", "13Q.", "14Q.", "15Q.", "16Q."))) %>%
    mutate(across(-c(1,2, 23), ~if_else(date %in% target_dates, NA_real_, .)))  
  
  observables.with.dates$date <- as.POSIXct(observables.with.dates$date)
  
  if(indic.cycle){
    observables.with.dates <- observables.with.dates  %>%
      inner_join(DATA.US %>% dplyr::select(date,US.log.cpi.cycle), by="date") %>%
      inner_join(DATA.G.US %>% dplyr::select(date,US.log.gdp.cycle), by="date") 
  }
  
}

# Store initial observables
observables.with.dates.all <- observables.with.dates

if(!indic.4th){
  observables.with.dates <- observables.with.dates %>% 
    dplyr::select(-contains(c("k3rd", "k4th")))
}

# Select horizons:
H.US <- matrix(c(12, 15, 18, 21, 24, 60, 120),7,1) # number of columns = number of areas
H.US.G <- matrix(c(12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48),13,1) # number of columns = number of areas


#### remove when PDS data available
H.US <- H.US[1:min(5,length(c(H.US)))]
select.inflation.types.US <- array(0,c(length(c(H.US)),1,4))
select.inflation.types.US[1:min(5,length(c(H.US))),1,3] <- 11 

if(indic.4th){
  select.inflation.types.US[1:min(5,length(c(H.US))),1,3] <- 1111  
}

H.US.G <- H.US
select.inflation.types.US.G <- array(0,c(length(c(H.US.G)),1,4))
select.inflation.types.US.G[1:length(c(H.US.G)),1,3] <- 11 

if(indic.4th){
  select.inflation.types.US.G[1:length(c(H.US.G)),1,3] <- 1111  
}

observables <- observables.with.dates %>% dplyr::select(-date)

# Extract the frequency
freq.count.all <- observables.with.dates %>%
  mutate(year =format(as.Date(date, format="%d/%m/%Y"),"%Y")) %>%
  group_by(year) %>% count
freq <- max(freq.count.all$n)

# ==============================================================================
# Specify Matrix of horizons:
max.H.nb <- max(length(H.US),length(H.US.G))
H <- matrix(NaN,max.H.nb,2)
H[1:length(H.US),1] <- H.US
H[1:length(H.US.G),2] <- H.US.G

# Specify the horizon of interest taking into account the frequency
#H <- cbind(H.US, H.US.G)/(12/freq)
H <- H/(12/freq)

# Specify select.inflation.types array
select.inflation.types <- array(NaN,c(max.H.nb,2,4))
select.inflation.types[1:length(H.US),1,] <- select.inflation.types.US
select.inflation.types[1:length(H.US.G),2,] <- select.inflation.types.US.G


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
area <- "US"
indic.observed <- "q.o.q"
var.type <- "infl.gdp"


## Calculate delta
delta.aux <- matrix(1,m,r)
aux <- apply(delta.aux,1,function(x){sum(x^2)}) #squared element of the column

## Estimate stdv for observables (delta)
stdv.measur <- make.stdv.measure(observables.with.dates,select.inflation.types,r)

## Delete 4Q data
if(indic.no.4Q){
columns_to_replace <- c("4Q.beta.pe","4Q.beta.var","4Q.beta.k3rd","4Q.beta.k4th")
observables[, grepl(paste(columns_to_replace, collapse = "|"), colnames(observables))] <- NA
}

## Delete HP filter data
if(!indic.cycle.use){
  observables[,dim(observables)[2] - c(1,0)] <- NA
  observables.with.dates[,dim(observables)[2] - c(1,0)+1] <- NA
}

## Delete higher order moments
if(!indic.3rd.use){
  columns_to_replace <- c(".3rd")
  #columns_to_replace <- c(".var",".3rd", ".4th")
  observables[, grepl(paste(columns_to_replace, collapse = "|"), colnames(observables))] <- NA
  observables.with.dates[, grepl(paste(columns_to_replace, collapse = "|"), colnames(observables.with.dates))] <- NA
  
}

## Delete higher order moments
if(!indic.4th.use){
  columns_to_replace <- c(".4th")
  #columns_to_replace <- c(".3rd", ".4th")
  #columns_to_replace <- c(".var",".3rd", ".4th")
  observables[, grepl(paste(columns_to_replace, collapse = "|"), colnames(observables))] <- NA
  observables.with.dates[, grepl(paste(columns_to_replace, collapse = "|"), colnames(observables.with.dates))] <- NA
  
}

# Convert in matrix form for simplicity to plot
observables <- as.matrix(observables)

if(indic.estimate){
  
  if (estimation_start == "saved") {
    Model.initial <- readRDS(file="results/US/US.Model.gdp.trend.cycle.4.5.quarterly.best.corr.gdp.with.3rd.4th.RData", refhook = NULL)
    
    model <- list(pi.bar = Model.initial$pi.bar,
                  delta.t = Model.initial$delta.t,
                  delta.c = Model.initial$delta.c,
                  Phi.Y.r = Model.initial$Phi.Y.r,
                  Theta = Model.initial$Theta,
                  Gamma.Y0.r = Model.initial$Gamma.Y0.r,
                  Gamma.Y1.r =  Model.initial$Gamma.Y1.r,
                  nu = Model.initial$nu, # 1st parameter of the non centered gamma process (AGP(nu,phi,mu))
                  phi = Model.initial$phi,
                  mu = Model.initial$mu, # 3rd parameter of the AGP
                  sigma.av = stdv.measur$sigma.av,
                  sigma.var = stdv.measur$sigma.var #, sigma.k3rd = stdv.measur$sigma.k3rd, sigma.k4th = stdv.measur$sigma.k4th
                  
    )
  } else if (estimation_start == "generic") {
    model <- list(pi.bar = matrix(c(mean(observables.with.dates$US.infl, na.rm = TRUE),
                                    mean(observables.with.dates$US.growth, na.rm = TRUE)), 1, r),
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
  
  
  if(indic.4th){
    model$sigma.av[,,3][,2] = model$sigma.av[,,3][,2]/2
    model$sigma.k3rd = stdv.measur$sigma.k3rd/3 #before no change => to fit correctly skewness
    model$sigma.k4th = stdv.measur$sigma.k4th*2#5 #before no change => no need to be more precise for kurtosis
    
  }
  
  
  # Replace element theta[1,2] and theta[2,4] by 0
  model$Theta[1,2] <- 0
  model$Theta[2,4] <- 0
  
  
  # Create vector with estimates to consider
  estimated.Model <- Make.thetas.indicator.trend.cycle.model(model, delta.t.s=TRUE,  delta.c.s=TRUE, Phi.Y.r.s=TRUE, Gamma.Y0.r.s=TRUE,
                                                             Gamma.Y1.r.s=TRUE, nu.s=TRUE, phi.s=TRUE, Theta.s =TRUE)
  
  
  # Estimate only delta.t and delta.c for the second column (GDP growth)
  #step 2 and final
  estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[c(1,(m+2):(2*m))],2] <- 0
  
  
  # We want delta.c1 to be negative and delta.c2 to be positive => put indicator transformation = 3 and 2
  estimated.Model[grep("delta.c", rownames(as.data.frame(estimated.Model)))[c(1:2,(m+1):(m+2))],3] <- c(3,2,2,2)
  
  
  #We want the autoregressive parameters associated to the cycle to be between 0 and 0.99.
  # Phi.Y.r[1:4] have an indic.trans = 4
  nbr.to.excl <- 0
  estimated.Model[grep("Phi.Y.r", rownames(as.data.frame(estimated.Model)))[seq(1,m*m-nbr.to.excl*(m+1),by=(m+1))],3] <- 4
  
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
  estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]
  estimated.Model.list$Theta[2,4] <- - estimated.Model.list$Theta[2,3]
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type, all_mat)
  
  ## Add other information to "estimated.Model"
  estimated.Model.list$areas <- area
  
  # Store all
  Model <- estimated.Model.list
  
  
  # ============================================
  
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
  for(j in 1:3){
    
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
  all.thetas.final.list$Theta[1,2] <- - all.thetas.final.list$Theta[1,1]
  all.thetas.final.list$Theta[2,4] <- - all.thetas.final.list$Theta[2,3]
  
  ## Compute the parameters of the Model
  all.thetas.final.list <- make.parameters.model(all.thetas.final.list, var.type, all_mat)
  
  Model.final <- all.thetas.final.list
  
} else{
  
  #Model.final <- readRDS(file="results/US/US.Model.gdp.trend.cycle.4.5.quarterly.best.corr.gdp.with.3rd.4th.RData", refhook = NULL)
  if(!indic.model.var.only){
    if(indic.4th.use){
      Model.final <- readRDS(file="results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.4th.errors.no.4Q.best.final.new.RData", refhook = NULL)
    } else{
      Model.final <- readRDS(file="results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.errors.no.4Q.best.RData", refhook = NULL)
    }
  } else{
    Model.final <- readRDS(file="results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.no.3rd.4th.no.4Q.RData", refhook = NULL)
   
  }
  
}

# With initial errors
# saveRDS(Model.final, file="results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.4th.errors.no.4Q.best.final.new.RData")
# ?saveRDS
# USE THIS !!! (new fit, correction) for 4-5 joint thetas uniques: Model.final <- readRDS(file="results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.4th.errors.no.4Q.best.final.new.RData", refhook = NULL)
#3109

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


# Save Results in Model.final
Model.final$KF.load.final <- KF.load.final
Model.final$KF.result3 <- KF.result3
Model.final$KF.result4 <- KF.result4

# Save Model.final with KF.results to have updated values stored
if(save_estimation_results && !indic.model.var.only){
  if(indic.4th.use){
    saveRDS(Model.final, file="results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.4th.errors.no.4Q.best.final.new.RData")
  } else{
    saveRDS(Model.final, file="results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.with.3rd.only.errors.no.4Q.best.RData")
  }
} else if(save_estimation_results){
  saveRDS(Model.final, file="results/US/US.Model.infl.gdp.trend.cycle.4.5.quarterly.joint.no.3rd.4th.no.4Q.RData")
}


# Check Results using cpp kalman filter

KF.load <- KF.load.final

# Create a list called StateSpace to lunch Rcpp KF.
StateSpace <- list()
Y_t <- KF.load$Y
StateSpace$nu_t <- KF.load$X%*%t(KF.load$all.parameters$mu)
StateSpace$H <- KF.load$all.parameters$F
StateSpace$N <- KF.load$all.parameters$sigma# Qfunction(KF.load$all.parameters$sigma,KF.load$X)#
StateSpace$mu_t <- KF.load$X%*%KF.load$all.parameters$A
StateSpace$G <- t(KF.load$all.parameters$H)
StateSpace$M <- KF.load$all.parameters$delta #diag(KF.load$all.parameters$delta)
StateSpace$Sigma_0 <- KF.load$P.00
StateSpace$rho_0  <- KF.load$xi.00

test <- KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)

# Store initial observables
observables <- observables.with.dates.all %>% dplyr::select(-date)
observables <- as.matrix(observables)
