# ============================================
#        * * * * PAPER INFLATION - PHD * * * *
# ============================================
#           
# --------------------------------------------
#           Adrien Jean-Paul TSCHOPP 
# ============================================
#         0 - FUNCTION - 0
# --------------------------------------------
#                 * ALL *
# ============================================

# ============================================
#      A.  KALMAN FILTER FUNCTION
# ============================================

#install.packages("doParallel")
#install.packages("mise")
#install.packages('matrixcalc')
#install.packages('mFilter')
#install.packages("roll")
library(Matrix)
library(doParallel)
library(parallel)
library(foreach)
library(doSNOW)
library(lubridate) 
library(readxl)
library(dplyr)
library(zoo)
library(matrixcalc)
library(tidyverse)
library(numDeriv)
library(MASS)
library(bitops)
library(RCurl)
library(optimx)
library(sandwich)   # install.packages("zoo")
library(gdata)
library(compiler)
library(stringr)
library(ggplot2)
library(moments)
library(mFilter)
library(roll)
library(Rcpp)
library(rlang)

options(scipen=999)
#options(digits = 20)
Sys.setenv(LANG = "en")

Kalman.filter <- function(all.parameters, Y, X, xi.00, P.00, S="default", indic.pos.z=0){
  # Note: the data matrices to be specified in the function contain all 
  #       the available data but the calculations are then done using one 
  #       time period at a time. 
  #       So don't worry if the specified dimensions are strange
  
  # Equations of the KF are the following (equation for a given t):
  ## Y = A'X + H'xi + W (W is a white noise of dim n x 1) 
  ## xi = mu + F lag(xi) + V (V is a white noise of dim r x 1)  
  ### R = E[WW'] = delta%*%t(delta)  sigma%*%t(sigma)
  ### Q = E[VV'] = sigma%*%t(sigma)
  ### S = E[WV'] = Cov(W,V)
  
  #### WARNING: IN THIS MODEL; WE ADAPT THE KALMAN FILTER TO INTEGRATE A 
  #### CONSTANT IN THE TRANSITION EQUATION 
  ####
  #### IF YOU WANT TO ADD A CONSTANT FOR THE DYNAMIC OF THE UNOBSERVED,
  #### YOU CAN ALSO CONSIDER THE CASE WITH ONE ADDITIONAL UNOBSERVED THAT IS 1.
  
  # Inputs:
  # T = number of time period for observed variables
  # n = number of observed variables 
  # k = number of deterministic component
  # r = number of unobserved “state” variables.
  # Y is of dimension T x n and is a vector of observed variables
  # X is of dimension T x k and is a a vector or deterministic components (constants, trends, seasonal components,...)
  # xi.00 is of dimension r x 1 and is an initial vector of unobserved variables (initialize with mean)
  # P.00 is of dimension r x r and is an initial matrix of variance-covariance.
  
  # Parameters :
  # F is r x r
  # Q is r x r
  # A is k x n
  # H is r x n
  # R is n x n
  # S is n x r (matrix of zero by default)
  # mu is r x 1 (vector of zero by default)
  r <- dim(xi.00)[1]
  T <- dim(Y)[1]
  n <- dim(Y)[2]
  
  # Extract matrices from the list "all.parameters
  mu <- all.parameters$mu
  F <- all.parameters$F
  sigma <- all.parameters$sigma
  A <- all.parameters$A
  H <- all.parameters$H
  #delta <- all.parameters$delta #comment for new delta
  delta <- diag(all.parameters$delta[1,]) #add for new delta
  
  # Initialize Q and R for t=0
  Q.t <- Qfunction(sigma,xi.00) #Qfunction defined below (transition eq.)
  list.Q.t <- list(Q.t)
  R.t <- delta%*%t(delta)
  list.R.t <- list(R.t)
  
  # Dealing with the default case:
  ## For S:
  if (S[1]=="default"){ # No correlation between W and V
    S <- matrix(0,n,r)
  } else { # Correlation between W and V
    S <- S
  }
  
  # Outputs:
  # xi.tt is r x 1
  # matrix.xi.tt is of dimension T x r
  matrix.xi.tt <- matrix(NA,T,r)
  matrix.xi.tt.aux <- matrix(NA,T,r)
  matrix.xi.ttm1 <- matrix(NA,T+1,r)
  list.xi.tt <- list(NA)
  # P.tt is r x r
  # matrix.P.tt  is of dimension T x (r^2)
  matrix.P.tt <- matrix(NA,T,r*r)
  matrix.P.ttm1 <- matrix(NA,T+1,r*r)
  list.P.tt <- list(NA)
  
  # Intermediate variables of T elements (filled with NA)
  list.xi.ttm1 <- list(NA) 
  list.P.ttm1 <- list(NA)
  
  # Initialize xi.1|0 and P.1|0 for the loop (Equations 1 and 3)
  #(avoids having the initial condition in the final results)
  list.xi.ttm1[[1]] <- mu + F %*% xi.00 # = xi.1|0
  matrix.xi.ttm1[1,] <- t(list.xi.ttm1[[1]])
  list.P.ttm1[[1]] <- F %*% P.00 %*% t(F) + Q.t  # = P.1|0
  matrix.P.ttm1[1,]<-t(vec(list.P.ttm1[[1]]))
  
  # Initialize gains to zero (K.t)
  list.K.t <- list(NA) 
  list.h.t <- list(NA) 
  
  # Initialize log likelihood to zero
  log.lhd <- 0
  log.lhd.vec <- NULL
  penalty.neg.t <- 0 #before NULL
  penalty.neg.vec <- NULL
  penalty.neg.factor <- 100000
  
  # Loop that calculate the filter
  for(t in 1:T){
    
    #2 Forecast of y_t: y.t|t-1 = A'X.t + H'xi.t|t-1
    y.t <- t(A) %*% matrix(X[t,]) + t(H) %*% list.xi.ttm1[[t]]
    
    #4 h.t = H'P.t|t-1 H + R + 2H'S'
    h.t <- t(H) %*% list.P.ttm1[[t]] %*% H + R.t + 2*t(H)%*%t(S)
    list.h.t[[t]] <- h.t
    
    #5 Gain equation: K.t = [P.t|t-1*H + S']h_t^-1
    #K.t <- (list.P.ttm1[[t]] %*% H + t(S)) %*% solve(h.t)
    #list.K.t[[t]] <- K.t
    
    #6 Forecast error: eta.t = y.t - y.t|t-1
    eta.t <- matrix(Y[t,]) - y.t
    
    # Updating step for missing variables 
    indic_notNaN <- !is.na(eta.t)
    
    # Condition when no observation and when some observations are missing or not.
    if (sum(indic_notNaN)==0) { # There is no observation on this date.
      
      #Define n.star to 0 as no observed value
      n.star <- 0
      
      #4 h.t = H'P.t|t-1H + R + 2H'S'
      # No need to redefine it in this case.
      
      #5 Gain equation: K.t = [P.t|t-1*H + S']h_t^-1
      # Use equation 5 to export K.t results.
      K.t.all <- (list.P.ttm1[[t]] %*% H + t(S)) %*% ginv(h.t, tol=1e-30)
      K.t.all[,!indic_notNaN] <- NA
      list.K.t[[t]] <- K.t.all
      
      #7 Final prediction of the latent variables:
      #  xi.t|t =  Fxi.t|t-1 - K.t*eta.t
      list.xi.tt[[t]] <- list.xi.ttm1[[t]]
      matrix.xi.tt[t,]<-t(list.xi.tt[[t]])
      matrix.xi.tt.aux[t,]<-t(list.xi.tt[[t]])
      
      #8 Final prediction for the variance of latent variables
      # P.t|t = P.t|t-1 - K.t(H'P.t|t-1 + S)
      list.P.tt[[t]] <- list.P.ttm1[[t]] 
      matrix.P.tt[t,]<-t(vec(list.P.tt[[t]]))
      
      # log-likelihood update: 
      # Past log.lik + -0.5*n*log(2*pi) - 0.5*log(det(h.t)) - 0.5*eta.t' h.t^-1 eta.t
      # Note: to calculate the log likelihood we do it incrementally (for each t).
      # First part of the log likelihood : -0.5*n*log(2*pi)
      # Log.lik = -0.5*n*T*log(2*pi) - 0.5*sum(det(h.t)) - 0.5*sum(eta.t' h.t eta.t)
      # n.star = length(eta.t.star)
      logl.aux <- -0.5*n.star*log(2*pi)
      
      # Second part of the leg likelihood (no update as h.t and eta.t.star are empty)
      log.lhd <- log.lhd + logl.aux 
      log.lhd.vec <- rbind(log.lhd.vec,
                           logl.aux )
      
    } else { # There is at least one observation on this date.
      
      # Adjust H, R, S, n and eta.t for missing variables (some y_t could be missing)    
      n.star <- sum(indic_notNaN)
      H.star <- matrix(H[,indic_notNaN], r, n.star)
      R.star <- R.t[indic_notNaN, indic_notNaN]
      S.star <- matrix(S[indic_notNaN,], n.star , r)
      eta.t.star <- matrix(eta.t[indic_notNaN], n.star , 1)
      
      # 4 h.t = H'P.t|t-1H + R + 2H'S'
      # h.t.star <- t(H.star) %*% list.P.ttm1[[t]] %*% H.star + R.star + 2*t(H.star)%*%t(S.star)
      h.t.star <- h.t[indic_notNaN, indic_notNaN]
      
      # 5 Gain equation: K.t = [P.t|t-1*H + S']h_t^-1
      K.t <- (list.P.ttm1[[t]] %*% H.star + t(S.star)) %*% 
        ginv(h.t.star, tol=1e-30)
      K.t.all <- (list.P.ttm1[[t]] %*% H + t(S)) %*% ginv(h.t, tol=1e-30)
      # K.t.all[,!indic_notNaN] <- NA
      list.K.t[[t]] <- K.t.all
      
      # 7 Final prediction of the latent variables:
      #  xi.t|t =  Fxi.t|t-1 - K.t*eta.t
      list.xi.tt[[t]] <- list.xi.ttm1[[t]] + K.t %*% eta.t.star
      
      # Adjustment we have to make to the filter pertains to the fact that factors
      # z_t are non negative. For this purpose, after each updating step of the 
      # algorithm, negative entries in the z_t estimate are replaced by 0.
      x.tt.aux <- list.xi.tt[[t]]
      if(sum(indic.pos.z==1)>0){
        penalty.neg.t <- sum(x.tt.aux[indic.pos.z==1,][x.tt.aux[indic.pos.z==1,]<0]^2)*penalty.neg.factor
        penalty.neg.vec <- rbind(penalty.neg.vec, penalty.neg.t)
        list.xi.tt[[t]][indic.pos.z==1,] = pmax(list.xi.tt[[t]][indic.pos.z==1,],0)
      }
      matrix.xi.tt[t,]<-t(list.xi.tt[[t]])
      matrix.xi.tt.aux[t,]<-t(x.tt.aux)
      
      #8 Final prediction for the variance of latent variables
      # P.t|t = P.t|t-1 - K.t(H'P.t|t-1 + S)
      list.P.tt[[t]] <- list.P.ttm1[[t]] - K.t %*% (t(H.star) %*% list.P.ttm1[[t]])
      matrix.P.tt[t,]<-t(vec(list.P.tt[[t]]))
      
      # log-likelihood update: 
      # Past log.lik + -0.5*n*log(2*pi) - 0.5*log(det(h.t)) - 0.5*eta.t' h.t^-1 eta.t
      # Note: to calculate the log likelihood we do it incrementally (for each t).
      # First part of the log likelihood : -0.5*n*log(2*pi)
      # Log.lik = -0.5*n*T*log(2*pi) - 0.5*sum(det(h.t)) - 0.5*sum(eta.t' h.t eta.t)
      # n.star = length(eta.t.star)
      logl.aux <- -0.5*n.star*log(2*pi)
      
      # Second part of the leg likelihood
      #log.lhd <- log.lhd + logl.aux  - (1/2*log(det(h.t))) -
      #  (1/2*t(eta.t.star)%*%ginv(h.t.star)%*%(eta.t.star))
      
      if(length(c(h.t.star))==1){
        det.h.t <- h.t.star
      }else{
        det.h.t <- det(h.t.star)
      }
      
      if(is.nan(det.h.t) || is.na(det.h.t) || is.infinite(det.h.t) || det.h.t <= 0){
        log.lhd <- log.lhd + logl.aux  - 70000000 - penalty.neg.t -
          (1/2*t(eta.t.star)%*%solve(h.t.star, tol=1e-30)%*%(eta.t.star))
        log.lhd.vec <- rbind(log.lhd.vec, 
                             logl.aux  - 70000000 - penalty.neg.t -
                               (1/2*t(eta.t.star)%*%solve(h.t.star, tol=1e-30)%*%(eta.t.star)))
      } else
      {
        log.lhd <- log.lhd + logl.aux  - (1/2*log(det.h.t)) - penalty.neg.t -
          (1/2*t(eta.t.star)%*%solve(h.t.star, tol=1e-30)%*%(eta.t.star))
        log.lhd.vec <- rbind(log.lhd.vec,
                             logl.aux  - (1/2*log(det.h.t)) - penalty.neg.t -
                               (1/2*t(eta.t.star)%*%solve(h.t.star, tol=1e-30)%*%(eta.t.star)))
      }
    }
    
    #1  Forecast latent variables for next iteration
    # xi.t|t-1 =  Fxi.t-1|t-1
    list.xi.ttm1[[t+1]] <- mu + F %*% list.xi.tt[[t]]
    
    # Adjustment we have to make to the filter pertains to the fact that factors
    # z_t are non negative. For this purpose, after each updating step of the 
    # algorithm, negative entries in the z_t estimate are replaced by 0.
    if(sum(indic.pos.z==1)>0){
      x.ttm1.aux <- list.xi.ttm1[[t+1]] 
      #penalty.neg.t <- sum(x.ttm1.aux[indic.pos.z==1,][x.ttm1.aux[indic.pos.z==1,]<0]^2)*penalty.neg.factor + penalty.neg.t
      list.xi.ttm1[[t+1]][indic.pos.z==1,] = pmax(list.xi.ttm1[[t+1]][indic.pos.z==1,],0)
      
    }
    
    # Update Q and R for the next iteration
    Q.t <- Qfunction(sigma,list.xi.tt[[t]])
    list.Q.t[[t+1]] <- Q.t
    
    if(t<T){delta <- diag(all.parameters$delta[t+1,])} #add for new delta/ comment for old delta
    R.t <- delta%*%t(delta)
    list.R.t[[t+1]] <- R.t
    
    #3 Forecast observed variables for next iteration
    #  P.t|t-1 = F P.t-1|t-1 F'  + Q
    list.P.ttm1[[t+1]] <- F %*% list.P.tt[[t]] %*% t(F) + Q.t
    
    # Prepare xi.ttm1 and P.ttm1 for extraction
    matrix.xi.ttm1[t+1,]<-t(list.xi.ttm1[[t+1]])
    matrix.P.ttm1[t+1,]<-t(vec(list.P.ttm1[[t+1]]))
    
  } # End of the loop
  
  # Compute fitted values for observables
  fitted.obs <- X %*% A + matrix.xi.tt%*%H # Fitted observables
  
  # NOTE: P.tt return the var-covariance matrix for all T each column being 
  #       the ith element (by column) of the variance-covariance matrix (vec form).
  return(list("xi.tt"=matrix.xi.tt,
              "xi.tt.aux"=matrix.xi.tt.aux,
              "P.tt"=matrix.P.tt,
              "P.tt-1"=matrix.P.ttm1,
              "P.tt-1.list"=list.P.ttm1,
              "xi.tt-1"=matrix.xi.ttm1,
              "xi.tt-1.list"=list.xi.ttm1,
              "K.t"=list.K.t,
              "h.t"=list.h.t,
              "Q.t"=list.Q.t,
              "R.t"=list.R.t,
              "log.lik"=c(log.lhd),
              "log.lik.vec"=c(log.lhd.vec),
              "penalty.neg.vec"=penalty.neg.vec,
              "fitted.obs"=fitted.obs
  ))
}

# Function that computes the Variance-Covariance Matrix of Transition Equation (Q)
# This function is useful in order to compute if we want to define a time-varying 
# Variance-covariance matrix, that is, a matrix that changes at each time 
# iteration of the Kalman Filter. 
Qfunction <- function(Model, X){  
  # Model: list that contains all the information. Model could also be simply Sigma.
  # X: Vector of latent variable at time t-1
  
  # In that case, additional Transition Equations have been introduced
  # Then Model is a list that contains the Sigma as well as indices.of.addit.var
  # This will add heteroskedasticity
  if(class(Model)=="list"){
    
    # Identify z latent variables
    z <- X[(Model$n+1):(Model$n+Model$q)]
    
    # Compute covariance matrix of transition-eq innovations
    ## Sigma_11 = Theta*diag(Gamma.z0 + Gamma.z1'z_t-1)*Theta' + diag(Gamma.z0 + Gamma.z1'(mu_z+ z_t-1))
    A <- Model$Theta %*% diag(c(Model$Gamma.z0 + t(Model$Gamma.z1)%*%z), Model$q) %*% t(Model$Theta) + 
      diag(c(Model$Gamma.Y0 + t(Model$Gamma.Y1)%*%(Model$mu.z + Model$Phi.z %*%z)))
    
    ## Sigma_22 = diag(Gamma.z0 + Gamma.z1'z_t-1)
    B <- diag(c(Model$Gamma.z0 + t(Model$Gamma.z1)%*%z), Model$q)
    
    ## Sigma_12 = Theta*Sigma_22
    C <- Model$Theta %*% B
    Q1 <- cbind(A,C)
    Q2 <- cbind(t(C),B)
    Q  <- rbind(Q1,Q2)
    
    # Normal case when no heteroskedasticity  
  }else{
    sigma <- Model
    Q <- sigma%*%t(sigma)
  }
  
  return(Q)
}

Rfunction <- function(M,rho,aux){
  # compute covariance matrix of measurement errors
  # indices.of.addit.var indicate to which basic measurement eq. relate the addiitonal ones.
  if(class(M)=="list"){
    # in that case, additional measurement equations have been introduced
    # Then M is a list that contains the sigma as well as indices.of.addit.var
    # step 1: finds connections between additonal var. and spf variables:
    r <- dim(select.inflation.types)[2] # number of areas
    nb.basic.eq <- r + sum(select.inflation.types==1,na.rm=TRUE) +
      2*sum(select.inflation.types==11,na.rm=TRUE)
    nb.addit.var <- length(M$sigmas) - nb.basic.eq # this should correspond to the length of indices.of.addit.var
    indices.of.addit.var <- M$indices.of.addit.var
    M2.augment <- c(M$sigmas^2)
    M2.augment[(nb.basic.eq+1):(nb.basic.eq+nb.addit.var)] <- 
      M2.augment[(nb.basic.eq+1):(nb.basic.eq+nb.addit.var)] +
      M2.augment[indices.of.addit.var]
    R <- diag(c(M2.augment))
    R[indices.of.addit.var,(nb.basic.eq+1):(nb.basic.eq+nb.addit.var)] <-
      diag(c(M2.augment[indices.of.addit.var]))
    R[(nb.basic.eq+1):(nb.basic.eq+nb.addit.var),indices.of.addit.var] <-
      diag(c(M2.augment[indices.of.addit.var]))
  }else{
    R <- diag(c(M^2))
  }
  return(R)
}

# ============================================
#     B. FUNCTION CREATING THE LOADINGS
#        a_1, b_1, a_h, b_h, alpha_1, beta_1, 
#        alpha_h, beta_h, A_h, B_h, a_z, b_z 
# ============================================

# Create the S_p matrix
## S_p = sum_{i=1}^P[e_i^P %x% e_i^P]e_i^P'
make.S.p <- function(p){
  Id <- diag(p)
  S.p <- matrix(0,p*p,p)
  for(i in 1:p){
    S.p <- S.p + (Id[,i]%x%Id[,i])%*%t(Id[,i])
  }
  return(S.p)
}

# Create the S_p_tilde matrix for third cumulants
make.S.p.tilde <- function(p){
  Id <- diag(p)
  S.p <- matrix(0,p^3,p)
  for(i in 1:p){
    S.p <- S.p + (Id[,i]%x%Id[,i]%x%Id[,i])%*%t(Id[,i])
  }
  return(S.p)
}

# Create the S_n_q_bar matrix for third cumulant
make.S.n.q.bar <- function(p,q){
  Id <- diag(p*q)
  Id_bis <- diag(q)
  S.p <- matrix(0,(p*q)^2,q^2)
  for(i in 1:(p*q)){
    
    j <- floor((i-1)/q)
    S.p <- S.p + (Id[,i]%x%Id[,i])%*%t(Id_bis[,i-q*j]%x%Id_bis[,i-q*j])
    
  }  
  return(S.p)
}

# Create matrix Λ_n defined in "A Quadratic Kalman Filter" 
make.S.n.q <- function(n,q){
  Id_n <- diag(n)
  Id_q <- diag(q)
  S.p <- matrix(0,(n*q),(n*q))
  for(n.ind in 1:n){
    for (q.ind in 1:q) {
      
      S.p <- S.p + (Id_n[,n.ind]%*%t(Id_q[,q.ind]))%x%(Id_q[,q.ind]%*%t(Id_n[,n.ind]))
      
    }
  }  
  return(S.p)
}

# Create matrix that allows to pass from vec(X%x%Y'%x%Z) to vec(Xx%Z'%x%Y)
# such that S%*%vec(X%x%Y'%x%Z) = vec(Xx%Z'%x%Y)
make.S.x.z.y <- function(i,j,k){
  # i: length of X, 
  # j: length of Y, 
  # k: length of Z
  Id_i <- diag(i)
  Id_j <- diag(j)
  Id_k <- diag(k)
  S.p <- matrix(0,(i*j*k),(i*j*k))
  for(i.ind in 1:i){
    for (j.ind in 1:j) {
      for(k.ind in 1:k){
        
        S.p <- S.p + (Id_k[,k.ind]%x%Id_i[,i.ind]%x%Id_j[,j.ind])%*%t(Id_j[,j.ind]%x%Id_i[,i.ind]%x%Id_k[,k.ind])
        
      }
    }
  }  
  return(S.p)
}

# Create matrix that allows to pass from vec(X%x%Y'%x%Z) to vec(Zx%Y'%x%X)
# such that S%*%vec(X%x%Y'%x%Z) = vec(Zx%Y'%x%X)
make.S.z.y.x <- function(i,j,k){
  # i: length of X, 
  # j: length of Y, 
  # k: length of Z
  Id_i <- diag(i)
  Id_j <- diag(j)
  Id_k <- diag(k)
  S.p <- matrix(0,(i*j*k),(i*j*k))
  for(i.ind in 1:i){
    for (j.ind in 1:j) {
      for(k.ind in 1:k){
        
        S.p <- S.p + (Id_j[,j.ind]%x%Id_k[,k.ind]%x%Id_i[,i.ind])%*%t(Id_j[,j.ind]%x%Id_i[,i.ind]%x%Id_k[,k.ind])
        
      }
    }
  }  
  return(S.p)
}

# Create matrix that allows to pass from vec(X%x%Y'%x%Z) to vec(Yx%Z'%x%X)
# such that S%*%vec(X%x%Y'%x%Z) = vec(Y%x%Z'%x%X)
make.S.y.z.x <- function(i,j,k){
  # i: length of X, 
  # j: length of Y, 
  # k: length of Z
  Id_i <- diag(i)
  Id_j <- diag(j)
  Id_k <- diag(k)
  S.p <- matrix(0,(i*j*k),(i*j*k))
  for(i.ind in 1:i){
    for (j.ind in 1:j) {
      for(k.ind in 1:k){
        
        S.p <- S.p + (Id_k[,k.ind]%x%Id_j[,j.ind]%x%Id_i[,i.ind])%*%t(Id_j[,j.ind]%x%Id_i[,i.ind]%x%Id_k[,k.ind])
      }
    }
  }  
  return(S.p)
}

# Create matrix that allows to pass from vec(X%x%Y'%x%Z) to vec(Zx%X'%x%Y)
# such that S%*%vec(X%x%Y'%x%Z) = vec(Z%x%X'%x%Y)
make.S.z.x.y <- function(i,j,k){
  # i: length of X, 
  # j: length of Y, 
  # k: length of Z
  Id_i <- diag(i)
  Id_j <- diag(j)
  Id_k <- diag(k)
  S.p <- matrix(0,(i*j*k),(i*j*k))
  for(i.ind in 1:i){
    for (j.ind in 1:j) {
      for(k.ind in 1:k){
        
        S.p <- S.p + (Id_i[,i.ind]%x%Id_k[,k.ind]%x%Id_j[,j.ind])%*%t(Id_j[,j.ind]%x%Id_i[,i.ind]%x%Id_k[,k.ind])
      }
    }
  }  
  return(S.p)
}

# Create matrix that allows to pass from vec(X%x%Y'%x%Z) to vec(Yx%X'%x%Z)
# such that S%*%vec(X%x%Y'%x%Z) = vec(Y%x%X'%x%Z)
make.S.y.x.z <- function(i,j,k){
  # i: length of X, 
  # j: length of Y, 
  # k: length of Z
  Id_i <- diag(i)
  Id_j <- diag(j)
  Id_k <- diag(k)
  S.p <- matrix(0,(i*j*k),(i*j*k))
  for(i.ind in 1:i){
    for (j.ind in 1:j) {
      for(k.ind in 1:k){
        
        S.p <- S.p + (Id_i[,i.ind]%x%Id_j[,j.ind]%x%Id_k[,k.ind])%*%t(Id_j[,j.ind]%x%Id_i[,i.ind]%x%Id_k[,k.ind])
      }
    }
  }  
  return(S.p)
}


# Create A1, A2, B1, B2, C1, C2, D1, D2 that allow to compute u_{3,t}(X_{t+1})
# based on u_{3,t}(Y_{t+1}), u_{3,t}(z_{t+1}), K_t(Y_{t+1}, Y_{t+1},z_{t+1})
# and  K_t(Y_{t+1}, z_{t+1},z_{t+1})
make.matrices.third.cum <- function(n,q){
  
  In <- diag(n)
  Iq <- diag(q)
  
  ### Build A1 and A2 for K(Y,Y,Y)
  A1 <- NULL
  for (i in 1:n){
    
    A1.n <- cbind(matrix(0,n,(n*(i-1))), In,matrix(0,n,(n-i)*n))
    A1.q <- matrix(0, q, n^2)
    A1 <- rbind(A1, A1.n, A1.q)
    
  }
  A1 <- rbind(A1, matrix(0, (n+q)*q, n^2))
  A2 <- cbind(In, matrix(0, n, q))
  
  ### Build B1 and B2 for K(z,z,z)
  B1 <- matrix(0, (n+q)*n, q^2)
  for (i in 1:q){
    
    B1.n <- matrix(0, n, q^2) 
    B1.q <- cbind(matrix(0,q,(q*(i-1))), Iq,matrix(0,q,(q-i)*q))
    B1 <- rbind(B1, B1.n, B1.q)
    
  }
  B2 <- cbind(matrix(0, q, n), Iq)
  
  ### Build C1 and C2 for K(y,y,z)
  C1 <- NULL
  for (i in 1:n){
    
    C1.n <- matrix(0, n, n*q) 
    C1.q <- cbind(matrix(0,q,(q*(i-1))), Iq,matrix(0,q,(n-i)*q))
    C1 <- rbind(C1, C1.n, C1.q)
    
  }
  C1 <- rbind(C1, matrix(0, (n+q)*q, n*q))
  C2 <- A2
  
  ### Build D1 and D2 for K(y,z,y)
  D1 <- A1
  D2 <- B2
  
  ### Build E1 and E2 for K(z,y,y)
  E1 <- matrix(0, (n+q)*n, q*n)
  for (i in 1:q){
    
    E1.n <- cbind(matrix(0,n,(n*(i-1))), In,matrix(0,n,(q-i)*n))
    E1.q <- matrix(0, q, n*q)
    E1 <- rbind(E1, E1.n, E1.q)
    
  }
  E2 <- A2
  
  ### Build F1 and F2 for K(z,z,y)
  F1 <- E1 
  F2 <- B2
  
  ### Build G1 and G2 for K(z,y,z)
  G1 <- B1 
  G2 <- A2
  
  ### Build H1 and H2 for K(y,z,z)
  H1 <- C1 
  H2 <- B2
  
  return(list("A1"=A1, "A2"=A2, "B1"=B1, "B2"=B2, "C1"=C1, "C2"=C2,
              "D1"=D1, "D2"=D2, "E1"=E1, "E2"=E2, "F1"=F1, "F2"=F2,
              "G1"=G1, "G2"=G2, "H1"=H1, "H2"=H2))
}

# Create the gamma_tilde_y_1 matrix
Make.semi.diag.matrix <- function(matrix){
  p <- dim(matrix)[2]
  Id <- diag(p)
  S.p <- NULL
  for(i in 1:p){
    S.p <- rbind(S.p, matrix%*%matrix((Id[,i]%x%Id[,i]),p,p))
  }
  return(S.p)
}

# Create the Lambda_0 and Lambda_1 matrices
make.Lambda <- function(model, mu.z, Phi.z, Gamma.z0, Gamma.z1){
  
  n <- dim(model$Phi.Y)[1] # number of Y factors
  q <- length(c(model$nu))
  S.n <- make.S.p(n)
  S.q <- make.S.p(q)
  Pi.0 <- matrix(0,n,n+q) 
  Pi.1 <- cbind(matrix(0,q,n), diag(q))
  Pi <- rbind(Pi.0,Pi.1)
  
  
  S0 <- diag((n+q)^2)
  
  AUX_S1 <- matrix(0,n+q,n+q)
  AUX_S1[1:n,1:n] <- 1
  S_1 <- S0[,c(AUX_S1)==1]
  
  AUX_S2 <- matrix(0,n+q,n+q)
  AUX_S2[(n+1):(n+q),1:n] <- 1
  S_2 <- S0[,c(AUX_S2)==1]
  
  AUX_S3 <- matrix(0,n+q,n+q)
  AUX_S3[1:n,(n+1):(n+q)] <- 1
  S_3 <- S0[,c(AUX_S3)==1]
  
  AUX_S4 <- matrix(0,n+q,n+q)
  AUX_S4[(n+1):(n+q),(n+1):(n+q)] <- 1
  S_4 <- S0[,c(AUX_S4)==1]
  
  Lambda.0 <- S_1%*%((model$Theta %x% model$Theta)%*%S.q%*%Gamma.z0 + 
                       S.n%*%(model$Gamma.Y0 + t(model$Gamma.Y1)%*%mu.z)) +
    S_2%*%(model$Theta %x% diag(q))%*%S.q%*%Gamma.z0 +
    S_3%*%(diag(q) %x% model$Theta)%*%S.q%*%Gamma.z0 +
    S_4%*%S.q%*%Gamma.z0
  
  Lambda.1 <- (S_1%*%((model$Theta %x% model$Theta)%*%S.q%*%t(Gamma.z1) + 
                        S.n%*%t(model$Gamma.Y1)%*%Phi.z) +
                 S_2%*%(model$Theta %x% diag(q))%*%S.q%*%t(Gamma.z1) +
                 S_3%*%(diag(q) %x% model$Theta)%*%S.q%*%t(Gamma.z1) +
                 S_4%*%S.q%*%t(Gamma.z1))#%*%Pi.1
  
  Lambda.1 <- cbind(matrix(0,(n+q)^2,n),Lambda.1)
  
  
  return(list(Lambda.0=Lambda.0, Lambda.1=Lambda.1))
  
}

# Create the M0 and M1 matrices
Make.M0.M1 <- function(model, mu.z, Phi.z, Gamma.z0, Gamma.z1, all_mat){
  
  n <- dim(model$Phi.Y)[1] # number of Y factors
  q <- length(c(model$nu))
  
  M0 <- Matrix(all_mat$t_A2_x_A1, sparse=T)%*%Matrix(model$Theta%x%(model$Theta%x%model$Theta)%*%all_mat$S.p.tilde.q%*%(2*model$nu*model$mu^3) + 
                                                       #(all_mat$S.z.y.x.n.n.n + all_mat$S.x.z.y.n.n.n + diag(Model$n^3))%*%all_mat$S.n2.n%*%(Model$Theta%x%(all_mat$S.p.n%*%t(Model$Gamma.Y1))%*%all_mat$S.p.q%*%Model$Gamma.z0)) +
                                                       all_mat$M0.M1.mu.3.Y.aux%*%(model$Theta%x%(all_mat$S.p.n%*%t(model$Gamma.Y1))%*%all_mat$S.p.q%*%Gamma.z0), sparse=T) +
    Matrix(all_mat$t_B2_x_B1, sparse=T)%*%Matrix(all_mat$S.p.tilde.q%*%(2*model$nu*model$mu^3),sparse=T) +
    #(t(all_mat$C2)%x%all_mat$C1 + (t(all_mat$D2)%x%all_mat$D1)%*%all_mat$S.x.z.y.n.n.q + (t(all_mat$E2)%x%all_mat$E1)%*%all_mat$S.z.y.x.n.n.q)%*%(Model$Theta%x%(Model$Theta%x%diag(Model$q))%*%all_mat$S.p.tilde.q%*%(2*Model$nu*Model$mu^3) + 
    Matrix(all_mat$M0.M1.K.Y.Y.z.aux,sparse=T)%*%Matrix(model$Theta%x%(model$Theta%x%diag(q))%*%all_mat$S.p.tilde.q%*%(2*model$nu*model$mu^3) + 
                                                          all_mat$S.n2.q%*%(diag(q)%x%(all_mat$S.p.n%*%t(model$Gamma.Y1))%*%all_mat$S.p.q%*%Gamma.z0), sparse = T) +
    #((t(all_mat$F2)%x%all_mat$F1)%*%all_mat$S.z.y.x.n.q.q + (t(all_mat$G2)%x%all_mat$G1)%*%all_mat$S.y.x.z.n.q.q + (t(all_mat$H2)%x%all_mat$H1))%*%(diag(Model$q)%x%(Model$Theta%x%diag(Model$q))%*%all_mat$S.p.tilde.q%*%(2*Model$nu*Model$mu^3))
    Matrix(all_mat$M0.M1.K.Y.z.z.aux,sparse=T)%*%Matrix(diag(q)%x%(model$Theta%x%diag(q))%*%all_mat$S.p.tilde.q%*%(2*model$nu*model$mu^3), sparse=T)
  
  M1 <- Matrix(all_mat$t_A2_x_A1, sparse=T)%*%Matrix(model$Theta%x%(model$Theta%x%model$Theta)%*%all_mat$S.p.tilde.q%*%(6*diag(c(model$mu)^3)*model$phi) + 
                                                       #(all_mat$S.z.y.x.n.n.n + all_mat$S.x.z.y.n.n.n + diag(Model$n^3))%*%all_mat$S.n2.n%*%(Model$Theta%x%(all_mat$S.p.n%*%t(Model$Gamma.Y1))%*%all_mat$S.p.q%*%Model$Gamma.z1)) +
                                                       all_mat$M0.M1.mu.3.Y.aux%*%(model$Theta%x%(all_mat$S.p.n%*%t(model$Gamma.Y1))%*%all_mat$S.p.q%*%Gamma.z1), sparse=T) +
    Matrix(all_mat$t_B2_x_B1, sparse=T)%*%Matrix(all_mat$S.p.tilde.q%*%(6*diag(c(model$mu)^3)*model$phi),sparse=T) +
    #(t(all_mat$C2)%x%all_mat$C1 + (t(all_mat$D2)%x%all_mat$D1)%*%all_mat$S.x.z.y.n.n.q + (t(all_mat$E2)%x%all_mat$E1)%*%all_mat$S.z.y.x.n.n.q)%*%(Model$Theta%x%(Model$Theta%x%diag(Model$q))%*%all_mat$S.p.tilde.q%*%(6*diag(c(Model$mu)^3)*Model$phi) + 
    Matrix(all_mat$M0.M1.K.Y.Y.z.aux, sparse=T)%*%Matrix(model$Theta%x%(model$Theta%x%diag(q))%*%all_mat$S.p.tilde.q%*%(6*diag(c(model$mu)^3)*model$phi) + 
                                                           all_mat$S.n2.q%*%(diag(q)%x%(all_mat$S.p.n%*%t(model$Gamma.Y1))%*%all_mat$S.p.q%*%Gamma.z1),sparse = T) +
    #((t(all_mat$F2)%x%all_mat$F1)%*%all_mat$S.z.y.x.n.q.q + (t(all_mat$G2)%x%all_mat$G1)%*%all_mat$S.y.x.z.n.q.q + (t(all_mat$H2)%x%all_mat$H1))%*%(diag(Model$q)%x%(Model$Theta%x%diag(Model$q))%*%all_mat$S.p.tilde.q%*%(6*diag(c(Model$mu)^3)*Model$phi))
    Matrix(all_mat$M0.M1.K.Y.z.z.aux, sparse=T)%*%Matrix(diag(q)%x%(model$Theta%x%diag(q))%*%all_mat$S.p.tilde.q%*%(6*diag(c(model$mu)^3)*model$phi), sparse=T)
  
  M1 <- cbind(matrix(0,(n+q)^3,n),M1)
  
  return(list("M0"=M0, "M1"=M1))
  
}

Make.M0.M1.fast <- function(model, mu.z, Phi.z, Gamma.z0, Gamma.z1, all_mat){
  
  n <- dim(model$Phi.Y)[1] # number of Y factors
  q <- length(c(model$nu))
  
  # Register a parallel backend with the number of cores you want to use
  registerDoParallel(cores = 7)
  
  # Define the matrices to multiply
  A <- Matrix(all_mat$t_A2_x_A1, sparse = TRUE)
  B <- Matrix(model$Theta%x%(model$Theta%x%model$Theta)%*%all_mat$S.p.tilde.q%*%(2*model$nu*model$mu^3) + 
                all_mat$M0.M1.mu.3.Y.aux%*%(model$Theta%x%(all_mat$S.p.n%*%t(model$Gamma.Y1))%*%all_mat$S.p.q%*%Gamma.z0), sparse = TRUE)
  C <- Matrix(all_mat$t_B2_x_B1, sparse=T)
  D <- Matrix(all_mat$S.p.tilde.q%*%(2*model$nu*model$mu^3),sparse=T)  
  E <- Matrix(all_mat$M0.M1.K.Y.Y.z.aux,sparse=T)
  F <- Matrix(model$Theta%x%(model$Theta%x%diag(q))%*%all_mat$S.p.tilde.q%*%(2*model$nu*model$mu^3) + 
                all_mat$S.n2.q%*%(diag(q)%x%(all_mat$S.p.n%*%t(model$Gamma.Y1))%*%all_mat$S.p.q%*%Gamma.z0), sparse = T) 
  G <- Matrix(all_mat$M0.M1.K.Y.z.z.aux,sparse=T)
  H <- Matrix(diag(q)%x%(model$Theta%x%diag(q))%*%all_mat$S.p.tilde.q%*%(2*model$nu*model$mu^3), sparse=T)
  
  I <- Matrix(model$Theta%x%(model$Theta%x%model$Theta)%*%all_mat$S.p.tilde.q%*%(6*diag(c(model$mu)^3)*model$phi) + 
                all_mat$M0.M1.mu.3.Y.aux%*%(model$Theta%x%(all_mat$S.p.n%*%t(model$Gamma.Y1))%*%all_mat$S.p.q%*%Gamma.z1), sparse=T)
  J <- Matrix(all_mat$S.p.tilde.q%*%(6*diag(c(model$mu)^3)*model$phi),sparse=T)
  K <- Matrix(model$Theta%x%(model$Theta%x%diag(q))%*%all_mat$S.p.tilde.q%*%(6*diag(c(model$mu)^3)*model$phi) + 
                all_mat$S.n2.q%*%(diag(q)%x%(all_mat$S.p.n%*%t(model$Gamma.Y1))%*%all_mat$S.p.q%*%Gamma.z1),sparse = T)
  L <- Matrix(diag(q)%x%(model$Theta%x%diag(q))%*%all_mat$S.p.tilde.q%*%(6*diag(c(model$mu)^3)*model$phi), sparse=T)
  
  
  # Multiply the matrices in parallel using foreach
  M0 <- foreach(i = 1:4, .combine = "+") %dopar% {
    if (i == 1) {
      A %*% B
    } else if (i == 2) {
      C %*% D
    } else if (i == 3) {
      E %*% F
    } else {
      G %*% H
    }
  }
  
  M1 <- foreach(i = 1:4, .combine = "+") %dopar% {
    if (i == 1) {
      A %*% I
    } else if (i == 2) {
      C %*% J
    } else if (i == 3) {
      E %*% K
    } else {
      G %*% L
    }
  }
  
  M1 <- cbind(matrix(0,(n+q)^3,n),M1)
  
  return(list("M0"=M0, "M1"=M1))
  
}

# Create the a_1
## a_1= gamma*mu_x
a.1 <- function(gamma,Model){
  return(sum(gamma*Model$mu.X))
}

# Create the b_1
## a_1= phi.x'*gamma
b.1 <- function(gamma,Model){
  return(t(Model$Phi.X)%*%gamma)
}

# Create the alpha_1
alpha.1 <- function(gamma,Model,S.n,S.q){
  gamma.Y <- gamma[1:Model$n]
  gamma.z <- gamma[(Model$n+1):(Model$n+Model$q)]
  res <- matrix(gamma.Y %x% gamma.Y,nrow=1) %*% 
    ((Model$Theta %x% Model$Theta) %*% S.q %*% Model$Gamma.z0 +
       S.n %*% Model$Gamma.Y0 + S.n %*% t(Model$Gamma.Y1) %*% Model$mu.z) + 
    matrix(gamma.z %x% gamma.z,nrow=1) %*% 
    S.q %*% Model$Gamma.z0 +
    2*matrix(gamma.z %x% gamma.Y,nrow=1) %*% 
    (diag(Model$q) %x% Model$Theta) %*% S.q %*% Model$Gamma.z0
  return(res)
}

# Create the beta_1
beta.1 <- function(gamma,Model,S.n,S.q){
  gamma.Y <- gamma[1:Model$n]
  gamma.z <- gamma[(Model$n+1):(Model$n+Model$q)]
  res <- t(matrix(gamma.Y %x% gamma.Y,nrow=1) %*% 
             ((Model$Theta %x% Model$Theta) %*% S.q %*% t(Model$Gamma.z1) +
                S.n %*% t(Model$Gamma.Y1) %*% Model$Phi.z) +
             matrix(gamma.z %x% gamma.z,nrow=1) %*% 
             S.q %*% t(Model$Gamma.z1) +
             2*matrix(gamma.z %x% gamma.Y,nrow=1) %*% 
             (diag(Model$q) %x% Model$Theta) %*% S.q %*% t(Model$Gamma.z1)
  )
  res <- rbind(matrix(0,Model$n,1),res)
  return(res)
}

# Create the alpha_1_dot.dot
alpha.dot.dot.1 <- function(gamma,Model,S.q){
  gamma.Y <- gamma[1:Model$n]
  gamma.z <- gamma[(Model$n+1):(Model$n+Model$q)]
  res <- 2*matrix((gamma.Y%*%Model$Theta + gamma.z)^3,nrow=1) %*% 
    (Model$nu*Model$mu^3) + 3*(t(Model$Gamma.Y1%*%(gamma.Y^2))%x%
                                 matrix((gamma.Y%*%Model$Theta + gamma.z),nrow=1))%*%
    S.q %*% Model$Gamma.z0
  return(res)
}

# Create the beta_1_dot.dot
beta.dot.dot.1 <- function(gamma,Model,S.q){
  gamma.Y <- gamma[1:Model$n]
  gamma.z <- gamma[(Model$n+1):(Model$n+Model$q)] 
  res <- t(6*matrix((gamma.Y%*%Model$Theta + gamma.z)^3,nrow=1) %*% 
             diag(c(Model$mu^3))%*%Model$phi  + 3*(t(Model$Gamma.Y1%*%(gamma.Y^2))%x%
                                                     matrix((gamma.Y%*%Model$Theta + gamma.z),nrow=1))%*%
             S.q %*% Model$Gamma.z1)
  res <- rbind(matrix(0,Model$n,1),res)
  return(res)
}

# Create the alpha_1_dot.dot.dot
alpha.dot.dot.dot.1 <- function(gamma,Model,S.q, S.q.tilde){
  gamma.Y <- gamma[1:Model$n]
  gamma.z <- gamma[(Model$n+1):(Model$n+Model$q)]
  res <- 6*matrix((gamma.Y%*%Model$Theta + gamma.z)^4,nrow=1) %*% 
    (Model$nu*Model$mu^4)  + 3*(t(Model$Gamma.Y1%*%(gamma.Y^2))%x%
                                  t(Model$Gamma.Y1%*%(gamma.Y^2)))%*%S.q %*% Model$Gamma.z0 +
    12*matrix((gamma.Y%*%Model$Theta + gamma.z),nrow=1)%x%(t(Model$Gamma.Y1%*%(gamma.Y^2))%x%
                                                             matrix((gamma.Y%*%Model$Theta + gamma.z),nrow=1))%*%S.q.tilde%*%(Model$nu*Model$mu^3)
  
  return(res)
}

# Create the beta_1_dot.dot.dot
beta.dot.dot.dot.1 <- function(gamma,Model,S.q, S.q.tilde){
  gamma.Y <- gamma[1:Model$n]
  gamma.z <- gamma[(Model$n+1):(Model$n+Model$q)] 
  res <- t(24*matrix((gamma.Y%*%Model$Theta + gamma.z)^4,nrow=1) %*% 
             diag(c(Model$mu^4))%*%Model$phi + 3*(t(Model$Gamma.Y1%*%(gamma.Y^2))%x%
                                                    t(Model$Gamma.Y1%*%(gamma.Y^2)))%*%S.q %*% t(Model$Gamma.z1) +
             36*matrix((gamma.Y%*%Model$Theta + gamma.z),nrow=1)%x%(t(Model$Gamma.Y1%*%(gamma.Y^2))%x%
                                                                      matrix((gamma.Y%*%Model$Theta + gamma.z),nrow=1))%*%S.q.tilde%*%(diag(c(Model$mu^4))%*%Model$phi))
  
  res <- rbind(matrix(0,Model$n,1),res)
  return(res)
}

# Compute the a_h, b_h, alpha_h, beta_h, alpha.dot_h and beta.dot_h for all 4 cases
compute.loadings <- function(Gamma,Model,H,S.n,S.q,S.q.tilde,indic.5y.in.Xy=0){
  # Here, loadings are computed for three cases:
  # 1. Inflation is average inflation between t and t+h, and inflation between 
  #   t-12 and t is of the form Gamma'X_t. k = h/12.
  #   In this case, we have b_h = 1/k delta (Phi.X^(k-1)*12 + ... + Phi.X^h-12 + Phi.X^h)
  # 2. Inflation between t+h-12 and t+h (simpler).
  #   In this case, we have b_h = delta*Phi.X^h
  # 3. "US.SPF-style": average of Price index for a given year over average Price 
  #   index for the preceding year (in quarter).
  #   In this case, we have b_h = 1/4*delta(Phi.X^h-9 + Phi.X^h-6 + Phi.X^h-3 + Phi.X^h)
  # 4. "KOF/CH SPF-style": average of Price index for a given year over average 
  #   Price index for the preceding year.
  #   In this case, we have b_h = 1/12*delta(Phi.X^h-11 + Phi.X^h-10 + ... + Phi.X^h)
  #
  # indic.5y.in.Xy is equal to one if, for horizons higher
  #     than 5y, what is given (for case-1 inflation) are 5y-in-Xy inflation rates
  #
  # Gamma: parameters that multiply the inflation equation: Gamma is (n+q)x1
  # pi_t = pi_bar + delta*Y_t or pi_t = pi_bar + Gamma*X_t with Gamma = [delta , 0_(n+q)xr ]
  # H: horizons of the considered surveys, ex. c(12, 15, 18, 21, 60, 63, 66, 69)
  
  
  # Create matrix for the loadings
  ## 1st case
  a.h.1st.case.mat <- matrix(NA,max(H),1)
  b.h.1st.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.h.1st.case.mat <- matrix(NA,max(H),1)
  beta.h.1st.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.dot.dot.h.1st.case.mat <- matrix(NA,max(H),1)
  beta.dot.dot.h.1st.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.dot.dot.dot.h.1st.case.mat <- matrix(NA,max(H),1)
  beta.dot.dot.dot.h.1st.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  
  ## 2nd case 
  a.h.2nd.case.mat <- matrix(NA,max(H),1)
  b.h.2nd.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.h.2nd.case.mat <- matrix(NA,max(H),1)
  beta.h.2nd.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.dot.dot.h.2nd.case.mat <- matrix(NA,max(H),1)
  beta.dot.dot.h.2nd.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.dot.dot.dot.h.2nd.case.mat <- matrix(NA,max(H),1)
  beta.dot.dot.dot.h.2nd.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  
  ## 3rd case
  a.h.3rd.case.mat <- matrix(NA,max(H),1)
  b.h.3rd.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.h.3rd.case.mat <- matrix(NA,max(H),1)
  beta.h.3rd.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.dot.dot.h.3rd.case.mat <- matrix(NA,max(H),1)
  beta.dot.dot.h.3rd.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.dot.dot.dot.h.3rd.case.mat <- matrix(NA,max(H),1)
  beta.dot.dot.dot.h.3rd.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  
  ## 4th case
  a.h.4th.case.mat <- matrix(NA,max(H),1)
  b.h.4th.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.h.4th.case.mat <- matrix(NA,max(H),1)
  beta.h.4th.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.dot.dot.h.4th.case.mat <- matrix(NA,max(H),1)
  beta.dot.dot.h.4th.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  alpha.dot.dot.dot.h.4th.case.mat <- matrix(NA,max(H),1)
  beta.dot.dot.dot.h.4th.case.mat <-matrix(NA,max(H),Model$n+Model$q)
  
  # Create the Gamma
  ## 2nd case
  Gamma.2nd.case <- matrix(0,max(H),Model$n+Model$q)
  Gamma.2nd.case[max(H),] <- c(Gamma)
  
  ## 1st case
  Gamma.1st.case <- Gamma.2nd.case
  #aux.index <- seq(max(H),1,by=-12)
  aux.index <- seq(max(H),1,by=-freq)
  Gamma.1st.case[aux.index,] <- matrix(1,nrow=length(aux.index),1) %*% matrix(Gamma,nrow=1)
  
  ## 3rd case
  Gamma.3rd.case <- matrix(0,max(H),Model$n+Model$q)
  #aux.index <- seq(max(H),max(H)-11,by=-3)
  aux.index <- seq(max(H), max(H)-(freq-1),by=-freq/4)
  Gamma.3rd.case[aux.index,] <- 1/4*matrix(1,4,1) %*% matrix(Gamma,nrow=1)
  
  ## 4th case
  Gamma.4th.case <- matrix(0,max(H),Model$n+Model$q)
  #aux.index <- seq(max(H),max(H)-11,by=-1)
  aux.index <- seq(max(H), max(H)-(freq-1),by=-1)
  #Gamma.4th.case[aux.index,] <- 1/12*matrix(1,12,1) %*% matrix(Gamma,nrow=1)
  Gamma.4th.case[aux.index,] <- 1/freq*matrix(1,freq,1) %*% matrix(Gamma,nrow=1)
  
  # Modify the loadings for h > 60 (only 1st case)
  if(indic.5y.in.Xy==1){# Then, for horizons X>60mths, we consider 5y-in-Xmths rates
    if(dim(Gamma.1st.case)[1]>60){
      Gamma.1st.case[1:(dim(Gamma.1st.case)[1]-60),] <-
        0*Gamma.1st.case[1:(dim(Gamma.1st.case)[1]-60),]
    }
  }
  
  # Maximum observed maturity of expected inflation
  maxH <- max(H) # This should also correspond to the number of rows of Gamma
  nb.horizons <- length(c(H))
  a     <- array(NaN,c(1,nb.horizons,4)) # four layers for the four cases
  alpha <- array(NaN,c(1,nb.horizons,4)) # four layers for the four cases
  alpha.dot.dot <- array(NaN,c(1,nb.horizons,4)) # four layers for the four cases
  alpha.dot.dot.dot <- array(NaN,c(1,nb.horizons,4)) # four layers for the four cases
  b     <- array(NaN,c(Model$n+Model$q,nb.horizons,4)) # four layers for the four cases
  beta  <- array(NaN,c(Model$n+Model$q,nb.horizons,4)) # four layers for the four cases
  beta.dot.dot  <- array(NaN,c(Model$n+Model$q,nb.horizons,4)) # four layers for the four cases
  beta.dot.dot.dot  <- array(NaN,c(Model$n+Model$q,nb.horizons,4)) # four layers for the four cases
  
  count <- 0 # will count the number of requested maturities already treated
  for(h in 1:maxH){
    
    if(h==1){
      # Create the a_1
      ## 1st case
      a.h.1st.case <- a.1(c(Gamma.1st.case[maxH,]),Model)
      b.h.1st.case <- b.1(c(Gamma.1st.case[maxH,]),Model)
      alpha.h.1st.case <- alpha.1(c(Gamma.1st.case[maxH,]),Model,S.n,S.q)
      beta.h.1st.case  <- beta.1(c(Gamma.1st.case[maxH,]),Model,S.n,S.q)
      alpha.dot.dot.h.1st.case <- alpha.dot.dot.1(c(Gamma.1st.case[maxH,]),Model,S.q)
      beta.dot.dot.h.1st.case  <- beta.dot.dot.1(c(Gamma.1st.case[maxH,]),Model,S.q)
      alpha.dot.dot.dot.h.1st.case <- alpha.dot.dot.dot.1(c(Gamma.1st.case[maxH,]),Model,S.q, S.q.tilde)
      beta.dot.dot.dot.h.1st.case  <- beta.dot.dot.dot.1(c(Gamma.1st.case[maxH,]),Model,S.q, S.q.tilde)
      
      ## 2nd case
      a.h.2nd.case <- a.h.1st.case
      b.h.2nd.case <- b.h.1st.case
      alpha.h.2nd.case <- alpha.h.1st.case
      beta.h.2nd.case  <- beta.h.1st.case
      alpha.dot.dot.h.2nd.case <- alpha.dot.dot.h.1st.case
      beta.dot.dot.h.2nd.case  <- beta.dot.dot.h.1st.case
      alpha.dot.dot.dot.h.2nd.case <- alpha.dot.dot.dot.h.1st.case
      beta.dot.dot.dot.h.2nd.case  <- beta.dot.dot.dot.h.1st.case
      
      ## 3rd case
      a.h.3rd.case <- a.1(c(Gamma.3rd.case[maxH,]),Model)
      b.h.3rd.case <- b.1(c(Gamma.3rd.case[maxH,]),Model)
      alpha.h.3rd.case <- alpha.1(c(Gamma.3rd.case[maxH,]),Model,S.n,S.q)
      beta.h.3rd.case  <- beta.1(c(Gamma.3rd.case[maxH,]),Model,S.n,S.q)
      alpha.dot.dot.h.3rd.case <- alpha.dot.dot.1(c(Gamma.3rd.case[maxH,]),Model,S.q)
      beta.dot.dot.h.3rd.case  <- beta.dot.dot.1(c(Gamma.3rd.case[maxH,]),Model,S.q)
      alpha.dot.dot.dot.h.3rd.case <- alpha.dot.dot.dot.1(c(Gamma.3rd.case[maxH,]),Model,S.q, S.q.tilde)
      beta.dot.dot.dot.h.3rd.case  <- beta.dot.dot.dot.1(c(Gamma.3rd.case[maxH,]),Model,S.q, S.q.tilde)
      
      ## 4th case
      a.h.4th.case <- a.1(c(Gamma.4th.case[maxH,]),Model)
      b.h.4th.case <- b.1(c(Gamma.4th.case[maxH,]),Model)
      alpha.h.4th.case <- alpha.1(c(Gamma.4th.case[maxH,]),Model,S.n,S.q)
      beta.h.4th.case  <- beta.1(c(Gamma.4th.case[maxH,]),Model,S.n,S.q)
      alpha.dot.dot.h.4th.case <- alpha.dot.dot.1(c(Gamma.4th.case[maxH,]),Model,S.q)
      beta.dot.dot.h.4th.case  <- beta.dot.dot.1(c(Gamma.4th.case[maxH,]),Model,S.q)
      alpha.dot.dot.dot.h.4th.case <- alpha.dot.dot.dot.1(c(Gamma.4th.case[maxH,]),Model,S.q, S.q.tilde)
      beta.dot.dot.dot.h.4th.case  <- beta.dot.dot.dot.1(c(Gamma.4th.case[maxH,]),Model,S.q, S.q.tilde)
      
    }else{
      # Create the a_h for h > 1
      ## 1st case
      ### Alphas
      a.h.1st.case <- a.h.1st.case + a.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model)
      alpha.h.1st.case <- alpha.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model,S.n,S.q) +
        alpha.h.1st.case + a.1(beta.h.1st.case,Model)
      alpha.dot.dot.h.1st.case <- alpha.dot.dot.h.1st.case + a.1(beta.dot.dot.h.1st.case,Model) +
        alpha.dot.dot.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model,S.q) +
        3*t((beta.h.1st.case) %x% (c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case))%*%Model$Lambda.0
      alpha.dot.dot.dot.h.1st.case <- alpha.dot.dot.dot.h.1st.case + a.1(beta.dot.dot.dot.h.1st.case,Model) +
        alpha.dot.dot.dot.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model,S.q, S.q.tilde) +
        4*t((beta.dot.dot.h.1st.case) %x% (c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case))%*%Model$Lambda.0 +
        3*alpha.1(beta.h.1st.case,Model,S.n, S.q) + 
        6*t((c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case) %x% (beta.h.1st.case %x% (c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case)))%*%Model$M.0
      ### Betas
      beta.dot.dot.dot.h.1st.case <-  b.1(beta.dot.dot.dot.h.1st.case,Model) + 
        beta.dot.dot.dot.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model,S.q,S.q.tilde) +
        4*t(Model$Lambda.1)%*%((beta.dot.dot.h.1st.case) %x% (c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case)) +
        3*beta.1(beta.h.1st.case,Model,S.n, S.q) + 
        6*t(Model$M.1)%*%((c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case) %x% (beta.h.1st.case %x% (c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case)))
      beta.dot.dot.h.1st.case <-  b.1(beta.dot.dot.h.1st.case,Model) + 
        beta.dot.dot.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model,S.q) +
        3*t(Model$Lambda.1)%*%((beta.h.1st.case) %x% (c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case))
      beta.h.1st.case <- beta.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model,S.n,S.q) +
        b.1(beta.h.1st.case,Model)
      b.h.1st.case <- b.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model)
      
      
      ## 2nd case
      ### Alphas
      a.h.2nd.case <- a.h.2nd.case + a.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model)
      alpha.h.2nd.case <- alpha.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model,S.n,S.q) +
        alpha.h.2nd.case + a.1(beta.h.2nd.case,Model)
      alpha.dot.dot.h.2nd.case <- alpha.dot.dot.h.2nd.case + a.1(beta.dot.dot.h.2nd.case,Model) +
        alpha.dot.dot.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model,S.q) +
        3*(t(beta.h.2nd.case) %x% t(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case))%*%Model$Lambda.0
      alpha.dot.dot.dot.h.2nd.case <- alpha.dot.dot.dot.h.2nd.case + a.1(beta.dot.dot.dot.h.2nd.case,Model) +
        alpha.dot.dot.dot.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model,S.q,S.q.tilde) +
        4*t((beta.dot.dot.h.2nd.case) %x% (c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case))%*%Model$Lambda.0 +
        3*alpha.1(beta.h.2nd.case,Model,S.n, S.q) + 
        6*t((c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case) %x% (beta.h.2nd.case %x% (c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case)))%*%Model$M.0
      ### Betas
      beta.dot.dot.dot.h.2nd.case <-  b.1(beta.dot.dot.dot.h.2nd.case,Model) + 
        beta.dot.dot.dot.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model,S.q,S.q.tilde) +
        4*t(Model$Lambda.1)%*%((beta.dot.dot.h.2nd.case) %x% (c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case)) +
        3*beta.1(beta.h.2nd.case,Model,S.n, S.q) + 
        6*t(Model$M.1)%*%((c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case) %x% (beta.h.2nd.case %x% (c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case)))
      beta.dot.dot.h.2nd.case <-  b.1(beta.dot.dot.h.2nd.case,Model) + 
        beta.dot.dot.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model,S.q) +
        3*t(Model$Lambda.1)%*%((beta.h.2nd.case) %x% (c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case))
      beta.h.2nd.case <- beta.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model,S.n,S.q) +
        b.1(beta.h.2nd.case,Model)
      b.h.2nd.case <- b.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model)
      
      
      ## 3rd case
      ### Alphas
      a.h.3rd.case <- a.h.3rd.case + a.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model)
      alpha.h.3rd.case <- alpha.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model,S.n,S.q) +
        alpha.h.3rd.case + a.1(beta.h.3rd.case,Model)
      alpha.dot.dot.h.3rd.case <- alpha.dot.dot.h.3rd.case + a.1(beta.dot.dot.h.3rd.case,Model) +
        alpha.dot.dot.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model,S.q) +
        3*(t(beta.h.3rd.case) %x% t(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case))%*%Model$Lambda.0
      alpha.dot.dot.dot.h.3rd.case <- alpha.dot.dot.dot.h.3rd.case + a.1(beta.dot.dot.dot.h.3rd.case,Model) +
        alpha.dot.dot.dot.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model,S.q,S.q.tilde) +
        4*t((beta.dot.dot.h.3rd.case) %x% (c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case))%*%Model$Lambda.0 +
        3*alpha.1(beta.h.3rd.case,Model,S.n, S.q) + 
        6*t((c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case) %x% (beta.h.3rd.case %x% (c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case)))%*%Model$M.0
      ### Betas
      beta.dot.dot.dot.h.3rd.case <-  b.1(beta.dot.dot.dot.h.3rd.case,Model) + 
        beta.dot.dot.dot.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model,S.q,S.q.tilde) +
        4*t(Model$Lambda.1)%*%((beta.dot.dot.h.3rd.case) %x% (c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case)) +
        3*beta.1(beta.h.3rd.case,Model,S.n, S.q) + 
        6*t(Model$M.1)%*%((c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case) %x% (beta.h.3rd.case %x% (c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case)))
      beta.dot.dot.h.3rd.case <-  b.1(beta.dot.dot.h.3rd.case,Model) + 
        beta.dot.dot.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model,S.q) +
        3*t(Model$Lambda.1)%*%((beta.h.3rd.case) %x% (c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case))
      beta.h.3rd.case <- beta.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model,S.n,S.q) +
        b.1(beta.h.3rd.case,Model)
      b.h.3rd.case <- b.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model)
      
      
      ## 4th case
      ### Alphas
      a.h.4th.case <- a.h.4th.case + a.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model)
      alpha.h.4th.case <- alpha.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model,S.n,S.q) +
        alpha.h.4th.case + a.1(beta.h.4th.case,Model)
      alpha.dot.dot.h.4th.case <- alpha.dot.dot.h.4th.case + a.1(beta.dot.dot.h.4th.case,Model) +
        alpha.dot.dot.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model,S.q) +
        3*(t(beta.h.4th.case) %x% t(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case))%*%Model$Lambda.0
      alpha.dot.dot.dot.h.4th.case <- alpha.dot.dot.dot.h.4th.case + a.1(beta.dot.dot.dot.h.4th.case,Model) +
        alpha.dot.dot.dot.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model,S.q,S.q.tilde) +
        4*t((beta.dot.dot.h.4th.case) %x% (c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case))%*%Model$Lambda.0 +
        3*alpha.1(beta.h.4th.case,Model,S.n, S.q) + 
        6*t((c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case) %x% (beta.h.4th.case %x% (c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case)))%*%Model$M.0
      ### Betas
      beta.dot.dot.dot.h.4th.case <-  b.1(beta.dot.dot.dot.h.4th.case,Model) + 
        beta.dot.dot.dot.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model,S.q,S.q.tilde) +
        4*t(Model$Lambda.1)%*%((beta.dot.dot.h.4th.case) %x% (c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case)) +
        3*beta.1(beta.h.4th.case,Model,S.n, S.q) + 
        6*t(Model$M.1)%*%((c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case) %x% (beta.h.4th.case %x% (c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case)))
      beta.dot.dot.h.4th.case <-  b.1(beta.dot.dot.h.4th.case,Model) + 
        beta.dot.dot.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model,S.q) +
        3*t(Model$Lambda.1)%*%((beta.h.4th.case) %x% (c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case))
      beta.h.4th.case <- beta.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model,S.n,S.q) +
        b.1(beta.h.4th.case,Model)
      b.h.4th.case <- b.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model)
      
      
    }
    
    if(sum(h==H)>0){
      # in that case, the loadings are expected in the output
      count <- count + 1
      #deflator <- (1-indic.5y.in.Xy) * (1+as.integer((h-1)/12))+
      #  indic.5y.in.Xy * min(5,1+as.integer((h-1)/12))
      deflator <- (1-indic.5y.in.Xy) * (1+as.integer((h-1)/freq))+
        indic.5y.in.Xy * min(5,1+as.integer((h-1)/freq))
      
      ## 1st case
      a[1,count,1] <- a.h.1st.case/deflator
      b[,count,1]  <- b.h.1st.case/deflator
      alpha[1,count,1] <- alpha.h.1st.case/deflator^2
      beta[,count,1]   <- beta.h.1st.case/deflator^2
      alpha.dot.dot[1,count,1] <- alpha.dot.dot.h.1st.case/deflator^3
      beta.dot.dot[,count,1]   <- beta.dot.dot.h.1st.case/deflator^3
      alpha.dot.dot.dot[1,count,1] <- alpha.dot.dot.dot.h.1st.case/deflator^4
      beta.dot.dot.dot[,count,1]   <- beta.dot.dot.dot.h.1st.case/deflator^4
      
      ## 2nd case
      a[1,count,2] <- a.h.2nd.case
      b[,count,2]  <- b.h.2nd.case
      alpha[1,count,2] <- alpha.h.2nd.case
      beta[,count,2]   <- beta.h.2nd.case
      alpha.dot.dot[1,count,2] <- alpha.dot.dot.h.2nd.case
      beta.dot.dot[,count,2]   <- beta.dot.dot.h.2nd.case
      alpha.dot.dot.dot[1,count,2] <- alpha.dot.dot.dot.h.2nd.case
      beta.dot.dot.dot[,count,2]   <- beta.dot.dot.dot.h.2nd.case
      
      ## 3rd case
      a[1,count,3] <- a.h.3rd.case
      b[,count,3]  <- b.h.3rd.case
      alpha[1,count,3] <- alpha.h.3rd.case
      beta[,count,3]   <- beta.h.3rd.case
      alpha.dot.dot[1,count,3] <- alpha.dot.dot.h.3rd.case
      beta.dot.dot[,count,3]   <- beta.dot.dot.h.3rd.case
      alpha.dot.dot.dot[1,count,3] <- alpha.dot.dot.dot.h.3rd.case
      beta.dot.dot.dot[,count,3]   <- beta.dot.dot.dot.h.3rd.case
      
      ## 4th case
      a[1,count,4] <- a.h.4th.case
      b[,count,4]  <- b.h.4th.case
      alpha[1,count,4] <- alpha.h.4th.case
      beta[,count,4]   <- beta.h.4th.case
      alpha.dot.dot[1,count,4] <- alpha.dot.dot.h.4th.case
      beta.dot.dot[,count,4]   <- beta.dot.dot.h.4th.case
      alpha.dot.dot.dot[1,count,4] <- alpha.dot.dot.dot.h.4th.case
      beta.dot.dot.dot[,count,4]   <- beta.dot.dot.dot.h.4th.case
    }
  }
  return(list(a = a, b = b, alpha = alpha, beta = beta, alpha.dot.dot = alpha.dot.dot, 
              beta.dot.dot = beta.dot.dot, alpha.dot.dot.dot = alpha.dot.dot.dot, 
              beta.dot.dot.dot = beta.dot.dot.dot, Gamma.1st.case=Gamma.1st.case, 
              Gamma.2nd.case=Gamma.2nd.case, Gamma.3rd.case=Gamma.3rd.case,
              Gamma.4th.case=Gamma.4th.case))
}

#compute.loadings.c <- cmpfun(compute.loadings)

# Create a_z, linked to gamma process
a_z <- function(Model,u_z){
  # u_z is of dimension Model$q x k or can be a scalar
  # k is the number of u's at which the function is evaluated
  # The output is of dimension 1 x k
  
  if(class(u_z)[1]=="matrix"){
    k <- dim(u_z)[2]
    mat.ones <- matrix(1,Model$q,k)
    mat.mu <- matrix(Model$mu,Model$q,k)
  }else{
    mat.ones <- 1
    mat.mu <- Model$mu
  }
  return(- t(Model$nu) %*% log(mat.ones - mat.mu*u_z))
}

# Create b_z, linked to gamma process
b_z <- function(Model,u_z){
  # u_z is of dimension Model$q x k or can be a scalar
  # k is the number of u's at which the function is evaluated
  # The output is of dimension Model$q x k
  
  if(class(u_z)[1]=="matrix"){
    k <- dim(u_z)[2]
    mat.ones <- matrix(1,Model$q,k)
    mat.mu <- matrix(Model$mu,Model$q,k)
  }else{
    mat.ones <- 1
    mat.mu <- Model$mu
  }
  return(Model$phi %*% (mat.mu*u_z/(mat.ones - mat.mu*u_z)))
}

# Create A_1
A_1 <- function(Model,u_Y,u_z){
  # u_Y is of dimension Model$n x k
  # u_z is of dimension Model$q x k
  # k is the number of u's at which the function is evaluated
  z.bar <- solve(diag(Model$q)-Model$Phi.z)%*%Model$mu.z
  return(
    a_z(Model,u_z + t(Model$Theta)%*%u_Y + .5 * Model$Gamma.Y1 %*% u_Y^2) -
      t(t(u_Y) %*% Model$Theta %*% z.bar) + .5 * t(Model$Gamma.Y0) %*% u_Y^2
  )
}

# Create B_1
B_1 <- function(Model,u_Y,u_z){
  # u_Y is of dimension Model$n x k
  # u_z is of dimension Model$q x k
  # k is the number of u's at which the function is evaluated
  return(rbind(
    t(Model$Phi.Y) %*% u_Y,
    b_z(Model,u_z + t(Model$Theta)%*%u_Y + .5 * Model$Gamma.Y1 %*% u_Y^2)
  ))
}



# Compute the a_h, b_h, alpha_h, beta_h, A_h and B_h for all 4 cases
compute.loadings.with.AB <- function(Gamma,Model,H,S.n,S.q,indic.5y.in.Xy=0){
  # Here, loadings are computed for three cases:
  # 1. Inflation is average inflation between t and t+h, and inflation between 
  #   t-12 and t is of the form Gamma'X_t. k = h/12.
  #   In this case, we have b_h = 1/k delta (Phi.X^(k-1)*12 + ... + Phi.X^h-12 + Phi.X^h)
  # 2. Inflation between t+h-12 and t+h (simpler).
  #   In this case, we have b_h = delta*Phi.X^h
  # 3. "US.SPF-style": average of Price index for a given year over average Price 
  #   index for the preceding year (in quarter).
  #   In this case, we have b_h = 1/4*delta(Phi.X^h-9 + Phi.X^h-6 + Phi.X^h-3 + Phi.X^h)
  # 4. "KOF/CH SPF-style": average of Price index for a given year over average 
  #   Price index for the preceding year.
  #   In this case, we have b_h = 1/4*delta(Phi.X^h-11 + Phi.X^h-10 + ... + Phi.X^h)
  #
  # indic.5y.in.Xy is equal to one if, for horizons higher
  #     than 5y, what is given (for case-1 inflation) are 5y-in-Xy inflation rates
  #
  # Gamma: parameters that multiply the inflation equation: Gamma is (n+q)x1
  # pi_t = pi_bar + delta*Y_t or pi_t = pi_bar + Gamma*X_t with Gamma = [delta , 0_(n+q)xr ]
  
  
  # Create the Gamma
  ## 2nd case
  Gamma.2nd.case <- matrix(0,max(H),Model$n+Model$q)
  Gamma.2nd.case[max(H),] <- c(Gamma)
  
  ## 1st case
  Gamma.1st.case <- Gamma.2nd.case
  aux.index <- seq(max(H),1,by=-12)
  Gamma.1st.case[aux.index,] <- matrix(1,nrow=length(aux.index),1) %*% matrix(Gamma,nrow=1)
  
  ## 3rd case
  Gamma.3rd.case <- matrix(0,max(H),Model$n+Model$q)
  aux.index <- seq(max(H),max(H)-11,by=-3)
  Gamma.3rd.case[aux.index,] <- 1/4*matrix(1,4,1) %*% matrix(Gamma,nrow=1)
  
  ## 4th case
  Gamma.4th.case <- matrix(0,max(H),Model$n+Model$q)
  aux.index <- seq(max(H),max(H)-11,by=-1)
  Gamma.4th.case[aux.index,] <- 1/12*matrix(1,12,1) %*% matrix(Gamma,nrow=1)
  
  # Modify the loadings for h > 60 (only 1st case)
  if(indic.5y.in.Xy==1){# Then, for horizons X>60mths, we consider 5y-in-Xmths rates
    if(dim(Gamma.1st.case)[1]>60){
      Gamma.1st.case[1:(dim(Gamma.1st.case)[1]-60),] <-
        0*Gamma.1st.case[1:(dim(Gamma.1st.case)[1]-60),]
    }
  }
  
  # Maximum observed maturity of expected inflation
  maxH <- max(H) # This should also correspond to the number of rows of Gamma
  nb.horizons <- length(c(H))
  a     <- array(NaN,c(1,nb.horizons,4)) # four layers for the four cases
  alpha <- array(NaN,c(1,nb.horizons,4)) # four layers for the four cases
  A     <- array(NaN,c(1,nb.horizons,4)) # four layers for the four cases
  b     <- array(NaN,c(Model$n+Model$q,nb.horizons,4)) # four layers for the four cases
  beta  <- array(NaN,c(Model$n+Model$q,nb.horizons,4)) # four layers for the four cases
  B  <- array(NaN,c(Model$n+Model$q,nb.horizons,4)) # four layers for the four cases
  
  
  count <- 0 # will count the number of requested maturities already treated
  for(h in 1:maxH){
    
    if(h==1){
      # Create the a_1, b_1, ...
      ## 1st case
      a.h.1st.case <- a.1(c(Gamma.1st.case[maxH,]),Model)
      b.h.1st.case <- b.1(c(Gamma.1st.case[maxH,]),Model)
      A.h.1st.case <- A_1(Model,c(Gamma.1st.case[maxH,1:Model$n]),
                          c(Gamma.1st.case[maxH,(Model$n+1):(Model$n+Model$q)]))
      alpha.h.1st.case <- alpha.1(c(Gamma.1st.case[maxH,]),Model,S.n,S.q)
      beta.h.1st.case  <- beta.1(c(Gamma.1st.case[maxH,]),Model,S.n,S.q)
      B.h.1st.case <- B_1(Model,c(Gamma.1st.case[maxH,1:Model$n]),
                          c(Gamma.1st.case[maxH,(Model$n+1):(Model$n+Model$q)]))
      
      
      ## 2nd case
      a.h.2nd.case <- a.h.1st.case
      b.h.2nd.case <- b.h.1st.case
      alpha.h.2nd.case <- alpha.h.1st.case
      beta.h.2nd.case  <- beta.h.1st.case
      A.h.2nd.case <- A.h.1st.case
      B.h.2nd.case <- B.h.1st.case
      
      ## 3rd case
      a.h.3rd.case <- a.1(c(Gamma.3rd.case[maxH,]),Model)
      b.h.3rd.case <- b.1(c(Gamma.3rd.case[maxH,]),Model)
      alpha.h.3rd.case <- alpha.1(c(Gamma.3rd.case[maxH,]),Model,S.n,S.q)
      beta.h.3rd.case  <- beta.1(c(Gamma.3rd.case[maxH,]),Model,S.n,S.q)
      A.h.3rd.case <- A_1(Model,c(Gamma.3rd.case[maxH,1:Model$n]),
                          c(Gamma.3rd.case[maxH,(Model$n+1):(Model$n+Model$q)]))
      B.h.3rd.case <- B_1(Model,c(Gamma.3rd.case[maxH,1:Model$n]),
                          c(Gamma.3rd.case[maxH,(Model$n+1):(Model$n+Model$q)]))
      
      
      
      ## 4th case
      a.h.4th.case <- a.1(c(Gamma.4th.case[maxH,]),Model)
      b.h.4th.case <- b.1(c(Gamma.4th.case[maxH,]),Model)
      alpha.h.4th.case <- alpha.1(c(Gamma.4th.case[maxH,]),Model,S.n,S.q)
      beta.h.4th.case  <- beta.1(c(Gamma.4th.case[maxH,]),Model,S.n,S.q)
      A.h.4th.case <- A_1(Model,c(Gamma.4th.case[maxH,1:Model$n]),
                          c(Gamma.4th.case[maxH,(Model$n+1):(Model$n+Model$q)]))
      B.h.4th.case <- B_1(Model,c(Gamma.4th.case[maxH,1:Model$n]),
                          c(Gamma.4th.case[maxH,(Model$n+1):(Model$n+Model$q)]))
      
      
      
    }else{
      # Create the a_h for h > 1
      ## 1st case
      a.h.1st.case <- a.h.1st.case + a.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model)
      b.h.1st.case <- b.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model)
      alpha.h.1st.case <- alpha.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model,S.n,S.q) +
        alpha.h.1st.case + a.1(beta.h.1st.case,Model)
      beta.h.1st.case <- beta.1(c(Gamma.1st.case[maxH-h+1,]) + b.h.1st.case,Model,S.n,S.q) +
        b.1(beta.h.1st.case,Model)
      A.h.1st.case <- A.h.1st.case + A_1(Model,
                                         c(Gamma.1st.case[maxH-h+1,1:Model$n]) + B.h.1st.case[1:(Model$n),],
                                         c(Gamma.1st.case[maxH-h+1,(Model$n+1):(Model$n+Model$q)]) + 
                                           B.h.1st.case[(Model$n+1):(Model$n+Model$q),]
      )
      B.h.1st.case <- B_1(Model,
                          c(Gamma.1st.case[maxH-h+1,1:Model$n]) + B.h.1st.case[1:(Model$n),],
                          c(Gamma.1st.case[maxH-h+1,(Model$n+1):(Model$n+Model$q)]) + 
                            B.h.1st.case[(Model$n+1):(Model$n+Model$q),]
      )
      
      ## 2nd case
      a.h.2nd.case <- a.h.2nd.case + a.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model)
      b.h.2nd.case <- b.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model)
      alpha.h.2nd.case <- alpha.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model,S.n,S.q) +
        alpha.h.2nd.case + a.1(beta.h.2nd.case,Model)
      beta.h.2nd.case <- beta.1(c(Gamma.2nd.case[maxH-h+1,]) + b.h.2nd.case,Model,S.n,S.q) +
        b.1(beta.h.2nd.case,Model)
      A.h.2nd.case <- A.h.2nd.case + A_1(Model,
                                         c(Gamma.2nd.case[maxH-h+1,1:Model$n]) + B.h.2nd.case[1:(Model$n),],
                                         c(Gamma.2nd.case[maxH-h+1,(Model$n+1):(Model$n+Model$q)]) + 
                                           B.h.2nd.case[(Model$n+1):(Model$n+Model$q),]
      )
      B.h.2nd.case <- B_1(Model,
                          c(Gamma.2nd.case[maxH-h+1,1:Model$n]) + B.h.2nd.case[1:(Model$n),],
                          c(Gamma.2nd.case[maxH-h+1,(Model$n+1):(Model$n+Model$q)]) +
                            B.h.2nd.case[(Model$n+1):(Model$n+Model$q),]
      )
      
      ## 3rd case
      a.h.3rd.case <- a.h.3rd.case + a.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model)
      b.h.3rd.case <- b.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model)
      alpha.h.3rd.case <- alpha.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model,S.n,S.q) +
        alpha.h.3rd.case + a.1(beta.h.3rd.case,Model)
      beta.h.3rd.case <- beta.1(c(Gamma.3rd.case[maxH-h+1,]) + b.h.3rd.case,Model,S.n,S.q) +
        b.1(beta.h.3rd.case,Model)
      A.h.3rd.case <- A.h.3rd.case + A_1(Model,
                                         c(Gamma.3rd.case[maxH-h+1,1:Model$n]) + B.h.3rd.case[1:(Model$n),],
                                         c(Gamma.3rd.case[maxH-h+1,(Model$n+1):(Model$n+Model$q)]) +
                                           B.h.3rd.case[(Model$n+1):(Model$n+Model$q),]
      )
      B.h.3rd.case <- B_1(Model,
                          c(Gamma.3rd.case[maxH-h+1,1:Model$n]) + B.h.3rd.case[1:(Model$n),],
                          c(Gamma.3rd.case[maxH-h+1,(Model$n+1):(Model$n+Model$q)]) +
                            B.h.3rd.case[(Model$n+1):(Model$n+Model$q),]
      )
      
      ## 4th case
      a.h.4th.case <- a.h.4th.case + a.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model)
      b.h.4th.case <- b.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model)
      alpha.h.4th.case <- alpha.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model,S.n,S.q) +
        alpha.h.4th.case + a.1(beta.h.4th.case,Model)
      beta.h.4th.case <- beta.1(c(Gamma.4th.case[maxH-h+1,]) + b.h.4th.case,Model,S.n,S.q) +
        b.1(beta.h.4th.case,Model)
      A.h.4th.case <- A.h.4th.case + A_1(Model,
                                         c(Gamma.4th.case[maxH-h+1,1:Model$n]) + B.h.4th.case[1:(Model$n),],
                                         c(Gamma.4th.case[maxH-h+1,(Model$n+1):(Model$n+Model$q)]) + 
                                           B.h.4th.case[(Model$n+1):(Model$n+Model$q),]
      )
      B.h.4th.case <- B_1(Model,
                          c(Gamma.4th.case[maxH-h+1,1:Model$n]) + B.h.4th.case[1:(Model$n),],
                          c(Gamma.4th.case[maxH-h+1,(Model$n+1):(Model$n+Model$q)]) +
                            B.h.4th.case[(Model$n+1):(Model$n+Model$q),]
      )
    }
    
    
    if(sum(h==H)>0){
      # in that case, the loadings are expected in the output
      count <- count + 1
      deflator <- (1-indic.5y.in.Xy) * (1+as.integer((h-1)/12))+
        indic.5y.in.Xy * min(5,1+as.integer((h-1)/12))
      
      ## 1st case
      a[1,count,1] <- a.h.1st.case/deflator
      b[,count,1]  <- b.h.1st.case/deflator
      alpha[1,count,1] <- alpha.h.1st.case/deflator^2
      beta[,count,1]   <- beta.h.1st.case/deflator^2
      A[1,count,1] <- A.h.1st.case/deflator
      B[,count,1]  <- B.h.1st.case/deflator
      
      ## 2nd case
      a[1,count,2] <- a.h.2nd.case
      b[,count,2]  <- b.h.2nd.case
      alpha[1,count,2] <- alpha.h.2nd.case
      beta[,count,2]   <- beta.h.2nd.case
      A[1,count,2] <- A.h.2nd.case
      B[,count,2]  <- B.h.2nd.case
      
      ## 3rd case
      a[1,count,3] <- a.h.3rd.case
      b[,count,3]  <- b.h.3rd.case
      alpha[1,count,3] <- alpha.h.3rd.case
      beta[,count,3]   <- beta.h.3rd.case
      A[1,count,3] <- A.h.3rd.case
      B[,count,3]  <- B.h.3rd.case
      
      ## 4th case
      a[1,count,4] <- a.h.4th.case
      b[,count,4]  <- b.h.4th.case
      alpha[1,count,4] <- alpha.h.4th.case
      beta[,count,4]   <- beta.h.4th.case
      A[1,count,4] <- A.h.4th.case
      B[,count,4]  <- B.h.4th.case
    }
  }
  return(list(a = a, b = b, alpha = alpha, beta = beta, A = A, B = B))
}


compute_AB <- function(Model,u_Y,u_z,HH,indic.average=0){
  # U_Y is of dimension Model$n x k
  # U_z is of dimension Model$q x k
  # k is the number of u's at which the function is evaluated
  # HH is a vector indicating the horizons, ex. seq(12,120,by=12)
  # If indic.average = 1, one considers E(exp(u'[X{+HH}+X{+HH-12}+X{+HH-24}+...]))
  k <- dim(u_Y)[2]
  nb.H <- length(HH)
  A <- array(NaN,c(1,k,nb.H))
  B <- array(NaN,c(Model$n+Model$q,k,nb.H))
  count <- 0
  for(h in 1:max(HH)){
    if(indic.average==2){# In that case, third case of inflation (US SPF type)
      if(h==1){
        A.h <- A_1(Model,u_Y,u_z)
        B.h <- B_1(Model,u_Y,u_z)
      }else{
        #((h-1)/3 == as.integer((h-1)/3))&(h<12)
        if(((h-1)/(freq/4) == as.integer((h-1)/(freq/4)))&(h<=freq)){
          u.1.Y <- u_Y
          u.1.z <- u_z
        }else{
          u.1.Y <- 0*u_Y
          u.1.z <- 0*u_z
        }
        A.h <- A.h_1 + A_1(Model,
                           u.1.Y + B.h_1[1:(Model$n),],
                           u.1.z + B.h_1[(Model$n+1):(Model$n+Model$q),])
        B.h <- B_1(Model,
                   u.1.Y + B.h_1[1:(Model$n),],
                   u.1.z + B.h_1[(Model$n+1):(Model$n+Model$q),])
      }
    }else{
      if(h==1){
        A.h <- A_1(Model,u_Y,u_z)
        B.h <- B_1(Model,u_Y,u_z)
      }else{
        #(indic.average==1)&((h-1)/12 == as.integer((h-1)/12))
        if((indic.average==1)&((h-1)/freq == as.integer((h-1)/freq))){
          # =====================================
          # =====================================
          # =====================================
          if((Indic.5y.in.Xy==1)&(h>61)){
            u.1.Y <- 0*u_Y
            u.1.z <- 0*u_z
          }else{
            u.1.Y <- u_Y
            u.1.z <- u_z
          }
          # =====================================
          # =====================================
          # =====================================
        }else{
          u.1.Y <- 0*u_Y
          u.1.z <- 0*u_z
        }
        A.h <- A.h_1 + A_1(Model,
                           u.1.Y + B.h_1[1:(Model$n),],
                           u.1.z + B.h_1[(Model$n+1):(Model$n+Model$q),])
        B.h <- B_1(Model,
                   u.1.Y + B.h_1[1:(Model$n),],
                   u.1.z + B.h_1[(Model$n+1):(Model$n+Model$q),])
      }
    }
    A.h_1 <- A.h
    B.h_1 <- B.h
    if(sum(h==HH)>0){
      count <- count + 1
      A[,,count] <- c(A.h)
      B[,,count] <- c(B.h)
    }
  }
  return(list(A=A,B=B))
}



#### WORK ON THAT
compute.derivative.for.cumulants <- function(u,Gamma,Model,H,S.n,S.q,indic.5y.in.Xy=0){
  
  # This function computes the multi-horizon Laplace transform of the different inflation rates.
  epsilon <- 10^-5
  
  u_Gamma <- u*Gamma
  
  AB.list <- compute.loadings.with.AB(u_Gamma,Model,H,S.n,S.q,indic.5y.in.Xy)
  
  
  u_Gamma.1epsilon <- (u+epsilon)*Gamma
  u_Gamma.2epsilon <- (u+2*epsilon)*Gamma
  u_Gamma.3epsilon <- (u+3*epsilon)*Gamma
  
  u_Gamma.1_epsilon <- (u-epsilon)*Gamma
  u_Gamma.2_epsilon <- (u-2*epsilon)*Gamma
  u_Gamma.3_epsilon <- (u-3*epsilon)*Gamma
  
  AB.list.1prime <- compute.loadings.with.AB(u_Gamma.1epsilon,Model,H,S.n,S.q,indic.5y.in.Xy)
  AB.list.2prime <- compute.loadings.with.AB(u_Gamma.2epsilon,Model,H,S.n,S.q,indic.5y.in.Xy)
  AB.list.3prime <- compute.loadings.with.AB(u_Gamma.3epsilon,Model,H,S.n,S.q,indic.5y.in.Xy)
  
  AB.list.1_prime <- compute.loadings.with.AB(u_Gamma.1_epsilon,Model,H,S.n,S.q,indic.5y.in.Xy)
  AB.list.2_prime <- compute.loadings.with.AB(u_Gamma.2_epsilon,Model,H,S.n,S.q,indic.5y.in.Xy)
  AB.list.3_prime <- compute.loadings.with.AB(u_Gamma.3_epsilon,Model,H,S.n,S.q,indic.5y.in.Xy)
  
  
  f.1prime.A <- 1/epsilon*(AB.list.1prime$A - AB.list$A)
  f.1prime.B <- 1/epsilon*(AB.list.1prime$B - AB.list$B)
  f.2prime.A  <- 1/(epsilon^2)*(AB.list.2prime$A - 2*AB.list.1prime$A + AB.list$A)
  f.2prime.B  <- 1/(epsilon^2)*(AB.list.2prime$B - 2*AB.list.1prime$B + AB.list$B)
  f.3prime.A <- 1/(epsilon^3)*(AB.list.3prime$A - 3*AB.list.2prime$A + 3*AB.list.1prime$A -
                                 AB.list$A)
  f.3prime.B <- 1/(epsilon^3)*(AB.list.3prime$B - 3*AB.list.2prime$B + 3*AB.list.1prime$B 
                               - AB.list$B)
  
  f.2_prime.A  <- 1/(epsilon^2)*(AB.list$A - 2*AB.list.1_prime$A + AB.list.2_prime$A)
  f.2_prime.B  <- 1/(epsilon^2)*(AB.list$B - 2*AB.list.1_prime$B + AB.list.2_prime$B)
  
  f.2prime.A  <- 1/2*(f.2prime.A + f.2_prime.A)
  f.2prime.B  <- 1/2*(f.2prime.B + f.2_prime.B)
  
  
  return(list("a.h"=f.1prime.A, 
              "b.h"=f.1prime.B, 
              "alpha.h"=f.2prime.A,
              "beta.h"=f.2prime.B, 
              "a_3.h"=f.3prime.A, 
              "b_3.h"=f.3prime.B))
}

third.cumulant.z.t <- function(alpha,phi,nu,mu,z0){
  # This function generate the cumulant for the z_t one period ahead. 
  
  # phi, nu and mu are the parameters of the distribution, should be in matrix form.
  # z0 is the initial value of the simulation, should be a matrix.
  
  # Ensure z0 and the rest are matrice and not scalars
  z0 <- as.matrix(z0)
  phi <- as.matrix(phi) 
  nu <-  as.matrix(nu)
  mu <- as.matrix(mu)
  
  # Third derivative of the a.z evaluated in w=0.
  alpha.dot.dot <- 2*t(alpha^3)%*%(nu*(mu^3))
  
  # Third derivative of the b.z evaluated in w=0.
  beta.dot.dot <- 6*t(alpha^3)%*%diag(c(mu)^3)%*%t(phi)
  
  z <- alpha.dot.dot + beta.dot.dot%*%as.matrix(z0[1,])
  
  return(z)
  
}

###### NEW !!!!!
sim.ARG.process <- function(nb.sim,phi,nu,mu,z0){
  # This function simulate an ARG process for one period ahead, n times. 
  # For n period ahead replace by z_1 <- z[,,i-1]
  ## z0 is a Txq vector.
  ## nu is a qx1 vector
  ## mu is a qx1 vector
  ## phi is a qxq matrix
  ### T = nbr of observations by sample (z_t)
  ### q = nbr of different z_t; nbr of variables
  ### nbr of simulation = number of sample to generate
  
  # nb.sim is the number of simulation that are performed
  # phi, nu and mu are the parameters of the distribution, should be in matrix form.
  # z0 is the initial value of the simulation, should be a matrix.
  # Notations: see our Appendix
  
  q <- dim(z0)[2] #dim(z0)[2] # dimension of the vector to simulate
  T <- dim(z0)[1] #dim(z0)[1] # nbr of observations by sample
  z    <- array(0,c(T,q,nb.sim+1))
  z[,,1] <- z0
  # i: nbr of simulation => loop for each vector
  # j: loop for each z_{i,t}
  
  for (i in 2:(nb.sim+1)){
    z_1   <- z[,,i-1]
    for(j in 1:q){
      Z_poisson <- rpois(T,as.matrix(z_1) %*% phi[,j])
      z[,j,i] <- rgamma(T, shape = nu[j] + Z_poisson, scale = mu[j])
    }
  }
  return(z)
}

simul.model <- function(Model, Y.0, z.0,horizon){
  # Y.0 and z.0 are of dimension T x n and T x q, respectively
  T <- dim(Y.0)[1]
  X <- array(NaN,c(T,Model$n+Model$q,horizon+1))
  X[,1:Model$n,1] <- Y.0
  z <- rVARG(horizon,Model,z.0)
  S.n <- make.S.p(Model$n)
  S.q <- make.S.p(Model$q)
  moments <- unc.moments(Model,S.n,S.q)
  X[,(Model$n+1):(Model$n+Model$q),] <- z
  for(t in 2:(horizon+1)){
    Y_1 <- as.matrix(X[,1:Model$n,t-1])
    Y <- Y_1 %*% t(Model$Phi.Y) + (z[,,t]-matrix(1,T,1)%*%c(moments$z.mom.1)) %*% t(Model$Theta) +
      matrix(rnorm(Model$n*T),T,Model$n) *
      sqrt(matrix(1,T,1)%*%c(Model$Gamma.Y0) + as.matrix(z[,,t]) %*% Model$Gamma.Y1)
    X[,1:Model$n,t] <- Y
  }
  X <- X[,,1:dim(X)[3]]
  return(X)
}

rVARG <- function(nb.sim,Model,z0){
  # z0 is the initial value of the simulation
  # Notations: see our Appendix
  q <- dim(z0)[2] # dimension of the vector to simulate
  T <- dim(z0)[1]
  Ones <- matrix(1,q,1)
  z    <- array(0,c(T,q,nb.sim+1))
  z[,,1] <- z0
  for (i in 2:(nb.sim+1)){
    z_1   <- z[,,i-1]
    for(j in 1:q){
      Z_poisson <- rpois(T,as.matrix(z_1) %*% Model$phi[,j])
      #       if(is.na(Z_poisson)){
      #         print(c(i,z_1))
      #         stop()}
      z[,j,i] <- rgamma(T,shape = Model$nu[j] + Z_poisson,scale = Model$mu[j])
    }
  }
  return(z)
}

compute.phi.pi <- function(Model,u_pi,HH,indic.average=0){
  # This function computes the multi-horizon Laplace transform of the different inflation rates.
  # u_pi is of dimension (r x k)
  # k is the number of points where we want to compute the Laplance transform.
  r <- Model$r
  k <- dim(u_pi)[2]
  
  u_Y <- Model$delta %*% u_pi
  u_z <- matrix(0,Model$q,k)
  
  AB.list <- compute_AB(Model,u_Y,u_z,HH,indic.average)
  
  return(list(A=AB.list$A,B=AB.list$B))
}


# ============================================
#     C. FUNCTION CREATING THE STANDARD
#         DEVIATION OF THE OBSERVABLES
# ============================================

make.stdv.measure <- function(observables.with.dates,select.inflation.types,r){
  # observables.with.dates = all observables data with date as first column
  # select.inflation.types = array with inflation types
  # r = number of area or key.var to be consider
  
  # CALILBRATION
  ## Calibrate measurement errors, the standard deviations of the measurement equations 
  ## that is, to the components of vectors σ_avg and σ_var. These standard 
  ## deviations are calibrated in a preliminary step. The approach of Green and 
  ## Silverman (1994) is used: Smoothing spline to the raw survey-based expectations 
  ## and variances (SPF_t and VSPF_t ). The standard deviations of the differences 
  ## between the survey-based series and their smoothed counterparts. 
  ## The standard errors of the measurement equations are set to these values.
  
  # To plot nicely the graph
  par(mfrow=c(3,3))
  par(plt=c(.15,.9,.2,.9))
  vec.stdv <- NULL
  
  # Loop to plot and derive the spline line.
  for(i in 2:dim(observables.with.dates)[2]){
    
    # Extract series of interest with the date as first column.
    series <- observables.with.dates[,c(1,i)]
    series.with.no.na <- as.data.frame(na.omit(series))
    series.with.no.na[series.with.no.na[,2] > mean(series.with.no.na[,2], na.rm=T) + 2*sd(series.with.no.na[,2], na.rm=T),2] <- NA
    series.with.no.na <- as.data.frame(na.omit(series.with.no.na))
    series.with.no.na[series.with.no.na[,2] < mean(series.with.no.na[,2], na.rm=T) - 2*sd(series.with.no.na[,2], na.rm=T),2] <- NA
    series.with.no.na <- as.data.frame(na.omit(series.with.no.na))
    title <- colnames(series)[2]
    
    # Condition that calculate the spline line only if we have more than 10 obs.
    # Otherwise, we just compute the sd of the series.
    if(dim(series.with.no.na)[1]>10){
      
      # Smooth data and calculate standard deviation based on that to calibrate
      # Note: the lower df, the smoother
      spl <- smooth.spline(series.with.no.na[,1],series.with.no.na[,2],df=max(10,dim(series.with.no.na)[1]/7))
      
      # Plot the data 
      plot(series.with.no.na[,1], series.with.no.na[,2], type="l", xlab = "", main=title)
      lines(spl, col="red")
      vec.stdv <- c(vec.stdv, sd(residuals(spl)))
    }else if(dim(series.with.no.na)[1]>0){
      plot(series.with.no.na[,1], series.with.no.na[,2], type="l",  xlab = "", main=title)
      vec.stdv <- c(vec.stdv, sd(series.with.no.na[,2]))
    } else {
      vec.stdv <- c(vec.stdv, NaN)
    }
  }
  
  # Create vector of NaN
  sigma.av <- NaN * select.inflation.types
  sigma.var <- NaN * select.inflation.types
  sigma.k3rd <- NaN * select.inflation.types #New
  sigma.k4th <- NaN * select.inflation.types #New
  
  count4 <- 0
  for(area.var in 1:r){
    for(infl.type in 1:4){
      if(infl.type>=2){
        count0 <- count4 + 0 
      }else{
        count0 <- count4 + 1 # position of the inflation series
      }
      # First type of inflation:
      ## Calculate the number of point estimate for area.var i and inflation type j.
      nb.pe  <- sum(select.inflation.types[,area.var,infl.type]>=1,na.rm=TRUE) # number of measurement eq on point estimates
      count1 <- count0 + nb.pe
      ## Calculate the number of variance estimate for area.var i and inflation type j.
      nb.var <- sum(select.inflation.types[,area.var,infl.type]>=11,na.rm=TRUE) # number of measurement eq on variances
      count2 <- count1 + nb.var
      nb.k3rd <- sum(select.inflation.types[,area.var,infl.type]>=111,na.rm=TRUE) # number of measurement eq on 3rd cumulants
      count3 <- count2 + nb.k3rd # New
      nb.k4th <- sum(select.inflation.types[,area.var,infl.type]>=1111,na.rm=TRUE) # number of measurement eq on 4th cumulants
      count4 <- count3 + nb.k4th # New
      ## Print counts
      print(c(count0,count1,count2,count3,count4))
      # For area.var==1
      
      # Replace sigma.av with the appropriate variance calculate for the point estimate
      sigma.av[(!is.na(select.inflation.types[,area.var,infl.type])&select.inflation.types[,area.var,infl.type]>=1),area.var,infl.type] <-
        vec.stdv[(count0+1):(count1)]
      # Replace sigma.var with the appropriate variance calculate for the point estimate
      sigma.var[(!is.na(select.inflation.types[,area.var,infl.type])&select.inflation.types[,area.var,infl.type]>=11),area.var,infl.type] <-
        vec.stdv[(count1+1):(count2)]
      # Replace sigma.var with the appropriate variance calculate for the point estimate
      sigma.k3rd[(!is.na(select.inflation.types[,area.var,infl.type])&select.inflation.types[,area.var,infl.type]>=111),area.var,infl.type] <-
        vec.stdv[(count2+1):(count3)]
      # Replace sigma.var with the appropriate variance calculate for the point estimate
      sigma.k4th[(!is.na(select.inflation.types[,area.var,infl.type])&select.inflation.types[,area.var,infl.type]>=1111),area.var,infl.type] <-
        vec.stdv[(count3+1):(count4)]
      
    }
  }
  return(list(sigma.pi=vec.stdv[1],sigma.av=sigma.av,sigma.var=sigma.var,
              sigma.k3rd=sigma.k3rd,sigma.k4th=sigma.k4th))
}


# ============================================
#     D. FUNCTION CREATING THE PARAMETERS
#                 OF THE MODEL
# ============================================

make.parameters.model <- function(model, var.type, all_mat=NULL){
  # model: contains the initial parameters of the model
  # This function use the initial parameter of the model to compute the parameters
  # associated with the q latent z_t.
  
  # Extract dimensions
  n <- dim(model$Theta)[1] # number of Y factors
  m <- dim(model[[2]])[1]
  q <- length(c(model$nu)) # number of z factors
  r <- length(c(model$pi.bar))
  
  # Compute the matrices of interest
  model$Phi.Y <- make.parameters.trend.cycle.model(model, var.type, n, m, q, r)$Phi.Y
  model$Gamma.Y1 <- make.parameters.trend.cycle.model(model, var.type, n, m, q, r)$Gamma.Y1
  model$Gamma.Y0 <- make.parameters.trend.cycle.model(model, var.type, n, m, q, r)$Gamma.Y0
  model$delta <- make.parameters.trend.cycle.model(model, var.type, n, m, q, r)$delta
  #Remove these two line to model annual inflation and GDP 
  ### Don't forget to remove, delta.q in the prepare.KF.model function
  
  if(indic.observed =="y.o.y" | indic.observed==""){
    model$delta.q <- model$delta
  } else if(indic.observed == "q.o.q"){
    
    model$delta.q <- cbind(model$delta[,1],model$delta[,2]
                           #model$delta[,2]*4
    )
    model$delta.q[(2*m+1):n,] <- 0 
    
  }
  
  #Compute all the parameters to the z_t
  ## mu_z = mu * nu (* product element by element)
  mu.z <- matrix(model$mu,q,1) * matrix(model$nu,q,1)
  
  ## Phi_z
  Phi.z <- (matrix(model$mu,q,1)%*%matrix(1,1,q)) * t(model$phi)
  ## Gamma_z0
  Gamma.z0 <- matrix(model$mu,q,1) * matrix(model$mu,q,1) * matrix(model$nu,q,1)
  ## Gamma_z1
  Gamma.z1 <- 2 * t(((matrix(model$mu,q,1) * matrix(model$mu,q,1))%*%matrix(1,1,q)) *
                      t(model$phi))
  ## mu.X
  mu.X <- rbind(matrix(- model$Theta %*% Phi.z %*% solve(diag(q) - Phi.z) %*% mu.z,ncol=1),
                matrix(mu.z,ncol=1))
  ## Phi.X
  Phi.X <- cbind(model$Phi.Y,model$Theta %*% Phi.z)
  Phi.X <- rbind(Phi.X,cbind(matrix(0,q,n),Phi.z))
  
  Lambda.0 <- make.Lambda(model, mu.z, Phi.z, Gamma.z0, Gamma.z1)$Lambda.0
  Lambda.1 <- make.Lambda(model, mu.z, Phi.z, Gamma.z0, Gamma.z1)$Lambda.1
  
  if(is.null(all_mat)){
    M.0 <- matrix(0,(n+q)^3,1)
    M.1 <- matrix(0,(n+q)^3,n+q)
    
  } else{
    #M.0.1 <- Make.M0.M1(model, mu.z, Phi.z, Gamma.z0, Gamma.z1, all_mat)
    M.0.1 <- Make.M0.M1.fast(model, mu.z, Phi.z, Gamma.z0, Gamma.z1, all_mat)
    M.0 <- as.matrix(M.0.1$M0)
    M.1 <- as.matrix(M.0.1$M1)}
  
  
  return(merge.list(model,list(r=r,n=n,m=m,q=q,mu.z=mu.z,Phi.z=Phi.z,Gamma.z0=Gamma.z0,Gamma.z1=Gamma.z1,
                               mu.X=mu.X,Phi.X=Phi.X, Lambda.0=Lambda.0, Lambda.1=Lambda.1, M.0 = M.0, M.1=M.1)))
}


make.parameters.trend.cycle.model <- function(model, var.type="gdp", n, m, q, r){
  # model: contains the initial parameters of the model
  # This function use the initial parameter of the model to compute the parameters
  # in the trend cycle model.
  # note that freq is equal to: 
  ## var.type=="gdp": freq = (n-2*m + 2)/2
  ## var.type=="infl": freq = n-2*m+1
  ## var.type=="infl": freq = (n-2*m + 3)/3
  
  if(n > m){
    
    # Calculate Gamma.Y0 and Gamma.Y1
    Gamma.Y0.big <- rbind(model$Gamma.Y0.r,as.matrix(rep(0,n-m)))
    Gamma.Y1.big <- cbind(model$Gamma.Y1.r,matrix(0,q,n-m))
    
    # Calculate Phi.Y 
    if(r == 1){
      Phi.Y.1 <- cbind(model$Phi.Y.r, matrix(0,m,n-m))
      Phi.Y.2 <- cbind(diag(m), matrix(0,m,n-m))
      Phi.Y.3 <- cbind(t(model$delta.t) + t(model$delta.c), -t(model$delta.c), matrix(0,1,n-2*m))
      Phi.Y.4 <- cbind(matrix(0,n-2*m-1,2*m), diag(n-2*m-1), matrix(0,n-2*m-1,1))
      
      Phi.Y.big <- rbind(Phi.Y.1 , Phi.Y.2, Phi.Y.3, Phi.Y.4)}
    else if(r == 2){
      # Compute freq
      #freq <- (n-2*m + 3)/3
      freq <- (n-2*m)/2 +1
      # Compute nbr of lag for inflation
      n.infl <- freq-1
      
      # Compute nbr of lag for gdp  
      #n.gdp <- 2*freq-2
      n.gdp <- freq-1
      
      Phi.Y.1 <- cbind(model$Phi.Y.r, matrix(0,m,n-m))
      Phi.Y.2 <- cbind(diag(m), matrix(0,m,n-m))
      Phi.Y.3 <- cbind(t(model$delta.t[,1]) + t(model$delta.c[,1]), -t(model$delta.c[,1]), matrix(0,1,n-2*m))
      Phi.Y.4 <- cbind(matrix(0,n.infl-1,2*m), diag(n.infl-1), matrix(0,n.infl-1,n -2*m - (n.infl-1)))
      Phi.Y.5 <- cbind(t(model$delta.t[,2]) + t(model$delta.c[,2]), -t(model$delta.c[,2]), matrix(0,1,n-2*m))
      Phi.Y.6 <- cbind(matrix(0,n.gdp-1,2*m + n.infl), diag(n.gdp-1), matrix(0,n.gdp-1, 1))
      
      Phi.Y.big <- rbind(Phi.Y.1 , Phi.Y.2, Phi.Y.3, Phi.Y.4, Phi.Y.5, Phi.Y.6)
    } else {
      stop('error: r > 2 !!! you should define r between 1 and 2')
    }
    
    # Calculate delta
    if(var.type=="gdp"){
      delta.big <- 1/((n-2*m+2)/2)*rbind(model$delta.t + model$delta.c, -model$delta.c, 
                                         matrix(c(1+seq(1,(n-2*m)/2), seq((n-2*m)/2,1)),n-2*m,1))
    } else if(var.type=="infl"){
      delta.big <- rbind(model$delta.t + model$delta.c, -model$delta.c, 
                         matrix(rep(1,(n-2*m)),n-2*m,1))
      #matrix(c(rep(1,(n-2*m)/2), rep(0,(n-2*m)/2)),n-2*m,1))
    } else if(var.type=="infl.gdp"){
      delta.infl <- rbind(matrix(model$delta.t[,1] + model$delta.c[,1],m,1), -matrix(model$delta.c[,1],m,1), 
                          matrix(c(rep(1,n.infl), rep(0,n.gdp)),n-2*m,1))
      # delta.gdp <- 1/((n-2*m+3)/3)*rbind(matrix(model$delta.t[,2] + model$delta.c[,2],m,1), -matrix(model$delta.c[,2],m,1), 
      #                     matrix(c(rep(0,n.infl), 1+seq(1,(n-2*m)/3), seq((n-2*m)/3,1)),n-2*m,1))
      delta.gdp <- rbind(matrix(model$delta.t[,2] + model$delta.c[,2],m,1), -matrix(model$delta.c[,2],m,1), 
                         matrix(c(rep(0,n.infl), rep(1,n.infl)),n-2*m,1))
      delta.big <- cbind(delta.infl, delta.gdp)
    } else {
      stop('error: the "var.type" defined does not exist !!! you should define var.type="gdp", var.type="infl" or var.type="infl.gdp".')
    }
    
  } else if(n == m){
    Phi.Y.big <- model$Phi.Y
    Gamma.Y0.big <- c(model$Gamma.Y0)
    Gamma.Y1.big <- model$Gamma.Y1
    delta.big <- model$delta 
    
    
  } else {
    stop("error: n < m or m doesn't exist !!! you should define n = m or n > m.")
  }
  
  return(list(Phi.Y=Phi.Y.big, Gamma.Y0=Gamma.Y0.big, Gamma.Y1=Gamma.Y1.big,
              delta=delta.big))  
}

# ============================================
#     E. FUNCTION CREATING UNCONDITIONAL
#                   MOMENT
# ============================================

# Function that computes the unconditional moments of z (first and second).
unc.moments <- function(Model,S.n.input=NaN,S.q.input=NaN){
  if(is.na(S.n.input[1])){
    S.n <- make.S.p(Model$n)
    S.q <- make.S.p(Model$q)
  }else{
    S.n <- S.n.input
    S.q <- S.q.input
  }
  # Expectation 
  z.mom.1 <- solve(diag(Model$q)-Model$Phi.z)%*%Model$mu.z
  
  # Variance
  z.vec.mom.2 <- solve(diag((Model$q)^2)-Model$Phi.z%x%Model$Phi.z)%*%S.q %*%(Model$Gamma.z0+t(Model$Gamma.z1)%*%z.mom.1)
  z.mom.2 <- matrix(z.vec.mom.2, Model$q, Model$q)
  return(list(z.mom.1 = z.mom.1, z.mom.2 = z.mom.2))
}

unc.moments.X <- function(Model,S.n.input=NaN,S.q.input=NaN){
  if(is.na(S.n.input[1])){
    S.n <- make.S.p(Model$n)
    S.q <- make.S.p(Model$q)
  }else{
    S.n <- S.n.input
    S.q <- S.q.input
  }
  # mu_z
  unc.mom.z <- unc.moments(Model,S.n,S.q) # Computes the unconditional moments of z
  # mu_x
  
  mu.y <- -Model$Theta%*%Model$Phi.z%*%solve(diag(Model$q) - Model$Phi.z)%*%Model$mu.z
  X.0   <- c(mu.y,unc.mom.z$z.mom.1) # unconditional mean of the state vector
  aux.Sigma0 <- Qfunction(Model,X.0) # Computes the unconditional value of the conditional covariance matrix
  #Sigma.0.vec <- ginv(diag((Model$n+Model$q)^2) - Model$Phi.X%x%Model$Phi.X) %*% c(aux.Sigma0)
  Sigma.0.vec <- Matrix::solve(diag((Model$n+Model$q)^2) - Model$Phi.X%x%Model$Phi.X, tol=10^(-100)) %*% c(aux.Sigma0)
  Sigma.0 <- matrix(Sigma.0.vec,Model$n+Model$q,Model$n+Model$q)
  return(list(X.mom.1 = X.0, X.mom.2 = Sigma.0))
}

# ============================================
#     F. FUNCTION CREATING THE LIST OF
#           PARAMATERS TO ESTIMATE
# ============================================


Make.thetas.indicator.trend.cycle.model <- function(model, pi.bar.s=FALSE, delta.t.s=FALSE,
                                                    delta.c.s=FALSE, Phi.Y.r.s=FALSE, Theta.s=FALSE, 
                                                    Gamma.Y0.r.s=FALSE,Gamma.Y1.r.s=FALSE, nu.s=FALSE, 
                                                    phi.s=FALSE, mu.s=FALSE, sigma.av.s=FALSE, 
                                                    sigma.var.s=FALSE, sigma.k3rd.s=FALSE, sigma.k4th.s=FALSE){
  
  # Note: All the parameters should be specified in the function
  #       The function creates a matrix with the 3 columns:
  #         - The first one returns all the parameters with their original values.
  #         - The second one returns a vector with a 1 for the parameters to be estimated (0 if not).
  #         - The third one returns a vector with a 1 if the estimated parameters should be between 0 and 1,
  #           a 2 if the parameter should only be positive and 3 if the parameter should be negative
  #           (0 otherwise).
  #       By default the function will only consider F to be estimated.
  
  pi.bar <- model$pi.bar
  delta.t <- model$delta.t
  delta.c <- model$delta.c
  Phi.Y.r <- model$Phi.Y.r
  Theta <- model$Theta
  Gamma.Y0.r <- model$Gamma.Y0.r
  Gamma.Y1.r <- model$Gamma.Y1.r
  nu <- model$nu
  phi <- model$phi
  mu <- model$mu
  sigma.av <- model$sigma.av[!is.na(model$sigma.av)]
  sigma.var <- model$sigma.var[!is.na(model$sigma.var)]
  sigma.k3rd <- model$sigma.k3rd[!is.na(model$sigma.k3rd)]
  sigma.k4th <- model$sigma.k4th[!is.na(model$sigma.k4th)]
  
  # Conditions for pi.bar
  if(pi.bar.s == TRUE){
    
    pi.bar.i <- rep(1,length(pi.bar))
    pi.bar.t <- rep(0, length(pi.bar))
    
  } else {pi.bar.i <- rep(0, length(pi.bar))
  pi.bar.t <- rep(0, length(pi.bar))}
  
  # Conditions for delta.t
  if(delta.t.s == TRUE){
    
    delta.t.i <- rep(1,length(delta.t))
    delta.t.i[which(delta.t == 0)] <- 0 #Remove if wants to estimate 0.
    delta.t.t <- rep(0, length(delta.t)) 
    #delta.t.t[which(delta.t < 0)] <- 3
    
  } else {delta.t.i <- rep(0, length(delta.t))
  delta.t.t <- rep(0, length(delta.t))}
  
  # Conditions for delta.c
  if(delta.c.s == TRUE){
    
    delta.c.i <- rep(1,length(delta.c))
    delta.c.i[which(delta.c == 0)] <- 0 #Remove if wants to estimate 0.
    delta.c.t <- rep(0, length(delta.c))
    #delta.c.t[which(delta.c < 0)] <- 3
    
  } else {delta.c.i <- rep(0, length(delta.c))
  delta.c.t <- rep(0, length(delta.c))}
  
  
  # Conditions for Phi.Y.r
  if(Phi.Y.r.s == TRUE){
    
    Phi.Y.r.i <- rep(1, length(Phi.Y.r))
    Phi.Y.r.i[which(Phi.Y.r == 0)] <- 0
    Phi.Y.r.t <- rep(0, length(Phi.Y.r))
    Phi.Y.r.diag <- Phi.Y.r
    diag(Phi.Y.r.diag) <- "diag"
    Phi.Y.r.t[which(Phi.Y.r.diag == "diag")] <- 1
    
  } else {Phi.Y.r.i <- rep(0, length(Phi.Y.r))
  Phi.Y.r.t <- rep(0, length(Phi.Y.r))}
  
  # Conditions for Theta
  if(Theta.s == TRUE){
    
    Theta.i <- rep(1,length(Theta))
    Theta.i[which(Theta == 0)] <- 0 #Remove if wants to estimate 0.
    Theta.t <- rep(0, length(Theta))
    Theta.t[which(Theta > 0)] <- 2
    Theta.t[which(Theta < 0)] <- 3
    
  } else {Theta.i <- rep(0, length(Theta))
  Theta.t <- rep(0, length(Theta))}
  
  # Conditions for Gamma.Y0.r
  if(Gamma.Y0.r.s == TRUE){
    
    Gamma.Y0.r.i <- rep(1,length(Gamma.Y0.r))
    Gamma.Y0.r.i[which(Gamma.Y0.r == 0)] <- 0 #Remove if wants to estimate 0.
    Gamma.Y0.r.t <- rep(2, length(Gamma.Y0.r))
    
  } else {Gamma.Y0.r.i <- rep(0, length(Gamma.Y0.r))
  Gamma.Y0.r.t <- rep(0, length(Gamma.Y0.r))}
  
  # Conditions for Gamma.Y1.r
  if(Gamma.Y1.r.s == TRUE){
    
    Gamma.Y1.r.i <- rep(1,length(Gamma.Y1.r))
    Gamma.Y1.r.i[which(Gamma.Y1.r == 0)] <- 0 #Remove if wants to estimate 0.
    Gamma.Y1.r.t <- rep(2, length(Gamma.Y1.r))
    
  } else {Gamma.Y1.r.i <- rep(0, length(Gamma.Y1.r))
  Gamma.Y1.r.t <- rep(0, length(Gamma.Y1.r))}
  
  # Conditions for nu
  if(nu.s == TRUE){
    
    nu.i <- rep(1,length(nu))
    nu.i[which(nu == 0)] <- 0 #Remove if wants to estimate 0.
    nu.t <- rep(2, length(nu))
    
  } else {nu.i <- rep(0, length(nu))
  nu.t <- rep(0, length(nu))}
  
  # Conditions for phi
  if(phi.s == TRUE){
    
    phi.i <- rep(1, length(phi))
    phi.i[which(phi == 0)] <- 0
    phi.t <- rep(2, length(phi))
    phi.diag <- phi
    diag(phi.diag) <- "diag"
    #phi.t[which(phi.diag == "diag")] <- 1
    phi.t[which(phi.diag == "diag")] <- 4
    
  } else {phi.i <- rep(0, length(phi))
  phi.t <- rep(0, length(phi))}
  
  # Conditions for mu
  if(mu.s == TRUE){
    
    mu.i <- rep(1,length(mu))
    mu.t <- rep(0, length(mu))
    
  } else {mu.i <- rep(0, length(mu))
  mu.t <- rep(0, length(mu))}
  
  # Conditions for sigma.av
  if(sigma.av.s == TRUE){
    
    sigma.av.i <- rep(1, length(sigma.av))
    sigma.av.i[which(sigma.av == 0)] <- 0
    sigma.av.t <- rep(2, length(sigma.av))
    
  } else {sigma.av.i <- rep(0, length(sigma.av))
  sigma.av.t <- rep(0, length(sigma.av))}
  
  # Conditions for sigma.var
  if(sigma.var.s == TRUE){
    
    sigma.var.i <- rep(1, length(sigma.var))
    sigma.var.i[which(sigma.var == 0)] <- 0
    sigma.var.t <- rep(2, length(sigma.var))
    
  } else {sigma.var.i <- rep(0, length(sigma.var))
  sigma.var.t <- rep(0, length(sigma.var))}
  
  # Conditions for sigma.k3rd
  if(sigma.k3rd.s == TRUE){
    
    sigma.k3rd.i <- rep(1, length(sigma.k3rd))
    sigma.k3rd.i[which(sigma.k3rd == 0)] <- 0
    sigma.k3rd.t <- rep(2, length(sigma.k3rd))
    
  } else {sigma.k3rd.i <- rep(0, length(sigma.k3rd))
  sigma.k3rd.t <- rep(0, length(sigma.k3rd))}
  
  # Conditions for sigma.k4th
  if(sigma.k4th.s == TRUE){
    
    sigma.k4th.i <- rep(1, length(sigma.k4th))
    sigma.k4th.i[which(sigma.k4th == 0)] <- 0
    sigma.k4th.t <- rep(2, length(sigma.k4th))
    
  } else {sigma.k4th.i <- rep(0, length(sigma.k4th))
  sigma.k4th.t <- rep(0, length(sigma.k4th))}
  
  # Create vector indicating if the parameter should be estimated 
  # and vector indicating if the parameter should be transformed.
  vector.indicator.estimate <- c(pi.bar.i, delta.t.i, delta.c.i, Phi.Y.r.i, Theta.i, Gamma.Y0.r.i,
                                 Gamma.Y1.r.i, nu.i, phi.i, mu.i, sigma.av.i, sigma.var.i,
                                 sigma.k3rd.i, sigma.k4th.i)
  vector.indicator.transformation <- c(pi.bar.t, delta.t.t,  delta.c.t, Phi.Y.r.t, Theta.t, Gamma.Y0.r.t,
                                       Gamma.Y1.r.t, nu.t, phi.t, mu.t, sigma.av.t, sigma.var.t,
                                       sigma.k3rd.t, sigma.k4th.t)
  
  # Create list with all the parameters in matrix form
  all.thetas <- model
  
  # Create the vector of parameters
  all.thetas.vec <- na.omit(unlist(all.thetas))
  
  return(cbind(all.thetas.vec, vector.indicator.estimate,vector.indicator.transformation))
}

Model.to.model <- function(Model){
  
  model <- list(pi.bar = Model$pi.bar,
                delta.t = Model$delta.t,
                delta.c = Model$delta.c,
                Phi.Y.r = Model$Phi.Y.r,
                Theta = Model$Theta,
                Gamma.Y0.r = Model$Gamma.Y0.r, #0
                Gamma.Y1.r =  Model$Gamma.Y1.r,
                nu = Model$nu, # 1st parameter of the non centered gamma process (AGP(nu,phi,mu))
                phi = Model$phi,
                mu = Model$mu, # 3rd parameter of the AGP
                sigma.av = Model$sigma.av,
                sigma.var = Model$sigma.var,
                sigma.k3rd = Model$sigma.k3rd,
                sigma.k4th = Model$sigma.k4th
  )
  
  return(model)
  
}

#### FOR INITIAL MODEL
Make.thetas.indicator.model <- function(model, pi.bar.s=FALSE, delta.s=FALSE, 
                                        Phi.Y.s=FALSE, Theta.s=FALSE, Gamma.Y0.s=FALSE,
                                        Gamma.Y1.s=FALSE, nu.s=FALSE, phi.s=FALSE, 
                                        mu.s=FALSE, sigma.av.s=FALSE, sigma.var.s =FALSE,
                                        sigma.k3rd.s=FALSE){
  
  # Note: All the parameters should be specified in the function
  #       The function creates a matrix with the 3 columns:
  #         - The first one returns all the parameters with their original values.
  #         - The second one returns a vector with a 1 for the parameters to be estimated (0 if not).
  #         - The third one returns a vector with a 1 if the estimated parameters should be between 0 and 1
  #           And a 2 if the parameter should only be positive (0 otherwise).
  #       By default the function will only consider F to be estimated.
  
  pi.bar <- model$pi.bar
  delta <- model$delta
  Phi.Y <- model$Phi.Y
  Theta <- model$Theta
  Gamma.Y0 <- model$Gamma.Y0
  Gamma.Y1 <- model$Gamma.Y1
  nu <- model$nu
  phi <- model$phi
  mu <- model$mu
  sigma.av <- model$sigma.av[!is.na(model$sigma.av)]
  sigma.var <- model$sigma.var[!is.na(model$sigma.var)]
  sigma.k3rd <- model$sigma.k3rd[!is.na(model$sigma.k3rd)]
  
  # Conditions for pi.bar
  if(pi.bar.s == TRUE){
    
    pi.bar.i <- rep(1,length(pi.bar))
    pi.bar.t <- rep(0, length(pi.bar))
    
  } else {pi.bar.i <- rep(0, length(pi.bar))
  pi.bar.t <- rep(0, length(pi.bar))}
  
  # Conditions for delta
  if(delta.s == TRUE){
    
    delta.i <- rep(1,length(delta))
    delta.t <- rep(0, length(delta))
    
  } else {delta.i <- rep(0, length(delta))
  delta.t <- rep(0, length(delta))}
  
  # Conditions for Phi.Y
  if(Phi.Y.s == TRUE){
    
    Phi.Y.i <- rep(1, length(Phi.Y))
    Phi.Y.i[which(Phi.Y == 0)] <- 0
    Phi.Y.t <- rep(0, length(Phi.Y))
    Phi.Y.diag <- Phi.Y
    diag(Phi.Y.diag) <- "diag"
    Phi.Y.t[which(Phi.Y.diag == "diag")] <- 1
    
  } else {Phi.Y.i <- rep(0, length(Phi.Y))
  Phi.Y.t <- rep(0, length(Phi.Y))}
  
  # Conditions for Theta
  if(Theta.s == TRUE){
    
    Theta.i <- rep(1,length(Theta))
    Theta.i[which(Theta == 0)] <- 0 #Remove if wants to estimate 0.
    Theta.t <- rep(0, length(Theta))
    
  } else {Theta.i <- rep(0, length(Theta))
  Theta.t <- rep(0, length(Theta))}
  
  # Conditions for Gamma.Y0
  if(Gamma.Y0.s == TRUE){
    
    Gamma.Y0.i <- rep(1,length(Gamma.Y0))
    Gamma.Y0.i[which(Gamma.Y0 == 0)] <- 0 #Remove if wants to estimate 0.
    Gamma.Y0.t <- rep(2, length(Gamma.Y0))
    
  } else {Gamma.Y0.i <- rep(0, length(Gamma.Y0))
  Gamma.Y0.t <- rep(0, length(Gamma.Y0))}
  
  # Conditions for Gamma.Y1
  if(Gamma.Y1.s == TRUE){
    
    Gamma.Y1.i <- rep(1,length(Gamma.Y1))
    Gamma.Y1.i[which(Gamma.Y1 == 0)] <- 0 #Remove if wants to estimate 0.
    Gamma.Y1.t <- rep(2, length(Gamma.Y1))
    
  } else {Gamma.Y1.i <- rep(0, length(Gamma.Y1))
  Gamma.Y1.t <- rep(0, length(Gamma.Y1))}
  
  # Conditions for nu
  if(nu.s == TRUE){
    
    nu.i <- rep(1,length(nu))
    nu.i[which(nu == 0)] <- 0 #Remove if wants to estimate 0.
    nu.t <- rep(2, length(nu))
    
  } else {nu.i <- rep(0, length(nu))
  nu.t <- rep(0, length(nu))}
  
  # Conditions for phi
  if(phi.s == TRUE){
    
    phi.i <- rep(1, length(phi))
    phi.i[which(phi == 0)] <- 0
    phi.t <- rep(0, length(phi))
    phi.diag <- phi
    diag(phi.diag) <- "diag"
    phi.t[which(phi.diag == "diag")] <- 1
    
  } else {phi.i <- rep(0, length(phi))
  phi.t <- rep(0, length(phi))}
  
  # Conditions for mu
  if(mu.s == TRUE){
    
    mu.i <- rep(1,length(mu))
    mu.t <- rep(0, length(mu))
    
  } else {mu.i <- rep(0, length(mu))
  mu.t <- rep(0, length(mu))}
  
  # Conditions for sigma.av
  if(sigma.av.s == TRUE){
    
    sigma.av.i <- rep(1, length(sigma.av))
    sigma.av.i[which(sigma.av == 0)] <- 0
    sigma.av.t <- rep(2, length(sigma.av))
    
  } else {sigma.av.i <- rep(0, length(sigma.av))
  sigma.av.t <- rep(0, length(sigma.av))}
  
  # Conditions for sigma.var
  if(sigma.var.s == TRUE){
    
    sigma.var.i <- rep(1, length(sigma.var))
    sigma.var.i[which(sigma.var == 0)] <- 0
    sigma.var.t <- rep(2, length(sigma.var))
    
  } else {sigma.var.i <- rep(0, length(sigma.var))
  sigma.var.t <- rep(0, length(sigma.var))}
  
  # Conditions for sigma.k3rd
  if(sigma.k3rd.s == TRUE){
    
    sigma.k3rd.i <- rep(1, length(sigma.k3rd))
    sigma.k3rd.i[which(sigma.k3rd == 0)] <- 0
    sigma.k3rd.t <- rep(2, length(sigma.k3rd))
    
  } else {sigma.k3rd.i <- rep(0, length(sigma.k3rd))
  sigma.k3rd.t <- rep(0, length(sigma.k3rd))}
  
  # Create vector indicating if the parameter should be estimated 
  # and vector indicating if the parameter should be transformed.
  vector.indicator.estimate <- c(pi.bar.i, delta.i, Phi.Y.i, Theta.i, Gamma.Y0.i,
                                 Gamma.Y1.i, nu.i, phi.i, mu.i, sigma.av.i, sigma.var.i,
                                 sigma.k3rd.i)
  vector.indicator.transformation <- c(pi.bar.t, delta.t, Phi.Y.t, Theta.t, Gamma.Y0.t,
                                       Gamma.Y1.t, nu.t, phi.t, mu.t, sigma.av.t, sigma.var.t,
                                       sigma.k3rd.t)
  
  # Create list with all the parameters in matrix form
  all.thetas <- model
  
  # Create the vector of parameters
  all.thetas.vec <- na.omit(unlist(all.thetas))
  
  return(cbind(all.thetas.vec, vector.indicator.estimate,vector.indicator.transformation))
}


### OLD FUNCTION
Make.thetas.indicator <- function(all.parameters, mu.s=FALSE, F.s=TRUE, 
                                  sigma.s=FALSE, A.s=FALSE, H.s=FALSE, delta.s=FALSE){
  
  # Note: All the parameters should be specified in the function
  #       The function creates a matrix with the 3 columns:
  #         - The first one returns all the parameters with their original values.
  #         - The second one returns a vector with a 1 for the parameters to be estimated (0 if not).
  #         - The third one returns a vector with a 1 if the estimated parameters should be between 0 and 1,
  #           a 2 if the parameter should only be positive and a 3 if negative (0 otherwise). 
  #       By default the function will only consider F to be estimated.
  
  # Conditions for mu
  if(mu.s == TRUE){
    
    mu.i <- rep(1,length(mu))
    #mu.i[which(mu == 0)] <- 0
    mu.t <- rep(0, length(mu))
    
  } else {mu.i <- rep(0, length(mu))
  mu.t <- rep(0, length(mu))}
  
  # Conditions for F
  if(F.s == TRUE){
    
    F.i <- rep(1, length(F))
    F.i[which(F == 0)] <- 0
    F.t <- rep(0, length(F))
    F.diag <- F
    diag(F.diag) <- "diag"
    F.t[which(F.diag == "diag")] <- 1
    
  } else {F.i <- rep(0, length(F))
  F.t <- rep(0, length(F))}
  
  # Conditions for sigma
  if(sigma.s == TRUE){
    
    sigma.i <- rep(1, length(sigma))
    sigma.i[which(sigma == 0)] <- 0
    sigma.t <- rep(0, length(sigma))
    sigma.diag <- sigma
    diag(sigma.diag) <- "diag"
    sigma.t[which(sigma.diag == "diag")] <- 2
    
  } else {sigma.i <- rep(0, length(sigma))
  sigma.t <- rep(0, length(sigma))}
  
  # Conditions for A
  if(A.s == TRUE){
    
    A.i <- rep(1, length(A))
    #A.i[which(A == 0)] <- 0
    A.t <- rep(0, length(H))
    
  } else {A.i <- rep(0, length(A))
  A.t <- rep(0, length(A))}
  
  # Conditions for H
  if(H.s == TRUE){
    
    H.i <- rep(1, length(H))
    H.i[which(H == 0)] <- 0
    H.t <- rep(0, length(H))
    
  } else {H.i <- rep(0, length(H))
  H.t <- rep(0, length(H))}
  
  # Conditions for delta
  if(delta.s == TRUE){
    
    delta.i <- rep(1, length(delta))
    delta.i[which(delta == 0)] <- 0
    delta.t <- rep(0, length(delta))
    delta.diag <- delta
    diag(delta.diag) <- "diag"
    delta.t[which(delta.diag == "diag")] <- 2
    
  } else {delta.i <- rep(0, length(delta))
  delta.t <- rep(0, length(delta))}
  
  # Create vector indicating if the parameter should be estimated 
  # and vector indicating if the parameter should be transformed.
  vector.indicator.estimate <-c(mu.i, F.i, sigma.i, A.i, H.i, delta.i)
  vector.indicator.transformation <-c(mu.t, F.t, sigma.t, A.t, H.t, delta.t)
  
  # Create list with all the parameters in matrix form
  all.thetas <- all.parameters
  
  # Create the vector of parameters
  all.thetas.vec <- unlist(all.thetas)
  
  return(cbind(all.thetas.vec, vector.indicator.estimate,vector.indicator.transformation))
}


# ============================================
#     G. FUNCTION RETRIEVING INITIAL
#               PARAMATERS 
# ============================================


## FOR TREND CYCLE MODEL
# Function to retrieve the the initial matrices of parameters. 
Retrieve.initial.par.trend.cycle.model <- function(all.thetas.initial, n, m, q, r, nbr.horizon.max){
  
  # Note: all.thetas.initial is matrix created with the function Make.thetas.indicator. 
  #       The function gives back the initial matrices of parameters
  pi.bar <- matrix(all.thetas.initial[grep("pi.bar", rownames(as.data.frame(all.thetas.initial))) ,1],1 , r)
  delta.t <-  matrix(all.thetas.initial[grep("delta.t", rownames(as.data.frame(all.thetas.initial))) ,1], m, r)
  delta.c <-  matrix(all.thetas.initial[grep("delta.c", rownames(as.data.frame(all.thetas.initial))) ,1], m, r)
  Phi.Y.r <- matrix(all.thetas.initial[grep("Phi.Y.r", rownames(as.data.frame(all.thetas.initial))) ,1], m, m)
  Theta <- matrix(all.thetas.initial[grep("Theta", rownames(as.data.frame(all.thetas.initial))) ,1], n, q)
  Gamma.Y0.r <- matrix(all.thetas.initial[grep("Gamma.Y0.r", rownames(as.data.frame(all.thetas.initial))) ,1], m, 1)
  Gamma.Y1.r <- matrix(all.thetas.initial[grep("Gamma.Y1.r", rownames(as.data.frame(all.thetas.initial))) ,1], q, m)
  nu <- matrix(all.thetas.initial[grep("nu", rownames(as.data.frame(all.thetas.initial))) ,1], q, 1)
  phi <- matrix(all.thetas.initial[grep("phi", rownames(as.data.frame(all.thetas.initial))) ,1], q, q)
  mu <- matrix(all.thetas.initial[grep("mu", rownames(as.data.frame(all.thetas.initial))) ,1], q, 1)
  sigma.av <- array(NaN,c(nbr.horizon.max,r,4))
  position.sigma.av <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                       str_sub(rownames(as.data.frame(all.thetas.initial)), start= -2)[grep("sigma.av", rownames(as.data.frame(all.thetas.initial)))]))
  sigma.av[position.sigma.av] <- all.thetas.initial[grep("sigma.av", rownames(as.data.frame(all.thetas.initial))) ,1]
  sigma.var <- array(NaN,c(nbr.horizon.max,r,4))
  position.sigma.var <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                        str_sub(rownames(as.data.frame(all.thetas.initial)), start= -2)[grep("sigma.var", rownames(as.data.frame(all.thetas.initial)))]))
  sigma.var[position.sigma.var] <- all.thetas.initial[grep("sigma.var", rownames(as.data.frame(all.thetas.initial))) ,1]
  sigma.k3rd <- array(NaN,c(nbr.horizon.max,r,4))
  position.sigma.k3rd <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                         str_sub(rownames(as.data.frame(all.thetas.initial)), start= -2)[grep("sigma.k3rd", rownames(as.data.frame(all.thetas.initial)))]))
  sigma.k3rd[position.sigma.k3rd] <- all.thetas.initial[grep("sigma.k3rd", rownames(as.data.frame(all.thetas.initial))) ,1]
  sigma.k4th <- array(NaN,c(nbr.horizon.max,r,4))
  position.sigma.k4th <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                         str_sub(rownames(as.data.frame(all.thetas.initial)), start= -2)[grep("sigma.k4th", rownames(as.data.frame(all.thetas.initial)))]))
  sigma.k4th[position.sigma.k4th] <- all.thetas.initial[grep("sigma.k4th", rownames(as.data.frame(all.thetas.initial))) ,1]
  
  
  return(list("pi.bar" = pi.bar,
              "delta.t" = delta.t,
              "delta.c" = delta.c,
              "Phi.Y.r" = Phi.Y.r,
              "Theta" = Theta,
              "Gamma.Y0.r" = Gamma.Y0.r,
              "Gamma.Y1.r" = Gamma.Y1.r,
              "nu" = nu, # 1st parameter of the non centered gamma process (AGP(nu,phi,mu))
              "phi" = phi, # 2nd parameter of the AGP
              "mu" = mu, # 3rd parameter of the AGP
              "sigma.av" = sigma.av,
              "sigma.var" = sigma.var,
              "sigma.k3rd" = sigma.k3rd,
              "sigma.k4th" = sigma.k4th
  ))
}


## FOR INITIAL MODEL
# Function to retrieve the the initial matrices of parameters. 
Retrieve.Initial.Par.Model <- function(all.thetas.initial, n, q, r, nbr.horizon.max){
  
  # Note: all.thetas.initial is matrix created with the function Make.thetas.indicator. 
  #       The function gives back the initial matrices of parameters
  pi.bar <- matrix(all.thetas.initial[grep("pi.bar", rownames(as.data.frame(all.thetas.initial))) ,1],1 , r)
  delta <-  matrix(all.thetas.initial[grep("delta", rownames(as.data.frame(all.thetas.initial))) ,1], n, r)
  Phi.Y <- matrix(all.thetas.initial[grep("Phi.Y", rownames(as.data.frame(all.thetas.initial))) ,1], n, n)
  Theta <- matrix(all.thetas.initial[grep("Theta", rownames(as.data.frame(all.thetas.initial))) ,1], n, q)
  Gamma.Y0 <- matrix(all.thetas.initial[grep("Gamma.Y0", rownames(as.data.frame(all.thetas.initial))) ,1], n, 1)
  Gamma.Y1 <- matrix(all.thetas.initial[grep("Gamma.Y1", rownames(as.data.frame(all.thetas.initial))) ,1], q, n)
  nu <- matrix(all.thetas.initial[grep("nu", rownames(as.data.frame(all.thetas.initial))) ,1], q, 1)
  phi <- matrix(all.thetas.initial[grep("phi", rownames(as.data.frame(all.thetas.initial))) ,1], q, q)
  mu <- matrix(all.thetas.initial[grep("mu", rownames(as.data.frame(all.thetas.initial))) ,1], q, 1)
  sigma.av <- array(NaN,c(nbr.horizon.max,r,4))
  position.sigma.av <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                       str_sub(rownames(as.data.frame(all.thetas.initial)), start= -2)[grep("sigma.av", rownames(as.data.frame(all.thetas.initial)))]))
  sigma.av[position.sigma.av] <- all.thetas.initial[grep("sigma.av", rownames(as.data.frame(all.thetas.initial))) ,1]
  sigma.var <- array(NaN,c(nbr.horizon.max,r,4))
  position.sigma.var <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                        str_sub(rownames(as.data.frame(all.thetas.initial)), start= -2)[grep("sigma.var", rownames(as.data.frame(all.thetas.initial)))]))
  sigma.var[position.sigma.var] <- all.thetas.initial[grep("sigma.var", rownames(as.data.frame(all.thetas.initial))) ,1]
  sigma.k3rd <- array(NaN,c(nbr.horizon.max,r,4))
  position.sigma.k3rd <- as.numeric(gsub(".*?([0-9]+).*", "\\1", 
                                         str_sub(rownames(as.data.frame(all.thetas.initial)), start= -2)[grep("sigma.k3rd", rownames(as.data.frame(all.thetas.initial)))]))
  sigma.k3rd[position.sigma.k3rd] <- all.thetas.initial[grep("sigma.k3rd", rownames(as.data.frame(all.thetas.initial))) ,1]
  
  
  return(list("pi.bar" = pi.bar,
              "delta" = delta,
              "Phi.Y" = Phi.Y,
              "Theta" = Theta,
              "Gamma.Y0" = Gamma.Y0,
              "Gamma.Y1" = Gamma.Y1,
              "nu" = nu, # 1st parameter of the non centered gamma process (AGP(nu,phi,mu))
              "phi" = phi, # 2nd parameter of the AGP
              "mu" = mu, # 3rd parameter of the AGP
              "sigma.av" = sigma.av,
              "sigma.var" = sigma.var,
              "sigma.k3rd" = sigma.k3rd
  ))
}

# OLD ...Function to retrieve the the initial matrices of parameters. 
Retrieve.Initial.Par <- function(all.thetas.initial, r, n, k){
  
  # Note: all.thetas.initial is matrix created with the function Make.thetas.indicator. 
  #       The function gives back the initial matrices of parameters
  mu <- matrix(all.thetas.initial[grep("mu", rownames(as.data.frame(all.thetas.initial))) ,1], r, 1)
  F <-  matrix(all.thetas.initial[grep("F", rownames(as.data.frame(all.thetas.initial))) ,1], r, r)
  sigma <- matrix(all.thetas.initial[grep("sigma", rownames(as.data.frame(all.thetas.initial))) ,1], r, r)
  A <- matrix(all.thetas.initial[grep("A", rownames(as.data.frame(all.thetas.initial))) ,1], k, n)
  H <- matrix(all.thetas.initial[grep("H", rownames(as.data.frame(all.thetas.initial))) ,1], r, n)
  delta <- matrix(all.thetas.initial[grep("delta", rownames(as.data.frame(all.thetas.initial))) ,1], n, n)
  
  return(list("mu"=mu, 
              "F"=F, 
              "sigma"=sigma,
              "A"=A, 
              "H"=H, 
              "delta"=delta))
}



# ============================================
#     H. MAPPING FUNCTION AND ITS INVERSE  
# ============================================

# "Mapping" function f
# Note: To ensure that matrix F provides us with a covariance stationary process,
# we need to ensure that F is triangular (upper or lower) and that the diagonal
# elements are between 0 and 1.
# This function ensure that the diagonal elements are between 0 and 1.
Mapping.function <- function(thetas, ind.trans) { 
  
  # ind.trans: - 1 if the estimated parameters should be between 0 and 1,
  #            - 2 if the parameter should only be positive,
  #            - 3 if the parameter should only be negative,
  #            - 4 if the parameter should be be between 0 and 0.99
  #            - 0 otherwise. 
  
  # For thetas with indicator = 1; transform with exp()/(1+exp()) 
  thetas[which(ind.trans == 1)] <- exp(thetas[which(ind.trans == 1)])/(1+exp(thetas[which(ind.trans == 1)]))
  
  # For thetas with indicator = 2; transform with exp()
  thetas[which(ind.trans == 2)] <- exp(thetas[which(ind.trans == 2)])
  
  # For thetas with indicator = 3; transform with -exp()
  thetas[which(ind.trans == 3)] <- -exp(thetas[which(ind.trans == 3)])
  
  # For thetas with indicator = 4; transform with 0.99*exp()/(1+exp()) 
  thetas[which(ind.trans == 4)] <- 0.99*exp(thetas[which(ind.trans == 4)])/(1+exp(thetas[which(ind.trans == 4)]))
  
  return(thetas)
}

# Inverse "Mapping" function f
Mapping.function.inverse <- function(thetas, ind.trans) {
  
  epsilon <- 1e-30 
  
  # For thetas with indicator = 1; transform with exp()/(1+exp()) 
  if(any(1-thetas[which(ind.trans == 1)] == 0)){
    thetas[which(ind.trans == 1)] <- log(thetas[which(ind.trans == 1)]/(1-thetas[which(ind.trans == 1)] + epsilon))
  } else{
    thetas[which(ind.trans == 1)] <- log(thetas[which(ind.trans == 1)]/(1-thetas[which(ind.trans == 1)]))
  }
  
  
  # For thetas with indicator = 2; transform with exp()
  thetas[which(ind.trans == 2)] <- log(thetas[which(ind.trans == 2)])
  
  # For thetas with indicator = 3; transform with -exp()
  thetas[which(ind.trans == 3)] <- log(-1*thetas[which(ind.trans == 3)])
  
  # For thetas with indicator = 4; transform with 0.99*exp()/(1+exp()) 
  if(any(0.99*1-thetas[which(ind.trans == 4)] <= 0)){
    thetas[which(ind.trans == 4)] <- log(thetas[which(ind.trans == 4)]/(0.99*1-round(thetas[which(ind.trans == 4)],10) + epsilon))
  } else{
    thetas[which(ind.trans == 4)] <- log(thetas[which(ind.trans == 4)]/(0.99*1-thetas[which(ind.trans == 4)]))
  }
  
  return(thetas)
} 


# Recover the parameters' standard deviation
f.SE.KF.MC <- function(optim_sol, ind.trans){
  
  J <- jacobian(Mapping.function, optim_sol$par, ind.trans=ind.trans)
  
  if(is.singular.matrix(optim_sol$hessian) == TRUE){
    mat_var_bar <- matrix(NaN, dim(optim_sol$hessian)[1], dim(optim_sol$hessian)[2])
  } else {mat_var_bar <- solve(optim_sol$hessian, tol=1e-100)} 
  
  #mat_var_bar <- solve(optim_sol$hessian, tol=1e-09)
  mat_var <- J%*%mat_var_bar%*%t(J)
  SE <- sqrt(diag(mat_var))
  
  return(SE)
}


# Recover the parameters' standard deviation - More general formula
f.SE.KF <- function(optim_par, hessian, ind.trans){
  
  J <- jacobian(Mapping.function, optim_par, ind.trans=ind.trans)
  
  if(is.singular.matrix(hessian) == TRUE){
    mat_var_bar <- matrix(NaN, dim(hessian)[1], dim(hessian)[2])
  } else {mat_var_bar <- solve(hessian, tol=1e-100)} 
  
  #mat_var_bar <- solve(optim_sol$hessian, tol=1e-09)
  mat_var <- J%*%mat_var_bar%*%t(J)
  SE <- sqrt(diag(mat_var))
  
  return(SE)
}


# ============================================
#     I. FUNCTION PREPARING KF PARAMETERS
#                 AND VARIABLES
# ============================================

# Function to conditionally bind
conditional_cbind <- function(mat1, mat2) {
  if (nrow(mat2) == 0 && ncol(mat2) == 0) {
    return(mat1)  # Return the original matrix if it's 0x0
  } else {
    return(cbind(mat1, mat2))  # Combine otherwise
  }
}

# Function that prepare the parameters and the variables for the KF (model)
# It computes inter-alia the loadings.
prepare.KF.model <- function(Model, observables){
  # Model : is a list that contains all the info
  # observables : observables data
  
  # Convert observables in matrix form (to apply KF)
  observables <- as.matrix(observables)
  
  # Fill NA values with the first and last available value
  observables.full <- na.locf(observables)
  observables.full <- na.locf(observables.full, fromLast = TRUE)
  
  # Store dimensions
  T <- dim(observables)[1] # number of dates
  m <- dim(observables)[2] # number of observed variables
  r <- length(c(Model$pi.bar)) # number of areas or key variables
  k <- length(Model$w) # number of additional measurement equations, if any.
  n <- Model$n
  q <- Model$q
  
  # Dynamics of latent variables (i.e. transition equations):
  mu <- matrix(Model$mu.X, n+q,1)
  F    <- Model$Phi.X
  sigma <- NULL
  
  # Measurement equations:
  A <- NULL #matrix(0, 1, m, byrow = TRUE)
  X <- matrix(1,T,1) # assume it is a constant
  HH <- NULL
  delta <- NULL
  
  # Vector of delta (inflation eq) with 0 for the observables of the variances
  Gamma <- rbind(Model$delta,matrix(0,Model$q,r))
  
  S.n <- make.S.p(n)
  S.q <- make.S.p(q)
  S.q.tilde <- make.S.p.tilde(q)
  
  # Prepare measurement SPF equations (level and variance, ie n+q per area or key variables)
  ## Preparation of the matrices for the Kalman Filter (done country by country)
  for(area.var in 1:r){
    # Compute  the different matrices
    ## Add country by country
    A <- cbind(A, Model$pi.bar[area.var])
    
    ## Add component related to realized inflation
    ### Condition if you we considered q.o.q inflation or y.o.y inflation.
    HH <- rbind(HH,matrix(c(Model$delta.q[,area.var],rep(0,Model$q)),1,n+q)) # realized inflation
    
    ## Get the Horizon of interest
    H.aux <- H[,area.var] # get horizons for area or key variables r
    H.aux <- H.aux[!is.na(H.aux)]
    
    # Compute the loadings in order to then fill matrix HH
    #loadings <- compute.loadings.c(c(Gamma[,area.var]),Model,H.aux,S.n,S.q,S.q.tilde,indic.5y.in.Xy = Indic.5y.in.Xy)
    loadings <- compute_loadings(as.matrix(c(Gamma[,area.var])),Model,H.aux,S.n,S.q,S.q.tilde,indic_5y_in_Xy = Indic.5y.in.Xy, freq=freq)
    
    
    # Compute A
    A <- cbind(A,
               # First type of inflation:
               t(matrix(4*Model$pi.bar[area.var] + c(loadings$a[1,which(select.inflation.types[,area.var,1]>=1),1]),sum(select.inflation.types[,area.var,1]>=1,na.rm=TRUE),1)), # corresponds to point estimates
               t(matrix(c(loadings$alpha[1,which(select.inflation.types[,area.var,1]>=11),1]),sum(select.inflation.types[,area.var,1]>=11,na.rm=TRUE),1)), # corresponds to variances
               t(matrix(c(loadings$alpha.dot.dot[1,which(select.inflation.types[,area.var,1]>=111),1]),sum(select.inflation.types[,area.var,1]>=111,na.rm=TRUE),1)), # corresponds to 3rd order cumulants
               t(matrix(c(loadings$alpha.dot.dot.dot[1,which(select.inflation.types[,area.var,1]>=1111),1]),sum(select.inflation.types[,area.var,1]>=1111,na.rm=TRUE),1)), # corresponds to 4th order cumulants
               # Second type of inflation:
               t(matrix(4*Model$pi.bar[area.var] + c(loadings$a[1,which(select.inflation.types[,area.var,2]>=1),2]),sum(select.inflation.types[,area.var,2]>=1,na.rm=TRUE),1)), # corresponds to point estimates
               t(matrix(c(loadings$alpha[1,which(select.inflation.types[,area.var,2]>=11),2]),sum(select.inflation.types[,area.var,2]>=11,na.rm=TRUE),1)), # corresponds to variances
               t(matrix(c(loadings$alpha.dot.dot[1,which(select.inflation.types[,area.var,2]>=111),2]),sum(select.inflation.types[,area.var,2]>=111,na.rm=TRUE),1)), # corresponds to 3rd order cumulants
               t(matrix(c(loadings$alpha.dot.dot.dot[1,which(select.inflation.types[,area.var,2]>=1111),2]),sum(select.inflation.types[,area.var,2]>=1111,na.rm=TRUE),1)), # corresponds to 4th order cumulants
               # Third type of inflation:
               t(matrix(4*Model$pi.bar[area.var] + c(loadings$a[1,which(select.inflation.types[,area.var,3]>=1),3]),sum(select.inflation.types[,area.var,3]>=1,na.rm=TRUE),1)), # corresponds to point estimates
               t(matrix(c(loadings$alpha[1,which(select.inflation.types[,area.var,3]>=11),3]),sum(select.inflation.types[,area.var,3]>=11,na.rm=TRUE),1)), # corresponds to variances
               t(matrix(c(loadings$alpha.dot.dot[1,which(select.inflation.types[,area.var,3]>=111),3]),sum(select.inflation.types[,area.var,3]>=111,na.rm=TRUE),1)), # corresponds to 3rd order cumulants
               t(matrix(c(loadings$alpha.dot.dot.dot[1,which(select.inflation.types[,area.var,3]>=1111),3]),sum(select.inflation.types[,area.var,3]>=1111,na.rm=TRUE),1)), # corresponds to 4th order cumulants
               # Fourth type of inflation:
               t(matrix(4*Model$pi.bar[area.var] + c(loadings$a[1,which(select.inflation.types[,area.var,4]>=1),4]),sum(select.inflation.types[,area.var,4]>=1,na.rm=TRUE),1)), # corresponds to point estimates
               t(matrix(c(loadings$alpha[1,which(select.inflation.types[,area.var,4]>=11),4]),sum(select.inflation.types[,area.var,4]>=11,na.rm=TRUE),1)), # corresponds to variances
               t(matrix(c(loadings$alpha.dot.dot[1,which(select.inflation.types[,area.var,4]>=111),4]),sum(select.inflation.types[,area.var,4]>=111,na.rm=TRUE),1)), # corresponds to 3rd order cumulants
               t(matrix(c(loadings$alpha.dot.dot.dot[1,which(select.inflation.types[,area.var,4]>=1111),4]),sum(select.inflation.types[,area.var,4]>=1111,na.rm=TRUE),1)) # corresponds to 4th order cumulants
    )
    
    # Compute HH
    HH <- rbind(HH,
                # First type of inflation:
                t(loadings$b[,which(select.inflation.types[,area.var,1]>=1),1]), # first type of inflation, point estimates
                t(loadings$beta[,which(select.inflation.types[,area.var,1]>=11),1]), # variances
                t(loadings$beta.dot.dot[,which(select.inflation.types[,area.var,1]>=111),1]), # 3rd order cumulants
                t(loadings$beta.dot.dot.dot[,which(select.inflation.types[,area.var,1]>=1111),1]), # 4th order cumulants
                # Second type of inflation:
                t(loadings$b[,which(select.inflation.types[,area.var,2]>=1),2]), # Second type of inflation, point estimates
                t(loadings$beta[,which(select.inflation.types[,area.var,2]>=11),2]), # variances
                t(loadings$beta.dot.dot[,which(select.inflation.types[,area.var,2]>=111),2]), # 3rd order cumulants
                t(loadings$beta.dot.dot.dot[,which(select.inflation.types[,area.var,2]>=1111),2]), # 4th order cumulants
                # Third type of inflation:
                t(loadings$b[,which(select.inflation.types[,area.var,3]>=1),3]), # Third type of inflation, point estimates
                t(loadings$beta[,which(select.inflation.types[,area.var,3]>=11),3]), # variances
                t(loadings$beta.dot.dot[,which(select.inflation.types[,area.var,3]>=111),3]), # 3rd order cumulants
                t(loadings$beta.dot.dot.dot[,which(select.inflation.types[,area.var,3]>=1111),3]), # 4th order cumulants
                # Fourth type of inflation:
                t(loadings$b[,which(select.inflation.types[,area.var,4]>=1),4]), # Fourth type of inflation, point estimates
                t(loadings$beta[,which(select.inflation.types[,area.var,4]>=11),4]), # variances
                t(loadings$beta.dot.dot[,which(select.inflation.types[,area.var,4]>=111),4]), # 3rd order cumulants
                t(loadings$beta.dot.dot.dot[,which(select.inflation.types[,area.var,4]>=1111),4]) # 4th order cumulants
    )
    
    # Compute delta
    ## Augment vector of measurement-eq standard deviations:
    # delta <- c(delta,
    #            #stdv.measur$sigma.pi,
    #            if(area.var==1){0.15}else{.2}, # realized inflation .1 before, 0.2 for US
    #            # First type of inflation:
    #            Model$sigma.av[!is.na(Model$sigma.av[,area.var,1]),area.var,1], # point estimates equations
    #            Model$sigma.var[!is.na(Model$sigma.var[,area.var,1]),area.var,1], # variances equations
    #            Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,1]),area.var,1], # 3rd order cumulants equations
    #            Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,1]),area.var,1], # 4th order cumulants equations
    #            # Second type of inflation:
    #            Model$sigma.av[!is.na(Model$sigma.av[,area.var,2]),area.var,2], # point estimates equations
    #            Model$sigma.var[!is.na(Model$sigma.var[,area.var,2]),area.var,2], # variances equations
    #            Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,2]),area.var,2], # 3rd order cumulants equations
    #            Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,2]),area.var,2], # 4th order cumulants equations
    #            # Third type of inflation:
    #            Model$sigma.av[!is.na(Model$sigma.av[,area.var,3]),area.var,3], # point estimates equations
    #            Model$sigma.var[!is.na(Model$sigma.var[,area.var,3]),area.var,3], # variances equations
    #            Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3], #  3rd order cumulants equations
    #            Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,3]),area.var,3], # 4th order cumulants equations
    #            # Fourth type of inflation:
    #            Model$sigma.av[!is.na(Model$sigma.av[,area.var,4]),area.var,4], # point estimates equations
    #            Model$sigma.var[!is.na(Model$sigma.var[,area.var,4]),area.var,4], # variances equations
    #            Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,4]),area.var,4], # 3rd order cumulants equations
    #            Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,4]),area.var,4] # 4th order cumulants equations
    # ) #comment for new delta
    
    # add for new delta
    #delta <- cbind(delta,rep(if(area.var==1){0.15}else{.2},T))
    
    if(area=="US"){
      delta <- cbind(delta,rep(if(area.var==1){0.15}else{.2},T))
    } else{
      delta <- cbind(delta,rep(if(area.var==1){0.1}else{.25},T))#0.15, 0.35
    }
    
    # First type of inflation:
    delta <- conditional_cbind(delta,matrix(rep(Model$sigma.av[!is.na(Model$sigma.av[,area.var,1]),area.var,1], each=T), ncol= length(Model$sigma.av[!is.na(Model$sigma.av[,area.var,1]),area.var,1]))) # point estimates equations
    delta <- conditional_cbind(delta,matrix(rep(Model$sigma.var[!is.na(Model$sigma.var[,area.var,1]),area.var,1], each=T), ncol= length(Model$sigma.var[!is.na(Model$sigma.var[,area.var,1]),area.var,1]))) # variances equations
    #delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,1]),area.var,1], each=T), ncol= length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,1]),area.var,1]))) # 3rd order cumulants equations
    #delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,1]),area.var,1], each=T), ncol= length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,1]),area.var,1]))) # 4th order cumulants equations

    start.indic.k3rd <- length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var-1,1]),area.var-1,1])+1 # Calculate the start of the range
    end.indic.k3rd <-  length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,1]),area.var,1])*area.var            # Calculate the end of the range
    alpha.k3rd <- Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,1]),area.var,1] / colMeans(as.matrix(observables[,grep("var", colnames(observables))[start.indic.k3rd:end.indic.k3rd]]^(3/2)), na.rm=T)
    if(is.na(sum(alpha.k3rd))){
      delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,1]),area.var,1], each=T), ncol= length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,1]),area.var,1]))) # 3rd order cumulants equations
    } else{
      delta <- conditional_cbind(delta,observables.full[,grep("var", colnames(observables.full))[start.indic.k3rd:end.indic.k3rd]]^(3/2) * matrix(rep(alpha.k3rd, each=T), ncol= length(alpha.k3rd))) # 3rd order cumulants equations
    }
    
    start.indic.k4th <- length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var-1,1]),area.var-1,1])+1 # Calculate the start of the range
    end.indic.k4th <-  length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,1]),area.var,1])*area.var            # Calculate the end of the range
    alpha.k4th <- Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,1]),area.var,1] / colMeans(as.matrix(observables[,grep("var", colnames(observables))[start.indic.k4th:end.indic.k4th]]^(2)), na.rm=T)
    if(is.na(sum(alpha.k4th))){
      delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,1]),area.var,1], each=T), ncol= length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,1]),area.var,1]))) # 4th order cumulants equations
    } else{
      delta <- conditional_cbind(delta,observables.full[,grep("var", colnames(observables.full))[start.indic.k4th:end.indic.k4th]]^(2) * matrix(rep(alpha.k4th, each=T), ncol= length(alpha.k4th))) # 4th order cumulants equations
    }
    
    # Second type of inflation:
    delta <- conditional_cbind(delta,matrix(rep(Model$sigma.av[!is.na(Model$sigma.av[,area.var,2]),area.var,2], each=T), ncol= length(Model$sigma.av[!is.na(Model$sigma.av[,area.var,2]),area.var,2]))) # point estimates equations
    delta <- conditional_cbind(delta,matrix(rep(Model$sigma.var[!is.na(Model$sigma.var[,area.var,2]),area.var,2], each=T), ncol= length(Model$sigma.var[!is.na(Model$sigma.var[,area.var,2]),area.var,2]))) # variances equations
    #delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,2]),area.var,2], each=T), ncol= length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,2]),area.var,2]))) # 3rd order cumulants equations
    #delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,2]),area.var,2], each=T), ncol= length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,2]),area.var,2]))) # 4th order cumulants equations

    start.indic.k3rd <- length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var-1,2]),area.var-1,2])+1 # Calculate the start of the range
    end.indic.k3rd <-  length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,2]),area.var,2])*area.var            # Calculate the end of the range
    alpha.k3rd <- Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,2]),area.var,2] / colMeans(as.matrix(observables[,grep("var", colnames(observables))[start.indic.k3rd:end.indic.k3rd]]^(3/2)), na.rm=T)
    if(is.na(sum(alpha.k3rd))){
      delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,2]),area.var,2], each=T), ncol= length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,2]),area.var,2]))) # 3rd order cumulants equations
    } else{
      delta <- conditional_cbind(delta,observables.full[,grep("var", colnames(observables.full))[start.indic.k3rd:end.indic.k3rd]]^(3/2) * matrix(rep(alpha.k3rd, each=T), ncol= length(alpha.k3rd))) # 3rd order cumulants equations
    }
    
    start.indic.k4th <- length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var-1,2]),area.var-1,2])+1 # Calculate the start of the range
    end.indic.k4th <-  length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,2]),area.var,2])*area.var            # Calculate the end of the range
    alpha.k4th <- Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,2]),area.var,2] / colMeans(as.matrix(observables[,grep("var", colnames(observables))[start.indic.k4th:end.indic.k4th]]^(2)), na.rm=T)
    if(is.na(sum(alpha.k4th))){
      delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,2]),area.var,2], each=T), ncol= length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,2]),area.var,2]))) # 4th order cumulants equations
    } else{
      delta <- conditional_cbind(delta,observables.full[,grep("var", colnames(observables.full))[start.indic.k4th:end.indic.k4th]]^(2) * matrix(rep(alpha.k4th, each=T), ncol= length(alpha.k4th))) # 4th order cumulants equations
    }
    
    # Third type of inflation:
    delta <- conditional_cbind(delta,matrix(rep(Model$sigma.av[!is.na(Model$sigma.av[,area.var,3]),area.var,3], each=T), ncol= length(Model$sigma.av[!is.na(Model$sigma.av[,area.var,3]),area.var,3]))) # point estimates equations
    delta <- conditional_cbind(delta,matrix(rep(Model$sigma.var[!is.na(Model$sigma.var[,area.var,3]),area.var,3], each=T), ncol= length(Model$sigma.var[!is.na(Model$sigma.var[,area.var,3]),area.var,3]))) # variances equations
    #delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3], each=T), ncol= length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3]))) # 3rd order cumulants equations
    #delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,3]),area.var,3], each=T), ncol= length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,3]),area.var,3]))) # 4th order cumulants equations
    
    start.indic.k3rd <- length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var-1,3]),area.var-1,3])+1 # Calculate the start of the range
    end.indic.k3rd <-  length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3])*area.var            # Calculate the end of the range
    alpha.k3rd <- Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3] / colMeans(as.matrix(observables[,grep("var", colnames(observables))[start.indic.k3rd:end.indic.k3rd]]^(3/2)), na.rm=T)
    if(is.na(sum(alpha.k3rd))){
      delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3], each=T), ncol= length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3]))) # 3rd order cumulants equations
    } else{
      delta <- conditional_cbind(delta,observables.full[,grep("var", colnames(observables.full))[start.indic.k3rd:end.indic.k3rd]]^(3/2) * matrix(rep(alpha.k3rd, each=T), ncol= length(alpha.k3rd))) # 3rd order cumulants equations
    }

    start.indic.k4th <- length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var-1,3]),area.var-1,3])+1 # Calculate the start of the range
    end.indic.k4th <-  length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,3]),area.var,3])*area.var            # Calculate the end of the range
    alpha.k4th <- Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,3]),area.var,3] / colMeans(as.matrix(observables[,grep("var", colnames(observables))[start.indic.k4th:end.indic.k4th]]^(2)), na.rm=T)
    if(is.na(sum(alpha.k4th))){
      delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,3]),area.var,3], each=T), ncol= length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,3]),area.var,3]))) # 4th order cumulants equations
    } else{
      delta <- conditional_cbind(delta,observables.full[,grep("var", colnames(observables.full))[start.indic.k4th:end.indic.k4th]]^(2) * matrix(rep(alpha.k4th, each=T), ncol= length(alpha.k4th))) # 4th order cumulants equations
    }
    
      
    #observables.full[,grep("k3rd", colnames(observables.full))[(1:length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3]))*area.var]] * matrix(rep(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3], each=T), ncol= length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,3]),area.var,3]))

    # Fourth type of inflation:
    delta <- conditional_cbind(delta,matrix(rep(Model$sigma.av[!is.na(Model$sigma.av[,area.var,4]),area.var,4], each=T), ncol= length(Model$sigma.av[!is.na(Model$sigma.av[,area.var,4]),area.var,4]))) # point estimates equations
    delta <- conditional_cbind(delta,matrix(rep(Model$sigma.var[!is.na(Model$sigma.var[,area.var,4]),area.var,4], each=T), ncol= length(Model$sigma.var[!is.na(Model$sigma.var[,area.var,4]),area.var,1]))) # variances equations
    #delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,4]),area.var,4], each=T), ncol= length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,4]),area.var,4]))) # 3rd order cumulants equations
    #delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,4]),area.var,4], each=T), ncol= length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,4]),area.var,4]))) # 4th order cumulants equations

    start.indic.k3rd <- length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var-1,4]),area.var-1,4])+1 # Calculate the start of the range
    end.indic.k3rd <-  length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,4]),area.var,4])*area.var            # Calculate the end of the range
    alpha.k3rd <- Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,4]),area.var,4] / colMeans(as.matrix(observables[,grep("var", colnames(observables))[start.indic.k3rd:end.indic.k3rd]]^(3/2)), na.rm=T)
    if(is.na(sum(alpha.k3rd))){
      delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,4]),area.var,4], each=T), ncol= length(Model$sigma.k3rd[!is.na(Model$sigma.k3rd[,area.var,4]),area.var,4]))) # 3rd order cumulants equations
    } else{
      delta <- conditional_cbind(delta,observables.full[,grep("var", colnames(observables.full))[start.indic.k3rd:end.indic.k3rd]]^(3/2) * matrix(rep(alpha.k3rd, each=T), ncol= length(alpha.k3rd))) # 3rd order cumulants equations
    }
    
    start.indic.k4th <- length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var-1,4]),area.var-1,4])+1 # Calculate the start of the range
    end.indic.k4th <-  length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,4]),area.var,4])*area.var            # Calculate the end of the range
    alpha.k4th <- Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,4]),area.var,4] / colMeans(as.matrix(observables[,grep("var", colnames(observables))[start.indic.k4th:end.indic.k4th]]^(2)), na.rm=T)
    if(is.na(sum(alpha.k4th))){
      delta <- conditional_cbind(delta,matrix(rep(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,4]),area.var,4], each=T), ncol= length(Model$sigma.k4th[!is.na(Model$sigma.k4th[,area.var,4]),area.var,4]))) # 4th order cumulants equations
    } else{
      delta <- conditional_cbind(delta,observables.full[,grep("var", colnames(observables.full))[start.indic.k4th:end.indic.k4th]]^(2) * matrix(rep(alpha.k4th, each=T), ncol= length(alpha.k4th))) # 4th order cumulants equations
    }
    
    }
  
  # Add two observables for cycle
  if(indic.cycle=="TRUE"){
    
    A <- cbind(A, 0, 0)
    HH <- rbind(HH, matrix(c(Model$delta.c[,area.var-1],rep(0,Model$n +Model$q-Model$m)),1,n+q),
                matrix(c(Model$delta.c[,area.var],rep(0,Model$n +Model$q-Model$m)),1,n+q))
    #delta <- c(delta, .5, .5) #1 #comment for new delta
    delta <- cbind(delta, rep(.5,T), rep(.5,T)) #add for new delta
    
    # A <- cbind(A, 0)
    # HH <- rbind(HH,matrix(c(Model$delta.c[,area.var],rep(0,Model$n +Model$q-Model$m)),1,n+q))
    # delta <- c(delta, .25) #1
    
  }
  
  # Indicator Position
  indic_pos <- c(rep(0,n),rep(1,q)) # Indicates that the z must be positive (will be used in the modified Kalman filter)
  
  # B- Prepare measurement for additional factors (there are k of them)
  # ====================================
  if(k>0){
    A <- cbind(A,matrix(1,1,1)%*%matrix(Model$w,nrow=1))
    HH <- rbind(HH,Model$W)
    delta <- list(sigmas = c(delta,Model$sigma.addit),
                  indices.of.addit.var = Model$indices.of.addit.var)
  }
  
  # Compute unconditional moment of X
  unc.mom.X <- unc.moments.X(Model,S.n,S.q)
  ## Expectation
  xi.00 <- as.matrix(unc.mom.X$X.mom.1)
  ## Variance
  P.00 <- unc.mom.X$X.mom.2
  
  all.parameters <- list( "mu"=mu, 
                          "F"=F, 
                          "sigma"=Model, 
                          "A"=A, 
                          "H"=t(HH), 
                          #"delta"=diag(delta) # comment for new delta
                          "delta"=delta # add for new delta
                          )
  
  Y <- observables
  
  return(list("all.parameters" = all.parameters,
              "Y"=Y,
              "X"=X,
              "xi.00" = xi.00 ,
              "P.00" = P.00,
              "indic_pos" = indic_pos
  ))
  
}

# ============================================
#      J.  LOG LIKELIHOOD FUNCTION 
#   (ALL THE PARAMETERS CAN BE ESTIMATED)
# ============================================

fit.log.lik.trend.cycle.joint.model.5.7 <- function(thetas, estimated.Model, n, m.Y, q, r, nbr.horizon.max, var.type, all_mat) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using the mapping function 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  
  estimated.Model.list <- Retrieve.initial.par.trend.cycle.model(all.thetas.inter, n, m.Y, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by -theta[1,1] and theta[2,4] by -theta[2,3]
  estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]
  estimated.Model.list$Theta[2,4] <- - estimated.Model.list$Theta[2,3]
  
  # Replace element nu[2,1] by nu[1,1], nu[4,1] by nu[3,1] and nu[6,1] by nu[5,1]
  estimated.Model.list$nu[2,1] <- estimated.Model.list$nu[1,1]
  estimated.Model.list$nu[4,1] <- estimated.Model.list$nu[3,1]
  estimated.Model.list$nu[6,1] <- estimated.Model.list$nu[5,1]
  
  # Replace element phi[2,2] by phi[1,1], phi[4,4] by phi[3,3] and phi[5,5] by phi[4,4]
  estimated.Model.list$phi[2,2] <- estimated.Model.list$phi[1,1]
  estimated.Model.list$phi[4,4] <- estimated.Model.list$phi[3,3]
  estimated.Model.list$phi[6,6] <- estimated.Model.list$phi[5,5]
  
  # Replace element phi[5,2] by phi[5,1] and phi[6,4] by phi[6,3]
  estimated.Model.list$phi[5,2] <- estimated.Model.list$phi[5,1]
  estimated.Model.list$phi[6,4] <- estimated.Model.list$phi[6,3]
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type, all_mat)
  Model <- estimated.Model.list
  
  # Generate matrices used by KF.
  KF.load <- prepare.KF.model(Model, observables)
  
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
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  # return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00, S="default",
  #                       indic.pos.z=KF.load$indic_pos)$log.lik)
  # return(-Kalman_filter_cpp(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00,
  #                       indic_pos_z=KF.load$indic_pos)$log.lik)
  return(-KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$loglik)
} 


fit.log.lik.trend.cycle.joint.model <- function(thetas, estimated.Model, n, m.Y, q, r, nbr.horizon.max, var.type, all_mat) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using the mapping function 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  
  estimated.Model.list <- Retrieve.initial.par.trend.cycle.model(all.thetas.inter, n, m.Y, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by -theta[1,1] and theta[2,4] by -theta[2,3]
  estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]
  estimated.Model.list$Theta[2,4] <- - estimated.Model.list$Theta[2,3]
  
  # Replace element nu[2,1] by nu[1,1] and nu[4,1] by nu[3,1]
  estimated.Model.list$nu[2,1] <- estimated.Model.list$nu[1,1]
  estimated.Model.list$nu[4,1] <- estimated.Model.list$nu[3,1]
  
  # Replace element phi[2,2] by phi[1,1] and phi[4,4] by phi[3,3]
  estimated.Model.list$phi[2,2] <- estimated.Model.list$phi[1,1]
  estimated.Model.list$phi[4,4] <- estimated.Model.list$phi[3,3]
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type, all_mat)
  Model <- estimated.Model.list
  
  # Generate matrices used by KF.
  KF.load <- prepare.KF.model(Model, observables)
  
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
  
  penalty.sign <- sum(abs(sign(Model$delta.c[,1]) - sign(Model$delta.t[,1])) > 1)*10000
  
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  # return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00, S="default",
  #                       indic.pos.z=KF.load$indic_pos)$log.lik)
  # return(-Kalman_filter_cpp(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00,
  #                       indic_pos_z=KF.load$indic_pos)$log.lik)
  return(-KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$loglik + penalty.sign)
} 



fit.log.lik.trend.cycle.joint.model.4.2 <- function(thetas, estimated.Model, n, m.Y, q, r, nbr.horizon.max, var.type, all_mat) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using the mapping function 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  
  estimated.Model.list <- Retrieve.initial.par.trend.cycle.model(all.thetas.inter, n, m.Y, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by -theta[1,1] and theta[2,4] by -theta[2,3]
  #estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]
  #estimated.Model.list$Theta[2,4] <- - estimated.Model.list$Theta[2,3]
  
  # Replace element nu[2,1] by nu[1,1] and nu[4,1] by nu[3,1]
  #estimated.Model.list$nu[2,1] <- estimated.Model.list$nu[1,1]
  #estimated.Model.list$nu[4,1] <- estimated.Model.list$nu[3,1]
  
  # Replace element phi[2,2] by phi[1,1] and phi[4,4] by phi[3,3]
  #estimated.Model.list$phi[2,2] <- estimated.Model.list$phi[1,1]
  #estimated.Model.list$phi[4,4] <- estimated.Model.list$phi[3,3]
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type, all_mat)
  Model <- estimated.Model.list
  
  # Generate matrices used by KF.
  KF.load <- prepare.KF.model(Model, observables)
  
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
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  # return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00, S="default",
  #                       indic.pos.z=KF.load$indic_pos)$log.lik)
  # return(-Kalman_filter_cpp(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00,
  #                       indic_pos_z=KF.load$indic_pos)$log.lik)
  
  penalty.sign <- sum((sign(Model$delta.c[,1])!=sign(Model$delta.t[,1])))*10000 + sum((sign(Model$delta.c[,2])!=sign(Model$delta.t[,2])))*10000
  
  
  # Penalty on the cycle => should be around -0.5;0.5 otherwise penalty
  if(indic.cycle.use=="FALSE"){
    
    abs.diff.expectation <- abs(rep(0, 2) - rowMeans(KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$Obs.updated)[dim(observables.with.dates)[2] - c(1,0)-1]) # Penalty in order that the mean of computed model based 4th cum is close to zero
    
    penalty.cycle.expectation <- ifelse(abs.diff.expectation >= 0.5, abs.diff.expectation^2*10, 0)
  } else{
    penalty.cycle.expectation <- 0
  }
  
  
  ## Penlaty on 4th moment
  
  if(indic.4th=="TRUE"){
    ### extract a.h and b.h from 4th cum
    a.h <- as.matrix(KF.load$all.parameters$A[grep("k4th", colnames(observables))])
    b.h <- KF.load$all.parameters$H[,grep("k4th", colnames(observables))]
    
    ### compute observed vs modelled mean and variance of the 4th cum
    moments.4th.cum.model <- compute.two.first.moments.loadings(Model, a.h, b.h)
    expectation.4th.observables <- rep(0, length(moments.4th.cum.model$Cum.mom.1)) #apply(observables.with.dates[,grep("k4th", colnames(observables.with.dates))], 2, mean, na.rm=TRUE)
    var.4th.observables <- apply(observables.with.dates[,grep("k4th", colnames(observables.with.dates))], 2, var, na.rm=TRUE)
    
    ### compute absolute differences between observed and modelled mean and variance of the 4th cum
    abs.diff.expectation <- abs(expectation.4th.observables - moments.4th.cum.model$Cum.mom.1) # Penalty in order that the unconditional expectation of model based 4th cum is close to zero
    abs.diff.expectation.2nd <- abs(expectation.4th.observables - rowMeans(KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$Obs.updated)[(grep("k4th", colnames(observables.with.dates))-1)]) # Penalty in order that the mean of computed model based 4th cum is close to zero
    abs.diff.var <- abs(var.4th.observables - moments.4th.cum.model$Cum.mom.2)
    
    ### penalty 
    penalty.4th.cum.expectation <- ifelse(abs.diff.expectation >= 0.5, abs.diff.expectation^2*10, 0)
    penalty.4th.cum.expectation.2nd <- ifelse(abs.diff.expectation.2nd >= 0.5, abs.diff.expectation.2nd^2*10, 0)
    penalty.4th.cum.var <- ifelse(abs.diff.var >= var.4th.observables*1/3, abs.diff.var^2*10, 0)
  } else{
    penalty.4th.cum.expectation <- 0
    penalty.4th.cum.expectation.2nd <- 0
    penalty.4th.cum.var <- 0
  }
  
  
  return(-KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$loglik + penalty.sign + sum(penalty.4th.cum.expectation) + sum(penalty.4th.cum.expectation.2nd) + sum(penalty.4th.cum.var) + sum(penalty.cycle.expectation))
} 


fit.log.lik.trend.cycle.joint.model.3.1 <- function(thetas, estimated.Model, n, m.Y, q, r, nbr.horizon.max, var.type, all_mat) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using the mapping function 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  
  estimated.Model.list <- Retrieve.initial.par.trend.cycle.model(all.thetas.inter, n, m.Y, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by -theta[1,1] and theta[2,4] by -theta[2,3]
  
  if(q>1){
    estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]
    if(q>3){
      estimated.Model.list$Theta[2,4] <- - estimated.Model.list$Theta[2,3]
    }
  }
  
  # Replace element nu[2,1] by nu[1,1] and nu[4,1] by nu[3,1]
  #estimated.Model.list$nu[2,1] <- estimated.Model.list$nu[1,1]
  #estimated.Model.list$nu[4,1] <- estimated.Model.list$nu[3,1]
  
  # Replace element phi[2,2] by phi[1,1] and phi[4,4] by phi[3,3]
  #estimated.Model.list$phi[2,2] <- estimated.Model.list$phi[1,1]
  #estimated.Model.list$phi[4,4] <- estimated.Model.list$phi[3,3]
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type, all_mat)
  Model <- estimated.Model.list
  
  # Generate matrices used by KF.
  KF.load <- prepare.KF.model(Model, observables)
  
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
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  # return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00, S="default",
  #                       indic.pos.z=KF.load$indic_pos)$log.lik)
  # return(-Kalman_filter_cpp(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00,
  #                       indic_pos_z=KF.load$indic_pos)$log.lik)
  
  #penalty.sign <- sum((sign(Model$delta.c[,1])!=sign(Model$delta.t[,1])))*10000 + sum((sign(Model$delta.c[,2])!=sign(Model$delta.t[,2])))*10000
  penalty.sign <- sum(abs(sign(Model$delta.c[,1]) - sign(Model$delta.t[,1])) > 1)*10000 + sum(abs(sign(Model$delta.c[,2]) - sign(Model$delta.t[,2])) > 1)*10000
  
  penalty.size <- (1+Model$delta.c[2:(m.Y+1)]^2)/abs(Model$delta.c[2:(m.Y+1)])
  
  # Penalty on the cycle => should be around -0.5;0.5 otherwise penalty
  if(indic.cycle.use=="FALSE"){
    
    abs.diff.expectation <- abs(rep(0, 2) - rowMeans(KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$Obs.updated)[dim(observables.with.dates)[2] - c(1,0)-1]) # Penalty in order that the mean of computed model based 4th cum is close to zero
    
    penalty.cycle.expectation <- ifelse(abs.diff.expectation >= 0.5, abs.diff.expectation^2*100, 0)
  } else{
    penalty.cycle.expectation <- 0
  }
  
  
  ## Penlaty on 4th moment
  
  # if(indic.4th=="TRUE"){
  #   ### extract a.h and b.h from 4th cum
  #   a.h <- as.matrix(KF.load$all.parameters$A[grep("k4th", colnames(observables))])
  #   b.h <- KF.load$all.parameters$H[,grep("k4th", colnames(observables))]
  #   
  #   ### compute observed vs modelled mean and variance of the 4th cum
  #   moments.4th.cum.model <- compute.two.first.moments.loadings(Model, a.h, b.h)
  #   expectation.4th.observables <- apply(observables.with.dates[,grep("k4th", colnames(observables.with.dates))], 2, mean, na.rm=TRUE)
  #   #var.4th.observables <- apply(observables.with.dates[,grep("k4th", colnames(observables.with.dates))], 2, var, na.rm=TRUE)
  #   
  #   ### compute absolute differences between observed and modelled mean and variance of the 4th cum
  #   abs.diff.expectation <- abs(expectation.4th.observables - moments.4th.cum.model$Cum.mom.1) # Penalty in order that the unconditional expectation of model based 4th cum is close to zero
  #   abs.diff.expectation.2nd <- abs(expectation.4th.observables - rowMeans(KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$Obs.updated)[(grep("k4th", colnames(observables.with.dates))-1)]) # Penalty in order that the mean of computed model based 4th cum is close to zero
  #   #abs.diff.var <- abs(var.4th.observables - moments.4th.cum.model$Cum.mom.2)
  #   
  #   ### penalty 
  #   penalty.4th.cum.expectation <- ifelse(abs.diff.expectation >= 1, abs.diff.expectation^2*10, 0)
  #   penalty.4th.cum.expectation.2nd <- ifelse(abs.diff.expectation.2nd >= 1, abs.diff.expectation.2nd^2*10, 0)
  #   #penalty.4th.cum.var <- ifelse(abs.diff.var >= var.4th.observables*1/3, abs.diff.var^2*10, 0)
  # } else{
  #   penalty.4th.cum.expectation <- 0
  #   penalty.4th.cum.expectation.2nd <- 0
  #   #penalty.4th.cum.var <- 0
  # }
  
  res <- -KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$loglik + penalty.sign +sum(penalty.cycle.expectation)
  
  # if(!is.numeric(res)){
  #   res <- 1000000000
  # }
  # 
  # # To avoid extreme negative value undesirable
  if(res < -1000000){
    res <- abs(res)
  }
  
  return(res)
} 


fit.log.lik.trend.cycle.joint.model.4.3 <- function(thetas, estimated.Model, n, m.Y, q, r, nbr.horizon.max, var.type, all_mat) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using the mapping function 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  
  estimated.Model.list <- Retrieve.initial.par.trend.cycle.model(all.thetas.inter, n, m.Y, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by -theta[1,1] and theta[2,4] by -theta[2,3]
  estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]
  #estimated.Model.list$Theta[2,4] <- - estimated.Model.list$Theta[2,3]
  
  # Replace element nu[2,1] by nu[1,1] and nu[4,1] by nu[3,1]
  estimated.Model.list$nu[2,1] <- estimated.Model.list$nu[1,1]
  #estimated.Model.list$nu[4,1] <- estimated.Model.list$nu[3,1]
  
  # Replace element phi[2,2] by phi[1,1] and phi[4,4] by phi[3,3]
  estimated.Model.list$phi[2,2] <- estimated.Model.list$phi[1,1]
  #estimated.Model.list$phi[4,4] <- estimated.Model.list$phi[3,3]
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type, all_mat)
  Model <- estimated.Model.list
  
  # Generate matrices used by KF.
  KF.load <- prepare.KF.model(Model, observables)
  
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
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  # return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00, S="default",
  #                       indic.pos.z=KF.load$indic_pos)$log.lik)
  # return(-Kalman_filter_cpp(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00,
  #                       indic_pos_z=KF.load$indic_pos)$log.lik)
  
  penalty.sign <- sum(abs(sign(Model$delta.c[,1]) - sign(Model$delta.t[,1])) > 1)*10000 + sum(abs(sign(Model$delta.c[,2]) - sign(Model$delta.t[,2])) > 1)*10000
  
  return(-KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$loglik + penalty.sign)
} 

fit.log.lik.trend.cycle.joint.model.5.5 <- function(thetas, estimated.Model, n, m.Y, q, r, nbr.horizon.max, var.type, all_mat) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using the mapping function 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  
  estimated.Model.list <- Retrieve.initial.par.trend.cycle.model(all.thetas.inter, n, m.Y, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by -theta[1,1] and theta[2,4] by -theta[2,3]
  estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]
  estimated.Model.list$Theta[2,4] <- - estimated.Model.list$Theta[2,3]
  
  # Replace element nu[2,1] by nu[1,1] and nu[4,1] by nu[3,1]
  estimated.Model.list$nu[2,1] <- estimated.Model.list$nu[1,1]
  estimated.Model.list$nu[4,1] <- estimated.Model.list$nu[3,1]
  
  # Replace element phi[2,2] by phi[1,1] and phi[4,4] by phi[3,3]
  estimated.Model.list$phi[2,2] <- estimated.Model.list$phi[1,1]
  estimated.Model.list$phi[4,4] <- estimated.Model.list$phi[3,3]
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type, all_mat)
  Model <- estimated.Model.list
  
  # Generate matrices used by KF.
  KF.load <- prepare.KF.model(Model, observables)
  
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
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  # return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00, S="default",
  #                       indic.pos.z=KF.load$indic_pos)$log.lik)
  # return(-Kalman_filter_cpp(KF.load$all.parameters, KF.load$Y, KF.load$X,
  #                       KF.load$xi.00, KF.load$P.00,
  #                       indic_pos_z=KF.load$indic_pos)$log.lik)
  
  penalty.sign <- sum((sign(Model$delta.c[,1])!=sign(Model$delta.t[,1])))*10000
  
  return(-KF_filter_cpp(Y_t, StateSpace, indic_pos_z=KF.load$indic_pos)$loglik + penalty.sign)
} 

fit.log.lik.trend.cycle.model <- function(thetas, estimated.Model, n, m.Y, q, r, nbr.horizon.max, var.type) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using the mapping function 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  
  estimated.Model.list <- Retrieve.initial.par.trend.cycle.model(all.thetas.inter, n, m.Y, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by theta[1,1]
  estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]
  
  # Replace element nu[2,1] by nu[1,1]
  estimated.Model.list$nu[2,1] <- estimated.Model.list$nu[1,1]
  
  # Replace element phi[2,2] by phi[1,1]
  estimated.Model.list$phi[2,2] <- estimated.Model.list$phi[1,1]
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type)
  
  Model <- estimated.Model.list
  
  KF.load <- prepare.KF.model(Model, observables)
  # Generate matrices used by KF.
  
  # KF.res <- Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X, 
  #                           KF.load$xi.00, KF.load$P.00, S="default", 
  #                           indic.pos.z=KF.load$indic_pos)
  # 
  # if(mean(KF.res$fitted.obs[,7])> 2){
  #   
  #   KF.res$log.lik <- KF.res$log.lik-10000000
  #   
  # }
  # 
  # return(-KF.res$log.lik)
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X,
                        KF.load$xi.00, KF.load$P.00, S="default",
                        indic.pos.z=KF.load$indic_pos)$log.lik)
} 



fit.log.lik.trend.cycle.model.general <- function(thetas, estimated.Model, n, m.Y, q, r, nbr.horizon.max, var.type) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using the mapping function 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  estimated.Model.list <- Retrieve.initial.par.trend.cycle.model(all.thetas.inter, n, m.Y, q, r, nbr.horizon.max)
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list, var.type)
  
  Model <- estimated.Model.list
  
  # Generate matrices used by KF.
  KF.load <- prepare.KF.model(Model, observables)
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X,
                        KF.load$xi.00, KF.load$P.00, S="default",
                        indic.pos.z=KF.load$indic_pos)$log.lik)
} 

fit.log.lik.model.special <- function(thetas, estimated.Model, n, q, r, nbr.horizon.max) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using the mapping function 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  
  # Replace element theta[1,2] by theta[1,1]
  #all.thetas.inter[grep("Theta4", rownames(as.data.frame(all.thetas.inter))) ,1] <- -all.thetas.inter[grep("Theta1", rownames(as.data.frame(all.thetas.inter))) ,1]
  
  estimated.Model.list <- Retrieve.Initial.Par.Model(all.thetas.inter, n, q, r, nbr.horizon.max)
  
  # Replace element theta[1,2] by theta[1,1]
  estimated.Model.list$Theta[1,2] <- - estimated.Model.list$Theta[1,1]
  
  # Replace element nu[2,1] by nu[1,1]
  estimated.Model.list$nu[2,1] <- estimated.Model.list$nu[1,1]
  
  # Replace element phi[2,2] by phi[1,1]
  estimated.Model.list$phi[2,2] <- estimated.Model.list$phi[1,1]
  
  # Replace element phi[3,2] by phi[3,1]
  estimated.Model.list$phi[3,2] <- estimated.Model.list$phi[3,1]
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list)
  
  Model <- estimated.Model.list
  
  KF.load <- prepare.KF.model(Model, observables)
  # Generate matrices used by KF.
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X, 
                        KF.load$xi.00, KF.load$P.00, S="default", 
                        indic.pos.z=KF.load$indic_pos)$log.lik)
} 



fit.log.lik.model <- function(thetas, estimated.Model, n, q, r, nbr.horizon.max) { 
  # thetas: specify in column.
  # estimated.Model: a vector with all the coefficients and the indicators for
  #                  the transformation and the estimation.
  
  all.thetas.inter <- estimated.Model
  indicator.transformation <- estimated.Model[,3][ which(!estimated.Model[,2] == 0)]
  
  # Create the vector of parameters with thetas using 
  all.thetas.inter[rownames(as.data.frame(thetas)),1] <- Mapping.function(thetas, indicator.transformation)
  
  estimated.Model.list <- Retrieve.Initial.Par.Model(all.thetas.inter, n, q, r, nbr.horizon.max)
  
  ## Compute the parameters of the Model
  estimated.Model.list <- make.parameters.model(estimated.Model.list)
  
  Model <- estimated.Model.list
  
  KF.load <- prepare.KF.model(Model, observables)
  # Generate matrices used by KF.
  
  
  # Return the log-likelihood only (but perform the whole Kalman Filter).
  return(-Kalman.filter(KF.load$all.parameters, KF.load$Y, KF.load$X, 
                        KF.load$xi.00, KF.load$P.00, S="default", 
                        indic.pos.z=KF.load$indic_pos)$log.lik)
} 

# ============================================
#      K.  KALMAN SMOOTHER FUNCTION
# ============================================

Kalman.smoother <- function(all.parameters, Y, X, xi.00, P.00, S="default", indic.pos.z=0){
  
  # Inputs:
  # T = number of time period for observed variables/ horizon of simulation.
  # n = number of observed variables (Y).
  # k = number of deterministic component (X).
  # r = number of unobserved “state” variables.
  r <- dim(xi.00)[1]
  T <- dim(Y)[1]
  n <- dim(Y)[2]
  
  # Extract matrices from the list "all.parameters
  mu <- all.parameters$mu
  F <- all.parameters$F
  sigma <- all.parameters$sigma
  A <- all.parameters$A
  H <- all.parameters$H
  delta <- all.parameters$delta
  
  # Compute Kalman Filter. This will allow to have xi.tT and P.tT 
  result.KF <- Kalman.filter(all.parameters, Y, X, xi.00, P.00, S, indic.pos.z)
  
  # Extract results from the Kalman Filter.
  ## We extract the list as it is easier to work with in the loop
  matrix.xi.tt <- result.KF$xi.tt
  matrix.P.tt <- result.KF$P.tt
  list.P.ttm1 <- result.KF$`P.tt-1.list`
  list.xi.ttm1 <- result.KF$`xi.tt-1.list`
  
  # Initialize P.tT and xi.tT with T elements (filled with NA).
  matrix.P.tT <- matrix(NA,T,r*r)
  matrix.xi.tT <- matrix(NA,T,r)
  
  # Initialize P.tT and xi.tT with the last computed observation by the filter.
  matrix.P.tT[T,] <- matrix.P.tt[T,]
  matrix.xi.tT[T,] <- matrix.xi.tt[T,]
  
  # Loop to calculate the Kalman smoother. 
  # The loop start with period t-1 and compute the filter recursively
  for(t in (T-1):1) {
    
    ## Calculate J_t (adjustment term)
    #J_t <- matrix(matrix.P.tt[t,],r,r)%*%t(F)%*%solve(list.P.ttm1[[t+1]])
    J_t <- matrix(matrix.P.tt[t,],r,r)%*%t(F)%*%ginv(list.P.ttm1[[t+1]], tol=1e-30)
    
    ## Compute the smoothed data and their errors for each t.
    Xi_tT <- matrix(matrix.xi.tt[t,],r,1) + J_t %*% (matrix(matrix.xi.tT[t+1,],r,1)-
                                                       list.xi.ttm1[[t+1]])
    
    # Adjustment we have to make to the filter pertains to the fact that factors
    # z_t are non negative. For this purpose, after each updating step of the 
    # algorithm, negative entries in the z_t estimate are replaced by 0.
    if(sum(indic.pos.z==1)>0){
      Xi_tT[indic.pos.z==1,] = pmax(Xi_tT[indic.pos.z==1,],0)
    }
    
    P_tT<- matrix(matrix.P.tt[t,],r,r) + J_t%*%(matrix(matrix.P.tT[t+1,],r,r)-
                                                  list.P.ttm1[[t+1]])%*%t(J_t)
    
    ## Store all the smoothed data and their errors in a matrix.
    matrix.xi.tT[t,] <- t(Xi_tT)
    matrix.P.tT[t,] <- t(vec(P_tT))
    
  }
  
  # Compute fitted values for observables
  fitted.obs.smoothed <- X %*% A + matrix.xi.tT%*%H # Fitted observables
  
  # Return the final elements
  return(list("xi.tt"=matrix.xi.tt,
              "P.tt"=matrix.P.tt,
              "P.tt-1"=result.KF$`P.tt-1`,
              "P.tt-1.list"=list.P.ttm1,
              "xi.tt-1"=result.KF$`x.tt-1`,
              "xi.tt-1.list"=list.xi.ttm1,
              "xi.tT"=matrix.xi.tT,
              "P.tT"=matrix.P.tT,
              "K.t"=result.KF$list.K.t,
              "Q.t" <- result.KF$list.Q.t,
              "R.t" <- result.KF$list.R.t,
              "log.lik"= result.KF$log.lik,
              "log.lik.vec"= result.KF$log.lik.vec,
              "fitted.obs"=result.KF$fitted.obs,
              "fitted.obs.smoothed"=fitted.obs.smoothed
  ))
  
}


# ============================================
# L.  CREATE ESTIMATED SPFs (EXPECTATION AND
#     VARIANCE) FOR ALL HORIZONS SPECIFIED
#             
# ============================================

Make.SPF.Mean.Var <- function(Model, xi.tt, inflation.types=2){
  # This function Create the term structure of the first and the second moment
  # Model: list that contains all the estimates/parameters of interest
  # xi.tt: latent variables (can be smoothed or not) in vector form
  #        ideally it should be the output of the KF function
  # inflation.types: specify the definition of inflation that should be considered
  #                  in order to compute the SPF.
  #     1. Inflation is average inflation between t and t+h, and inflation
  #        between t-12 and t is of the form Gamma'X_t. k = h/12.
  #        In this case, we have b_h = 1/k delta (Phi.X^(k-1)*12 + ... + Phi.X^h-12 + Phi.X^h).
  #     2. Inflation between t+h-12 and t+h (simpler).
  #        In this case, we have b_h = delta*Phi.X^h
  #        USE " BY DEFAULT TO COMPUTE THE SPFs (IT IS EASIEST WAY TO COMPARE BETWEEN COUNTRIES)
  #     3. "US.SPF-style": average of Price index for a given year over average Price 
  #         index for the preceding year (in quarter).
  #         In this case, we have b_h = 1/4*delta(Phi.X^h-9 + Phi.X^h-6 + Phi.X^h-3 + Phi.X^h)
  #     4. "KOF/CH SPF-style": average of Price index for a given year over average 
  #         Price index for the preceding year.
  #         In this case, we have b_h = 1/4*delta(Phi.X^h-11 + Phi.X^h-10 + ... + Phi.X^h)
  
  
  # Define horizon we want to display
  H <- 12:120
  
  # Extract dimensions of interest
  T <- dim(xi.tt)[1]
  r <- length(c(Model$pi.bar)) # number of areas or key variables
  k <- length(Model$w) # number of additional measurement equations, if any.
  n <- Model$n
  q <- Model$q
  
  ## Get the Horizon of interest
  H.aux <- H[!is.na(H)]
  
  # Vector of delta (inflation eq) with 0 for the observables of the variances
  Gamma <- rbind(Model$delta,matrix(0,Model$q,r))
  
  # S.n and S.q
  S.n <- make.S.p(n)
  S.q <- make.S.p(q)
  
  # Compute the loadings in order to then fill matrix HH
  loadings <- compute.loadings.c(c(Gamma),Model,H.aux,S.n,S.q,indic.5y.in.Xy = Indic.5y.in.Xy)
  
  A1 <- matrix(Model$pi.bar[1] + c(loadings$a[1,,inflation.types]),length(H.aux),T)
  
  A2 <- matrix(c(loadings$alpha[1,,inflation.types]),length(H.aux),T)
  
  A3 <- matrix(c(loadings$alpha.dot.dot[1,,inflation.types]),length(H.aux),T)
  
  H1 <- t(loadings$b[,,inflation.types])
  
  H2 <- t(loadings$beta[,,inflation.types])
  
  H3 <- t(loadings$beta.dot.dot[,,inflation.types])
  
  # each column represents one maturity for all time period available (line)
  SPF.mean <- A1 + H1%*%t(xi.tt)
  SPF.var <- A2 + H2%*%t(xi.tt)
  SPF.k3rd <- A3 + H3%*%t(xi.tt)
  
  # Column represent the maturity and the line the time period 
  return(list(
    "SPF.mean" = t(SPF.mean),
    "SPF.var" = t(SPF.var),
    "SPF.k3rd" = t(SPF.k3rd)
  ))
}


# ============================================
# M.  BETA DISTRIBUTION FUNCTION
#             
# ============================================

# Compute the incomplete beta function, i.e. Beta(x;a,b)
# Beta[x;a,b]= F(x)*B(a,b)
# Notes : - pbeta() is a function that provide the ability for generating probability 
#         density values, cumulative probability density values and moment about 
#         zero values for the Beta Distribution bounded between [0,1]. 
#         - beta is a function that compute the beta function value for the given parameters.
incomplete.beta <- function(x,a,b){
  return(pbeta(x,a,b)*beta(a,b))
}

# Compute the CDF of the beta distribution called "The regularized incomplete beta function"
# With parameters to center and reduce the distribution (c,d).
# 4 parameters beta dist to adjust support (location and scale); 
# c: min of dist; d: max of dist
# F(x)=Beta[(x-c)/(d-c);a,b]/B(a,b)
cdf.incompl.beta <- function(x,param,max.c,min.d){  
  #max.c: minimum value for c, based on the upper limit of the smallest bins with observations
  #min.c: minimum value for d, based the upper limit of the smallest bins with observations
  a <- 1+abs(param[1])
  b <- 1+abs(param[2])
  c <- -abs(param[3]) + max.c
  d <- abs(param[4]) + min.d
  return(
    incomplete.beta((x-c)/(d-c),a,b)/beta(a,b)
  )
}


cdf.incompl.beta.PE <- function(x,param,max.c,min.d, PE=NA){  
  #max.c: minimum value for c, based on the upper limit of the smallest bins with observations
  #min.c: minimum value for d, based the upper limit of the smallest bins with observations
  
  # Convert the parameter if we have a PE that we want to perfectly fitted
  c <- -abs(param[3]) + max.c
  d <- abs(param[4]) + min.d
  
  if(is.na(PE)){# When we have no PE, the four parameters can fluctuate
    a <- 1+abs(param[1])
    b <- 1+abs(param[2])
  } else{# When we want to fit PE, a or b should be fixed depending on the ratio (PE-c)/(d-PE)
    if((PE-c)/(d-PE)>1){# Ratio > 1 => fix a
      b <- 1+abs(param[2])
      a <- b*(PE-c)/(d-PE)  
    } else{# Ratio < 1 => fix b
      a <- 1+abs(param[1])
      b <- a*(d-PE)/(PE-c)
    }
  }
  
  return(
    incomplete.beta((x-c)/(d-c),a,b)/beta(a,b)
  )
}

# CDF with 4 parameters but in vectorial form
cdf.incompl.beta.vectorial <- function(x,param,max.c,min.d){  #Appendix section 7.4.1; F(x)=Beta[(x-c)/(d-c);a,b]/B(a,b)
  a <- 1+abs(param[,1])
  b <- 1+abs(param[,2])
  c <- -abs(param[,3]) + max.c
  d <- abs(param[,4]) + min.d
  return(
    incomplete.beta((x-c)/(d-c),a,b)/beta(a,b)
  )
}

# Compute the PDF of the beta distribution
# With parameters to center and reduce the distribution (c,d).
# 4 parameters beta dist to adjust support (location and scale); 
# c: min of dist; d: max of dist
# pdf=[F(x2)-F(x1)]/[x2-x1]
pdf.incompl.beta2 <- function(x,param,max.c,min.d, increment=0.05){  
  x2 <- x - increment
  return(
    (cdf.incompl.beta(x,param,max.c,min.d)-cdf.incompl.beta(x2,param,max.c,min.d))/(x-x2)
  )
}

# Compute the PDF of the beta distribution
# With parameters to center and reduce the distribution (c,d).
# 4 parameters beta dist to adjust support (location and scale); 
# c: min of dist; d: max of dist
# f(x;a,b,c,d)=1/((d-c)*beta(a,b))*((x-c)/(d-c))^(a-1)*((d-x)/(d-c))^(b-1)
pdf.incompl.beta <- function(x,param,max.c,min.d){  
  a <- 1+abs(param[1])
  b <- 1+abs(param[2])
  c <- -abs(param[3]) + max.c
  d <- abs(param[4]) + min.d
  return(
    ((x-c)/(d-c))^(a-1)*((d-x)/(d-c))^(b-1)/((d-c)*beta(a,b))
  )
}

# Return the squared error between the observed cdf values and the model implied
# cdf values (function to fit the cdf)
dist.2.incompl.beta <- function(param,max.c,min.d,x,cdf.values){
  model.implied.cdf <- cdf.incompl.beta(x,param,max.c,min.d)
  return(sum((cdf.values - model.implied.cdf)^2))
}


dist.2.incompl.beta.PE <- function(param,max.c,min.d,x,cdf.values,PE=NA, var.min=0){
  #model.implied.cdf <- cdf.incompl.beta.PE(x,param,max.c,min.d, PE=NA)
  model.implied.cdf <- cdf.incompl.beta.PE(x,param,max.c,min.d, PE)
  moments <- moments.incompl.beta.PE(param,max.c,min.d,PE)
  penalty <- max(var.min-moments$Variance,0)*1
  
  return(sum((cdf.values - model.implied.cdf)^2) + penalty)
}

# Function that fits the observed cdf values with a beta distribution.
# The function estimates the parameters a,b,c,d to construct the beta distribution
fit.cdf <- function(x,cdf.values,param.0=c(5,5,-1,6)){
  # Look for the parameterization (a,b,c,d) that provides the best fit of the CDF,
  # The x vector contains the bins.
  # The cdf.values vector contains the cumulative probabilities (sum of all previous bins).
  # param.0 is a 4-dimensional vector containing the initial conditions for the parameter values.
  max.c <- x[min(which(round(cdf.values,3)>0))] - max(x- lag(x), na.rm = T)
  min.d <- x[max(which(round(cdf.values,3)<1))] + max(x- lag(x), na.rm = T)
  #max.c <- min(x) - max(x- lag(x), na.rm = T)
  #min.d <- max(x) + max(x- lag(x), na.rm = T)
  
  for(i in 1:4){
    res.optim <- optim(param.0,dist.2.incompl.beta,max.c=max.c, min.d=min.d,
                       x=x, cdf.values=cdf.values,
                       gr = NULL,
                       method="Nelder-Mead",
                       #method="BFGS",
                       #                     method="CG",
                       control=list(trace=FALSE,maxit=300))
    param.0 <- res.optim$par
    res.optim <- optim(param.0,dist.2.incompl.beta,max.c=max.c, min.d=min.d,
                       x=x, cdf.values=cdf.values,
                       gr = NULL,
                       #method="Nelder-Mead",
                       method="BFGS",
                       #                     method="CG",
                       control=list(trace=FALSE,maxit=10))
    param.0 <- res.optim$par
  }
  
  param <- param.0
  model.implied.cdf <- cdf.incompl.beta(x,param,max.c,min.d)
  return(list(
    param=param,model.implied.cdf=model.implied.cdf,
    cdf.values=cdf.values,x=x,max.c=max.c,min.d=min.d
  ))
}


fit.cdf.PE <- function(x,cdf.values,param.0=c(5,5,-1,6),PE=NA,max.c=NA, var.min=0){
  # Look for the parameterization (a,b,c,d) that provides the best fit of the CDF,
  # The x vector contains the bins.
  # The cdf.values vector contains the cumulative probabilities (sum of all previous bins).
  # param.0 is a 4-dimensional vector containing the initial conditions for the parameter values.
  if(is.na(max.c)){
    max.c <- x[min(which(round(cdf.values,3)>0))] - max(x- lag(x), na.rm = T)
  } 
  min.d <- x[max(which(round(cdf.values,3)<1))] + max(x- lag(x), na.rm = T)
  
  # Optimization 
  for(i in 1:4){
    res.optim <- optim(param.0,dist.2.incompl.beta.PE,max.c=max.c, min.d=min.d,
                       x=x, cdf.values=cdf.values, PE=PE, var.min=var.min,
                       gr = NULL,
                       method="Nelder-Mead",
                       #method="BFGS",
                       #                     method="CG",
                       control=list(trace=FALSE,maxit=300))
    param.0 <- res.optim$par
    # res.optim <- optim(param.0,dist.2.incompl.beta.PE,max.c=max.c, min.d=min.d,
    #                    x=x, cdf.values=cdf.values, PE=PE,
    #                    gr = NULL,
    #                    #method="Nelder-Mead",
    #                    method="BFGS",
    #                    #                     method="CG",
    #                    control=list(trace=FALSE,maxit=10))
    tryCatch({
      res.optim <- optim(param.0,dist.2.incompl.beta.PE,max.c=max.c, min.d=min.d,
                         x=x, cdf.values=cdf.values, PE=PE, var.min=var.min,
                         gr = NULL,
                         #method="Nelder-Mead",
                         method="BFGS",
                         #                     method="CG",
                         control=list(trace=FALSE,maxit=10))
    }, error = function(err) {
      # Handle the optimization error here
      cat("Optimization error occurred: ", conditionMessage(err), "at iteration", i , "\n")
      
      # Implement an alternative approach or provide default values
      res.optim <- optim(param.0,dist.2.incompl.beta.PE,max.c=max.c, min.d=min.d,
                         x=x, cdf.values=cdf.values, PE=PE, var.min=var.min,
                         gr = NULL,
                         method="Nelder-Mead",
                         #method="BFGS",
                         #                     method="CG",
                         control=list(trace=FALSE,maxit=300))
      
    })
    param.0 <- res.optim$par
  }
  
  # Convert the parameter if we have a PE that we want to perfectly fitted
  c <- -abs(param.0[3]) + max.c
  d <- abs(param.0[4]) + min.d
  
  if(is.na(PE)){# When we have no PE, the four parameters can fluctuate
    a <- 1+abs(param.0[1])
    b <- 1+abs(param.0[2])
  } else{# When we want to fit PE, a or b should be fixed depending on the ratio (PE-c)/(d-PE)
    if((PE-c)/(d-PE)>1){# Ratio > 1 => fix a
      b <- 1+abs(param.0[2])
      a <- b*(PE-c)/(d-PE)
      param.0[1] <- a-1
    } else{# Ratio < 1 => fix b
      a <- 1+abs(param.0[1])
      b <- a*(d-PE)/(PE-c)
      param.0[2] <- b-1
    }
  }
  
  param <- param.0
  model.implied.cdf <- cdf.incompl.beta(x,param,max.c,min.d)
  error.squared <-sum((cdf.values - model.implied.cdf)^2)
  return(list(
    param=param,model.implied.cdf=model.implied.cdf, error.squared=error.squared,
    cdf.values=cdf.values,x=x,max.c=max.c,min.d=min.d
  ))
}

# Function that computes the mean and the variance of the estimated beta distribution
# 4 parameters beta dist to adjust support (location and scale); 
# c: min of dist; d: max of dist
# For a standard beta distribution, the general formula for E(X) and V(X) are:
# E(X) = a/(a+b) ; V(X) = a*b/[(a+b)^2*(a+b+1)] 
moments.incompl.beta <- function(param,max.c,min.d){ 
  a <- 1+abs(param[1])
  b <- 1+abs(param[2])
  c <- -abs(param[3]) + max.c
  d <- abs(param[4]) + min.d
  Mean <- (a*d + b*c)/(a+b)
  Variance <- (d-c)^2*a*b/(a+b)^2/(a+b+1)
  Skewness <- 2*(b-a)*sqrt(a+b+1)/((a+b+2)*sqrt(a*b))
  Kurtosis <- 6*((a-b)^2*(a+b+1)-a*b*(a+b+2))/((a*b)*(a+b+2)*(a+b+3)) + 3
  return(list(
    Mean = Mean, Variance = Variance, Skewness = Skewness, Kurtosis= Kurtosis
  ))
}

moments.incompl.beta.PE <- function(param,max.c,min.d,PE=NA){ 
  
  # Convert the parameter if we have a PE that we want to perfectly fitted
  c <- -abs(param[3]) + max.c
  d <- abs(param[4]) + min.d
  
  if(is.na(PE)){# When we have no PE, the four parameters can fluctuate
    a <- 1+abs(param[1])
    b <- 1+abs(param[2])
  } else{# When we want to fit PE, a or b should be fixed depending on the ratio (PE-c)/(d-PE)
    if((PE-c)/(d-PE)>1){# Ratio > 1 => fix a
      b <- 1+abs(param[2])
      a <- b*(PE-c)/(d-PE)  
    } else{# Ratio < 1 => fix b
      a <- 1+abs(param[1])
      b <- a*(d-PE)/(PE-c)
    }
  }
  
  Mean <- (a*d + b*c)/(a+b)
  Variance <- (d-c)^2*a*b/(a+b)^2/(a+b+1)
  Skewness <- 2*(b-a)*sqrt(a+b+1)/((a+b+2)*sqrt(a*b))
  Kurtosis <- 6*((a-b)^2*(a+b+1)-a*b*(a+b+2))/((a*b)*(a+b+2)*(a+b+3)) + 3
  return(list(
    Mean = Mean, Variance = Variance, Skewness = Skewness, Kurtosis= Kurtosis
  ))
}


from.m.to.moment.of.interest <- function(m.1,m.2,m.3,m.4){
  
  Mean <- m.1
  k.1 <- Mean
  Variance <- m.2-m.1^2
  k.2 <- Variance 
  Skewness <- (m.3-3*m.1*Variance-m.1^3)/(sqrt(Variance))^3
  k.3 <- m.3 - 3*m.1*k.2 - m.1^3
  k.4 <- m.4 - 4*k.3*k.1 - 3*k.2^2 - 6*k.2*k.1^2 - k.1^4
  Kurtosis <- (k.4 + 3*(k.2)^2)/k.2^2
  
  return(list(
    Mean = Mean,
    Variance = Variance,
    Skewness = Skewness,
    Kurtosis = Kurtosis,
    k.1 = k.1,
    k.2 = k.2,
    k.3 = k.3,
    k.4 = k.4
  ))
  
}


make.cumulant.until.order.4 <- function(Mean,Variance,Skewness,Kurtosis){
  
  # Compute the moments
  m.1 <- Mean
  m.2 <- Variance + m.1^2
  m.3 <- (sqrt(m.2 - m.1^2))^3*Skewness + m.1^3 + 3*m.1*Variance
  #m.3 <- (sqrt(m.2 - m.1^2))^3*Skewness + m.1^3 #### FALSE 
  m.3.bis <- (sqrt(Variance))^3*Skewness + m.1^3 + 3*m.1*Variance
  #m.3.bis <- (sqrt(Variance))^3*Skewness + m.1^3 #### FALSE 
  
  # Compute the cumulants  
  k.1 <- Mean
  k.2 <- Variance 
  k.3 <- m.3 - 3*m.1*k.2 - m.1^3
  k.3.bis <- k.2^(3/2)*Skewness
  k.3.bis.bis <- m.3 - 3*m.1*m.2 + 2*m.1^3
  k.4 <- Kurtosis*k.2^2 - 3*(k.2)^2
  m.4 <- k.4 + 4*k.3*k.1 + 3*k.2^2 + 6*k.2*k.1^2 + k.1^4
  
  
  return(list(
    m.1 = m.1,
    m.2 = m.2,
    m.3 = m.3,
    m.3.bis = m.3.bis,
    m.4 =m.4,
    k.1 = k.1,
    k.2 = k.2,
    k.3 = k.3,
    k.3.bis = k.3.bis,
    k.3.bis.bis = k.3.bis.bis,
    k.4 = k.4
  ))
}

cum.3 <- function(data){
  # Function that computes the 3rd cumulant based on data
  
  # Compute the cumulants  
  k.3 <- var(data)^(3/2)*skewness(data)
  return(k.3)
  
}

#For your info, compute density versus Kernel (normal)
# dnorm(x2) = create pdf or pdf=[F(x2)-F(x1)]/[x2-x1]
# pnorm() = create cdf or 
# simul.norm <- rnorm(3000000)
# plot(seq(-5.5,5.5,by=0.01),dnorm(seq(-5.5,5.5,by=0.01)), main="pdf normal", type="l", ylab ="")
# hist(simul.norm, freq=FALSE, main="computed kernel density based on simulations")
# lines(density(simul.norm))

# Function that plots the distribution with the histogram of observed survey data
# data: non aggregated data 
# x: values at which the pdf should be evaluated
# param; max.c; min.d: parameters of beta dist.
plot.fit.survey.distribution <- function(data,x,param,max.c,min.d, mean.class, break2, xtitle="Inflation rate %"){
  
  proba <- matrix(as.numeric(data))/100
  nbr.obs.class <- round(proba*1000,0)
  #mean.class <- seq(-2.25, 3.75, 0.5)
  serie <- NULL
  for(i in 1:length(nbr.obs.class)){
    serie <- c(serie,rep(mean.class[i],nbr.obs.class[i]))
  }
  
  #break1 <- seq(-2.5, 4, 0.5)
  #break2 <- sort(c(seq(-2.5, 4, 0.5), seq(-1.55, 3.45, by =0.5)))
  pdf <- pdf.incompl.beta(x,param,max.c,min.d)
  #xlim = c(min(break2)-1,max(break2)+1)
  h <- hist(serie, breaks = break2, xlim = c(min(break2)-1,max(break2)+1), main="", border="grey",
            freq = FALSE, xlab=xtitle, ylim=c(0, max(pdf,proba/0.75, na.rm = T)+0.1), cex.lab=0.75, cex.axis=0.75, cex.main=0.75, cex.sub=0.75)
  box(col = "black") #ylim=c(0, max(pdf,proba/0.45, na.rm = T)+0.1)
  lines(x,pdf)
  
}

# ============================================
# M.bis  GAUSSIAN MIXTURE DISTRIBUTION FUNCTION
#             
# ============================================

# CDF of gaussian mixture distribution. There are five parameters.
# probability p, means of two gaussians, stdv of thwo gaussians
cdf.mixture <- function(x,param,min.bin,max.bin, min.sigma=0.25){
  # is the min of mu_1 and mu_2, b is the max of mu_1 and mu_2
  ## We remove min.sigma such that one mean is at least one sd 
  ## smaller or bigger than max and mean (mean not at border)
  a <- min.bin+min.sigma
  b <- max.bin-min.sigma
  if(param[1]>700){param[1]<-700}
  p <- 0.999*exp(param[1])/(1+exp(param[1]))
  #min.sigma <- 0.25
  max.sigma <- (max.bin-min.bin)/6
  #max.sigma <- (max.bin-min.bin)/2
  sigma1 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[2])) #0.25 + abs(param[2])
  sigma2 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[3])) #0.25 + abs(param[3])
  mu1 <- a+(b-a)/(1+exp(-param[4])) #param[4] #5*(2*exp(param[4])/(1+exp(param[4]))-1) 
  mu2 <- a+(b-a)/(1+exp(-param[5])) #param[5] #5*(2*exp(param[5])/(1+exp(param[5]))-1) 
  return( p*pnorm((x-mu1)/sigma1) + (1-p)*pnorm((x-mu2)/sigma2) )
}

# CDF of gaussian mixture distribution, fitting the point estimates
cdf.mixture.PE <- function(x,param,min.bin,max.bin,min.sigma=0.25,PE=NaN){
  a <- min.bin+min.sigma
  b <- max.bin-min.sigma
  if(param[1]>700){param[1] <- 700}
  p <- 0.999*exp(param[1])/(1+exp(param[1]))
  #min.sigma <- 0.25
  max.sigma <- (max.bin-min.bin)/6
  #max.sigma <- (max.bin-min.bin)/2
  sigma1 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[2])) #0.25 + abs(param[2])
  sigma2 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[3])) #0.25 + abs(param[3])
  mu1 <- a+(b-a)/(1+exp(-param[4])) #param[4] #5*(2*exp(param[4])/(1+exp(param[4]))-1) 
  
  if(is.na(PE)){# When we have no PE, the four parameters can fluctuate
    mu2 <- a+(b-a)/(1+exp(-param[5])) #param[5] #5*(2*exp(param[5])/(1+exp(param[5]))-1) 
  }else{# When we want to fit PE, mu_2 should be fixed
    mu2 <- (PE - mu1*p)/(1-p)
  }
  return( p*pnorm((x-mu1)/sigma1) + (1-p)*pnorm((x-mu2)/sigma2) )
}

# PDF of gaussian mixture distribution, fitting the point estimates
pdf.mixture.PE <- function(x,param,min.bin,max.bin,min.sigma=0.25,PE=NaN){
  a <- min.bin+min.sigma
  b <- max.bin-min.sigma
  if(param[1]>700){param[1]<-700}
  p <- 0.999*exp(param[1])/(1+exp(param[1]))
  #min.sigma <- 0.25
  max.sigma <- (max.bin-min.bin)/6
  #max.sigma <- (max.bin-min.bin)/2
  sigma1 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[2])) #0.25 + abs(param[2])
  sigma2 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[3])) #0.25 + abs(param[3])
  mu1 <- a+(b-a)/(1+exp(-param[4])) #param[4] #5*(2*exp(param[4])/(1+exp(param[4]))-1) 
  if(is.na(PE)){# When we have no PE, the four parameters can fluctuate
    mu2 <- a+(b-a)/(1+exp(-param[5])) #param[5] #5*(2*exp(param[5])/(1+exp(param[5]))-1) 
  }else{# When we want to fit PE, a or b should be fixed depending on the ratio (PE-c)/(d-PE)
    mu2 <- (PE - mu1*p)/(1-p)
  }
  return( p*dnorm((x-mu1)/sigma1)/sigma1 + (1-p)*dnorm((x-mu2)/sigma2)/sigma2 )
}

# PDF of gaussian mixture distribution
pdf.mixture <- function(x,param,min.bin,max.bin,min.sigma=0.25){
  a <- min.bin+min.sigma
  b <- max.bin-min.sigma
  if(param[1]>700){param[1]<-700}
  p <- 0.999*exp(param[1])/(1+exp(param[1]))
  #min.sigma <- 0.25
  max.sigma <- (max.bin-min.bin)/6
  #max.sigma <- (max.bin-min.bin)/2
  sigma1 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[2])) #0.25 + abs(param[2])
  sigma2 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[3])) #0.25 + abs(param[3])
  mu1 <- a+(b-a)/(1+exp(-param[4])) #param[4] #5*(2*exp(param[4])/(1+exp(param[4]))-1) 
  mu2 <- a+(b-a)/(1+exp(-param[5])) #param[5] #5*(2*exp(param[5])/(1+exp(param[5]))-1) #param[5]
  return( p*dnorm((x-mu1)/sigma1)/sigma1 + (1-p)*dnorm((x-mu2)/sigma2)/sigma2 )
}


find_peaks <- function(data) {
  peaks <- which(diff(sign(diff(data))) == -2) + 1
  return(length(peaks))
}

penalty_small_peak <- function(data) {
  min<-which(diff(sign(diff(data))) == 2) + 1
  max<-which(diff(sign(diff(data))) == -2) + 1
  return(min(data[max]-data[min]))
}

#Penalty on first or last bins if two much different
penalty.first.last.bins <- function(cdf.values,model.implied.cdf){
  
  return(abs(abs(cdf.values-model.implied.cdf)-0.01)*1)
  
}

# Return the squared error between the observed cdf values and the model implied
# cdf values (function to fit the cdf)
dist.2.mixture.PE <- function(param,x,cdf.values,min.bin,max.bin,min.sigma=0.25,PE=NA){
  model.implied.cdf <- cdf.mixture.PE(x,param,min.bin,max.bin,min.sigma,PE)
  x.bis <- seq(min(x-0.5),max(x+0.5),0.05)
  penalty <- (find_peaks(pdf.mixture.PE(x.bis,param,min.bin,max.bin,min.sigma,PE))-1)*0.1
  penalty2 <- 0
  if(penalty>0){penalty2 <- penalty_small_peak(pdf.mixture.PE(x.bis,param,min.bin,max.bin,min.sigma,PE))/2 }
  
  #penalty last or first bins
  l.bin <- length(model.implied.cdf)
  penalty3 <-0
  penalty4 <- 0
  if(abs(cdf.values[1] - model.implied.cdf[1]) > 0.01){
    penalty3 <- penalty.first.last.bins(cdf.values[1],model.implied.cdf[1])
  }
  if(abs(cdf.values[l.bin] - model.implied.cdf[l.bin]) > 0.01){
    penalty4 <- penalty.first.last.bins(cdf.values[l.bin],model.implied.cdf[l.bin])
  }
  
  return(sum((cdf.values - model.implied.cdf)^2)  + penalty2 + penalty3 + penalty4)
}


# Function that fits the observed cdf values with a gaussian mixture distribution.
# The function estimates five parameters to construct the gaussian mixture distribution
fit.cdf.mixture.PE <- function(x,cdf.values,param.0=c(0,1,1,2,3),min.bin,max.bin,min.sigma=0.25,PE=NA){
  # Look for the parameterization (a,b,c,d,e) that provides the best fit of the CDF,
  # The x vector contains the bins.
  # The cdf.values vector contains the cumulative probabilities (sum of all previous bins).
  # param.0 is a 5-dimensional vector containing the initial conditions for the parameter values.
  
  # Optimization 
  for(i in 1:8){
    res.optim <- optim(param.0,dist.2.mixture.PE,
                       x=x, cdf.values=cdf.values,
                       min.bin=min.bin,max.bin=max.bin,min.sigma=min.sigma,
                       PE=PE,gr = NULL,
                       method="Nelder-Mead",
                       control=list(trace=FALSE,maxit=500))
    param.0 <- res.optim$par
    
    tryCatch({
      res.optim <- optim(param.0,dist.2.mixture.PE,
                         x=x, cdf.values=cdf.values,
                         min.bin=min.bin,max.bin=max.bin,min.sigma=min.sigma,
                         PE=PE,gr = NULL,
                         method="BFGS",
                         control=list(trace=FALSE,maxit=20))
    }, error = function(err) {
      # Handle the optimization error here
      cat("Optimization error occurred: ", conditionMessage(err), "at iteration", i , "\n")
      
      # Implement an alternative approach or provide default values
      res.optim <- optim(param.0,dist.2.mixture.PE,
                         x=x, cdf.values=cdf.values, 
                         min.bin=min.bin,max.bin=max.bin,min.sigma=min.sigma,
                         PE=PE, gr = NULL,
                         method="Nelder-Mead",
                         control=list(trace=FALSE,maxit=300))
      
    })
    param.0 <- res.optim$par
  }
  a <- min.bin+min.sigma
  b <- max.bin-min.sigma
  if(param.0[1]>700){param.0[1]<-700}
  p <- 0.999*exp(param.0[1])/(1+exp(param.0[1]))
  #min.sigma <- 0.25
  max.sigma <- (max.bin-min.bin)/6
  #max.sigma <- (max.bin-min.bin)/2
  sigma1 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param.0[2])) #0.25 + abs(param[2])
  sigma2 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param.0[3])) #0.25 + abs(param[3])
  
  if(is.na(PE)){# When we have no PE, the five parameters can fluctuate
    mu1 <- a+(b-a)/(1+exp(-param.0[4])) #param[4] #5*(2*exp(param[4])/(1+exp(param[4]))-1) 
    mu2 <- a+(b-a)/(1+exp(-param.0[5])) #param[5] #5*(2*exp(param[5])/(1+exp(param[5]))-1) 
  }else{# When we want to fit PE, a or b should be fixed depending on the ratio (PE-c)/(d-PE)
    mu1 <- a+(b-a)/(1+exp(-param.0[4])) #param[4] #5*(2*exp(param[4])/(1+exp(param[4]))-1) 
    mu2 <- a+(b-a)/(1+exp(-param.0[5]))
    increment <- 0.00000
    if((PE-mu1*p)/(1-p)>=min.bin & (PE-mu1*p)/(1-p)<=max.bin){ #p<0.5
      mu1 <- a+(b-a)/(1+exp(-param.0[4])) #param[4] #5*(2*exp(param[4])/(1+exp(param[4]))-1) 
      mu2 <- (PE - mu1*p)/(1-p)
      param.0[5] <- -log((b-(mu2-increment))/((mu2+increment)-a)) #mu2 #log(z/(1-z)) #mu2
    } else{
      mu2 <- a+(b-a)/(1+exp(-param.0[5]))
      mu1 <- (PE - mu2*(1-p))/p
      param.0[4] <- -log((b-(mu1-increment))/((mu1+increment)-a)) #mu2 #log(z/(1-z)) #mu2
    }
  }
  
  param <- param.0
  model.implied.cdf <- cdf.mixture(x,param,min.bin,max.bin,min.sigma)
  error.squared <-sum((cdf.values - model.implied.cdf)^2)
  return(list(
    param=param,
    model.implied.cdf=model.implied.cdf,
    error.squared=error.squared,
    cdf.values=cdf.values,
    x=x))
}


# Function that computes the moments for gaussian mixture distribution
# param has dim 5
moments.mixture.PE <- function(param,min.bin,max.bin,min.sigma=0.25,PE=NA){ 
  a <- min.bin+min.sigma
  b <- max.bin-min.sigma
  if(param[1]>700){param[1]<-700}
  p <- 0.999*exp(param[1])/(1+exp(param[1]))
  #min.sigma <- 0.25
  max.sigma <- (max.bin-min.bin)/6
  #max.sigma <- (max.bin-min.bin)/2
  sigma1 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[2])) #0.25 + abs(param[2])
  sigma2 <- min.sigma+(max.sigma-min.sigma)/(1+exp(-param[3])) #0.25 + abs(param[3])
  
  mu1 <- a+(b-a)/(1+exp(-param[4])) #param[4] #5*(2*exp(param[4])/(1+exp(param[4]))-1) #param[4]
  if(is.na(PE)){# When we have no PE, the four parameters can fluctuate
    mu2 <- a+(b-a)/(1+exp(-param[5])) #param[5] #5*(2*exp(param[5])/(1+exp(param[5]))-1) #param[5]
  }else{# When we want to fit PE, a or b should be fixed depending on the ratio (PE-c)/(d-PE)
    mu2 <- (PE - mu1*p)/(1-p)
  }
  
  m.1 <- mu1*p + mu2*(1-p)
  m.2 <-  p*(mu1^2+sigma1^2) + (1-p)*(mu2^2+sigma2^2) 
  m.3 <- p*(mu1^3 + 3*mu1*sigma1^2) + (1-p)*(mu2^3 + 3*mu2*sigma2^2)
  m.4 <- p*(mu1^4 + 6*mu1^2*sigma1^2 + 3*sigma1^4) +
    (1-p)*(mu2^4 + 6*mu2^2*sigma2^2 + 3*sigma2^4)
  
  Mean <- from.m.to.moment.of.interest(m.1,m.2,m.3,m.4)$Mean
  Variance <- from.m.to.moment.of.interest(m.1,m.2,m.3,m.4)$Variance
  Skewness <- from.m.to.moment.of.interest(m.1,m.2,m.3,m.4)$Skewness
  Kurtosis <- from.m.to.moment.of.interest(m.1,m.2,m.3,m.4)$Kurtosis
  #Mean <- mu1*p + mu2*(1-p)
  #Variance <- p*sigma1^2 + (1-p)*sigma2^2 + p*(1-p)*(mu1 - mu2)^2
  #Skewness <- (3*(p-(1-p))*(mu1-mu2)^3/((sigma1^2 + sigma2^2 + (mu1-mu2)^2)^(3/2)))/sqrt(Variance)
  #Kurtosis <- p*3 + (1-p)*3 + 3*(p-(1-p))^2
  
  #Skewness <- p*(mu1^3 + 3*mu1*sigma1^2) + (1-p)*(mu2^3 + 3*mu2*sigma2^2)
  #Kurtosis <- p*(mu1^4 + 6*mu1^2*sigma1^2 + 3*sigma1^4) +
  #  (1-p)*(mu2^4 + 6*mu2^2*sigma2^2 + 3*sigma2^4)
  return(list(
    Mean = Mean, Variance = Variance,
    Skewness = Skewness, Kurtosis= Kurtosis,
    m.1=m.1, m.2=m.2, m.3=m.3, m.4=m.4))
}


plot.fit.survey.distribution.mixture <- function(data,x,param,
                                                 min.bin,
                                                 max.bin,
                                                 min.sigma=0.25,
                                                 mean.class,
                                                 break2,
                                                 xtitle="Inflation rate %",
                                                 display.mixture=TRUE,
                                                 cex=.75){
  
  proba <- matrix(as.numeric(data))/100
  nbr.obs.class <- round(proba*1000,0)
  serie <- NULL
  for(i in 1:length(nbr.obs.class)){
    serie <- c(serie,
               rep(mean.class[i],nbr.obs.class[i]))
  }
  
  pdf <- pdf.mixture(x,param,
                     min.bin,
                     max.bin,
                     min.sigma)
  h <- hist(serie,
            breaks = break2,
            xlim = c(min(break2)-1,max(break2)+1),
            main =  "",
            ylab = "",
            border = "grey",
            freq = FALSE,
            xlab = xtitle,
            ylim = c(0, max(pdf,proba, na.rm = T)+0.10),
            cex.lab  = cex,
            cex.axis = cex,
            cex.main = cex,
            cex.sub  = cex,
            las=1)
  
  box(col = "black") 
  
  if(display.mixture){
    lines(x,pdf,lwd=2,lty=3)
  }
  return(h)
}


################

Make.SPF.Mean.Var <- function(Model, xi.tt, inflation.types){
  # This function Create the term structure of the first and the second moment
  
  # Define horizon we want to display
  H <- 12:120
  
  # Extract dimensions of interest
  T <- dim(xi.tt)[1]
  r <- length(c(Model$pi.bar)) # number of areas or key variables
  k <- length(Model$w) # number of additional measurement equations, if any.
  n <- Model$n
  q <- Model$q
  
  ## Get the Horizon of interest
  H.aux <- H[!is.na(H)]
  
  # Vector of delta (inflation eq) with 0 for the observables of the variances
  Gamma <- rbind(Model$delta,matrix(0,Model$q,r))
  
  # S.n and S.q
  S.n <- make.S.p(n)
  S.q <- make.S.p(q)
  
  # Compute the loadings in order to then fill matrix HH
  loadings <- compute.loadings.c(c(Gamma),Model,H.aux,S.n,S.q,indic.5y.in.Xy = Indic.5y.in.Xy)
  
  A1 <- matrix(Model$pi.bar[1] + c(loadings$a[1,,inflation.types]),length(H.aux),T)
  
  A2 <- matrix(c(loadings$alpha[1,,inflation.types]),length(H.aux),T)
  
  H1 <- t(loadings$b[,,inflation.types])
  
  H2 <- t(loadings$beta[,,inflation.types])
  
  # each column represents one maturity for all time period available (line)
  SPF.mean <- A1 + H1%*%t(xi.tt)
  SPF.var <- A2 + H2%*%t(xi.tt)
  
  # Column represent the maturity and the line the time period 
  return(list(
    "SPF.mean" = t(SPF.mean),
    "SPF.var" = t(SPF.var)
  ))
}


# ============================================
# N.  MODEL IMPLIED DISTRIBUTION FUNCTIONS
#             
# ============================================

# Function that computes the model implied probability for all time.
# The function return the implied probability => one point of the distribution for T periods.
# This is done for one y => calculate the probability with the same y for T periods and HH horizons.
compute.proba <- function(Model,X,b,y,HH,step,max.v,
                          indic.average=0,
                          AB.list.external=NaN,
                          Indic.5y.in.Xy=1){ #used in the function below compute.distri.plus.stdv
  # Computes the proba that bX is below y at horizons given in HH, ie,
  # E(1{b' * X_{+HH} < y}|X)
  # based on Duffie, Pan and Singleton (2000, Econometrica), formula 2.12
  #
  # step and max.step are used to numerically evaluate the integral of the formula:
  # i.e. v = .000001, step, 2 step, ..., max.v
  #
  # X is a matrix of dimension T x (Model$n + Model$q)
  # 
  # If indic.average = 1, then one will consider:
  #    E(1{b'(X_{+HH} + X_{+HH-12} + X_{+HH-24} + ...) < (HH/12) * y}|X)
  # 
  
  #### STEP 1: INITIALIZE ALL ####
  
  Proba <- NULL
  T <- dim(X)[1]
  #v <- matrix(seq(.00000001,max.v,by=step),nrow=1)
  # New way of parameterzing v:
  max.x <- log(max.v)
  v <- matrix(exp(seq(-10,max.x,by=step)),nrow=1)
  dim.v <- length(v)
  
  ## In case the AB have not been computed
  if(!is.list(AB.list.external)){# recompute Laplace transform
    u_Y <- 1i * matrix(b[1:Model$n],ncol=1) %*% v  # our paper, section 7.3 Appendix i*v*gamma for Y
    u_z <- 1i * matrix(b[(Model$n+1):(Model$n+Model$q)],ncol=1) %*% v  # our paper, section 7.3 Appendix i*v*gamma for z
    AB.list <- compute_AB(Model,u_Y,u_z,HH,indic.average)
  }else{# use Laplace transform already computed
    AB.list <- AB.list.external
  }
  if(indic.average==1){
    # =========================
    # =========================
    Y <- y * as.integer((HH-Indic.5y.in.Xy*(HH-60)*(HH>60))/12) # because in that case (indic.average==1), compute_AB
    #                            considers E(exp[u'{X{+HH}+X{+HH-12}+X{+HH-24}...+}])
    # =========================
    # =========================
  }else if(indic.average==2){# in that case (indic.average==2), compute_AB
    #                            considers E(exp[u'{X{+HH}+X{+HH-3}+X{+HH-6}+X{+HH-9}}])
    Y <- rep(y*4,length(HH))
  }else{
    Y <- rep(y,length(HH))
  }
  
  #### STEP 2: CALCULATIONS ####
  
  # psi.B contains all the B%*%X_t. psi.B is a matrix of dimension (T*dim(v)) x nb.HH.
  # This implies that each column contains all the B%*%X_t in vector.
  # i.e as.vector(X %*% AB.list$B[,,1]) == psi.B[,1] (for one horizon)  
  psi.B <- apply(AB.list$B,3,function(layer){X %*% layer})
  
  # Initialize the count for the loop
  count <- 0
  
  # extend v to start at 0 and finish at max.v.
  # This is done to calculate the epsilon namely first = v[1] - 0.
  # The final dimension of epsilon will be 1 x dim(v).
  extended.v <- cbind(0,v,max.v)
  dim.ext.v <- dim.v + 2
  
  # Weights to compute the integral, i.e the epsilon in the Riemann sum.
  # integral(f(v)dv) = sum(epsilon f(v))
  # Compute weights as 0.5(epsilon(j+1)- epsilon(j)) + 0.5(epsilon(j+2)-epsilon(j+1)).
  weights <- matrix(.5*(extended.v[2:(dim.ext.v-1)]-extended.v[1:(dim.ext.v-2)])+
                      .5*(extended.v[3:dim.ext.v]-extended.v[2:(dim.ext.v-1)]),nrow=1)
  
  # Loop for each horizon considered
  for(h in HH){
    
    # Increment the count
    count <- count + 1
    
    # Compute psi.h= exp(A_h + B_h X_t). psi.h is of dimension T x dim(v).
    # Each column contain the info of one increment for the CDF.
    psi <- exp(
      matrix(1,T,1) %*% AB.list$A[,,count] + matrix(psi.B[,count],T,dim.v)
    )
    
    # Compute Im(psi.h(...)e^(-i*v*y)). numerator is of dimension T x dim(v).
    numerator <- Im(psi * (matrix(1,T,1) %*% exp(-1i*v*Y[count])))
    
    # Compute the integrand part (just divide by v). integrand is of dimension T x dim(v).
    integrand <- numerator / (matrix(1,T,1) %*% v)
    
    # Compute probability. To do so, we have to multiply by the weights (Riemann sum).
    Proba <- cbind(Proba, .5 - 1/pi * integrand %*% t(weights))
  }
  
  ## Return the probability => dimension is T x dim(HH).
  ## Where dim(HH) is the number of horizon considered and T the time periods.
  return(Proba)
}

# Function that computes the model implied distribution (pdf and CDF) for all time.
# The function uses the implied probability => one bins of the Riemann sum.
# WARNING : THIS IS DONE FOR ONLY ONE HORIZON.
# PDF and CDF obtained are of dimension T x (dim(vec.y)-2)
# Get rid of 2 y's to do the increment.
# PDF.x = vec.y without first and last observation. 
compute.distri.plus.stdv <- function(Model,X,b,HH,  # used in compute.analytical.distri
                                     min.bx,max.bx,step.distri,
                                     step.4.integral,max.v.4.integral,
                                     indic.average=0,
                                     Sigma.X=0,
                                     AB.list.external=NaN){
  # Computes the distribution of b'X
  # HH should be of dimension 1.
  # X is of dimension T x (Model$n + Model$q)
  # Sigma.X is of dimension T x ((Model$n + Model$q) x (Model$n + Model$q))
  # If Sigma.X = 0, no standard deviation for the distri function is computed.
  # min.bx (mas.bx) is the minimun (maximum) considered value of the random variable b'X
  # The distrinbution function is evaluated J points, where J = (max.bx - min.bx)/step.distri
  # The output is a list that contains:
  #     - the matrix (PDF) of evaluations of the distribution function (dimension T x J).
  #       This means that PDF[1,] is the pdf for the first period considered
  #     - the matrix (PDF.stdv) of std. dev, of the evaluations of the distribution function (dimension T x J)
  #       (this stdv takes into account the uncertainty of the Xs, as relfected in Sigma.X)
  #     - the vector (PDF.x) of points where the distribution is evaluated
  
  # vec.y is the vector of y, that is P(gamma_i X_(t+i) < y).
  vec.y <- seq(min.bx,max.bx,by=step.distri)
  CDF <- NULL
  
  # Loop for the different => different bins of the Riemann sum.
  for(iii in 1:length(vec.y)){
    y <- vec.y[iii]
    #print(y)
    
    # Function that computes the model implied probability for all time.
    # The function return the implied probability => one point of the distribution.
    # This is done for one y => calculate the probability with the same y for T periods and HH horizons.
    # The loop in the script make HH as a scalar => one horizon computed at a time.
    aux <- compute.proba(Model,X,b,y,HH,
                         step.4.integral,
                         max.v.4.integral,
                         indic.average,
                         AB.list.external,
                         Indic.5y.in.Xy)
    
    # Compute CDF:
    ## - each column contains P(gamma_i X_(t+i) < y) for a different y and all time period (for a given horizon (maturity)).
    ## - each line contains a distribution for a considered time and horizon (maturity).
    CDF <- cbind(CDF,aux)
  }
  
  # Compute the PDF:
  # [F(x2)-F(x1)]/[x2-x1]
  #### CHECK THAT !!!!!!!!
  PDF <- (CDF[,2:dim(CDF)[2]] - CDF[,1:(dim(CDF)[2]-1)])/(vec.y[2:(length(vec.y))]-vec.y[1:(length(vec.y)-1)])  # [F(x2)-F(x1)]/[x2-x1]
  PDF.x <- vec.y[1:(length(vec.y)-1)] + step.distri/2 + Model$pi.bar[area.var]*freq # add freq
  #PDF <- (CDF[,3:dim(CDF)[2]] - CDF[,1:(dim(CDF)[2]-2)])/2  # [F(x2)-F(x1)]/[x2-x1]
  #PDF.x <- vec.y[2:(length(vec.y)-1)]
  if(length(Sigma.X)==1){
    Grad.PDF <- NULL
    Grad.CDF <- NULL
    PDF.stdv <- NULL
    CDF.stdv <- NULL
  }else{
    # Computation of the stdv of the proba by the Delta method.
    
    # ================================
    # Step 1: Compute the gradient of the distri matrix wrt each of the X ->
    # The gradient is put in a 3d array of dimension T x J x (n + q)
    Grad.PDF <- array(NaN,c(dim(X)[1],length(PDF.x),Model$n+Model$q))
    Grad.CDF <- array(NaN,c(dim(X)[1],length(PDF.x)+2,Model$n+Model$q))
    
    # computation of the marginal standard deviations of the Xs
    # (will be used to calibrate the perturbation when computing the gradient of the proba wrt the Xs)
    stdv.X <- sqrt(diag(unc.moments.X(Model)$X.mom.2))
    
    for(dim.X in 1:(Model$n + Model$q)){
      epsilon <- stdv.X[dim.X]/100
      CDF.aux <- NULL
      for(iii in 1:length(vec.y)){
        y <- vec.y[iii]
        print(c(dim.X,y))
        X.aux <- X
        X.aux[,dim.X] <- X.aux[,dim.X] + epsilon
        aux <- compute.proba(Model,X.aux,b,y,HH,step.4.integral,
                             max.v.4.integral,
                             indic.average,
                             AB.list.external,
                             Indic.5y.in.Xy)
        CDF.aux <- cbind(CDF.aux,aux)
      }
      #### CHECK THAT !!!!!!!!
      Grad.PDF[,,dim.X] <- (
        (CDF.aux[,3:dim(CDF)[2]] - CDF.aux[,1:(dim(CDF)[2]-2)])/2
        - PDF
      )/epsilon
      Grad.CDF[,,dim.X] <- (CDF.aux - CDF)/epsilon
    }
    PDF.stdv <- NULL
    CDF.stdv <- NULL
    for(x in 1:dim(Grad.PDF)[2]){
      Grad.PDF.aux <- Grad.PDF[,x,]
      Grad.Grad.PDF <- t(apply(Grad.PDF.aux,1,function(x){x %x% x})) # each row is the outer-product of the rows of Grad
      aux <- apply(Grad.Grad.PDF * Sigma.X,1,sum)  # delta method: var[f(x)]=[f'(x)]^2 * sigma_X^2
      PDF.stdv <- cbind(PDF.stdv,sqrt(aux))
      Grad.CDF.aux <- Grad.CDF[,x+1,]
      Grad.Grad.CDF <- t(apply(Grad.CDF.aux,1,function(x){x %x% x})) # each row is the outer-product of the rows of Grad
      aux <- apply(Grad.Grad.CDF * Sigma.X,1,sum)
      CDF.stdv <- cbind(CDF.stdv,sqrt(aux))
    }
  } #CDF[,1:(dim(CDF)[2]-1)
  return(list(CDF=CDF[,1:(dim(CDF)[2]-1)],PDF=PDF,PDF.x=PDF.x,
              PDF.stdv=PDF.stdv,
              CDF.stdv=CDF.stdv,
              Grad.PDF=Grad.PDF,
              Grad.CDF=Grad.CDF))
}

# ============================================
# O.  MODEL IMPLIED CORRELATIONS
# ============================================


compute.lambda.h <- function(Model, h){
  # This compute compute the lambda.t for all horizon h.
  #Model: list with all the needed parameters
  #h : max horizon, we should compute the lambda.t
  
  #Create the array to store the results
  lambda.0.h <- array(NaN,c(dim(Model$Lambda.0)[1],dim(Model$Lambda.0)[2],h))
  lambda.1.h <- array(NaN,c(dim(Model$Lambda.1)[1],dim(Model$Lambda.1)[2],h))
  
  #Calculate Lambda.1
  lambda.0.h[,,1] <- Model$Lambda.0
  lambda.1.h[,,1] <- Model$Lambda.1
  
  #Loop to compute the all the lambda.h until h.
  for (i in 2:h) {
    
    lambda.0.h[,,i] <- lambda.0.h[,,i-1] + lambda.1.h[,,i-1]%*%Model$mu.X + (Model$Phi.X^(i-1) %x% 
                                                                               Model$Phi.X^(i-1))%*%Model$Lambda.0
    lambda.1.h[,,i] <- lambda.1.h[,,i-1]%*%Model$Phi.X + (Model$Phi.X^(i-1) %x% Model$Phi.X^(i-1))%*%Model$Lambda.1
    
  }
  
  return(list(lambda.0.h=lambda.0.h, lambda.1.h=lambda.1.h))
  
}



compute.var.cov.X.tT <- function(Model, X.tT, HH){
  #This function compute the variance covariance matrice for X.t+h in vector form
  #Model: list with all the needed parameters
  #X.tT: latent variables, X_t (of dimension Tx(n+q))
  #HH: horizon of interest to compute the variance-covariance matrice of X_t 
  
  h <- max(HH)
  HHH <- length(HH)
  all.lambda <- compute.lambda.h(Model,h)
  
  var.cov.vector <- array(NaN,c(dim(Model$Lambda.1)[1],dim(X.tT)[1],HHH))
  count <- 1
  
  for (i in HH) {
    
    var.cov.vector[,,count] <- (all.lambda$lambda.0.h[,,i] + all.lambda$lambda.1.h[,,i]%*%t(X.tT))
    count <- count + 1   
    
  }
  
  return(var.cov.vector)
  #var.cov.vector <- matrix(var.cov.vector,dim(var.cov.vector)[1]^(1/2), dim(var.cov.vector)[1]^(1/2))
  
}

# Compute the modelled expectation and variance using loading
compute.two.first.moments.loadings <- function(Model, a.h, b.h){
  
  n <- Model$n
  q <- Model$q
  S.n <- make.S.p(n)
  S.q <- make.S.p(q)
  # Compute unconditional moment of X
  unc.mom.X <- unc.moments.X(Model,S.n,S.q)
  ## Unconditional Expectation of X
  E.X <- as.matrix(unc.mom.X$X.mom.1)
  ## Unconditional Variance of X
  V.X <- unc.mom.X$X.mom.2
  
  ## Unconditional Expectation of cumulant
  E.cum <- c(a.h + t(b.h)%*%E.X)
  ## Unconditional Variance of cumulant
  V.cum <- diag(t(b.h)%*%V.X%*%b.h)
  
  return(list("Cum.mom.1"=E.cum,"Cum.mom.2"=V.cum))
  
}
