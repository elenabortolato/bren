
rm(list=ls())
require(mnormt)
require(gee)
require(parallel)
require(MLGdata) 


## function to estimates an equicorrelated normal model with regression parameters
## using various bias reduction methods .

# bren
# Input:
# y : (n*q x 1) vector of responses (long dataset format), n sample size, q variate response (replication over the same unit)
# x : (n*q x p) design matrix
# id : a vector of integers in [1,n] indicating the individual 
# method: "ML" maximum likelihood, 
#         "correction" for Bias correction,
#          "AS_mean" Bias reduction (Firth)
#          "AS_median" Median bias reduction (Kenne et al.)
#          "MPL_Jeffreys" Jeffreys-type penalization
# a: power for Jeffreys penalization, a=-0.5 default for Jeffreys penalization
# maxit: Newton Raphson iterations
# epsilon: tolerance for convergenceof Newton Raphson
# max_step_factor: length of the sequence of rescaling factors for updating the parameters (2^-1, 2^-2, ... ,2^-max_step_factor)
# slowit: 1 for Newton Raphson, <1 for smaller steps
# Output:
# estimate = beta coefficients, sigma2, rho
# sigma2 
# rho
# se: standard errors
# converged : flag for convergence
# iter : iterations to converge 
# grad : values of the score components at final iteration
# method : method used, as in Input
 
bren <- function(x, y, id, method = "ML", a = 0.5,
                  maxit = 1000, epsilon = 10^(-11), max_step_factor = 11, slowit = 1) {
  x <- as.matrix(x)
  p <- ncol(x)
  n <- length(table(id))
  q <- as.vector(table(id)[1])
  #initial estimate
  betas <-  solve(crossprod(x,x)) %*% (t(x) %*% y)
  #prediction
  mus <- drop(x %*% betas)
  #residuals
  residuals <- y - mus
  residualsmat <- matrix(residuals, nrow = n, byrow = T)
  globalmean <- mean(residuals)
  #within sum of squares
  SW <- apply(residualsmat, 1, function (x) sum((x - mean(x))^2))
  SSW <- sum(SW)
  subjectmeans <- apply(residualsmat, 1, mean)
  #between sum of squares
  SSB <-  sum((subjectmeans - globalmean)^2)
  #initial estimates of sigma and rho
  sigma2 <- (SSW + q * (SSB))/(n * q)
  rho <- (1 - q/(q - 1)*SSW/(SSW + q*SSB) )
  
  
  #update beta estimates with covariance structure
  varcov= varcov=sigma2*(matrix(rho,ncol=q,nrow=q)+diag(1-rho,q))
  betas  <-  (solve(t(x) %*% (diag(n)%x% solve(varcov) )%*% x)) %*% (t(x) %*% ( (diag(n)%x%  solve(varcov) )%*% y))
  mus <- drop(x %*% betas)
  residuals <- y - mus
  residualsmat <- matrix(residuals, nrow = n, byrow = T)
  globalmean <- mean(residuals)
  SW <- apply(residualsmat, 1, function (x) sum((x - mean(x))^2))
  SSW <- sum(SW)
  subjectmeans <- apply(residualsmat, 1, mean)
  SSB <-  sum((subjectmeans - globalmean)^2)
  sigma2 <- (SSW + q * (SSB))/(n * q) 
  rho <-  (1 - q/(q - 1)*SSW/(SSW + q*SSB) )
  if (rho <= -1/(q-1)) rho <- 0; flag = 1
  pars <- c(sigma2, rho) 
  
  
  ## score function of sigma2 and rho
  score <- function(pars, y, mus) {
    sigma2 <- pars[1]
    rho <- pars[2]
    omega_s2 <- matrix( -rho/((rho - 1) * (q * rho - rho + 1) * sigma2^2),ncol = q, nrow = q)
    diag(omega_s2) <-  (q * rho - 2 * rho + 1)/((rho - 1) * (q * rho - rho + 1) * sigma2^2)
    omega_rho <- matrix(-(q * rho^2 - rho^2 + 1)/((rho - 1)^2 * (q * rho - rho + 1)^2 * sigma2) ,
                        ncol = q, nrow = q)
    diag(omega_rho) <- ((q - 1) * rho * (q * rho - 2 * rho + 2))/((rho - 1)^2 * (q * rho - rho + 1)^2 * sigma2)
    #return two components
    c( ( -n*q/(2 * sigma2) -1/2 * sum(t(y - mus) %*% ( diag(n)%x%omega_s2) %*% (y - mus)) ),
       (-(n * (q - 1) * q * rho)/(2 * (rho - 1) * (q * rho - rho + 1))-1/2 * sum(t(y - mus) %*% ( diag(n)%x%omega_rho) %*% (y - mus)) )
    )
  }
  
  # expected info
  information <- function(pars, level = 0, inverse = FALSE) {
    sigma2 <- pars[1]
    rho <- pars[2]
    if (level == 0) {
      Iq <- diag(q)
      ones <- matrix(1, q, q)
      omega <- (1/(sigma2 * (1 - rho)))*(Iq - ones * rho/(1 + rho * (q - 1)))
      if (inverse) {
        return(solve(t(x) %*%  (diag(n)%x%omega) %*% x))
      }
      else {
        return((t(x) %*%  (diag(n)%x%omega) %*% x))
      }
    }
    if (level == 1) {
      #block matrix sigma-rho 
      ans <- matrix(NA, 2, 2)
      if (inverse) {
        ans[1,1] <- ( (2 * (q * rho^2 - rho^2 + 1) * sigma2^2)/(n * q) )
        ans[1,2] <- ans[2,1] <- ( -(2 * (rho - 1) * rho * (q * rho - rho + 1) * sigma2)/(n * q) ) 
        ans[2,2] <- ( (2 * (rho - 1)^2 * (q * rho - rho + 1)^2)/(n * (q - 1) * q) )
        return(ans)
      }
      else {
        ans[1,1] <- n*q/(2 * sigma2^2)
        ans[1,2] <- ans[2,1] <- ( (n * q * (q - 1) * rho)/(2 * (rho - 1) * (q * rho - rho + 1) * sigma2) )
        ans[2,2] <- ( n * (q - 1) * q * (q * rho^2 - rho^2 + 1)/(2 * (rho - 1)^2 * (q * rho - rho + 1)^2) )
        return(ans)  
      }
    }
  }
  
  # Adjustment functions
  # Jeffreys adjustment 
  as_Jeffreys_adjustment <-
    function(pars) {
      sigma2 <- pars[1]
      rho <- pars[2]
      w_rho_d = ((q - 1) * rho * (q * rho - 2 * rho + 2)) / ((rho - 1) ^ 2 *
                                                               (q * rho - rho + 1) ^ 2 * sigma2)
      w_rho_od = -(q * rho ^ 2 - rho ^ 2 + 1) / ((rho - 1) ^ 2 * (q * rho -
                                                                    rho + 1) ^ 2 * sigma2)
      Omega_rho = matrix(w_rho_od, ncol = q, nrow = q) - diag(w_rho_od - w_rho_d, ncol =
                                                                q, nrow = q)
      x.omega.rho.x =  (t(x) %*%  (diag(n) %x% Omega_rho) %*% x)
      d = det(x.omega.rho.x)
      #return two components
      a * c(-(2 + p) / sigma2, -(4 * q * rho - 4 * rho - 2 * q + 4) / ((rho - 1) * (q * rho - rho + 1)) -
              log(abs(d)))
    }
  # mean adjustment
  as_mean_adjustment <- function(pars) {
    sigma2 <- pars[1]
    rho <- pars[2]
    #off diagonal elements of Omega
    w0od = rho / ((rho - 1) * (q * rho - rho + 1) * sigma2)
    #diagonal elements of Omega
    w0d = -(q * rho - 2 * rho + 1) / ((rho - 1) * (q * rho - rho + 1) *
                                        sigma2)
    Omega = matrix(w0od, ncol = q, nrow = q) - diag(w0od - w0d, nrow = q, ncol =
                                                      q)
    
    #derivarives of Omega
    w_rho_d = ((q - 1) * rho * (q * rho - 2 * rho + 2)) / ((rho - 1) ^ 2 *
                                                             (q * rho - rho + 1) ^ 2 * sigma2)
    w_rho_od = -(q * rho ^ 2 - rho ^ 2 + 1) / ((rho - 1) ^ 2 * (q * rho -
                                                                  rho + 1) ^ 2 * sigma2)
    Omega_rho = matrix(w_rho_od, ncol = q, nrow = q) - diag(w_rho_od - w_rho_d, ncol =
                                                              q, nrow = q)
    x.omega.rho.x =  (t(x) %*%  (diag(n) %x% Omega_rho) %*% x)
    x.omega.x =  solve(t(x) %*%  (diag(n) %x% Omega) %*% x)
    xprodx = x.omega.rho.x %*% x.omega.x
    c(
      -(2 * q * rho ^ 2 - 2 * rho ^ 2 - p) / (2 * sigma2),
      #-(((q - 1) * (2 * q * rho^3 - 2 * rho^3 - p * rho + 2 * rho + p))/(2 * (rho - 1) * (q * rho - rho + 1)))-(q -
      1) * rho * (q * rho ^ 2 - rho ^ 2 + 1) / ((rho - 1) * (q * rho - rho + 1)) -
      0.5 * sum(diag(xprodx))
    
    )
  }
  
  ## median adjustment
  as_median_adjustment <- function(pars) {
    sigma2 <- pars[1]
    rho <- pars[2]
    
    as_mean_adjustment(pars) - c(
      -(
        3 * q ^ 2 * rho ^ 4 - 6 * q * rho ^ 4 + 3 * rho ^ 4 + 6 * q * rho ^ 2 -
          6 * rho ^ 2 - q * rho + 2 * rho + 1
      ) / (3 * (q * rho ^ 2 - rho ^ 2 + 1) * sigma2),-(
        3 * q ^ 3 * rho ^ 5 - 9 * q ^ 2 * rho ^ 5 + 9 * q * rho ^ 5 - 3 * rho ^
          5 + 9 * q ^ 2 * rho ^ 3 - 18 * q * rho ^ 3 + 9 * rho ^ 3 - 2 * q ^ 2 * rho ^
          2 + 6 * q * rho ^ 2
        - 4 * rho ^ 2 + 4 * q * rho - 4 * rho - q + 2
      ) / (3 * (rho - 1) * (q * rho - rho + 1) * (q * rho ^ 2 - rho ^ 2 + 1))
    )
  }
  
  ## Bias ###
  
  bias <- function(pars) {
    sigma2 <- pars[1]
    rho <- pars[2]
    c(-(p * (q * rho - rho + 1) * sigma2)/(n * q),
      ((rho - 1) * (2 * rho + p) * (q * rho - rho + 1))/(n * q)
    )
  }
  
  ##chose adjustment term function with method
  
  adjustment_function <- switch(method, 
                                AS_mean = as_mean_adjustment, 
                                AS_median = as_median_adjustment,
                                MPL_Jeffreys = as_Jeffreys_adjustment,
                                correction = function(pars, ...) 0,
                                ML = function(pars, ...) 0)
  
  inverse_info <- information(pars, level = 1, inverse = TRUE)
  adjusted_grad <- score(pars, y, mus) + adjustment_function(pars)
  step <- drop(inverse_info %*% adjusted_grad)
  
  if (maxit == 0) {
    iter <- 0
    failed <- FALSE
  }
  else 
  {
    for (iter in seq.int(maxit)) { #start iterations
      step_factor <- 0
      testhalf <- TRUE
      while (testhalf & step_factor < max_step_factor) {
        step_previous <- step
        pars <- pars + slowit * 2^(-step_factor) * step
        varcov= varcov=pars[1]*(matrix(pars[2],ncol=q,nrow=q)+diag(1-pars[2],q))
        betas_new <-  (solve(t(x) %*% (diag(n)%x% solve(varcov) )%*% x)) %*% (t(x) %*% ( (diag(n)%x%  solve(varcov) )%*% y))
        mus <- drop(x %*% betas_new)
        
        inverse_info <- information(pars, level = 1, inverse = TRUE)
        adjusted_grad <- score(pars, y = y, mus = mus) + adjustment_function(pars)
        step <- drop(inverse_info %*% adjusted_grad)
        
        
        
        if (step_factor == 0 & iter == 1) {
          testhalf <- TRUE
        }
        else {
          testhalf <- sum(abs(step), na.rm = TRUE) > sum(abs(step_previous),  na.rm = TRUE)
        }
        step_factor <- step_factor + 1
      }
      if (( all(abs( adjusted_grad ) < epsilon)) & all((betas_new-betas)<10^-11)) {
        break
      }
      betas=betas_new
    }
  }
  if(iter >= maxit) {
    convergence <- FALSE
    warning("optimization failed to converge")
  }
  else
  {
    convergence <- TRUE
  }
  if (method == "correction") 
  {
    pars <- pars - bias(pars)
  }
  
  inverse_info_betas <- information(pars, level = 0, inverse = TRUE)
  inverse_info_sigma2rho <- information(pars, level = 1, inverse = TRUE)
  se <- c(sqrt(diag(inverse_info_betas)),
          sqrt(diag(inverse_info_sigma2rho)))
  
  list(estimate = c(drop(betas),pars[1],pars[2]), sigma2 = pars[1], rho = pars[2], se = se, converged = convergence,
       iter =iter, grad = c(rep(0,p),adjusted_grad ), method = method)
  
}







############################################
#Functions to simulate/output
############################################
all_sim=function(q=q,p=p,n=n,pars=pars,x=xx, mis=FALSE, cor=NULL){
  
  beta=pars[1:(p-2)]; s2=pars[p-1]; rho=pars[p]
  varcov=s2*(matrix(rho,ncol=q,nrow=q)+diag(1-rho,q))
  id=rep(1:n,each=q)
  
  if(mis==TRUE)  {
    varcov=s2*cor
  }
  
  e=rmnorm(n=n, mean = c(rep(0,q)), varcov = varcov)
  e=as.vector(t(e))
  y=e+x%*%beta
  mle=bren( y = y, x=x , id=id)
  bc=bren( y = y, x=x , id=id,  method = "correction" )
  jeff_n05=bren( y = y, x=x , id=id, method = "MPL_Jeffreys",   a=-0.5)
  mbr=bren( y = y, x=x , id=id,  method = "AS_median" )
  br=bren(y = y, x=x , id=id, method = "AS_mean"  )
  rob=gee(y~ .,data=as.data.frame(x[,-1]), id=id,corstr = "exchangeable",
          maxiter = 1000,tol = 10^(-11),silent = T)
  
  # library(nlme)
  #lme <- lme(y ~ x[,2]+x[,3]+x[,4]+x[,5], random = ~ 1 | id)
  #                   data = Stroke1)
  return(list(samples=y,mle=mle,bc=bc, jeff_n05=jeff_n05, br=br, mbr=mbr,rob=rob, sim.set=pars))
} 


all_sim_id=function(q=q,p=p,n=n,pars=pars,x=xx, mis=FALSE, cor=NULL){
  
  beta=pars[1:(p-2)]; s2=pars[p-1]; rho=pars[p]
  varcov=s2*(matrix(rho,ncol=q,nrow=q)+diag(1-rho,q))
  id=rep(1:n,each=q)
  
  if(mis==TRUE)  {
    varcov=s2*cor
  }
  
  e=rmnorm(n=n, mean = c(rep(0,q)), varcov = varcov)
  e=as.vector(t(e))
  y=e+x%*%beta
  mle=bren( y = y, x=x , id=id)
  bc=bren( y = y, x=x , id=id,  method = "correction" )
  jeff_n05=bren( y = y, x=x , id=id, method = "MPL_Jeffreys",   a=-0.5)
  mbr=bren( y = y, x=x , id=id,  method = "AS_median" )
  br=bren(y = y, x=x , id=id, method = "AS_mean"  )
  rob=gee(y~ 1 , id=id,corstr = "exchangeable",
          maxiter = 1000,tol = 10^(-11),silent = T)
  
  # library(nlme)
  #lme <- lme(y ~ x[,2]+x[,3]+x[,4]+x[,5], random = ~ 1 | id)
  #                   data = Stroke1)
  return(list(samples=y,mle=mle,bc=bc, jeff_n05=jeff_n05, br=br, mbr=mbr,rob=rob, sim.set=pars))
} 
out_all=function(input,p=p, alpha=0.05){
  #converged
  
  sim.set=input[[1]]$sim.set 
  sim.set=sim.set[1:p]
  n=length(input)
  
  mle=matrix(unlist(lapply(input, function (x) x$mle$estimate   )),ncol=p, nrow=n, byrow=T)
  bc=matrix(unlist(lapply(input, function (x) x$bc$estimate   )),ncol=p, nrow=n, byrow=T)
  jeff_n05=matrix(unlist(lapply(input, function (x) x$jeff_n05$estimate  )),ncol=p, nrow=n,byrow=T)
  br=matrix(unlist(lapply(input, function (x) x$br$estimate   )),ncol=p, nrow=n, byrow=T)
  mbr=matrix(unlist(lapply(input, function (x) x$mbr$estimate   )),ncol=p, nrow=n, byrow=T)
  gee=matrix(unlist(lapply(input, function (x) c(x$rob$coefficient,x$rob$scale, min(x$rob$working.correlation)))),
             ncol=p, nrow=n, byrow=T)
  rob=matrix(unlist(lapply(input, function (x) c(x$rob$coefficient,x$rob$scale, min(x$rob$working.correlation)))),
             ncol=p, nrow=n, byrow=T)
  #appartenenza a spazio param 
  
  se_mle=matrix(unlist(lapply(input, function (x) x$mle$se  )),ncol=p, nrow=n,byrow=T)
  se_bc=matrix(unlist(lapply(input, function (x) x$bc$se  )),ncol=p, nrow=n,byrow=T)
  se_jeff_n05=matrix(unlist(lapply(input, function (x) x$jeff_n05$se   )) ,ncol=p, nrow=n,byrow=T)
  se_br=matrix(unlist(lapply(input, function (x) x$br$se  )),ncol=p, nrow=n,byrow=T)
  se_mbr=matrix(unlist(lapply(input, function (x) x$mbr$se  )),ncol=p, nrow=n,byrow=T)
  se_gee=matrix(unlist(lapply(input, function (x) c(sqrt(diag(x$rob$naive.variance)),NA,NA))),ncol=p, nrow=n,byrow=T)
  se_rob=matrix(unlist(lapply(input, function (x) c(sqrt(diag(x$rob$robust.variance)),NA,NA))),ncol=p, nrow=n,byrow=T)
  
  err_mle=sweep(mle,2,  sim.set, "-")
  err_bc=sweep(bc,2,  sim.set, "-")
  err_jeff_n05=sweep(jeff_n05,2,  sim.set, "-")
  err_br=sweep(br,2,  sim.set, "-")
  err_mbr=sweep(mbr,2,  sim.set, "-")
  err_gee=sweep(gee,2,  sim.set, "-")
  err_rob=sweep(rob,2,  sim.set, "-")
  
  bias_mle=colMeans(err_mle)  
  bias_bc=colMeans(err_bc)
  bias_jeff_n05=colMeans(err_jeff_n05)  
  bias_br=colMeans(err_br)    
  bias_mbr=colMeans(err_mbr)  
  bias_gee=colMeans(err_gee)  
  bias_rob=colMeans(err_rob)  
  
  bias=cbind(bias_mle, bias_bc, bias_jeff_n05,  bias_br,  bias_mbr, bias_gee, bias_rob)
  
  var_mle=apply(mle,2, function (x) var(x ))
  var_bc=apply(bc,2, function (x) var(x ))
  var_jeff_n05=apply(jeff_n05,2, function (x) var(x )) 
  var_br=apply(br,2, function (x) var(x ))
  var_mbr=apply(mbr,2, function (x) var(x ))
  var_gee=apply(gee,2, function (x) var(x ))
  var_rob=apply(rob,2, function (x) var(x ))
  
  perc_underest=100*cbind(colMeans((err_mle) < 0),
                          colMeans((err_bc) < 0), 
                          colMeans((err_jeff_n05) < 0),
                          colMeans((err_br) < 0),
                          colMeans((err_mbr) < 0),
                          colMeans((err_gee) < 0),
                          colMeans((err_rob) < 0))
  
  mse_mle=colMeans(err_mle^2)
  mse_bc=colMeans(err_bc^2)
  mse_jeff_n05=colMeans(err_jeff_n05^2)
  mse_br=colMeans(err_br^2)
  mse_mbr=colMeans(err_mbr^2)
  mse_gee=colMeans(err_gee^2)
  mse_rob=colMeans(err_rob^2)
  
  rel_bias=100*cbind( bias_mle/sim.set, 
                      bias_bc/sim.set,
                      bias_jeff_n05/sim.set ,
                      bias_br/sim.set,
                      bias_mbr/sim.set,
                      bias_gee/sim.set,
                      bias_rob/sim.set)
  
  #for(i in 1:p)   if(is.infinite(rel_bias[i])) rel_bias[i,]=bias[i,] #mu=0 or rho=0
  wald=100*cbind(colMeans(abs(err_mle)/se_mle < qnorm(1-alpha/2)),
                 colMeans(abs(err_bc)/se_bc < qnorm(1-alpha/2)),
                 colMeans(abs(err_jeff_n05)/se_jeff_n05 < qnorm(1-alpha/2)),
                 colMeans(abs(err_br)/se_br < qnorm(1-alpha/2)),
                 colMeans(abs(err_mbr)/se_mbr < qnorm(1-alpha/2)),
                 colMeans(abs(err_gee)/se_gee < qnorm(1-alpha/2)),
                 colMeans(abs(err_rob)/se_rob < qnorm(1-alpha/2)))                
  
  ibmse=100*cbind(bias_mle^2/(var_mle) ,
                  bias_bc^2/(var_bc) ,
                  bias_jeff_n05^2/(var_jeff_n05),
                  bias_br^2/(var_br), 
                  bias_mbr^2/(var_mbr) ,
                  bias_gee^2/(var_gee) ,
                  bias_rob^2/(var_rob) )                 
  
  
  out=rbind( bias, perc_underest, rel_bias, wald, ibmse)
  colnames(out)=c("MLE", "BC",   "J(-0.5)","BR", "MBR","gee", "rob")
  rownames(out)=c( "bias",rep(" ", p-1), "PU", rep(" ", p-1),   "RB", rep(" ", p-1),"WALD",rep(" ", p-1), "IBMSE",rep(" ", p-1))
  out
  
}
############################################
#Simulations
R=10^4
############################################

#Case 1 : investigation on q, sigma^2, rho, N=20,40, MU constant

n_=c(20,100)
q_=c(2,5)
rho_=c(-0.1,0.2,0.9)
sigma2_=c(1,5)

case1=NULL
for(i in 1:length(n_)){
  for(j in 1:length(q_)){
    for(k in 1:length(sigma2_)){
      for (l in 1:length(rho_)) {
        case1=rbind(case1, cbind(n_[i], q_[j], sigma2_[k], rho_[l]))
        
      }
    }
  }
}

case1
#R=10000

case1_ris=NULL
for(i in 1:nrow(case1)){
  
  rho=case1[i,4]
  s2=case1[i,3]
  beta=c(1)
  
  q=case1[i,2]
  n=case1[i,1]
  
  pars=c(beta,s2,rho)
  set=c(pars, q,n)
  
  p=3
  # Xmatrix  
  set.seed(123)
  x=cbind(c(rep(1,n)))
  x=matrix(apply(x, 1 ,function(r) rep(t(r),(q))),ncol=(p-2), byrow=T)
  
  est=rep(NA,R);risultati_all=list(NULL);ris_out=NULL
  
  system.time(
    if(0==0){
      seed0=(987654321)
      set.seed(seed0)
      cl <- makeCluster(getOption("cl.cores", 8))
      clusterExport(cl, c("breqn","all_sim_id","pars","x","p","q","n","rmnorm","gee"))
      risultati_all =clusterApply(cl,est,fun =   function (z)  
        all_sim_id(cor = NULL,p=p,q=q,n=n,x = x, pars=pars, mis=FALSE))
    })
  
  case1_ris[[i]]=out_all(risultati_all,p = 3)
  cat(i)
}

 
#saveRDS(case1_ris, "iid10000_brenq.RDS")



for (i in 1:nrow(case1)) {
  print(case1[i,])
  print(knitr::kable(case1_ris[[i]],digits=2))
  
}


plot_summary_iid=function(ris,p=3,listcase){
  
  l=length(ris);  listcase=as.data.frame(listcase)
  tmp=data=NULL;  ris2=list()
  
  for( i in 1:l){
    tmp=ris[[i]]
    ris2 = NULL
    for(j in 1:7){
      cc=cbind(tmp[,c(j )],listcase[i,] , colnames(tmp)[j], 
               c(rep("bias",3),
                 rep("PU",3),
                 rep("RB",3),
                 rep("WALD",3),
                 rep("IBMSE",3)))
      colnames(cc)=NA
      ris2  =rbind(ris2 , cc)
    }
    data=rbind(data, ris2 )
  }
  data=cbind(data, c(rep(c("mu", "s2","rho"),nrow(data)/3)))
  colnames(data)=c("val", "n", "q","s2", "rho", "Estimator", "stat","par")
  data$Estimator <- factor(data$Estimator, levels = c("MLE","gee", "rob","BC", "BR", "MBR", "J(-0.5)"))
  
  dataW=data[data$stat=="WALD",]
  dataW=dataW[dataW$par=="mu",]
  dataW20=dataW[dataW$n==20,]
  p2=ggplot(dataW20 , aes(  y=val,x=Estimator  )) +
    theme_minimal() +  ggtitle("WALD n=20")+     #geom_boxplot( ) +    # scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"))+
    labs(x = "  ",y = " ", labeller = label_parsed) + 
    #theme(legend.position = "none") +
    geom_hline( yintercept = 95, size=0.5, color="black")+
    #  geom_vline( xintercept=xref, size=0.5, color="black")+
    geom_point(aes(shape=factor(q)), size=1.2,   alpha=1,
               position = position_jitterdodge(dodge.width = 0.7,jitter.width = 0,jitter.height = 0)) 
  
  dataW40=dataW[dataW$n==40,]
  p3=ggplot(dataW40 , aes(  y=val,x=Estimator  )) +
    theme_minimal() +  ggtitle("WALD n=40")+
    #geom_boxplot( ) +    # scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"))+
    labs(x = "  ",y = " ", labeller = label_parsed) + 
    #theme(legend.position = "none") +
    geom_hline( yintercept = 95, size=0.5, color="black")+
    #  geom_vline( xintercept=xref, size=0.5, color="black")+
    geom_point(aes(shape=factor(q)), size=1.2,   alpha=1,
               position = position_jitterdodge(dodge.width = 0.7,jitter.width = 0,jitter.height = 0)) 
  
  grid.arrange(p2, p3, ncol=2)
}


library(ggplot2)
library("gridExtra")
#case1_ris=readRDS(file.choose()) #id1000.RDS
#iid=plot_summary_iid(case1_ris, listcase = case1)


###############################################################
#COVARIATES FIXED
#set2

n_=c(20,50)
q_=c( 5, 10)
rho_=c(0.9)
sigma2_=c(5)
mis_=c(F,T,T,T)

case2=NULL
for(i in 1:length(n_)){
  for(j in 1:length(q_)){
    for(k in 1:length(sigma2_)){
      for (l in 1:length(rho_)) {
        for (m in 1:length(mis_)) {
          case2=rbind(case2, cbind(n_[i], q_[j], sigma2_[k], rho_[l], mis_[m]))
        }
      }
    }
  }
}

case2

case2_ris=list(16)
cor_str=list(16)


set.seed(125)
cor_un=cor(matrix(runif(5*5),ncol=5))


cor_un10=cor(matrix(runif(20*20),ncol=10))
cor_cluster=matrix(c(1,-0.2,0,0,0,
                     -0.2,1,0,0,0,
                     0,0,1,0.7,0.7,
                     0,0,0.7,1,0.7,
                     0,0,0.7,0.7,1), ncol = 5)
cor_cluster
cor_cluster10=matrix(c(0.6,0.5, 0.5,0.6), ncol=2)%x%cor_cluster
rho=0.9
cor_ar=matrix( c(1, rho, rho^2, rho^3, rho^4, 
                 rho, 1, rho,rho^2,rho^3, 
                 rho^2,rho, 1, rho,rho^2, 
                 rho^3,rho^2,rho,1, rho,
                 rho^4, rho^3, rho^2,rho,1), ncol=5)
cor_ar10=matrix( c(1, rho, rho^2, rho^3, rho^4,rho^5,rho^6, rho^7,  rho^8, rho^9 ,
                   rho, 1, rho,rho^2,rho^3, rho^4,rho^5,rho^6, rho^7,  rho^8, 
                   rho^2,rho, 1, rho,rho^2, rho^3, rho^4,rho^5,rho^6, rho^7,
                   rho^3,rho^2,rho,1, rho,rho^2, rho^3, rho^4,rho^5,rho^6, 
                   rho^4, rho^3, rho^2,rho,1,rho, rho^2, rho^3, rho^4,rho^5,
                   rho^5, rho^4, rho^3, rho^2,rho,1, rho,rho^2, rho^3, rho^4,
                   rho^6, rho^5, rho^4, rho^3, rho^2,rho,1, rho,rho^2, rho^3,  
                   rho^7,rho^6, rho^5, rho^4, rho^3, rho^2,rho,1,rho, rho^2,   
                   rho^8,rho^7,rho^6, rho^5, rho^4, rho^3, rho^2,rho,1, rho, 
                   rho^9,rho^8,rho^7,rho^6, rho^5, rho^4, rho^3, rho^2,rho,1
), ncol=10)
cor_str[[1]]= cor_str[[9]]=diag(5) #exch
cor_str[[2]]=cor_str[[10]]=cor_un
cor_str[[3]]=cor_str[[11]]=cor_cluster
cor_str[[4]]=cor_str[[12]]=cor_ar


cor_str[[5]]=cor_str[[13]]=diag(10)
cor_str[[6]]=cor_str[[14]]=cor_un10
cor_str[[7]]=cor_str[[15]]=cor_cluster10
cor_str[[8]]=cor_str[[16]]=cor_ar10
R=10000
for(i in 1:16){
  
  
  rho=case2[i,4]
  s2=case2[i,3]
  beta=c(2, 0.3,-1, 3,-0.5)
  q=(case2[i,2])
  n=(case2[i,1])
  pars=c(beta,s2,rho)
  set=c(pars, q,n)
  
  mis=as.logical(case2[i,5])
  cor0=cor_str[[i]]
  p=7
  # Xmatrix  
  set.seed(123)
  x=cbind(c(rep(1,n)), runif(n, -10,10), rexp(n, 0.5) , rbinom(n,1, 0.5),rbinom(n,1, 0.2))
  x=matrix(apply(x, 1 ,function(r) rep(t(r),(q))),ncol=(p-2), byrow=T)
  x=x[,1:(p-2)]    
  est=rep(NA,R);risultati_all=list(NULL);ris_out=NULL
  
  
  if(0==0){
    seed0=(987654321)
    set.seed(seed0)
    cl <- makeCluster(getOption("cl.cores", 8))
    clusterExport(cl, c("breqn","all_sim","pars","x","p","q","n","rmnorm","gee", "cor0", "mis"))
    risultati_all =clusterApply(cl,est,fun =   function (z)  
      all_sim(cor = cor0,p=p,q=q,n=n,x = x, pars=pars, mis= mis))
  }
  
  case2_ris[[i]]=out_all(risultati_all,p = 7)
  cat(i)}

#saveRDS(case2_ris, "fixed10000_brenq.RDS")#q=5  #vector length
#saveRDS(case3_ris, "variable10000_brenq_1_10.RDS")#q=5  #vector length

#case2_ris=readRDS(file.choose())#"fixed1000.RDS"
#case3_ris=readRDS(file.choose())#"fixed1000.RDS"
 




table_summary=function(case2_ris){
  pars=6:7
  Estimator=rep(rep(c("MLE", "BC", "JEF",  "BR" ,"MBR"), each=length(pars)),4)
  coefficient=c("sigma2","rho")
  w=c(1,5,9,13)
  taball=list()
  for(i in 1:length(w)){
    tab=(case2_ris[[w[i]]][,c(1,2,3,4,5)])[c(pars+7,pars+14, pars+21, pars+28),]
    tab=tab[,c(5,3,4,2,1)]
    seq_s2=seq(from=1,to=8,by = 2)
    seq_rho=seq(from=1,to=8,by = 2)+1
    tab=rbind(tab[seq_s2,], tab[seq_rho,])
    tab=t(tab)
    tab=rbind(tab[,1:4], tab[,5:8])
    taball[[i]]=tab}
  taball
  
  rbind(cbind(taball[[1]],taball[[2]]), cbind(taball[[3]],taball[[4]]))
}
#xtable::xtable(table_summary(case2_ris))
#xtable::xtable(table_summary(case3_ris))



##################################################################
#some graphical representations  
plot_summary=function(ris,p=7, stat="bias", pars=NULL){
  
  Estimator=rep(rep(c("MLE", "BC", "JEF",  "BR" ,"MBR", "GEE", "ROB"), each=length(pars)),4)
  coefficient=c("beta 1","beta 2","beta 3","beta 4","beta 5")
  pars=1:5
  
  select_rows=function(h, t1 ,t2 ,t3 ,t4 ){
    t1h=t1[h:(h+p-1),]; t2h=t2[h:(h+p-1),]; t3h=t3[h:(h+p-1),]; t4h=t4[h:(h+p-1),]
    t1h=t1h[pars,];   t2h=t2h[pars,];     t3h=t3h[pars,];    t4h=t4h[pars,]
    t1h=as.vector(t1h);    t2h=as.vector(t2h);    t3h=as.vector(t3h);    t4h=as.vector(t4h)
    t1h=cbind(t1h, c(rep("exch", length(t1h)))) 
    t2h=cbind(t2h, c(rep("unstr", length(t2h))))
    t3h=cbind(t3h, c(rep("clus", length(t3h))))
    t4h=cbind(t4h, c(rep("AR1", length(t4h))))
    colnames(t1h)[ncol(t1h)]=colnames(t3h)[ncol(t3h)]=
      colnames(t4h)[ncol(t4h)]=colnames(t2h)[ncol(t2h)]="corstr"
    
    tab=rbind(t1h,t2h,t3h, t4h)
    tab
  }
  
  data_scenario=function( t1,t2,t3,t4){
    #WALD
    h0=1+3*p
    tab=select_rows(h=h0 ,t1,t2,t3,t4)
    dataw=data.frame(Estimator, tab, (coefficient ))
    dataw[,3]=as.factor(dataw[,3])
    
    colnames(dataw)[2]="val"
    dataw$val=as.numeric(dataw$val)
    dataw$Estimator <- factor(dataw$Estimator, levels = c( "MBR", "JEF", "BR", "BC", "MLE", "GEE", "ROB"))
    dataw$corstr <- factor(dataw$corstr, levels = c("exch","AR1", "clus","unstr"))
    dataw[,ncol(dataw)+1]="WALD"
    colnames(dataw[ncol(dataw)])="stat"
    
    #RB
    h0=1+2*p
    tab=select_rows(h=h0 ,t1,t2,t3,t4)
    datarb=data.frame(Estimator, tab, (coefficient ))
    datarb[,3]=as.factor(datarb[,3])
    
    colnames(datarb)[2]="val"
    datarb$val=as.numeric(datarb$val)
    datarb$Estimator <- factor(datarb$Estimator, levels = c( "MBR", "JEF", "BR", "BC", "MLE", "GEE", "ROB"))
    datarb$corstr <- factor(datarb$corstr, levels = c("exch","AR1", "clus","unstr"))
    datarb[,ncol(datarb)+1]="RB"
    colnames(datarb[ncol(datarb)])="stat"
    
    #PU
    h0=1+p
    tab=select_rows(h=h0, t1,t2,t3,t4)
    datapu=data.frame(Estimator, tab,(coefficient  ))
    datapu[,3]=as.factor(datapu[,3])
    
    colnames(datapu)[2]="val"
    datapu$val=as.numeric(datapu$val)
    datapu$Estimator <- factor(datapu$Estimator, levels = c( "MBR", "JEF", "BR", "BC", "MLE", "GEE", "ROB"))
    datapu$corstr <- factor(datapu$corstr, levels = c("exch","AR1", "clus","unstr"))
    datapu[,ncol(datapu)+1]="PU"
    colnames(datapu[ncol(datapu)])="stat"
    
    data=rbind(dataw, datarb, datapu)
    #leave=(data$corstr!="exch"&data$coefficient=="rho")
    #data=data[leave==FALSE,]
    colnames(data)[ncol(data)]="stat"
    data$stat <- factor(data$stat, levels = c("PU","RB", "WALD" ))
    data
  }
  
  t1=ris[[1]];  t2=ris[[2]]; t3=ris[[3]]; t4=ris[[4]]
  data1=data_scenario(t1,t2,t3,t4)
  data1[,ncol(data1)+1]=20
  data1[,ncol(data1)+1]="q=5"
  colnames(data1)[ncol(data1)-1]="n"
  colnames(data1)[ncol(data1)]="q"
  colnames(data1)[ncol(data1)-3]="coefficient"
  
  t1=ris[[5]];  t2=ris[[6]]; t3=ris[[7]]; t4=ris[[8]]
  data2=data_scenario(t1,t2,t3,t4)
  data2[,ncol(data2)+1]=20
  data2[,ncol(data2)+1]="q=10"
  colnames(data2)[ncol(data2)-1]="n"
  colnames(data2)[ncol(data2)]="q"
  colnames(data2)[ncol(data2)-3]="coefficient"
  
  
  t1=ris[[9]];  t2=ris[[10]]; t3=ris[[11]]; t4=ris[[12]]
  data3=data_scenario(t1,t2,t3,t4)
  data3[,ncol(data3)+1]=50
  data3[,ncol(data3)+1]="q=5"
  colnames(data3)[ncol(data3)-1]="n"
  colnames(data3)[ncol(data3)]="q"
  colnames(data3)[ncol(data3)-3]="coefficient"
  
  
  t1=ris[[13]];  t2=ris[[14]]; t3=ris[[15]]; t4=ris[[16]]
  data4=data_scenario(t1,t2,t3,t4)
  data4[,ncol(data4)+1]=50
  data4[,ncol(data4)+1]="q=10"
  
  colnames(data4)[ncol(data4)-1]="n"
  colnames(data4)[ncol(data4)]="q"
  colnames(data4)[ncol(data4)-3]="coefficient"
  
  plotscenario=function(data, nchoose){
    data=data[data$n==nchoose,]
    # lab=paste(paste("n=",unique(data[,ncol(data)-1]), collapse = ""), paste("q=",unique(data[,ncol(data)]), collapse = ""), collapse = ", ")
    hline_dat=data.frame(stat=c("PU","RB", "WALD"), threshold=c(50,0,95))
    hline_dat$stat=factor(hline_dat$stat,levels = c("PU","RB", "WALD"))
    data$ref=0
    data$ref[ data$stat=="WALD"]=95
    data$ref[ data$stat=="PU"]  =50
    data$q=factor(data$q, levels = c("q=5", "q=10"))
    #data$q2=sapply(qq, function (x) parse("q=", x))
    # print(data)
    p=ggplot(data, aes(x=Estimator, y=val  )) + 
      theme_bw() +
      scale_shape_manual(values = c(1:5))+
      theme(axis.text.x = element_text(angle = 90,size = 7, vjust = 0.05))+
      labs(y = paste("n =", nchoose),x = " ", labeller = label_parsed) + 
      
      geom_hline(  aes(yintercept=ref), col = "grey" )+
      scale_y_continuous( )+
      geom_point(aes(shape=coefficient, col=coefficient), size=0.9,   alpha=1,
                 position = position_jitterdodge(dodge.width = 0.0,jitter.width = 0,jitter.height = 0))+
      #geom_hline(data,ref)+
      theme(legend.position = "top")+
      #facet_grid(cols = vars(corstr), rows = vars(stat), as.table = T, scales = "free")
      facet_nested( (q) +stat~ corstr , labeller = label_value , scales = "free")
    p
  }
  #p1=plotscenario(data1) 
  #  p2=plotscenario(data2)
  #p3=plotscenario(data3)
  #p4=plotscenario(data4)
  
  dataall=rbind(data1,data2, data3,data4)
  
  library(ggh4x)
  p1 <- plotscenario(dataall, 20)
  p2 <- plotscenario(dataall, 50)
  #p + facet_nested(vs + cyl ~ am + gear, labeller = label_both)
  # install.packages("ggh4x")
  library(gridExtra)
  return(list(p1, p2))
  
  #plot_grid(p1, p2, p3,p4, labels=c("a", "b", "c","d"), ncol = 2, nrow = 2, align = "h")
  #grid.arrange(p1,p2,p3,p4, ncol=2)
}

par(mfrow=c(2,1))
plot_summary(case3_ristot,p = 7, pars=1:5)

#setEPS()
#postscript(file="p3.pdf", width=4, height=5)


pdf(file="p3.pdf", width=6, height=6 )
plot_summary(case3ristot, p=7, pars = 1:5)[[1]]
dev.off()
pdf(file="p4.pdf", width=6, height=6 )
plot_summary(case3ristot,p = 7, pars=1:5)[[2]]
dev.off()
#########################################################

#pdf(file="p3.pdf", width=6, height=6 )
plot_summary(case3_ris,p = 7, pars=1:5)[[1]]
#dev.off()
#pdf(file="p4.pdf", width=6, height=6 )
plot_summary(case3_ris,p = 7, pars=1:5)[[2]]
d#ev.off()




########################################
R=10000
#case 3 covariate time non-constant
#set1

n_=c(20,50)
q_=c( 5, 10)
rho_=c(0.9)
sigma2_=c(5)
mis_=c(F,T,T,T)

case3=NULL
for(i in 1:length(n_)){
  for(j in 1:length(q_)){
    for(k in 1:length(sigma2_)){
      for (l in 1:length(rho_)) {
        for (m in 1:length(mis_)) {
          case3=rbind(case3, cbind(n_[i], q_[j], sigma2_[k], rho_[l], mis_[m]))
        }
      }
    }
  }
}
case3
case3_ris=list()
R=10000
for(i in 11:16){
  
  rho=case3[i,4]
  s2=case3[i,3]
  beta=c(2, 0.3,-1, 3,-0.5)
  q=(case3[i,2])
  n=(case3[i,1])
  pars=c(beta,s2,rho)
  set=c(pars, q,n)
  
  mis=as.logical(case3[i,5])
  cor0=cor_str[[i]]
  p=7
  # Xmatrix  
  set.seed(123)
  x=cbind(c(rep(1,n)),  rbinom(n,1, 0.5),rbinom(n,1, 0.2))
  x=matrix(apply(x, 1 ,function(r) rep(t(r),(q))),ncol=(ncol(x)), byrow=T)
  
  x=cbind(x,  rnorm(n*q,0,10), rep(1:q,n) )
  est=rep(NA,R);risultati_all=list(NULL);ris_out=NULL   
  
  if(0==0){
    seed0=(987654321)
    set.seed(seed0)
    cl <- makeCluster(getOption("cl.cores", 8))
    clusterExport(cl, c("breqn","all_sim","pars","x","p","q","n","rmnorm","gee", "cor0", "mis"))
    risultati_all =clusterApply(cl,est,fun =   function (z)  
      all_sim(cor = cor0,p=p,q=q,n=n,x = x, pars=pars, mis= mis))
  }
  
  case3_ris[[i]]=out_all(risultati_all,p = 7)
  cat(i)}


#saveRDS(case3_ris, "variable10000.RDS")#q=5  #vector length

#case3_ris=readRDS(file.choose())#variable10000.RDS
#case3_ris2=readRDS(file.choose())#variable10000.RDS
 
plot_summary(case3_ris,p = 7, pars=1:5)


#setEPS()
#postscript(file="epileptic_plot_sim.eps", width=4, height=4.5)



library(ggpubr)
library(gridExtra)
plot_grid(p1, p2, p3,p4, labels=c("", "", "",""), ncol = 1, nrow = 4, align = "v")
grid.arrange(p1,p2,p3)

 


