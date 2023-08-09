
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



