# LIBRARIES ----------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(viridis)
library(xtable)


#######################################################################
##### pnormLS
#######################################################################
# Normal Distribution Function (LS)
pnormLS <- function(x, mu, sigma, ...) {pnorm(x, mu, sigma)}

#######################################################################
##### rnormLS
#######################################################################
# Random Value Generation for the Normal Distribution (LS)
rnormLS <- function(n, mu, sigma, ...) {rnorm(n, mu, sigma)}

#######################################################################
##### pexpLS
#######################################################################
# Shifted Exponential Distribution Function (LS)
pexpLS <- function(x, mu, sigma, ...) {pexp((x - (mu - sigma)) / sigma, rate = 1) * (x >= (mu - sigma))}

#######################################################################
##### rexpLS
#######################################################################
# Random Value Generation for the Shifted Exponential Distribution (LS)
rexpLS <- function(n, mu, sigma, ...) {(mu - sigma) + sigma * rexp(n, rate = 1)}

#######################################################################
##### punifLS
#######################################################################
# Location-Scale Uniform Distribution Function Derived from Standard Form
punifLS <- function(x, mu, sigma, ...){punif(x, mu - sigma*sqrt(3), mu + sigma*sqrt(3))}

#######################################################################
##### runifLS
#######################################################################
# Random Variate Generation for Location-Scale Uniform Distribution from Standard Form
runifLS <- function(n, mu, sigma, ...){mu + sigma*runif(n, -sqrt(3), sqrt(3))}

#######################################################################
##### plogisLS
#######################################################################
# Logistic Distribution Function (LS)
plogisLS <- function(x, mu, sigma, ...){plogis((x-mu)/sigma, 0, sqrt(3)/pi)}

#######################################################################
##### rlogisLS
#######################################################################
# Random Value Generation for the Logistic Distribution (LS)
rlogisLS <- function(n, mu, sigma, ...){(mu)+sigma*rlogis(n, 0, sqrt(3)/pi)}

#######################################################################
##### mle.logis
#######################################################################
# Retrieves location and scale MLE for the logistic distribution.
mle.logis <- function(data){
  log.lik <- function(par, data){
    mu <- par[1]; sigma <- par[2]
    if(sigma>0) return(sum(dlogis(data, par[1], par[2]* sqrt(3)/pi, log=TRUE)))
    else return(-Inf)  
  }
  mle.logis <- optim(par = c(mean(data), sd(data)), fn = log.lik, data = data, control = list(fnscale = -1),
                     method = "Nelder-Mead")
  return(mle.logis$par)
}


#######################################################################
##### plaplaceLS
#######################################################################
# Location-Scale Laplace Distribution Function Derived from Standard Form
plaplaceLS <- function(q, mu, sigma, ...) {
  b <- sigma / sqrt(2)
  p <- ifelse(q < mu,
              0.5 * exp((q - mu) / b),
              1 - 0.5 * exp(-(q - mu) / b))
  return(p)
}

#######################################################################
##### rlaplaceLS
#######################################################################
# Random Value Generation for the Laplace Distribution (LS)
rlaplaceLS <- function(n, mu, sigma, ...) {
  b <- sigma / sqrt(2)
  u <- runif(n, -0.5, 0.5)
  mu - b * sign(u) * log(1 - 2 * abs(u))
}

#######################################################################
##### mle.laplace
#######################################################################
# Retrieves location and scale MLE for the Laplace distribution.
mle.laplace<- function(x, start = NULL, method = "Nelder-Mead") {
  n <- length(x)
  
  # Log-likelihood
  loglik <- function(par, x) {
    mu <- par[1]
    sigma <- par[2]
    
    if (sigma <= 0 || !is.finite(mu) || !is.finite(sigma)) {
      return(1e10)  # penalization
    }
    
    sqrt2 <- sqrt(2)
    
    # Density
    ll <- ifelse(x < mu,
                 log(sqrt2 / (2 * sigma)) - sqrt2 * (-x + mu) / sigma,
                 log(sqrt2 / (2 * sigma)) - sqrt2 * (x - mu) / sigma)
    
    return(-sum(ll))  # Minimise negative log-likelihood
  }
  
  # Initial values
  if (is.null(start)) {
    mu0 <- median(x)
    sigma0 <- mean(abs(x - mu0))
    start <- c(mu0, sigma0)
  }
  
  # Optimization
  fit <- optim(
    par = start,
    fn = loglik,
    x = x,
    method = method,
    control = list(fnscale = 1)
  )
  
  # Results
  list(
    start = start,
    mle = setNames(fit$par, c("mu", "sigma")),
    logLik = -fit$value,
    convergence = fit$convergence,
    message = fit$message
  )
}

#######################################################################
##### distributions.list
#######################################################################
# List of available distributions with embedded lists of the corresponding MLE 
# expressions as functions.

distributions.list <- list(
  normLS = list(
    mu.mle = function(x, ...) mean(x), 
    sigma.mle = function(x, ...) sqrt(mean((x-mean(x))^2))
  ),
  expLS = list(
    mu.mle = function(x, ...) mean(x),
    sigma.mle = function(x, ...) mean(x)-min(x) 
  ),
  logisLS = list(
    mu.mle = function(x, ...) mle.logis(x)[1],
    sigma.mle = function (x, ...) mle.logis(x)[2]
  ),
  unifLS = list(
    mu.mle = function(x, ...) (1/2)*(min(x) + max(x)),
    sigma.mle = function(x, ...) sqrt((1/12)*(max(x) - min(x))^2)
  ),
  laplaceLS = list(
    mu.mle = function(x, ...) mle.laplace(x)$mle[1],
    sigma.mle = function(x, ...) mle.laplace(x)$mle[2]
  )
)