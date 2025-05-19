###### get.ad.statistic
########################################
# Computes the Anderson-Darling statistic (using estimators of choice: MLE or mean-sd) for a 
# given dataset (data) against a specified theoretical family distribution (dist).

# Output:
# - Returns the parameter based Anderson-Darling statistic as a numeric value.

get.ad.statistic <- function(data, dist = "norm",shape = NULL, shape1 = NULL, shape2 = NULL, df = NULL, method = NULL) {
  data <- sort(data[complete.cases(data)])
  n <- length(data)
  
  if (method == "MLE"){
    mu.MLE<<-as.double(get.mle.estimators(dist, data, shape = shape, shape1 = shape1, shape2 = shape2, df = df)$mu.mle)
    sigma.MLE<<-as.double(get.mle.estimators(dist, data, shape = shape, shape1 = shape1, shape2 = shape2, df = df)$sigma.mle)
  }
  else {
    mu.MLE<<-as.double(get.mu.sd.estimators(dist, data, shape = shape, shape1 = shape1, shape2 = shape2, df = df)$mu.mle) # .mle is just the name used in general.
    sigma.MLE<<-as.double(get.mu.sd.estimators(dist, data, shape = shape, shape1 = shape1, shape2 = shape2, df = df)$sigma.mle)  # .mle is just the name used in general.
  }
  data.sort<-sort(data)
  
  pdist <<- tryCatch({
    get(paste0("p", dist), mode = "function")
  }, error = function(e) stop("'dist' must be a valid function name"))
  
  # F0
  F0 <- pdist(data, mu.MLE, sigma.MLE, shape = shape, shape1 = shape1, shape2 = shape2, df = df)
  F0[F0 == 0] <- .Machine$double.eps
  F0[F0 == 1] <- 1 - .Machine$double.eps
  
  # GET AD STATISITIC
  i <- 1:n
  ad.stat <- -n - sum((2*i - 1) * (log(F0) + log(1 - rev(F0)))) / n
  
  return(ad.stat)
}


get.ad.statistic.add.distribution<- function(data, mu.mle, sigma.mle, pdist, shape = NULL, shape1 = NULL, shape2 = NULL, df = NULL, method = NULL){
  data <- sort(data[complete.cases(data)])
  n <- length(data)
  
  if (method == "MLE"){
    mu.MLE<<-as.double(mu.mle(data, shape = shape, shape1 = shape1, shape2 = shape2, df = df))
    sigma.MLE<<-as.double(sigma.mle(data, shape = shape, shape1 = shape1, shape2 = shape2, df = df))
  }
  else {
    get.mu.sd.estimators(dist,data)
    mu.MLE<<-as.double(get.mu.sd.estimators(dist, data, shape = shape, shape1 = shape1, shape2 = shape2, df = df)$mu.mle) # .mle is just the name used in general.
    sigma.MLE<<-as.double(get.mu.sd.estimators(dist, data, shape = shape, shape1 = shape1, shape2 = shape2, df = df)$sigma.mle)  # .mle is just the name used in general.
  }
  data.sort<-sort(data)
  F0 <- pdist(data, mu.MLE, sigma.MLE, shape = shape, shape1 = shape1, shape2 = shape2, df = df)
  F0[F0 == 0] <- .Machine$double.eps
  F0[F0 == 1] <- 1 - .Machine$double.eps
  
  i <- 1:n
  ad.stat <- -n - sum((2*i - 1) * (log(F0) + log(1 - rev(F0)))) / n
  
  return(ad.stat)
}

#########################################
###### ad.test.LoScF
########################################

# Performs the parameter based Anderson-Darling test to assess whether 
# a dataset (x) conforms to a specified location-scale family distribution (dist).

# Output:
# Returns a structured object of class "htest".
#########################################

ad.test.LoScF<-function(x,dist="norm", shape = NULL, shape1 = NULL, shape2 = NULL, df = NULL, alfa=0.05,simulate.distH0=TRUE, nr=10^4, add.distribution=FALSE, method = "MLE"){
  DNAME <- deparse(substitute(x))
  x <- sort(x[complete.cases(x)])
  n <- length(x)
  TIES <- FALSE
  if (length(unique(x)) < n) {
    warning("ties should not be present for the test")
    TIES <- TRUE
  }
  if(isFALSE(add.distribution)){
    
    rdist <- tryCatch({
      get(paste0("r", dist), mode = "function")
    }, error = function(e) stop("'dist' must be a valid function name")) # Made global so user can check the generator.
    
    if (isTRUE(simulate.distH0)){
      ad.statistic.dist <- replicate(nr,{
        iter.sample<-rdist(length(x),0,1, shape = shape, shape1 = shape1, shape2 = shape2, df = df)
        get.ad.statistic(iter.sample,dist=dist, shape = shape, shape1 = shape1, shape2 = shape2, df = df, method = method)
      })
      
    }
    else if (is.vector(simulate.distH0)) {
      ad.statistic.dist <- simulate.distH0
    } 
    else {
      stop("'simulate.distH0' must be TRUE or a numeric vector.")
    }
    
    ad.statistic<-get.ad.statistic(x,dist=dist, shape = shape, shape1 = shape1, shape2 = shape2, df = df, method = method)
  }
  else if(is.list(add.distribution) && !dist %in% names(distributions.list)){
    
    rdist_name <- paste0("r", dist)
    pdist_name <- paste0("p", dist)
    
    rdist <- if (rdist_name %in% names(add.distribution)) {
      if (is.function(add.distribution[[rdist_name]])) {
        add.distribution[[rdist_name]]
      } else {
        stop("The elements in the 'add.distribution' list must be functions.")
      }
    } else {
      stop(paste("The function", rdist_name, "is not found in the 'add.distribution' list.
                  Make sure the correct format is being used. The 'dist' parameter 
                  introduced must coincide with the function r<dist> given 
                  by the user in the add.distribution list."))
    }
    
    
    pdist <- if (pdist_name %in% names(add.distribution)) {
      if (is.function(add.distribution[[pdist_name]])) {
        add.distribution[[pdist_name]]
      } else {
        stop("The elements in the 'add.distribution' list must be functions.")
      }
    } else {
      stop(paste("The function", pdist_name, "is not found in the 'add.distribution' list.
                  Make sure the correct format is being used. The 'dist' parameter 
                  introduced must coincide with the function p<dist> given 
                  by the user in the add.distribution list."))
    }
    
    mu.mle <- if ("mu.mle" %in% names(add.distribution)) {
      if (is.function(add.distribution$mu.mle)) {
        add.distribution[["mu.mle"]]
      } else {
        stop("The elements in the 'add.distribution' list must be functions.")
      }
    } else {
      stop(paste("The function mu.MLE is not found in the 'add.distribution' list.
                  Make sure the correct format is being used. The 'mu.mle' parameter 
                  introduced by the user must be correctly named mu.mle in the 
                  'add.distribution' list."))
    }
    sigma.mle <- if ("sigma.mle" %in% names(add.distribution)) {
      if (is.function(add.distribution$sigma.mle)) {
        add.distribution[["sigma.mle"]]
      } else {
        stop("The elements in the 'add.distribution' list must be functions.")
      }
    } else {
      stop(paste("The function sigma.mle is not found in the 'add.distribution' list.
                  Make sure the correct format is being used. The 'sigma.mle' parameter 
                  introduced by the user must be correctly named sigma.mle in the 
                  'add.distribution' list."))
    }
    
    if (isTRUE(simulate.distH0)){
      ad.statistic.dist <- replicate(nr,{
        iter.sample<-rdist(length(x),0,1, shape, shape1, shape2, df)
        get.ad.statistic.add.distribution(iter.sample, mu.mle, sigma.mle, pdist, shape = shape, shape1 = shape1, shape2 = shape2, df = df, method = method)
      })
    }
    else if (is.vector(simulate.distH0)) {
      ad.statistic.dist <- simulate.distH0
    } 
    else {
      stop("'simulate.distH0' must be TRUE or a numeric vector.")
    }
    ad.statistic<-get.ad.statistic.add.distribution(x, mu.mle, sigma.mle, pdist, shape = shape, shape1 = shape1, shape2 = shape2, df = df, method = method)
  }
  else  {
    stop("'add.distribution' must be a list containing the Location-Scale
          MLE parameter general functions depending on the data.")
  }
  PVAL<-mean(ad.statistic.dist>=ad.statistic)
  
  structure(
    list(statistic = setNames(ad.statistic, "D"), p.value = PVAL, 
         method = "Anderson-Darling type based test using MLE or mean-sd estimators", data.name = DNAME, 
         alternative = "True distribution differs from the expected distribution", mu.mle = mu.MLE, sigma.mle = sigma.MLE
    ),
    class = "htest"
  )
  
}

#########################################
###### simulate.distH0.ad
########################################

# Simulates the parameter based Anderson-Darling test statistic distribution under H0 using estimated parameters of choice (MLE or mean-sd).


simulate.distH0.ad <- function(size = 100, mu = 0, sigma = 1, nr = 10^4, dist, shape = NULL, shape1 = NULL, shape2 = NULL, df = NULL, method = "MLE"){
  if(dist %in% names(distributions.list)){
    rdist_name <- paste0("r", dist)
    rdist <- tryCatch({
      get(paste0("r", dist), mode = "function")
    }, error = function(e) stop("'dist' must be a valid function name"))
    replicate(nr,{
      iter.sample<-eval(rdist(size,mu,sigma, shape = shape, shape1 = shape1, shape2 = shape2, df = df))
      get.ad.statistic(iter.sample,dist = dist, shape = shape, shape1 = shape1, shape2 = shape2, df = df, method = method)
    })
  }
  else return("Argument 'dist' not supported or not found in distributions.list")
}

#########################################
###### build.optional.args
########################################
build.optional.args <- function(params) {
  if (length(params) == 0) return("")
  paste(paste0(names(params), "=", params), collapse = ", ")
}


#########################################
###### potencia.ad
########################################
# Performs a power analysis for the parameter based AD test
potencia.ad <- function(n, str.distH0, str.dist.used,
                        params.H0 = list(), params.used = list(), method = "MLE") {
  set.seed(1)
  alpha <- 0.05
  
  args.H0 <- c(
    list(size = n, mu = 0, sigma = 1, nr = 10^4, dist = str.distH0, method = method),
    params.H0
  )
  distH0 <- do.call(simulate.distH0.ad, args.H0)
  
  rechazo <- replicate(10^4, {
    args.used <- c(
      list(n, 0, 1),
      params.used
    )
    rdist.name <- paste0("r", str.dist.used)
    rdist <- get(rdist.name, mode = "function")
    muestra <- do.call(rdist, args.used)
    
    pvalor <- ad.test.LoScF(muestra, str.distH0, simulate.distH0 = distH0, method = method)$p.value
    pvalor <= alpha
  })
  mean(rechazo)
}


# Distribution families definition
d1 <- c('laplace', 'unif', 'exp', 'logis', 'norm')
d2 <- c('laplace', 'unif', 'exp', 'logis', 'norm')

# Mapping -> real distribution
dist.map <- list(
  unif = 'unifLS',
  laplace = 'laplaceLS',
  exp = 'expLS',
  logis = 'logisLS',
  norm = 'normLS'
)

# Specific parameters
parametros <- list(
  laplace = list(),
  unif = list()
)