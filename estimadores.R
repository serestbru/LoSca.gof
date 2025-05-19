## CASUISTRY: MAXIMUM LIKELIHOOD ESTIMATION

get.mle.estimators <- function(dist_name, data, shape = NULL, shape1 = NULL, shape2 = NULL, df = NULL) {
  if (!dist_name %in% names(distributions.list)) {
    stop("Distribution cannot be recognised")
  }
  
  dist <- distributions.list[[dist_name]]
  
  mu.estimator <- dist$mu.mle(data, shape = shape, shape1 = shape1, shape2 = shape2, df = df)
  sigma.estimator <- if (!is.null(dist$sigma.mle)){
    dist$sigma.mle(data, shape = shape, shape1 = shape1, shape2 = shape2, df = df)
  }
  else{
    stop("Distribution cannot be recognised")
  }
  
  return(list(mu.mle = mu.estimator, sigma.mle = sigma.estimator))
}

## CASUISTRY: MEAN AND STANDARD DEVIATION.

get.mu.sd.estimators <- function(dist_name, data, shape = NULL, shape1 = NULL, shape2 = NULL, df = NULL) {
  if (!dist_name %in% names(distributions.list)) {
    stop("Distribution cannot be recognised")
  }
  
  dist <- distributions.list[[dist_name]]
  
  mu.estimator <- mean(data)
  sigma.estimator <- sd(data)
  
  return(list(mu.mle = mu.estimator, sigma.mle = sigma.estimator))
}