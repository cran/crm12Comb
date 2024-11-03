#' Title: Bayesian likelihood inference
#' @references: R package dfcrm

#' Hyperbolic tangent model with exponential prior
#' @param beta parameter.
#' @param beta_mean A number.
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return l -> likelihood function
#' @export
tanh_ExpPriorLikelihood <- function(beta, beta_mean, x, y){
  l <- beta_mean*exp(-beta_mean*beta)
  for (i in 1:length(x)){
    l <- l * ((((tanh(x[i])+1)/2)^beta)^y[i])*((1-(((tanh(x[i])+1)/2)^beta))^(1-y[i]))
  }
  return(l)
}

tanh_ExpPriorLikelihood_est <- function(beta, beta_mean, x, y){
  l <- beta * beta_mean*exp(-beta_mean*beta)
  for (i in 1:length(x)){
    l <- l * ((((tanh(x[i])+1)/2)^beta)^y[i])*((1-(((tanh(x[i])+1)/2)^beta))^(1-y[i]))
  }
  return(l)
}


#' empiric model with normal prior
#' @param beta parameter.
#' @param beta_mean A number.
#' @param beta_sd A number.
#' @param x A numeric vector. 
#' @param y A numeric vector. 
#'
#' @return l -> likelihood function
#' @export
empiric_NormalPriorLikelihood <- function(beta, beta_mean, beta_sd, x, y){
  l <- (1/sqrt(2*pi*beta_sd^2))*exp(-(beta-beta_mean)^2/(2*beta_sd^2))
  for (i in 1:length(x)){
    l <- l * ((x[i]^(exp(beta)))^y[i])*((1-x[i]^(exp(beta)))^(1-y[i]))
  }
  return(l)
}

empiric_NormalPriorLikelihood_est <- function(beta, beta_mean, beta_sd, x, y){
  l <- beta * (1/sqrt(2*pi*beta_sd^2))*exp(-(beta-beta_mean)^2/(2*beta_sd^2))
  for (i in 1:length(x)){
    l <- l * ((x[i]^(exp(beta)))^y[i])*((1-x[i]^(exp(beta)))^(1-y[i]))
  }
  return(l)
}

#' empiric model with gamma prior
#' @param beta parameter.
#' @param beta_shape A number.
#' @param beta_inverse_scale A number.
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return l -> likelihood function
#' @export
empiric_GammaPriorLikelihood <- function(beta, beta_shape, beta_inverse_scale, x, y){
  l <- ((beta^(beta_shape-1))*exp(-beta_inverse_scale*beta)*(beta_inverse_scale^beta_shape))/gamma(beta_shape)
  for (i in 1:length(x)){
    l <- l * ((x[i]^(beta))^y[i])*((1-x[i]^(beta))^(1-y[i]))
  }
  return(l)
}

empiric_GammaPriorLikelihood_est <- function(beta, beta_shape, beta_inverse_scale, x, y){
  l <- beta * ((beta^(beta_shape-1))*exp(-beta_inverse_scale*beta)*(beta_inverse_scale^beta_shape))/gamma(beta_shape)
  for (i in 1:length(x)){
    l <- l * ((x[i]^(beta))^y[i])*((1-x[i]^(beta))^(1-y[i]))
  }
  return(l)
}


#' one-parameter logistic model with normal prior
#' @param alpha1 parameter.
#' @param alpha1_mean A number.
#' @param alpha1_sd A number.
#' @param intcpt A number.
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return l -> likelihood function
#' @export
logisticOnePara_NormalPriorLikelihood <- function(alpha1, alpha1_mean, alpha1_sd, intcpt, x, y){
  l <- (1/sqrt(2*pi*alpha1_sd^2))*exp(-(alpha1-alpha1_mean)^2/(2*alpha1_sd^2))
  for (i in 1:length(x)){
    l <- l * ((1/(1 + exp(-intcpt-exp(alpha1)*x[i])))^y[i])*((1-1/(1 + exp(-intcpt-exp(alpha1)*x[i])))^(1-y[i]))
  }
  return(l)
}

logisticOnePara_NormalPriorLikelihood_est <- function(alpha1, alpha1_mean, alpha1_sd, intcpt, x, y){
  l <- alpha1 * (1/sqrt(2*pi*alpha1_sd^2))*exp(-(alpha1-alpha1_mean)^2/(2*alpha1_sd^2))
  for (i in 1:length(x)){
    l <- l * ((1/(1 + exp(-intcpt-exp(alpha1)*x[i])))^y[i])*((1-1/(1 + exp(-intcpt-exp(alpha1)*x[i])))^(1-y[i]))
  }
  return(l)
}

#' one-parameter logistic model with gamma prior
#' @param alpha1 parameter.
#' @param alpha1_shape A number.
#' @param alpha1_inverse_scale A number.
#' @param intcpt A number.
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return l -> likelihood function
#' @export
logisticOnePara_GammaPriorLikelihood <- function(alpha1, alpha1_shape, alpha1_inverse_scale, intcpt, x, y){
  l <- ((alpha1^(alpha1_shape-1))*exp(-alpha1_inverse_scale*alpha1)*(alpha1_inverse_scale^alpha1_shape))/gamma(alpha1_shape)
  for (i in 1:length(x)){
    l <- l * ((1/(1 + exp(-intcpt-alpha1*x[i])))^y[i])*((1-1/(1 + exp(-intcpt-alpha1*x[i])))^(1-y[i]))
  }
  return(l)
}

logisticOnePara_GammaPriorLikelihood_est <- function(alpha1, alpha1_shape, alpha1_inverse_scale, intcpt, x, y){
  l <- alpha1 * ((alpha1^(alpha1_shape-1))*exp(-alpha1_inverse_scale*alpha1)*(alpha1_inverse_scale^alpha1_shape))/gamma(alpha1_shape)
  for (i in 1:length(x)){
    l <- l * ((1/(1 + exp(-intcpt-alpha1*x[i])))^y[i])*((1-1/(1 + exp(-intcpt-alpha1*x[i])))^(1-y[i]))
  }
  return(l)
}


#' two-parameter logistic model with normal prior
#' @param alpha0 parameter.
#' @param alpha1 parameter.
#' @param alpha0_mean A number.
#' @param alpha0_sd A number.
#' @param alpha1_mean A number.
#' @param alpha1_sd A number.
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return l -> likelihood function
#' @export
logisticTwoPara_NormalPriorLikelihood <- function(alpha0, alpha1, alpha0_mean, alpha0_sd, alpha1_mean, alpha1_sd, x, y){
  l <- (1/sqrt(2*pi*alpha0_sd^2))*exp(-(alpha0-alpha0_mean)^2/(2*alpha0_sd^2)) * (1/sqrt(2*pi*alpha1_sd^2))*exp(-(alpha1-alpha1_mean)^2/(2*alpha1_sd^2))
  for (i in 1:length(x)){
    l <- l * ((1/(1 + exp(-alpha0-exp(alpha1)*x[i])))^y[i])*((1-1/(1 + exp(-alpha0-exp(alpha1)*x[i])))^(1-y[i]))
  }
  return(l)
}

#' two-parameter logistic model with gamma prior
#' @param alpha0 parameter.
#' @param alpha1 parameter.
#' @param alpha0_shape A number.
#' @param alpha0_inverse_scale A number.
#' @param alpha1_shape A number.
#' @param alpha1_inverse_scale A number.
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @return l -> likelihood function
#' @export
logisticTwoPara_GammaPriorLikelihood <- function(alpha0, alpha1, alpha0_shape, alpha0_inverse_scale, alpha1_shape, alpha1_inverse_scale, x, y){
  l <- (((alpha0^(alpha0_shape-1))*exp(-alpha0_inverse_scale*alpha0)*(alpha0_inverse_scale^alpha0_shape))/gamma(alpha0_shape)) * (((alpha1^(alpha1_shape-1))*exp(-alpha1_inverse_scale*alpha1)*(alpha1_inverse_scale^alpha1_shape))/gamma(alpha1_shape))
  for (i in 1:length(x)){
    l <- l * ((1/(1 + exp(-alpha0-alpha1*x[i])))^y[i])*((1-1/(1 + exp(-alpha0-alpha1*x[i])))^(1-y[i]))
  }
  return(l)
}


