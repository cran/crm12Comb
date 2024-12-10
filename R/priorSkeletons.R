# function priorSkeletons() is modified based on the function getprior() from the package dfcrm and
# reference paper:
# Lee SM, Cheung YK. Model calibration in the continual reassessment method. Clinical Trials 2009; 6(3): 227–238.
# that can further apply to more models and prior distributions

# input values:
# (1) beta_mean for hyperbolic tangent model with exponential prior
# (2) beta_mean for empiric model/one parameter logistic model with normal prior
# (3) beta_shape, beta_inverse_scale for empiric model/one parameter logistic model with gamma prior
# (4) alpha_mean, beta_mean for two parameter logistic model with normal prior for both alpha and beta
# (5) alpha_shape, alpha_inverse_scale, beta_shape, beta_inverse_scale for two parameter logistic model with normal prior for both alpha and beta


#' Title: Get the skeletons of toxicity and efficacy
#'
#' @param updelta A number for half-width of indifference intervals.
#' @param target A number for target DLT rate.
#' @param npos A number for the prior guess MTD position.
#' @param ndose A number for the number of testing doses.
#' @param model A string for indicating which model to use.
#' @param prior A string for indicating which prior distribution to use.
#' @param alpha_mean A number.
#' @param beta_mean A number.
#' @param a0 A number.
#' @param alpha_shape A number.
#' @param alpha_inverse_scale A number.
#' @param beta_shape A number.
#' @param beta_inverse_scale A number.
#'
#' @return val -> A numeric vector of skeleton
#' @export
#'
#' @references Modified based on the function getprior() from the package dfcrm
#'             and reference paper: Lee SM, Cheung YK. Model calibration in the continual reassessment method. Clinical Trials 2009; 6(3): 227–238.


priorSkeletons <- function (updelta, target, npos, ndose, model = "empiric", prior = "normal",
                            alpha_mean=NULL, beta_mean=0, a0 = 3,
                            alpha_shape=NULL, alpha_inverse_scale=NULL, beta_shape=NULL, beta_inverse_scale=NULL) {
  dosescaled <- rep(NA, ndose)

  if (model == "tanh") {
    if (prior == "exponential"){
      dosescaled[npos] <- atanh(2*target^(1/beta_mean)-1)
    }
    for (i in npos:2) {
      if (npos > 1) {
        dosescaled[i - 1] <- -1/2 * log((1+exp(-2*dosescaled[i]))^(log(target - updelta)/log(target + updelta))-1)
      }
    }
    for (i in npos:(ndose - 1)) {
      if (npos < ndose) {
        dosescaled[i + 1] <- -1/2 * log((1+exp(-2*dosescaled[i]))^(log(target + updelta)/log(target - updelta))-1)
      }
    }
    if (prior == "exponential"){
      val <- ((tanh(dosescaled)+1)/2)^beta_mean
    }
  }

  else if (model == "empiric") {
    if (prior == "normal"){
      dosescaled[npos] <- target^exp(-beta_mean)
    } else if (prior == "gamma"){
      dosescaled[npos] <- target^(beta_inverse_scale/beta_shape)
    }
    for (i in npos:2) {
      if (npos > 1) {
        dosescaled[i - 1] <- exp((log(target - updelta) * log(dosescaled[i])) / log(target + updelta))
      }
    }
    for (i in npos:(ndose - 1)) {
      if (npos < ndose) {
        dosescaled[i + 1] <- exp((log(target + updelta) * log(dosescaled[i])) / log(target - updelta))
      }
    }
    if (prior == "normal"){
      val <- dosescaled^exp(beta_mean)
    } else if (prior == "gamma"){
      val <- dosescaled^(beta_shape/beta_inverse_scale)
    }
  }

  else if (model == "logistic") {
    if (prior == "normal"){
      dosescaled[npos] <- (log(target/(1 - target)) - a0) / exp(beta_mean)
    } else if (prior == "gamma"){
      dosescaled[npos] <- (log(target/(1-target)) - a0)/(beta_shape/beta_inverse_scale)
    }
    for (i in npos:2) {
      if (npos > 1) {
        dosescaled[i - 1] <- (log((target - updelta)/(1 - target + updelta)) - a0) * dosescaled[i] / (log((target + updelta)/(1 - target - updelta)) - a0)
      }
    }
    for (i in npos:(ndose - 1)) {
      if (npos < ndose) {
        dosescaled[i + 1] <- (log((target + updelta)/(1 - target - updelta)) - a0) * dosescaled[i] / (log((target - updelta)/(1 - target + updelta)) - a0)
      }
    }
    if (prior == "normal"){
      val <- val <- 1/(1 + exp(-a0-exp(beta_mean)*dosescaled))
    } else if (prior == "gamma"){
      val <- val <- 1/(1 + exp(-a0-(beta_shape/beta_inverse_scale)*dosescaled))
    }
  }

  else if (model == "logistic2") {
    if (prior == "normal"){
      dosescaled[npos] <- (log(target/(1 - target)) - alpha_mean) / exp(beta_mean)
      for (i in npos:2) {
        if (npos > 1) {
          dosescaled[i - 1] <- (log((target - updelta)/(1 - target + updelta)) - alpha_mean) * dosescaled[i] / (log((target + updelta)/(1 - target - updelta)) - alpha_mean)
        }
      }
      for (i in npos:(ndose - 1)) {
        if (npos < ndose) {
          dosescaled[i + 1] <- (log((target + updelta)/(1 - target - updelta)) - alpha_mean) * dosescaled[i] / (log((target - updelta)/(1 - target + updelta)) - alpha_mean)
        }
      }
      val <- 1/(1+exp(-alpha_mean - exp(beta_mean)*dosescaled))

    } else if (prior == "gamma"){
      dosescaled[npos] <- (log(target/(1-target)) - (alpha_shape/alpha_inverse_scale))/(beta_shape/beta_inverse_scale)
      for (i in npos:2) {
        if (npos > 1) {
          dosescaled[i - 1] <- (log((target - updelta)/(1 - target + updelta)) - (alpha_shape/alpha_inverse_scale)) * dosescaled[i] / (log((target + updelta)/(1 - target - updelta)) - (alpha_shape/alpha_inverse_scale))
        }
      }
      for (i in npos:(ndose - 1)) {
        if (npos < ndose) {
          dosescaled[i + 1] <- (log((target + updelta)/(1 - target - updelta)) - (alpha_shape/alpha_inverse_scale)) * dosescaled[i] / (log((target - updelta)/(1 - target + updelta)) - (alpha_shape/alpha_inverse_scale))
        }
      }
      val <- 1/(1+exp(-alpha_shape/alpha_inverse_scale - (beta_shape/beta_inverse_scale)*dosescaled))
    }
  }

  return(val)
}
