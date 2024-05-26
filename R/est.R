#' Title: patient allocations
#'
#' @param pE_est A numeric vector.
#' @param seed_m A number.
#'
#' @return d -> index of next recommended combined dose level
#' @export
#'
maximization_phase <- function(pE_est, seed_m=NULL){
  set.seed(seed_m)
  d <- ifelse(length(which(pE_est==max(pE_est)))==1,
             which(pE_est==max(pE_est)),
             sample(which(pE_est==max(pE_est)),1)) # allocation to dose combination for next patient
  return(d)
}

#' Title: patient allocations
#'
#' @param pE_est A numeric vector.
#' @param seed_r A number.
#'
#' @return d -> index of next recommended combined dose level
#' @export
#'
randomization_phase <- function(pE_est, seed_r=NULL){
  R <- pE_est/sum(pE_est)
  set.seed(seed_r)
  ratio1 <- stats::rmultinom(1, 1, prob = R)
  d <- which(ratio1==1)
  return(d)
}


#' Title: toxicity estimation
#'
#' @param Dat A data frame.
#' @param I A number.
#' @param M A number.
#' @param M_prob A numeric vector.
#' @param DLT_skeleton A numeric list.
#' @param DLT_thresh A number.
#' @param model A string.
#' @param para_prior A string.
#' @param beta_mean A number.
#' @param beta_sd A number.
#' @param intcpt_lgst1 A number.
#' @param beta_shape A number.
#' @param beta_inverse_scale A number.
#' @param alpha_mean A number.
#' @param alpha_sd A number.
#' @param alpha_shape A number.
#' @param alpha_inverse_scale A number.
#' @param seed A number.
#'
#' @return list(AR=AR, M_prob=M_prob) -> AR for acceptable set, M_prob for posterior density of toxicity orderings
#' @export

toxicity_est <- function(Dat, I, M, M_prob, DLT_skeleton, DLT_thresh,
                        model, para_prior,
                        beta_mean, beta_sd, intcpt_lgst1, beta_shape, beta_inverse_scale,
                        alpha_mean, alpha_sd, alpha_shape, alpha_inverse_scale,
                        seed=NULL){

  currN <- nrow(Dat)

  if (currN == 0){
    if (is.null(M_prob)){
      M_prob <- rep(1/M, M) # equal prob for prior information of toxicity orderings
      set.seed(seed)
      mStar <- sample(c(1:M), 1)
    } else {
      mStar_candidate <- which(M_prob == max(M_prob))
      set.seed(seed)
      mStar <- ifelse(length(mStar_candidate) == 1, mStar_candidate, sample(mStar_candidate, 1))
    }
    DLT_skeleton_mStar <- DLT_skeleton[[mStar]]

    AR <- which(DLT_skeleton_mStar <= DLT_thresh) # acceptable set {doses}
    M_prob <- M_prob

  } else {
    tau_m <- M_prob

    temp <- rep(NA, M)
    omega <- rep(NA, M)
    beta_hat <- rep(NA, M)
    alpha_hat <- rep(NA, M)
    y <- Dat$DLT

    for (m in 1:M){
      DLT_skeleton_m <- DLT_skeleton[[m]]

      if (model == "tanh"){
        if (para_prior == "exponential"){
          dosescaled <- atanh(2*DLT_skeleton_m^(1/beta_mean)-1)
          x_m <- dosescaled[Dat$DoseLevel]
          temp[m] <- stats::integrate(tanh_ExpPriorLikelihood, 0, Inf, beta_mean, x_m, y)$value
          omega[m] <- tau_m[m]*stats::integrate(tanh_ExpPriorLikelihood, 0, Inf, beta_mean, x_m, y)$value
          beta_hat[m] <- stats::integrate(tanh_ExpPriorLikelihood_est, 0, Inf, beta_mean, x_m, y)$value/temp[m]
        }
      } else if (model == "empiric"){
        if (para_prior == "normal"){
          dosescaled <- DLT_skeleton_m^exp(-beta_mean)
          x_m <- dosescaled[Dat$DoseLevel]
          temp[m] <- stats::integrate(empiric_NormalPriorLikelihood, -Inf, Inf, beta_mean, beta_sd, x_m, y)$value
          omega[m] <- tau_m[m]*stats::integrate(empiric_NormalPriorLikelihood, -Inf, Inf, beta_mean, beta_sd, x_m, y)$value
          beta_hat[m] <- stats::integrate(empiric_NormalPriorLikelihood_est, -Inf, Inf, beta_mean, beta_sd, x_m, y)$value/temp[m]
        } else if (para_prior == "gamma"){
          dosescaled <- DLT_skeleton_m^(beta_inverse_scale/beta_shape)
          x_m <- dosescaled[Dat$DoseLevel]
          temp[m] <- stats::integrate(empiric_GammaPriorLikelihood, 0, Inf, beta_shape, beta_inverse_scale, x_m, y)$value
          omega[m] <- tau_m[m]*stats::integrate(empiric_GammaPriorLikelihood, 0, Inf, beta_shape, beta_inverse_scale, x_m, y)$value
          beta_hat[m] <- stats::integrate(empiric_GammaPriorLikelihood_est, 0, Inf, beta_shape, beta_inverse_scale, x_m, y)$value/temp[m]
        }
      } else if (model == "logistic"){
        if (para_prior == "normal"){
          dosescaled <- (log(DLT_skeleton_m/(1-DLT_skeleton_m)) - intcpt_lgst1)/exp(beta_mean)
          x_m <- dosescaled[Dat$DoseLevel]
          temp[m] <- stats::integrate(logisticOnePara_NormalPriorLikelihood, -Inf, Inf, beta_mean, beta_sd, intcpt_lgst1, x_m, y)$value
          omega[m] <- tau_m[m]*stats::integrate(logisticOnePara_NormalPriorLikelihood, -Inf, Inf, beta_mean, beta_sd, intcpt_lgst1, x_m, y)$value
          beta_hat[m] <- stats::integrate(logisticOnePara_NormalPriorLikelihood_est, -Inf, Inf, beta_mean, beta_sd, intcpt_lgst1, x_m, y)$value/temp[m]
        } else if (para_prior == "gamma"){
          dosescaled <- (log(DLT_skeleton_m/(1-DLT_skeleton_m)) - intcpt_lgst1)/(beta_shape/beta_inverse_scale)
          x_m <- dosescaled[Dat$DoseLevel]
          temp[m] <- stats::integrate(logisticOnePara_GammaPriorLikelihood, 0, Inf, beta_shape, beta_inverse_scale, intcpt_lgst1, x_m, y)$value
          omega[m] <- tau_m[m]*stats::integrate(logisticOnePara_GammaPriorLikelihood, 0, Inf, beta_shape, beta_inverse_scale, intcpt_lgst1, x_m, y)$value
          beta_hat[m] <- stats::integrate(logisticOnePara_GammaPriorLikelihood_est, 0, Inf, beta_shape, beta_inverse_scale, intcpt_lgst1, x_m, y)$value/temp[m]
        }
      } else if (model == "logistic2"){
        if (para_prior == "normal"){
          dosescaled <- (log(DLT_skeleton_m/(1-DLT_skeleton_m)) - alpha_mean)/exp(beta_mean)
          x_m <- dosescaled[Dat$DoseLevel]

          logisticTwoPara_NormalPriorLikelihood_int <- function(intcpt){
            stats::integrate(logisticTwoPara_NormalPriorLikelihood, lower=-Inf, upper=Inf, alpha0=intcpt, alpha_mean, alpha_sd, beta_mean, beta_sd, x_m, y, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          }
          logisticTwoPara_NormalPriorLikelihood_slope <- function(slope){
            stats::integrate(logisticTwoPara_NormalPriorLikelihood, lower=-Inf, upper=Inf, alpha1=slope, alpha_mean, alpha_sd, beta_mean, beta_sd, x_m, y, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          }

          temp[m] <- stats::integrate(Vectorize(logisticTwoPara_NormalPriorLikelihood_int), -Inf, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          omega[m] <- tau_m[m]*temp[m]
          alpha_hat[m] = stats::integrate(function(intcpt) {intcpt*logisticTwoPara_NormalPriorLikelihood_int(intcpt)}, -Inf, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value/temp[m]
          beta_hat[m] = stats::integrate(function(slope) {slope*logisticTwoPara_NormalPriorLikelihood_slope(slope)}, -Inf, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value/temp[m]

        } else if (para_prior == "gamma"){
          dosescaled <- (log(DLT_skeleton_m/(1-DLT_skeleton_m)) - (alpha_shape/alpha_inverse_scale))/(beta_shape/beta_inverse_scale)
          x_m <- dosescaled[Dat$DoseLevel]

          logisticTwoPara_GammaPriorLikelihood_int <- function(intcpt){
            stats::integrate(logisticTwoPara_GammaPriorLikelihood, lower=0, upper=Inf, alpha0=intcpt, alpha_shape, alpha_inverse_scale, beta_shape, beta_inverse_scale, x_m, y, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          }
          logisticTwoPara_GammaPriorLikelihood_slope <- function(slope){
            stats::integrate(logisticTwoPara_GammaPriorLikelihood, lower=0, upper=Inf, alpha1=slope, alpha_shape, alpha_inverse_scale, beta_shape, beta_inverse_scale, x_m, y, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          }

          temp[m] <- stats::integrate(Vectorize(logisticTwoPara_GammaPriorLikelihood_int), 0, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          omega[m] <- tau_m[m]*temp[m]
          alpha_hat[m] = stats::integrate(function(intcpt) {intcpt*logisticTwoPara_GammaPriorLikelihood_int(intcpt)}, 0, Inf, rel.tol = 1e-8, subdivisions = 1000)$value/temp[m]
          beta_hat[m] = stats::integrate(function(slope) {slope*logisticTwoPara_GammaPriorLikelihood_slope(slope)}, 0, Inf, rel.tol = 1e-8, subdivisions = 1000)$value/temp[m]
        }
      }
    }

    M_prob <- omega/sum(omega)
    mStar_candidate <- which(M_prob==max(M_prob))
    mStar <- ifelse(length(mStar_candidate) == 1, mStar_candidate, sample(mStar_candidate, 1))

    # get scaled skeletons
    DLT_skeleton_mStar <- DLT_skeleton[[mStar]]
    if (model == "tanh"){
      if (para_prior == "exponential"){
        DLT_skeleton_mStar1 <- atanh(2*DLT_skeleton_mStar^(1/beta_mean)-1)
      }
    } else if (model == "empiric"){
      if (para_prior == "normal"){
        DLT_skeleton_mStar1 <- DLT_skeleton_mStar^exp(-beta_mean)
      } else if (para_prior == "gamma"){
        DLT_skeleton_mStar1 <- DLT_skeleton_mStar^(beta_inverse_scale/beta_shape)
      }
    } else if (model == "logistic"){
      if (para_prior == "normal"){
        DLT_skeleton_mStar1 <- (log(DLT_skeleton_mStar/(1-DLT_skeleton_mStar)) - intcpt_lgst1)/exp(beta_mean)
      } else if (para_prior == "gamma"){
        DLT_skeleton_mStar1 <- (log(DLT_skeleton_mStar/(1-DLT_skeleton_mStar)) - intcpt_lgst1)/(beta_shape/beta_inverse_scale)
      }
    } else if (model == "logistic2"){
      if (para_prior == "normal"){
        DLT_skeleton_mStar1 <- (log(DLT_skeleton_mStar/(1-DLT_skeleton_mStar)) - alpha_mean)/exp(beta_mean)
      } else if (para_prior == "gamma"){
        DLT_skeleton_mStar1 <- (log(DLT_skeleton_mStar/(1-DLT_skeleton_mStar)) - (alpha_shape/alpha_inverse_scale))/(beta_shape/beta_inverse_scale)
      }
    }

    if (model == "tanh"){
      piT_hat <- ((tanh(DLT_skeleton_mStar1+1)/2)^beta_hat[mStar])
    } else if (model == "empiric"){
      piT_hat <- DLT_skeleton_mStar1^(exp(beta_hat[mStar]))
    } else if (model == "logistic"){
      piT_hat <- 1/(1 + exp(-intcpt_lgst1-exp(beta_hat[mStar])*DLT_skeleton_mStar1))
    } else if (model == "logistic2"){
      piT_hat <- 1/(1 + exp(-alpha_hat[mStar]-exp(beta_hat[mStar])*DLT_skeleton_mStar1))
    }

    AR <- which(piT_hat <= DLT_thresh)
  }
  return(list(AR=AR, M_prob=M_prob))
}


#' Title: efficacy estimation
#'
#' @param Dat A data frame.
#' @param AR A numeric vector.
#' @param I A number.
#' @param K A number.
#' @param K_prob A numeric vector or Null.
#' @param efficacy_skeleton A numeric list.
#' @param Nphase A number.
#' @param model A string.
#' @param para_prior A string.
#' @param theta_mean A number.
#' @param theta_sd A number.
#' @param theta_intcpt_lgst1 A number.
#' @param theta_shape A number.
#' @param theta_inverse_scale A number.
#' @param alphaT_mean A number.
#' @param alphaT_sd A number.
#' @param alphaT_shape A number.
#' @param alphaT_inverse_scale A number.
#' @param seed A number.
#' @param seed_rand A number.
#' @param seed_max A number.
#'
#' @return list(di=di, K_prob=K_prob) -> di for next recommended dose level, K_prob for posterior density of efficacy orderings
#' @export

efficacy_est <- function(Dat, AR, I, K, K_prob, efficacy_skeleton, Nphase,
                        model, para_prior,
                        theta_mean, theta_sd, theta_intcpt_lgst1, theta_shape, theta_inverse_scale,
                        alphaT_mean, alphaT_sd, alphaT_shape, alphaT_inverse_scale,
                        seed=NULL, seed_rand=NULL, seed_max=NULL){
  currN <- nrow(Dat)
  
  if (currN == 0){
    if (is.null(K_prob)){
      K_prob <- rep(1/K, K) # equal prob for prior information of toxicity orderings
      set.seed(seed)
      kStar <- sample(c(1:K), 1)
    } else {
      kStar_candidate <- which(K_prob == max(K_prob))
      set.seed(seed)
      kStar <- ifelse(length(kStar_candidate) == 1, kStar_candidate, sample(kStar_candidate, 1))
    }

    efficacy_skeleton_kStar <- efficacy_skeleton[[kStar]]

    # 1st patient allocation ratio
    if (length(AR)==0){
      di <- 1 # allocation to the lowest does for next patient
    } else{
      piE_hat <- efficacy_skeleton_kStar[AR]
      di <- randomization_phase(piE_hat, seed_rand) # x1=di
    }

  } else {
    psi_k <- K_prob

    temp1 <- rep(NA, K)
    omega1 <- rep(NA, K)
    theta_hat <- rep(NA, K)
    alpha_hat <- rep(NA, K)
    y <- Dat$ORR

    for (k in 1:K){
      efficacy_skeleton_k = efficacy_skeleton[[k]]

      if (model == "tanh"){
        if (para_prior == "exponential"){
          scaled1 <- atanh(2*efficacy_skeleton_k^(1/theta_mean)-1)
          x_k <- scaled1[Dat$DoseLevel]
          temp1[k] <- stats::integrate(tanh_ExpPriorLikelihood, 0, Inf, theta_mean, x_k, y)$value
          omega1[k] <- psi_k[k]*stats::integrate(tanh_ExpPriorLikelihood, 0, Inf, theta_mean, x_k, y)$value
          theta_hat[k] <- stats::integrate(tanh_ExpPriorLikelihood_est, 0, Inf, theta_mean, x_k, y)$value/temp1[k]
        }
      } else if (model == "empiric"){
        if (para_prior == "normal"){
          scaled1 <- efficacy_skeleton_k^exp(-theta_mean)
          x_k <- scaled1[Dat$DoseLevel]
          temp1[k] <- stats::integrate(empiric_NormalPriorLikelihood, -Inf, Inf, theta_mean, theta_sd, x_k, y, stop.on.error = FALSE)$value
          omega1[k] <- psi_k[k]*stats::integrate(empiric_NormalPriorLikelihood, -Inf, Inf, theta_mean, theta_sd, x_k, y, stop.on.error = FALSE)$value
          theta_hat[k] <- stats::integrate(empiric_NormalPriorLikelihood_est, -Inf, Inf, theta_mean, theta_sd, x_k, y, stop.on.error = FALSE)$value/temp1[k]
        } else if (para_prior == "gamma"){
          scaled1 <- efficacy_skeleton_k^(theta_inverse_scale/theta_shape)
          x_k <- scaled1[Dat$DoseLevel]
          temp1[k] <- stats::integrate(empiric_GammaPriorLikelihood, 0, Inf, theta_shape, theta_inverse_scale, x_k, y, stop.on.error = FALSE)$value
          omega1[k] <- psi_k[k]*stats::integrate(empiric_GammaPriorLikelihood, 0, Inf, theta_shape, theta_inverse_scale, x_k, y, stop.on.error = FALSE)$value
          theta_hat[k] <- stats::integrate(empiric_GammaPriorLikelihood_est, 0, Inf, theta_shape, theta_inverse_scale, x_k, y, stop.on.error = FALSE)$value/temp1[k]
        }
      } else if (model == "logistic"){
        if (para_prior == "normal"){
          scaled1 <- (log(efficacy_skeleton_k/(1-efficacy_skeleton_k)) - theta_intcpt_lgst1)/exp(theta_mean)
          x_k <- scaled1[Dat$DoseLevel]
          temp1[k] <- stats::integrate(logisticOnePara_NormalPriorLikelihood, -Inf, Inf, theta_mean, theta_sd, theta_intcpt_lgst1, x_k, y)$value
          omega1[k] <- psi_k[k]*stats::integrate(logisticOnePara_NormalPriorLikelihood, -Inf, Inf, theta_mean, theta_sd, theta_intcpt_lgst1, x_k, y)$value
          theta_hat[k] <- stats::integrate(logisticOnePara_NormalPriorLikelihood_est, -Inf, Inf, theta_mean, theta_sd, theta_intcpt_lgst1, x_k, y)$value/temp1[k]
        } else if (para_prior == "gamma"){
          scaled1 <- (log(efficacy_skeleton_k/(1-efficacy_skeleton_k)) - theta_intcpt_lgst1)/(theta_shape/theta_inverse_scale)
          x_k <- scaled1[Dat$DoseLevel]
          temp1[k] <- stats::integrate(logisticOnePara_GammaPriorLikelihood, 0, Inf, theta_shape, theta_inverse_scale, theta_intcpt_lgst1, x_k, y)$value
          omega1[k] <- psi_k[k]*stats::integrate(logisticOnePara_GammaPriorLikelihood, 0, Inf, theta_shape, theta_inverse_scale, theta_intcpt_lgst1, x_k, y)$value
          theta_hat[k] <- stats::integrate(logisticOnePara_GammaPriorLikelihood_est, 0, Inf, theta_shape, theta_inverse_scale, theta_intcpt_lgst1, x_k, y)$value/temp1[k]
        }
      } else if (model == "logistic2"){
        if (para_prior == "normal"){
          scaled1 <- (log(efficacy_skeleton_k/(1-efficacy_skeleton_k)) - alphaT_mean)/exp(theta_mean)
          x_k <- scaled1[Dat$DoseLevel]

          logisticTwoPara_NormalPriorLikelihood_int <- function(intcpt){
            stats::integrate(logisticTwoPara_NormalPriorLikelihood, lower=-Inf, upper=Inf, alpha0=intcpt, alphaT_mean, alphaT_sd, theta_mean, theta_sd, x_k, y, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          }
          logisticTwoPara_NormalPriorLikelihood_slope <- function(slope){
            stats::integrate(logisticTwoPara_NormalPriorLikelihood, lower=-Inf, upper=Inf, alpha1=slope, alphaT_mean, alphaT_sd, theta_mean, theta_sd, x_k, y, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          }

          temp1[k] <- stats::integrate(Vectorize(logisticTwoPara_NormalPriorLikelihood_int), -Inf, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          omega1[k] <- psi_k[k]*temp1[k]
          alpha_hat[k] = stats::integrate(function(intcpt) {intcpt*logisticTwoPara_NormalPriorLikelihood_int(intcpt)}, -Inf, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value/temp1[k]
          theta_hat[k] = stats::integrate(function(slope) {slope*logisticTwoPara_NormalPriorLikelihood_slope(slope)}, -Inf, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value/temp1[k]

        } else if (para_prior == "gamma"){
          scaled1 <- (log(efficacy_skeleton_k/(1-efficacy_skeleton_k)) - (alphaT_shape/alphaT_inverse_scale))/(theta_shape/theta_inverse_scale)
          x_k <- scaled1[Dat$DoseLevel]

          logisticTwoPara_GammaPriorLikelihood_int <- function(intcpt){
            stats::integrate(logisticTwoPara_GammaPriorLikelihood, lower=0, upper=Inf, alpha0=intcpt, alphaT_shape, alphaT_inverse_scale, theta_shape, theta_inverse_scale, x_k, y, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          }
          logisticTwoPara_GammaPriorLikelihood_slope <- function(slope){
            stats::integrate(logisticTwoPara_GammaPriorLikelihood, lower=0, upper=Inf, alpha1=slope, alphaT_shape, alphaT_inverse_scale, theta_shape, theta_inverse_scale, x_k, y, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          }

          temp1[k] <- stats::integrate(Vectorize(logisticTwoPara_GammaPriorLikelihood_int), 0, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value
          omega1[k] <- psi_k[k]*temp1[k]
          alpha_hat[k] = stats::integrate(function(intcpt) {intcpt*logisticTwoPara_GammaPriorLikelihood_int(intcpt)}, 0, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value/temp1[k]
          theta_hat[k] = stats::integrate(function(slope) {slope*logisticTwoPara_GammaPriorLikelihood_slope(slope)}, 0, Inf, stop.on.error = FALSE, rel.tol = 1e-8, subdivisions = 1000)$value/temp1[k]

        }
      }
    }

    K_prob <- omega1/sum(omega1)
    kStar_candidate <- which(K_prob==max(K_prob))
    kStar <- ifelse(length(kStar_candidate) == 1, kStar_candidate, sample(kStar_candidate, 1))

    if (length(AR) == 0){
      di = 1
    } else {
      efficacy_skeleton_AR <- efficacy_skeleton[[kStar]][AR]
      # get scaled skeletons
      if (model == "tanh"){
        if (para_prior == "exponential"){
          efficacy_skeleton_AR1 <- atanh(2*efficacy_skeleton_AR^(1/theta_mean)-1)
        }
      } else if (model == "empiric"){
        if (para_prior == "normal"){
          efficacy_skeleton_AR1 <- efficacy_skeleton_AR^exp(-theta_mean)
        } else if (para_prior == "gamma"){
          efficacy_skeleton_AR1 <- efficacy_skeleton_AR^(theta_inverse_scale/theta_shape)
        }
      } else if (model == "logistic"){
        if (para_prior == "normal"){
          efficacy_skeleton_AR1 <- (log(efficacy_skeleton_AR/(1-efficacy_skeleton_AR)) - theta_intcpt_lgst1)/exp(theta_mean)
        } else if (para_prior == "gamma"){
          efficacy_skeleton_AR1 <- (log(efficacy_skeleton_AR/(1-efficacy_skeleton_AR)) - theta_intcpt_lgst1)/(theta_shape/theta_inverse_scale)
        }
      } else if (model == "logistic2"){
        if (para_prior == "normal"){
          efficacy_skeleton_AR1 <- (log(efficacy_skeleton_AR/(1-efficacy_skeleton_AR)) - alphaT_mean)/exp(theta_mean)
        } else if (para_prior == "gamma"){
          efficacy_skeleton_AR1 <- (log(efficacy_skeleton_AR/(1-efficacy_skeleton_AR)) - (alphaT_shape/alphaT_inverse_scale))/(theta_shape/theta_inverse_scale)
        }
      }


      if (model == "tanh"){
        piE_hat <- ((tanh(efficacy_skeleton_AR1+1)/2)^theta_hat[kStar])
      } else if (model == "empiric"){
        piE_hat <- efficacy_skeleton_AR1^(exp(theta_hat[kStar]))
      } else if (model == "logistic"){
        piE_hat <- 1/(1 + exp(-theta_intcpt_lgst1-exp(theta_hat[kStar])*efficacy_skeleton_AR1))
      } else if (model == "logistic2"){
        piE_hat <- 1/(1 + exp(-alpha_hat[kStar]-exp(theta_hat[kStar])*efficacy_skeleton_AR1))
      }
      
      if (currN >= Nphase){
        di = maximization_phase(piE_hat, seed_m=seed_max)
      } else {
        di = randomization_phase(piE_hat, seed_r=seed_rand)
      }
    }
  }

  return(list(di=di, K_prob=K_prob))
}
