#' Title: Single simulation of phase I/II adaptive design for drug combinations based on CRM design
#'
#' @param nsim A number for number of simulations.
#' @param Nmax A number for maximum sample size of each trial.
#' @param DoseComb A numeric matrix for true toxicity and efficacy probabilities for each dose combination.
#' @param input_doseComb_forMat Either a numeric matrix or a numeric vector or a numeric list.
#' @param input_type_forMat A string.
#' @param input_Nphase A number for number of patients to determine phases.
#' @param input_DLT_skeleton A numeric vector.
#' @param input_efficacy_skeleton A numeric vector.
#' @param input_DLT_thresh A number.
#' @param input_efficacy_thresh A number.
#' @param input_M_prob A numeric vector.
#' @param input_K_prob A numeric vector.
#' @param input_cohortsize A number.
#' @param input_corr A number.
#' @param input_early_stopping_safety_thresh A number.
#' @param input_early_stopping_futility_thresh A number.
#' @param input_model A String.
#' @param input_para_prior A String.
#' @param input_beta_mean A number.
#' @param input_beta_sd A number.
#' @param input_intcpt_lgst1 A number.
#' @param input_beta_shape A number.
#' @param input_beta_inverse_scale A number.
#' @param input_theta_mean A number.
#' @param input_theta_sd A number.
#' @param input_theta_intcpt_lgst1 A number.
#' @param input_theta_shape A number.
#' @param input_theta_inverse_scale A number.
#' @param input_alpha_mean A number.
#' @param input_alpha_sd A number.
#' @param input_alpha_shape A number.
#' @param input_alpha_inverse_scale A number.
#' @param input_alphaT_mean A number.
#' @param input_alphaT_sd A number.
#' @param input_alphaT_shape A number.
#' @param input_alphaT_inverse_scale A number.
#'
#' @return list of operating characteristics
#' @export

SIM_phase_I_II <- function(nsim, Nmax, DoseComb, input_doseComb_forMat, input_type_forMat, input_Nphase,
                           input_DLT_skeleton, input_efficacy_skeleton,
                           input_DLT_thresh=0.3, input_efficacy_thresh=0.3,
                           input_M_prob=NULL, input_K_prob=NULL,
                           input_cohortsize=1, input_corr=0,
                           input_early_stopping_safety_thresh=0.33,
                           input_early_stopping_futility_thresh=0.2,
                           input_model="empiric", input_para_prior="normal",
                           input_beta_mean=0, input_beta_sd=1,
                           input_intcpt_lgst1=3,
                           input_beta_shape=1, input_beta_inverse_scale=1,
                           input_theta_mean=0, input_theta_sd=1,
                           input_theta_intcpt_lgst1=3,
                           input_theta_shape=1, input_theta_inverse_scale=1,
                           input_alpha_mean=3, input_alpha_sd=1,
                           input_alpha_shape=1, input_alpha_inverse_scale=1,
                           input_alphaT_mean=3, input_alphaT_sd=1,
                           input_alphaT_shape=1, input_alphaT_inverse_scale=1){

  ndoses <- length(DoseComb)

  orderings <- get_ordering(doseComb_forMat=input_doseComb_forMat, type_forMat=input_type_forMat)
  DLT_orderings <- lapply(orderings, function(or){input_DLT_skeleton[or]})
  ORR_orderings <- lapply(orderings, function(or){input_efficacy_skeleton[or]})
  
  nM <- length(DLT_orderings)
  nK <- length(ORR_orderings)

  Npatient <- rep(NA, nsim)
  O_DLT <- rep(NA, nsim)
  O_ORR <- rep(NA, nsim)
  ODC <- rep(NA, nsim)

  stop_futility <- 0
  stop_safety <- 0
  toxic <- 0
  target <- 0
  ineffective <- 0

  prop_ODC <- list()
  dn_check_all <- NULL
  datALL <- list()

  for (n in 1:nsim){
    currDat <- data.frame(DoseLevel=integer(), DLT=integer(), ORR=integer())
    stop <- 0
    j <- 0
    while (j <= Nmax){

      tox <- toxicity_est(Dat=currDat, I=ndoses, M=nM, M_prob=input_M_prob,
                          DLT_skeleton=DLT_orderings, DLT_thresh=input_DLT_thresh,
                          model=input_model, para_prior=input_para_prior,
                          beta_mean=input_beta_mean, beta_sd=input_beta_sd,
                          intcpt_lgst1=input_intcpt_lgst1,
                          beta_shape=input_beta_shape, beta_inverse_scale=input_beta_inverse_scale,
                          alpha_mean=input_alpha_mean, alpha_sd=input_alpha_sd,
                          alpha_shape=input_alpha_shape, alpha_inverse_scale=input_alpha_inverse_scale)
      input_M_prob <- tox$M_prob

      eff <- efficacy_est(Dat=currDat, AR=tox$AR, I=ndoses, K=nK, K_prob=input_K_prob,
                          efficacy_skeleton=ORR_orderings, Nphase=input_Nphase,
                          model=input_model, para_prior=input_para_prior,
                          theta_mean=input_theta_mean, theta_sd=input_theta_sd,
                          theta_intcpt_lgst1=input_theta_intcpt_lgst1,
                          theta_shape=input_theta_shape, theta_inverse_scale=input_theta_inverse_scale,
                          alphaT_mean=input_alphaT_mean, alphaT_sd=input_alphaT_sd,
                          alphaT_shape=input_alphaT_shape, alphaT_inverse_scale=input_alphaT_inverse_scale)
      input_K_prob <- eff$K_prob

      tempDat <- data.frame(rep(eff$di,input_cohortsize), rBin2Corr(cohortsize=input_cohortsize,
                                                                    pT=DoseComb[eff$di,1],
                                                                    pE=DoseComb[eff$di,2],
                                                                    psi=input_corr))
      names(tempDat) <- c("DoseLevel", "DLT", "ORR")
      currDat <- rbind(currDat, tempDat)

      j <- j + input_cohortsize

      if (j <= Nmax){
        # early stopping for safety
        lower_piT_d1 <- stats::binom.test(length(which(currDat$DLT == 1)), nrow(currDat),
                                          DoseComb[1,1], alternative="two.sided")$conf.int[1]
        if (lower_piT_d1 > input_early_stopping_safety_thresh){
          stop <- 1
          stop_safety <- stop_safety+1
          Npatient[n] <- nrow(currDat)
          O_DLT[n] <- ifelse(Npatient[n] == 0, 0, length(which(currDat$DLT == 1))/Npatient[n])
          O_ORR[n] <- ifelse(Npatient[n] == 0, 0, length(which(currDat$ORR == 1))/Npatient[n])
          ODC[n] <- NA
          prop_ODC[[n]] <- currDat$DoseLevel
          datALL[[n]] <- currDat
          break
        }

        # early stopping for futility for maximization phase only
        if (j >= input_Nphase) {
          upper_piE_dj = stats::binom.test(length(which(currDat$ORR == 1)), nrow(currDat),
                                           DoseComb[eff$di,2], alternative="two.sided")$conf.int[2]
          if (upper_piE_dj < input_early_stopping_futility_thresh){
            stop <- 1
            stop_futility <- stop_futility + 1
            Npatient[n] <- nrow(currDat)
            O_DLT[n] <- ifelse(Npatient[n] == 0, 0, length(which(currDat$DLT == 1))/Npatient[n])
            O_ORR[n] <- ifelse(Npatient[n] == 0, 0, length(which(currDat$ORR == 1))/Npatient[n])
            ODC[n] <- NA
            prop_ODC[[n]] <- currDat$DoseLevel
            datALL[[n]] <- currDat
            break
          }
        }

      } # end early stopping conditions
    } # end of trial

    if (stop == 0){
      temp_nrow <- nrow(currDat)
      currN <- temp_nrow-input_cohortsize
      Npatient[n] <- currN
      ODC[n] <- currDat$DoseLevel[temp_nrow]
      prop_ODC[[n]] <- currDat$DoseLevel[1:currN]
      O_DLT[n] = length(which(currDat$DLT[1:currN] == 1))/currN
      O_ORR[n] = length(which(currDat$ORR[1:currN] == 1))/currN
      datALL[[n]] <- currDat[1:currN, ]
    }

  } # end of simulation



  toxic <- NULL
  safe <- NULL
  target <- NULL

  prop_ODC1 <- NULL

  toxic <- which(DoseComb[,1] > input_DLT_thresh)
  target <- which(DoseComb[,1] <= input_DLT_thresh & DoseComb[,2] >= input_efficacy_thresh)
  safe <- which(DoseComb[,1] <= input_DLT_thresh & DoseComb[,2] < input_efficacy_thresh)

  prob_safe <- length(which(ODC %in% safe))/nsim
  prob_target <- length(which(ODC %in% target))/nsim
  prob_toxic <- length(which(ODC %in% toxic))/nsim
  mean_SS <- round(mean(Npatient), 3)
  mean_ODC <- round(mean(unlist(sapply(prop_ODC, function(x){length(which(x %in% target))/length(x)}))), 3)

  prob_stop_safety <- stop_safety/nsim
  prob_stop_futility <- stop_futility/nsim
  mean_DLT <- round(mean(O_DLT), 3)
  mean_ORR <- round(mean(O_ORR), 3)

  return(list(prob_safe=prob_safe, prob_target=prob_target, prob_toxic=prob_toxic,
              mean_SS=mean_SS, mean_ODC=mean_ODC,
              prob_stop_safety=prob_stop_safety, prob_stop_futility=prob_stop_futility,
              mean_DLT=mean_DLT, mean_ORR=mean_ORR,
              Npatient=Npatient, ODC=ODC,
              prop_ODC=prop_ODC, datALL=datALL))
}
