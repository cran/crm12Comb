#' Title: generate correlated bivariate binary variables
#'
#' @param cohortsize A number.
#' @param pT A number.
#' @param pE A number.
#' @param psi A number.
#' @param seed A number.
#'
#' @return t(out) -> matrix of bivariate binary outcomes
#' @export
#'
#' @references Murtaugh, P. A., & Fisher, L. D. (1990). Bivariate binary models of efficacy and toxicity in dose-ranging trials. Communications in Statistics-Theory and Methods, 19(6), 2003-2020.Chicago
#'             \url{https://www.tandfonline.com/doi/abs/10.1080/03610929008830305}
#' association parameter for efficacy and toxicity: psi, psi=0 means independent

rBin2Corr <- function(cohortsize, pT, pE, psi, seed=NULL){
  pi_11 <- pE*pT + pE*(1-pE)*pT*(1-pT)*(exp(psi)-1)/(exp(psi)+1)
  pi_00 <- (1-pE)*(1-pT) + pE*(1-pE)*pT*(1-pT)*(exp(psi)-1)/(exp(psi)+1)
  pi_01 <- pE*(1-pT) - pE*(1-pE)*pT*(1-pT)*(exp(psi)-1)/(exp(psi)+1)
  pi_10 <- (1-pE)*pT - pE*(1-pE)*pT*(1-pT)*(exp(psi)-1)/(exp(psi)+1)

  possibles <- matrix(c(1,1,0,0,0,1,1,0),nrow=4, byrow=TRUE)
  set.seed(seed)
  nxtclass <- stats::rmultinom(n=cohortsize, size=1, prob=c(pi_11,pi_00,pi_01,pi_10))
  out <- apply(nxtclass, 2, function(v){possibles[which(v==1),]})

  return(t(out))
}
