#' Title: plot results of each outcome from the given 6 scenarios
#' 4 settings with multiple inputs: need to fix three and plot each outcome vs. the remaining one by 6 scenarios
#' 1. N: maximum sample size -> 40, 50, 60
#' 2. nR: subset sample size -> 10, 20, 30
#' 3. Skeleton: two sets of skeletons for toxicity and efficacy
#' 4. corr: correlation between toxicity and efficacy binary outcomes -> 0, -2.049, 0.814
#'
#' @param dat A data frame.
#' @param outcome 
#' @param outname A string.
#' @param N A number or Null.
#' @param nR A number or Null.
#' @param Skeleton A number or Null.
#' @param corr A number or Null.
#'
#' @return plot (ggplot output)
#' @export

sample_plot <- function(dat, outcome, outname, N = NULL, nR = NULL, Skeleton = NULL, corr = NULL){
  Model <- NULL
  if (is.null(N)){
    subdat <- dat[which(dat$Skeleton == Skeleton & dat$corr == corr & dat$nR == nR),]
    p <- ggplot2::ggplot(subdat, ggplot2::aes_string(x="N", y=outcome, group="Model")) +
      ggplot2::geom_point(ggplot2::aes(color=Model)) +
      ggplot2::geom_line(ggplot2::aes(color=Model)) +
      ggplot2::scale_color_manual(values = c("#9900CC", "#0099CC", "#00CC66", "#FF6600")) +
      ggplot2::facet_wrap(~Scenario, labeller = ggplot2::labeller(Scenario = ~ paste("Scenario:", .), .multi_line = FALSE), scales = "free_y") +
      ggplot2::scale_x_continuous(breaks=unique(dat$N)) + 
      ggplot2::ylab(outname) +
      ggplot2::theme(legend.position = "bottom")
  } else if (is.null(nR)){
    subdat <- dat[which(dat$Skeleton == Skeleton & dat$corr == corr & dat$N == N),]
    p <- ggplot2::ggplot(subdat, ggplot2::aes_string(x="nR", y=outcome, group="Model")) +
      ggplot2::geom_point(ggplot2::aes(color=Model)) +
      ggplot2::geom_line(ggplot2::aes(color=Model)) +
      ggplot2::scale_color_manual(values = c("#9900CC", "#0099CC", "#00CC66", "#FF6600")) +
      ggplot2::facet_wrap(~Scenario, labeller = ggplot2::labeller(Scenario = ~ paste("Scenario:", .), .multi_line = FALSE), scales = "free_y") +
      ggplot2::scale_x_continuous(breaks=unique(dat$nR)) + 
      ggplot2::ylab(outname) +
      ggplot2::theme(legend.position = "bottom")
  } else if (is.null(corr)){
    subdat <- dat[which(dat$Skeleton == Skeleton & dat$nR == nR & dat$N == N),]
    p <- ggplot2::ggplot(subdat, ggplot2::aes_string(x="corr", y=outcome, group="Model")) +
      ggplot2::geom_point(ggplot2::aes(color=Model)) +
      ggplot2::geom_line(ggplot2::aes(color=Model)) +
      ggplot2::scale_color_manual(values = c("#9900CC", "#0099CC", "#00CC66", "#FF6600")) +
      ggplot2::facet_wrap(~Scenario, labeller = ggplot2::labeller(Scenario = ~ paste("Scenario:", .), .multi_line = FALSE), scales = "free_y") +
      ggplot2::scale_x_continuous(breaks=unique(dat$corr)) + 
      ggplot2::ylab(outname) +
      ggplot2::theme(legend.position = "bottom")
  } else if (is.null(Skeleton)){
    subdat <- dat[which(dat$corr == corr & dat$nR == nR & dat$N == N),]
    p <- ggplot2::ggplot(subdat, ggplot2::aes_string(x="Skeleton", y=outcome, group="Model")) +
      ggplot2::geom_point(ggplot2::aes(color=Model)) +
      ggplot2::geom_line(ggplot2::aes(color=Model)) +
      ggplot2::scale_color_manual(values = c("#9900CC", "#0099CC", "#00CC66", "#FF6600")) +
      ggplot2::facet_wrap(~Scenario, labeller = ggplot2::labeller(Scenario = ~ paste("Scenario:", .), .multi_line = FALSE), scales = "free_y") +
      ggplot2::scale_x_continuous(breaks=unique(dat$Skeleton)) + 
      ggplot2::ylab(outname) +
      ggplot2::theme(legend.position = "bottom")
  }
  return(p)
}