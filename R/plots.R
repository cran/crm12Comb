#' Title: Plot the patient enrollment sequentially
#'
#' @param data A data frame.
#'
#' @return ggplot plot output.
#' @export

enroll_patient_plot <- function(data){
  id <- DoseLevel <- NULL
  data$id <- as.numeric(rownames(data))

  p <-  ggplot2::ggplot() + 
    ggforce::geom_arc_bar(data=data[which(data$DLT==1), ], mapping=ggplot2::aes(x0 = id, y0 = DoseLevel, r0 = 0, r = 0.4, start=-pi, end=0, fill="Toxicity"), color=NA) +
    ggforce::geom_arc_bar(data=data[which(data$DLT==0), ], mapping=ggplot2::aes(x0 = id, y0 = DoseLevel, r0 = 0, r = 0.4, start=-pi, end=0, fill = "NonToxicity"), color=NA) +
    ggforce::geom_arc_bar(data=data[which(data$ORR==1), ], mapping=ggplot2::aes(x0 = id, y0 = DoseLevel, r0 = 0, r = 0.4, start=0, end=pi, fill="Efficacy"), color=NA) +
    ggforce::geom_arc_bar(data=data[which(data$ORR==0), ], mapping=ggplot2::aes(x0 = id, y0 = DoseLevel, r0 = 0, r = 0.4, start=0, end=pi, fill = "NonEfficacy"), color=NA) +
    ggplot2::scale_fill_manual(values = c(Toxicity="red", NonToxicity="gray40", Efficacy="green", NonEfficacy="gray70"),
                               labels = c(Toxicity="Toxicity", NonToxicity="Non-Toxicity", Efficacy="Efficacy", NonEfficacy="Non-Efficacy"),
                               limits = c("Toxicity", "NonToxicity", "Efficacy", "NonEfficacy")) +
    ggplot2::theme(legend.position="top") + 
    ggplot2::theme(legend.title=ggplot2::element_blank()) + 
    ggplot2::scale_x_continuous(breaks=data$id) +
    ggplot2::scale_y_continuous(breaks=seq(1, max(data$DoseLevel))) +
    ggplot2::xlab("Patients") + ggplot2::ylab("Dose Combination")

  return(p)
}

#' Title: Plot the patient allocation for single trial
#'
#' @param data A data frame.
#'
#' @return ggplot plot output.
#' @export
#' 
patient_allocation_plot <- function(data){
  count <- DoseLevel <- NULL
  p <- ggplot2::ggplot(data, ggplot2::aes(x=DoseLevel)) +
    ggplot2::geom_bar() + ggplot2::geom_text(stat="count", ggplot2::aes(label=ggplot2::after_stat(count)), vjust=-1) + 
    ggplot2::scale_x_continuous(breaks=seq(1, max(data$DoseLevel))) + ggplot2::xlab("Dose Combination") + ggplot2::ylab("Number of Patients")
  
  return(p)
}

#' Title: Plot the patient allocation for multiple trials
#'
#' @param SimsRes A list.
#'
#' @return ggplot plot output.
#' @export
#' 
ODC_plot <- function(SimsRes){
  count <- ODC <- NULL
  dat <- data.frame(ODC=unlist(SimsRes$ODC))
  n <- length(SimsRes)
  dat <- stats::na.omit(dat)
  
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = ODC)) +
    ggplot2::geom_bar() + ggplot2::geom_text(stat="count", ggplot2::aes(label=ggplot2::after_stat(count)), vjust=-1) + 
    ggplot2::scale_x_continuous(breaks=seq(1, max(dat$ODC))) + ggplot2::xlab("Dose Combination") + ggplot2::ylab("Number of Trials")
  
  return(p)
}

