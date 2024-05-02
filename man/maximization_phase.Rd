\name{maximization_phase}
\alias{maximization_phase}

\title{
Maximization phase
}
\description{
This function is used to perform maximization to select the dose level with maximum efficacy probability for next patient or cohort of patients allocation when the current sample size is greater than or equal to a pre-specified number.
}
\usage{
maximization_phase(pE_est)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pE_est}{A vector of estimated efficacy probability in the acceptable set.}
}

\details{
If several dose combinations have the same maximum estimated efficacy probability, then randomly select one dose level for next enrolled patient or cohort of patients.
}
\value{
A number is returned indicating the dose level for next patient or cohort of patients allocation.
}

\examples{
# Assume the estimated prbabilities for each dose combination in the acceptable set as:
p_est <- c(0.1, 0.2, 0.3, 0.4)
# Dose level for next enrolled patient or cohort of patients is:
d <- maximization_phase(p_est)
}
