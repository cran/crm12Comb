\name{randomization_phase}
\alias{randomization_phase}

\title{
Adaptive randomization
}
\description{
This function is used to perform adaptive randomization for next patient or cohort of patients allocation when the current sample size is less than a pre-specified number.
}
\usage{
randomization_phase(pE_est, seed_r=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pE_est}{A vector of estimated efficacy probability in the acceptable set.}
  \item{seed_r}{An integer for the seed to generate random numbers used in randomization phase, default is NULL.}
}

\details{
The dose combination for next patient or cohort of patients allocation is \eqn{d_i} with probability
\deqn{R_i = \frac{\hat{\pi}_E(d_i)}{\sum_i\hat{\pi}_E(d_i)}.}
}

\value{
A number is returned indicating the dose level for next patient or cohort of patients allocation.
}

\examples{
# Assume the estimated prbabilities for each dose combination in the acceptable set as:
p_est <- c(0.1, 0.2, 0.3, 0.4)
# Dose level for next enrolled patient or cohort of patients is:
d <- randomization_phase(p_est)
}
