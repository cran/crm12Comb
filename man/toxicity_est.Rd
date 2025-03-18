\name{toxicity_est}
\alias{toxicity_est}
\title{
Toxicity estimation
}
\description{
Estimate toxicity using Bayesian inference for each enrolled patient or cohort of patients given the current accumulated data.
}
\usage{
toxicity_est(Dat, I, M, M_prob, DLT_skeleton, DLT_thresh,
             model, para_prior,
             beta_mean, beta_sd, intcpt_lgst1, beta_shape, beta_inverse_scale,
             alpha_mean, alpha_sd, alpha_shape, alpha_inverse_scale, seed=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Dat}{A data frame for current data with three columns (DoseLevel, DLT and ORR).}
  \item{I}{Number of dose combinations.}
  \item{M}{Number of toxicity orderings.}
  \item{M_prob}{A vector of length \eqn{M} denoting prior probabilities of toxicity orderings (sum is 1).}
  \item{DLT_skeleton}{A list of vector with length \eqn{M} containing toxicity orderings.}
  \item{DLT_thresh}{DLT threshold to define acceptable set.}
  \item{model}{A character string to specify the model used, must be one of "empiric", "tanh", "logistic" or "logistic2".}
  \item{para_prior}{A character string to specify the prior distribution used for parameters, must be one of "normal" or "gamma" (when model is either "empiric" or "logistic" or "logistic2") or "exponential" (when model is "tanh").}
  \item{beta_mean}{The mean of parameter used when prior="exponential" or "normal", otherwise need to specify NULL.}
  \item{beta_sd}{The standard deviation of parameter used when prior="normal", otherwise need to specify NULL.}
  \item{intcpt_lgst1}{A constant value of intercept from a one-parameter logistic model only used when model="logistic" (suggested value is 3), otherwise need to specify NULL.}
  \item{beta_shape}{The shape parameter used when prior="gamma", otherwise need to specify NULL.}
  \item{beta_inverse_scale}{The scale parameter used when prior="gamma", otherwise need to specify NULL.}
  \item{alpha_mean}{The mean of intercept parameter of two-parameter logistic model only used when model="logistic2" and prior="normal", otherwise need to specify NULL.}
  \item{alpha_sd}{The standard deviation of intercept parameter of two-parameter logistic model only used when model="logistic2" and prior="normal", otherwise need to specify NULL.}
  \item{alpha_shape}{The shape parameter of intercept parameter from a two-parameter logistic model only used when model="logistic2" and prior="gamma", otherwise need to specify NULL.}
  \item{alpha_inverse_scale}{The scale parameter of intercept parameter from a two-parameter logistic model only used when model="logistic2" and prior="gamma", otherwise need to specify NULL.}
  \item{seed}{An integer for the seed to generate random numbers used for equal ordering prior probabilities, default is NULL.}
}
\details{
The toxicity estimation is based on Bayesian framework by calculating the likelihood function under each orderings, and select the ordering with maximum posterior probability. Then, estimated parameters can be obtained, which will be used for toxicity estimation based on the corresponding link function (specified in \eqn{model} statement).
}
\value{
A list is returned containing the following components:

\item{mStar}{An integer for toxicity ordering with maximum posterior density given the current accumulated data.}
\item{piT_hat}{A vector for estimated toxicity probabilities given the current accumulated data.}
\item{AR}{A vector for dose levels belonging to acceptable set (estimated DLT rate less than pre-specified DLT threshold).}
\item{M_prob}{A vector for posterier density of toxicity orderings and will be used as prior information for next enrolled patient or cohort of patients.}
}
\references{
Wages, N. A., & Conaway, M. R. (2014). Phase I/II adaptive design for drug combination oncology trials. Statistics in medicine, 33(12), 1990-2003. \doi{10.1002/sim.6097}
}

\seealso{
\code{\link{priorSkeletons}}, \code{\link{get_ordering}}
}
\examples{
# Generate a data including 3 columns: DoseLevel, DLT, ORR (DLT and ORR are binary outcomes)
currDat <- data.frame(sample(1:9, 6, replace=TRUE), rbinom(6, 1, 0.2), rbinom(6, 1, 0.5))
names(currDat) <- c("DoseLevel", "DLT", "ORR")

# Generate toxicity skeleton
DLT_skeleton_p <- priorSkeletons(updelta = 0.045, target = 0.3, npos= 5, ndose = 9,
                                 model = "logistic", prior = "normal", beta_mean = 0, a0 = 3)

# Obtain 6 complete orderings for toxicity skeleton
orderings <- get_ordering(doseComb_forMat=c(3,3), type_forMat="matrix")
DLT_skeleton_l <- lapply(orderings, function(or){DLT_skeleton_p[or]})

# estimate toxicity
tox <- toxicity_est(Dat=currDat, I=9, M=6, M_prob=rep(1/6, 6),
                    DLT_skeleton=DLT_skeleton_l, DLT_thresh=0.3,
                    model="logistic", para_prior="normal",
                    beta_mean=0, beta_sd=1, intcpt_lgst1=3,
                    beta_shape=NULL, beta_inverse_scale=NULL,
                    alpha_mean=NULL, alpha_sd=NULL,
                    alpha_shape=NULL, alpha_inverse_scale=NULL,
                    seed=42)
}

