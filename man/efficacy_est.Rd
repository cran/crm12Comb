\name{efficacy_est}
\alias{efficacy_est}
\title{
Efficacy estimation
}
\description{
Estimate efficacy using Bayesian inference for each enrolled patient or cohort of patients given the current accumulated data and toxicity estimation.
}
\usage{
efficacy_est(Dat, AR, I, K, K_prob, efficacy_skeleton, Nphase,
             model, para_prior,
             theta_mean, theta_sd, theta_intcpt_lgst1, theta_shape, theta_inverse_scale,
             alphaT_mean, alphaT_sd, alphaT_shape, alphaT_inverse_scale,
             seed=NULL, seed_rand=NULL, seed_max=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Dat}{A data frame for current data with three columns (DoseLevel, DLT and ORR).}
  \item{AR}{A vector for acceptable set.}
  \item{I}{Number of dose combinations.}
  \item{K}{Number of efficacy orderings.}
  \item{K_prob}{A vector of length \eqn{K} denoting prior probabilities of efficacy orderings (sum is 1).}
  \item{efficacy_skeleton}{A list of vector with length \eqn{K} containing efficacy orderings.}
  \item{Nphase}{Number of patients for determination of randomization phase (current number of patients less than Nphase) or maximization phase (current number of patients larger than or equal to Nphase).}
  \item{model}{A character string to specify the model used, must be one of "empiric", "tanh", "logistic" or "logistic2".}
  \item{para_prior}{A character string to specify the prior distribution used for parameters, must be one of "normal" or "gamma" (when model is either "empiric" or "logistic" or "logistic2") or "exponential" (when model is "tanh").}
  \item{theta_mean}{The mean of parameter used when prior="exponential" or "normal", otherwise need to specify NULL.}
  \item{theta_sd}{The standard deviation of parameter used when prior="normal", otherwise need to specify NULL.}
  \item{theta_intcpt_lgst1}{A constant value of intercept from a one-parameter logistic model only used when model="logistic" (suggested value is 3), otherwise need to specify NULL.}
  \item{theta_shape}{The shape parameter used when prior="gamma", otherwise need to specify NULL.}
  \item{theta_inverse_scale}{The scale parameter used when prior="gamma", otherwise need to specify NULL.}
  \item{alphaT_mean}{The mean of intercept parameter of two-parameter logistic model only used when model="logistic2" and prior="normal", otherwise need to specify NULL.}
  \item{alphaT_sd}{The standard deviation of intercept parameter of two-parameter logistic model only used when model="logistic2" and prior="normal", otherwise need to specify NULL.}
  \item{alphaT_shape}{The shape parameter of intercept parameter from a two-parameter logistic model only used when model="logistic2" and prior="gamma", otherwise need to specify NULL.}
  \item{alphaT_inverse_scale}{The scale parameter of intercept parameter from a two-parameter logistic model only used when model="logistic2" and prior="gamma", otherwise need to specify NULL.}
  \item{seed}{An integer for the seed to generate random numbers used for equal ordering prior probabilities, default is NULL.}
  \item{seed_rand}{An integer for the seed to generate random numbers used in randomization phase, default is NULL.}
  \item{seed_max}{An integer for the seed to generate random numbers used in maximization phase, default is NULL.}
}
\details{
The efficacy estimation is based on Bayesian framework by calculating the likelihood function under each orderings, and select the ordering with maximum posterior probability. Then, estimated parameters can be obtained, which will be used for efficacy estimation based on the corresponding link function (specified in \eqn{model} statement).
}
\value{
A list is returned containing the following components:

\item{kStar}{An integer for efficacy ordering with maximum posterior density given the current accumulated data.}
\item{piE_hat}{A vector for estimated efficacy probabilities from acceptable set given the current accumulated data.}
\item{di}{Number of dose level for next enrolled patient or cohort of patients}
\item{K_prob}{A vector for posterier density of efficacy orderings and will be used as prior information for next enrolled patient or cohort of patients.}
}
\references{
Wages, N. A., & Conaway, M. R. (2014). Phase I/II adaptive design for drug combination oncology trials. Statistics in medicine, 33(12), 1990-2003. \doi{10.1002/sim.6097}
}

\seealso{
\code{\link{priorSkeletons}}, \code{\link{get_ordering}}, \code{\link{toxicity_est}}, \code{\link{randomization_phase}}, \code{\link{maximization_phase}}
}
\examples{
### follow same steps as the example in function toxicity_est()
# Generate a data including 3 columns: DoseLevel, DLT, ORR (DLT and ORR are binary outcomes)
currDat <- data.frame(sample(1:9, 6, replace=TRUE), rbinom(6, 1, 0.2), rbinom(6, 1, 0.5))
names(currDat) <- c("DoseLevel", "DLT", "ORR")

# Generate toxicity and efficacy skeleton
DLT_skeleton_p <- priorSkeletons(updelta = 0.045, target = 0.3, npos= 5, ndose = 9,
                                 model = "logistic", prior = "normal", beta_mean = 0, a0 = 3)
eff_skeleton_p <- priorSkeletons(updelta = 0.045, target = 0.5, npos= 5, ndose = 9,
                                 model = "logistic", prior = "normal", beta_mean = 0, a0 = 3)

# Obtain 6 complete orderings for toxicity skeleton and efficacy skeleton
orderings <- get_ordering(doseComb_forMat=c(3,3), type_forMat="matrix")
DLT_skeleton_l <- lapply(orderings, function(or){DLT_skeleton_p[or]})
eff_skeleton_l <- lapply(orderings, function(or){eff_skeleton_p[or]})

# estimate toxicity
tox <- toxicity_est(Dat=currDat, I=9, M=6, M_prob=rep(1/6, 6),
                    DLT_skeleton=DLT_skeleton_l, DLT_thresh=0.33,
                    model="logistic", para_prior="normal",
                    beta_mean=0, beta_sd=1, intcpt_lgst1=3,
                    beta_shape=NULL, beta_inverse_scale=NULL,
                    alpha_mean=NULL, alpha_sd=NULL,
                    alpha_shape=NULL, alpha_inverse_scale=NULL,
                    seed=42)

### efficacy estimation
eff <- efficacy_est(Dat=currDat, AR=tox$AR, I=9, K=6, K_prob=rep(1/6, 6),
                    efficacy_skeleton=eff_skeleton_l, Nphas=20,
                    model="logistic", para_prior="normal",
                    theta_mean=0, theta_sd=1, theta_intcpt_lgst1=3,
                    theta_shape=NULL, theta_inverse_scale=NULL,
                    alphaT_mean=NULL, alphaT_sd=NULL,
                    alphaT_shape=NULL, alphaT_inverse_scale=NULL,
                    seed=1, seed_rand=2, seed_max=3)
}

