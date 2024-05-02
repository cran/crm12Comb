\name{priorSkeletons}
\alias{priorSkeletons}
\title{
Generate the skeletons of toxicity and efficacy
}
\description{
This function is used to generate skeletons of toxicity and efficacy. This is a modifed version based on \code{\link[dfcrm]{getprior}}, which keep the same procedure using empiric and one-parameter logistic models assumed normal priors with \eqn{mean=0} and further add multiple models with various prior distributions including hyperbolic tangent model with exponential prior, empiric/one-parameter logistic models with normal prior and self-input mean values as well as with gamma prior, and two-parameter logistic model with normal/gamma priors.
}
\usage{
priorSkeletons(updelta, target, npos, ndose,
               model = "empiric", prior = "normal",
               alpha_mean=NULL, beta_mean=0, a0 = 3,
               alpha_shape=NULL, alpha_inverse_scale=NULL,
               beta_shape=NULL, beta_inverse_scale=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{updelta}{The half-width of the indifference intervals.}
  \item{target}{The target DLT rate.}
  \item{npos}{The prior guess of the position of MTD.}
  \item{ndose}{The number of testing doses.}
  \item{model}{A character string to specify the model used. The default model is "empiric". Other models include hyperbolic tangent model specified by "tanh", one-parameter logistic model specified by "logistic", and two-parameter logistic model specified by "logistic2".}
  \item{prior}{A character sting to specify the prior distribution of parameter. The default prior is "normal" used together with the model="empiric". Other prior distributions include "exponential" when model="tanh", "gamma" when model="empiric", "normal" and "gamma" when model="logistic" and "logistic2".}
  \item{alpha_mean}{The mean of intercept parameter of two-parameter logistic model only used when model="logistic2" and prior="normal", otherwise will be ignored.}
  \item{beta_mean}{The mean of parameter used when prior="exponential" or "normal", otherwise will be ignored.}
  \item{a0}{A constant value of intercept from a one-parameter logistic model only used when model="logistic" with default value 3, otherwise will be ignored.}
  \item{alpha_shape}{The shape parameter of intercept parameter only used when model="logistic2" and prior="gamma", otherwise will be ignored.}
  \item{alpha_inverse_scale}{The scale parameter of intercept parameter only used when model="logistic2" and prior="gamma", otherwise will be ignored.}
  \item{beta_shape}{The shape parameter used when prior="gamma", otherwise will be ignored.}
  \item{beta_inverse_scale}{The scale parameter used when prior="gamma", otherwise will be ignored.}
}

\value{
A vector of length \eqn{ndose} is returned.
}
\references{
Lee, S. M., & Cheung, Y. K. (2009). Model calibration in the continual reassessment method. Clinical Trials, 6(3), 227-238. \doi{10.1177/1740774509105076}
}
\note{
The skeletons can be either specified by clinical researchers based on history information or directly generated based on this function given specific model and prior distribution.
}

\examples{
# generate skeleton based on empiric model with normal prior
prior <- priorSkeletons(updelta = 0.01, target = 0.25, npos= 5, ndose = 9, beta_mean = 0)

# generate skeleton based on one-parameter logistic model with normal prior
prior <- priorSkeletons(updelta = 0.01, target = 0.25, npos= 5, ndose = 9,
                        model = "logistic", prior = "normal", beta_mean = 0, a0 = 3)
}
