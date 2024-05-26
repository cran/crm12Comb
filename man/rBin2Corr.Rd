\name{rBin2Corr}
\alias{rBin2Corr}

\title{Generate correlated binary variables}

\description{
Generate correlated bivariate binary outcomes of toxicity and efficacy for a cohort number of patients.
}
\usage{
rBin2Corr(cohortsize, pT, pE, psi, seed=NULL)
}

\arguments{
  \item{cohortsize}{Number of patients in each cohort.}
  \item{pT}{Toxicity probability.}
  \item{pE}{Efficacy probability.}
  \item{psi}{Association parameter for efficacy and toxicity, where psi=0 means toxicity and efficacy is independent.}
  \item{seed}{An integer for the seed to generate random numbers, default is NULL.}
}
\details{
The formula for generating correlated binary variables is
\deqn{\pi_{i,j} = (\pi_E)^i(1-\pi_E)^{1-i}(\pi_T)^j(1-\pi_T)^{1-j} + (-1)^{i+j}\pi_E(1-\pi_E)\pi_T(1-\pi_T)\left(\frac{e^{\psi}-1}{e^{\psi}+1}\right),}
where \eqn{i, j = 0, 1}, so that four probabilities can be calculated for the possible combinations of (toxicity, efficacy) including \eqn{(1,1), (0,0), (0,1), (1,0)} given \eqn{\pi_T} and \eqn{\pi_E}. Multinomial distribution \code{\link[stats]{rmultinom}} is further used to generate bivariate binary outcomes (number equals to cohortsize) based on the four calculated probabilities.

}
\value{
Return a \eqn{cohortsize \times 2} matrix with columns corresponding to toxicity and efficacy, and rows for each observations of binary outcome with 0 for no toxicity (no efficacy) and 1 for toxicity (efficacy) at the first (second) column.
}
\references{
Murtaugh, P. A., & Fisher, L. D. (1990). Bivariate binary models of efficacy and toxicity in dose-ranging trials. Communications in Statistics-Theory and Methods, 19(6), 2003-2020. \doi{10.1080/03610929008830305}

Thall, P. F., & Cook, J. D. (2004). Dose‐finding based on efficacy–toxicity trade‐offs. Biometrics, 60(3), 684-693. \doi{10.1111/j.0006-341X.2004.00218.x}
}
%\author{
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
rBin2Corr(cohortsize = 1, pT = 0.2, pE = 0.4, psi = 0, seed=12)
}
