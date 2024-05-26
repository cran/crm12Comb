\name{sample_plot}
\alias{sample_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample plot for a given output results
}
\description{
This function is used to generate the plot of relationships between outcomes and total average sample size, number of patients for determination of randomization phase, skeleton, and Association parameter for efficacy and toxicity binary outcome by different scenarios and link functions.
}
\usage{
sample_plot(dat, outcome, outname, N = NULL, nR = NULL, Skeleton = NULL, corr = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Input data used for plot.}
  \item{outcome}{Column name for the outcome used in the plot.}
  \item{outname}{A string for the name of outcome.}
  \item{N}{Maximum sample size, if not fixed, use the default value NULL.}
  \item{nR}{Number of patients for determination of randomization phase, if not fixed, use the default value NULL.}
  \item{Skeleton}{Two skeletons with number 1 and 2, if not fixed, use the default value NULL.}
  \item{corr}{Association parameter for efficacy and toxicity, if not fixed, use the default value NULL.}
}

\details{
4 settings with multiple inputs: need to fix three and plot each outcome vs. the remaining one by 6 scenarios:

1. N: maximum sample size -> 40, 50, 60 \cr
2. nR: subset sample size -> 10, 20, 30 \cr
3. Skeleton: two sets of skeletons for toxicity and efficacy \cr
4. corr: correlation between toxicity and efficacy binary outcomes -> 0, -2.049, 0.814
}

\value{
Returns a ggplot object.
}

\examples{
# load the data stored in the crm12Comb package
file_path <- system.file("extdata", "examples_results.RData", package = "crm12Comb")
if (file_path == "") stop("Data file not found")
load(file_path)

# fix the number of patients for determination of randomization phase, skeleton, 
# and Association parameter for efficacy and toxicity binary outcome
# plot the relationship between 
# "Probability of ODC as target combinations" vs. total average sample size
sample_plot(examples_results, outcome = "prob_target", 
  outname = "Probability of ODC as target combinations", 
  N = NULL, nR = 20, Skeleton = 1, corr = 0)

}
