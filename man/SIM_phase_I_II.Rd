\name{SIM_phase_I_II}
\alias{SIM_phase_I_II}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Single simulation of phase I/II adaptive design for drug combinations based on CRM design
}
\description{
This function is to realize the simulations of whole process of combined phase I and phase II design for drug combinations with early stopping rules and flexible user-defined model or prior distribution, returning operating characters for performace evaluations and comparisons.
}
\usage{
SIM_phase_I_II(nsim, Nmax, DoseComb, input_doseComb_forMat, 
               input_type_forMat, input_Nphase,
               input_DLT_skeleton, input_efficacy_skeleton,
               input_DLT_thresh=0.3, input_efficacy_thresh=0.3,
               input_M_prob=NULL, input_K_prob=NULL,
               input_cohortsize=1, input_corr=0,
               input_early_stopping_safety_thresh=0.33,
               input_early_stopping_futility_thresh=0.2,
               input_model="empiric", input_para_prior="normal",
               input_beta_mean=0, input_beta_sd=1,
               input_intcpt_lgst1=3,
               input_beta_shape=1, input_beta_inverse_scale=1,
               input_theta_mean=0, input_theta_sd=1,
               input_theta_intcpt_lgst1=3,
               input_theta_shape=1, input_theta_inverse_scale=1,
               input_alpha_mean=3, input_alpha_sd=1,
               input_alpha_shape=1, input_alpha_inverse_scale=1,
               input_alphaT_mean=3, input_alphaT_sd=1,
               input_alphaT_shape=1, input_alphaT_inverse_scale=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{Number of simulation trials.}
  \item{Nmax}{Number of maximum enrolled patients for each trial.}
  \item{DoseComb}{A matrix with rows of length equaling number of dose combinations and columns of length equaling number of drugs for pre-defined true toxicity and true efficacy probabilities corresponding to each dose combinations.}
  \item{input_doseComb_forMat}{For 2 drugs, either a matrix with columns for all combinations and rows for two drugs when input_type_forMat is "comb", or a vector of length 2 indicating the number of levels of each drug when input_type_forMat is "matrix" should be input. For more than 2 drugs, a list of all possible orderings should be input when input_type_forMat is "self".}
  \item{input_type_forMat}{A character string to indicate the input type of dose combinations. The type_forMat is either "comb" for inputting all tested dose combinations for 2 drugs, "matrix" for entering number of levels of two drugs so that a matrix of combinations can be constructed, or "self" for more than 3 drugs with directly input of all possible orderings.}
  \item{input_Nphase}{Number of patients for determination of randomization phase (current number of patients less than Nphase) or maximization phase (current number of patients larger than or equal to Nphase).}
  \item{input_DLT_skeleton}{A vector with same length of dose combinations for toxicity skeleton for each dose combination.}
  \item{input_efficacy_skeleton}{A vector with same length of dose combinations for efficacy skeleton for each dose combination.}
  \item{input_DLT_thresh}{DLT threshold to define acceptable set.}
  \item{input_efficacy_thresh}{Efficacy threshold to define target combinations.}
  \item{input_M_prob}{A vector with same length of possible toxicity orderings denoting prior probabilities of toxicity orderings (sum is 1).}
  \item{input_K_prob}{A vector with same length of possible efficacy orderings denoting prior probabilities of efficacy orderings (sum is 1).}
  \item{input_cohortsize}{Number of patients in each cohort.}
  \item{input_corr}{Association parameter for efficacy and toxicity, where 0 means toxicity and efficacy is independent.}
  \item{input_early_stopping_safety_thresh}{Safety threshold for early stopping condition.}
  \item{input_early_stopping_futility_thresh}{Futility threshold for early stopping condition.}
  \item{input_model}{A character string to specify the model used, must be one of "empiric", "tanh", "logistic" or "logistic2".}
  \item{input_para_prior}{A character string to specify the prior distribution used for parameters, must be one of "normal" or "gamma" (when input_model is either "empiric" or "logistic" or "logistic2") or "exponential" (when input_model is "tanh").}
  \item{input_beta_mean}{The mean of parameter used when input_para_prior="exponential" or "normal", otherwise will be ignored.}
  \item{input_beta_sd}{The standard deviation of parameter used when input_para_prior="normal", otherwise will be ignored.}
  \item{input_intcpt_lgst1}{A constant value of intercept from a one-parameter logistic model only used when input_model="logistic" (suggested value is 3), otherwise will be ignored.}
  \item{input_beta_shape}{The shape parameter used when input_para_prior="gamma", otherwise will be ignored.}
  \item{input_beta_inverse_scale}{The scale parameter used when input_para_prior="gamma", otherwise will be ignoredL.}
  \item{input_theta_mean}{The mean of parameter used when input_para_prior="exponential" or "normal", otherwise will be ignored.}
  \item{input_theta_sd}{The standard deviation of parameter used when input_para_prior="normal", otherwise will be ignored.}
  \item{input_theta_intcpt_lgst1}{A constant value of intercept from a one-parameter logistic model only used when input_model="logistic" (default value is 3), otherwise will be ignored.}
  \item{input_theta_shape}{The shape parameter used when input_para_prior="gamma", otherwise will be ignored.}
  \item{input_theta_inverse_scale}{The scale parameter used when input_para_prior="gamma", otherwise will be ignored.}
  \item{input_alpha_mean}{The mean of intercept parameter of two-parameter logistic model only used when model="logistic2" and input_para_prior="normal", otherwise will be ignored.}
  \item{input_alpha_sd}{The standard deviation of intercept parameter of two-parameter logistic model only used when input_model="logistic2" and input_para_prior="normal", otherwise will be ignored.}
  \item{input_alpha_shape}{The shape parameter of intercept parameter from a two-parameter logistic model only used when input_model="logistic2" and input_para_prior="gamma", otherwise will be ignored.}
  \item{input_alpha_inverse_scale}{The scale parameter of intercept parameter from a two-parameter logistic model only used when input_model="logistic2" and input_para_prior="gamma", otherwise will be ignored.}
  \item{input_alphaT_mean}{The mean of intercept parameter of two-parameter logistic model only used when model="logistic2" and input_para_prior="normal", otherwise will be ignored.}
  \item{input_alphaT_sd}{The standard deviation of intercept parameter of two-parameter logistic model only used when input_model="logistic2" and input_para_prior="normal", otherwise will be ignored.}
  \item{input_alphaT_shape}{The shape parameter of intercept parameter from a two-parameter logistic model only used when input_model="logistic2" and input_para_prior="gamma", otherwise will be ignored.}
  \item{input_alphaT_inverse_scale}{The scale parameter of intercept parameter from a two-parameter logistic model only used when input_model="logistic2" and input_para_prior="gamma", otherwise will be ignored.}
}
\details{
This function is to realize the whole process of combined phase I and phase II adaptive design for drug combinations based on CRM amonng number of simulation trials. For each trial, basic steps include starting the trial for first patient or cohor of patients allocation, toxicity and efficacy estimations by current data, adaptive randomization or maximization phase for next patient or cohort of patients allocation, updating current data of patients enrollment, determining early stoppings for safety or fuility, and selecting optimal dose combination (ODC) after reaching maximum sample size or stopping early.
}
\value{
A list of operating characteristics is returned containing the following components:

\item{prob_safe}{Probability of simulation trials with ODC identified as safe/ineffective combinations.}
\item{prob_target}{Probability of simulation trials with ODC identified as target combinations.}
\item{prob_toxic}{Probability of simulation trials with ODC identified as toxic combinations.}
\item{mean_SS}{Average number of patients enrolled.}
\item{mean_ODC}{Average proportion of patients allocated to target ODS(s).}
\item{prob_stop_safety}{Probability of simulation trials stopping early for safety.}
\item{prob_stop_futility}{Probability of simulation trials stopping early for futility.}
\item{mean_DLT}{Average observed DLT rate.}
\item{mean_ORR}{Average observed response rate.}
\item{Npatient}{A vector of number of patients enrolled for all simulation trials.}
\item{ODC}{A vector of ODC for all simulation trials.}
\item{prop_ODC}{A list of vectors storing dose levels for each patients allocation.}
\item{datALL}{A list of data frames storing data of (doselevel, toxicity, efficacy) for each patient.}
}
\references{
Wages, N. A., & Conaway, M. R. (2014). Phase I/II adaptive design for drug combination oncology trials. Statistics in medicine, 33(12), 1990-2003. \doi{10.1002/sim.6097}
}


\seealso{
\code{\link{priorSkeletons}}, \code{\link{get_ordering}}, \code{\link{toxicity_est}}, \code{\link{efficacy_est}}, \code{\link{rBin2Corr}}, \code{\link[stats]{binom.test}}
}

\examples{
# input the scenario with pre-defined true toxicity and efficacy probabilities
# use the second scenario from the paper in the reference
scenario <- matrix(c(0.06, 0.05,
                     0.08, 0.10,
                     0.12, 0.20,
                     0.10, 0.10,
                     0.14, 0.25,
                     0.20, 0.40,
                     0.16, 0.20,
                     0.25, 0.35,
                     0.35, 0.50), ncol=2, byrow = TRUE)

# generate skeletons
DLT_skeleton <- priorSkeletons(updelta=0.045, target=0.3, npos=5, ndose=9,
                               model = "empiric", prior = "normal", beta_mean=0)
Efficacy_skeleton <- priorSkeletons(updelta=0.045, target=0.5, npos=5, ndose=9,
                               model = "empiric", prior = "normal", beta_mean=0)

# simulate 1 trial under the same model and prior distribution
simRes <- SIM_phase_I_II(nsim=1, Nmax=40, DoseComb=scenario, input_doseComb_forMat=c(3,3), 
                         input_type_forMat="matrix", input_Nphase=20,
                         input_DLT_skeleton=DLT_skeleton, input_efficacy_skeleton=Efficacy_skeleton,
                         input_DLT_thresh=0.3, input_efficacy_thresh=0.3,
                         input_cohortsize=1, input_corr=0,
                         input_early_stopping_safety_thresh=0.33,
                         input_early_stopping_futility_thresh=0.2,
                         input_model="empiric", input_para_prior="normal",
                         input_beta_mean=0, input_beta_sd=sqrt(1.34),
                         input_theta_mean=0, input_theta_sd=sqrt(1.34))
}

