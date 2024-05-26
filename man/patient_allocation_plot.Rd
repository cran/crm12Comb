\name{patient_allocation_plot}
\alias{patient_allocation_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot patient allocation for a single trial
}
\description{
This function is used to generate the plot of patient allocation by dose combinations of a single trial.
}
\usage{
patient_allocation_plot(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data frame with 3 columns DoseLevel, DLT for toxicity outcome, and ORR for efficacy outcome.}
}

\value{
Returns a ggplot object.
}

\examples{
# input the scenario with pre-defined true toxicity and efficacy probabilities
scenario <- matrix(c(0.02, 0.05,
                     0.04, 0.10,
                     0.08, 0.15,
                     0.12, 0.32,
                     0.06, 0.10,
                     0.10, 0.15,
                     0.14, 0.25,
                     0.20, 0.35,
                     0.12, 0.18,
                     0.16, 0.22,
                     0.22, 0.35,
                     0.25, 0.40,
                     0.20, 0.24,
                     0.24, 0.35,
                     0.35, 0.45,
                     0.40, 0.50), ncol=2, byrow = TRUE)

# generate skeletons
DLT_skeleton <- priorSkeletons(updelta=0.025, target=0.3, npos=10, ndose=16, 
                               model = "empiric", prior = "normal", beta_mean=0)
Efficacy_skeleton <- priorSkeletons(updelta=0.025, target=0.5, npos=10, ndose=16, 
                                    model = "empiric", prior = "normal", beta_mean=0)

# simulate 1 trial under the same model and prior distribution
simRes <- SIM_phase_I_II(nsim=1, Nmax=40, DoseComb=scenario, input_doseComb_forMat=c(4,4), 
                         input_type_forMat="matrix", input_Nphase=20,
                         input_DLT_skeleton=DLT_skeleton, 
                         input_efficacy_skeleton=Efficacy_skeleton,
                         input_DLT_thresh=0.3, input_efficacy_thresh=0.3,
                         input_cohortsize=1, input_corr=0,
                         input_early_stopping_safety_thresh=0.33,
                         input_early_stopping_futility_thresh=0.2,
                         input_model="empiric", input_para_prior="normal",
                         input_beta_mean=0, input_beta_sd=1,
                         input_theta_mean=0, input_theta_sd=1)

patient_allocation_plot(simRes$datALL[[1]])

}
