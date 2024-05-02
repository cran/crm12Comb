## ----setup, eval = FALSE------------------------------------------------------
#  install.packages("./crm12Comb_0.1.0.tar.gz", repos = NULL, type = "source")
#  library(crm12Comb)
#  help(package="crm12Comb")

## ---- eval=FALSE--------------------------------------------------------------
#  ### main function, please use this function to run simulations
#  ?SIM_phase_I_II
#  
#  ### helper functions need to run before simulations (pre-defined settings)
#  ?priorSkeletons
#  
#  ### helper functions included in the main function SIM_phase_I_II()
#  ?rBin2Corr
#  ?ger_ordering
#  ?toxicity_est
#  ?efficacy_est
#  ?randomization_phase
#  ?maximization_phase

## ----table1, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'--------
tabl1 <- "
|          |             | Scenario    |             |             |
|----------|-------------|-------------|-------------|-------------|
|          |  	         | Doses of A  |  	         |  	         |
|Doses of B|      1	     |      2	     |      3	     |      4	     |
|4	       | (0.20,0.24) | **(0.24,0.35)** | (0.35,0.45) | (0.40,0.50) |
|3	       | (0.12,0.18) | (0.16,0.22) | **(0.22,0.35)** | **(0.25,0.40)** |
|2	       | (0.06,0.10) | (0.10,0.15) | (0.14,0.25) | **(0.20,0.35)** |
|1	       | (0.02,0.05) | (0.04,0.10) | (0.08,0.15) | **(0.12,0.32)** |
"
cat(tabl1)

## ---- eval = TRUE-------------------------------------------------------------
library(crm12Comb)
# generate skeletons
DLT_skeleton <- priorSkeletons(updelta=0.025, target=0.3, npos=10, ndose=16, 
                               model = "empiric", prior = "normal", beta_mean=0)
print(paste0("DLT skeleton is: ", paste(round(DLT_skeleton,3), collapse=", ")))

Efficacy_skeleton <- priorSkeletons(updelta=0.025, target=0.5, npos=10, ndose=16, 
                                    model = "empiric", prior = "normal", beta_mean=0)
print(paste0("Efficacy skeleton is: ", paste(round(Efficacy_skeleton,3), collapse=", ")))

## ---- eval = TRUE-------------------------------------------------------------
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

# simulate 100 trials under the same model and prior distribution
simRes <- SIM_phase_I_II(nsim=100, Nmax=40, DoseComb=scenario, input_doseComb_forMat=c(4,4), 
                         input_type_forMat="matrix", input_Nphase=20,
                         input_DLT_skeleton=DLT_skeleton, 
                         input_efficacy_skeleton=Efficacy_skeleton,
                         input_DLT_thresh=0.3, input_efficacy_thresh=0.3,
                         input_cohortsize=1, input_corr=0,
                         input_early_stopping_safety_thresh=0.33,
                         input_early_stopping_futility_thresh=0.2,
                         input_model="empiric", input_para_prior="normal",
                         input_beta_mean=0, input_beta_sd=sqrt(1.34),
                         input_theta_mean=0, input_theta_sd=sqrt(1.34))
simRes[1:9]

## ---- eval = TRUE, fig.width=10, fig.height=5---------------------------------
# generate plots of patient enrollment of the first trial
enroll_patient_plot(simRes$datALL[[1]])

## ---- eval = TRUE, fig.width=8, fig.height=6----------------------------------
# generate plots of patient allocations by dose levels of the first trial
patient_allocation_plot(simRes$datALL[[1]])

## ---- eval = TRUE, fig.width=8, fig.height=6----------------------------------
# generate plots of ODC selections among all trials
ODC_plot(simRes)

## ----table2, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'--------
tabl2 <- "
|          |          |         | True  (toxicity, efficacy)  probabilities   ||
|----------|----------|------------------|------------------|------------------|
|          | Doses of |                  |     Doses of B   |                  |
| Scenario | A        |      1           |      2           |      3           |
| 1        | 3        | (0.08, 0.15)     | (0.10, 0.20)     | **(0.18, 0.40)** |
|          | 2        | (0.04, 0.10)     | (0.06, 0.16)     | (0.08, 0.20)     |
|          | 1        | (0.02, 0.05)     | (0.04, 0.10)     | (0.06, 0.15)     |
| 2        | 3        | (0.16, 0.20)     | **(0.25, 0.35)** | (0.35, 0.50)     |
|          | 2        | (0.10, 0.10)     | (0.14, 0.25)     | **(0.20, 0.40)** |
|          | 1        | (0.06, 0.05)     | (0.08, 0.10)     | (0.12, 0.20)     |
| 3        | 3        | **(0.24, 0.40)** | (0.33, 0.50)     | (0.40, 0.60)     |
|          | 2        | (0.16, 0.20)     | **(0.22, 0.40)** | (0.35, 0.50)     |
|          | 1        | (0.08, 0.10)     | (0.14, 0.25)     | **(0.20, 0.35)** |
| 4        | 3        | (0.33, 0.50)     | (0.40, 0.60)     | (0.55, 0.70)     |
|          | 2        | **(0.18, 0.35)** | **(0.25, 0.45)** | (0.42, 0.55)     |
|          | 1        | (0.12, 0.20)     | **(0.20, 0.40)** | (0.35, 0.50)     |
| 5        | 3        | (0.45, 0.55)     | (0.55, 0.65)     | (0.75, 0.75)     |
|          | 2        | **(0.20, 0.36)** | (0.35, 0.49)     | (0.40, 0.62)     |
|          | 1        | (0.15, 0.20)     | **(0.20, 0.35)** | **(0.25, 0.50)** |
| 6        | 3        | (0.65, 0.60)     | (0.80, 0.65)     | (0.85, 0.70)     |
|          | 2        | (0.55, 0.55)     | (0.70, 0.60)     | (0.75, 0.65)     |
|          | 1        | (0.50, 0.50)     | (0.55, 0.55)     | (0.65, 0.60)     |
"
cat(tabl2)

## ---- eval = FALSE------------------------------------------------------------
#  scenario1 <- matrix(c(0.02, 0.05,
#                        0.04, 0.10,
#                        0.06, 0.15,
#                        0.04, 0.10,
#                        0.06, 0.16,
#                        0.08, 0.20,
#                        0.08, 0.15,
#                        0.10, 0.20,
#                        0.18, 0.40), ncol=2, byrow = TRUE)

## ---- eval = FALSE------------------------------------------------------------
#  orderings <- function(DLT1, DLT2, ORR1, ORR2){
#    input_Nphase <- c(10, 20, 30)
#    input_corr <- c(0, -2.049, 0.814)
#    input_N <- c(40, 50, 60)
#    DLTs <- list(DLT1, DLT2)
#    ORRs <- list(ORR1, ORR2)
#  
#    conds <- list()
#    i <- 1
#    for (s in 1:2){
#      for (n in 1:3){
#        for (np in 1:3){
#          for (c in 1:3){
#            conds[[i]] <- list(DLT=DLTs[[s]], ORR=ORRs[[s]],
#                               N=input_N[n], Nphase=input_Nphase[np], corr=input_corr[c])
#            i <- i+1
#          }
#        }
#      }
#    }
#  
#    return(conds)
#  }

## ---- eval = FALSE------------------------------------------------------------
#  output <- data.frame(setting = character(), safe = double(),
#                       target = double(), toxic = double(),
#                       avgSS = double(), prop_ODC = double(),
#                       stop_safety = double(), stop_futility = double(),
#                       o_DLT = double(), o_ORR = double())
#  
#  colnames(output) <- c("Design Settings",
#  "Probability of recommending safe/ineffective combinations as ODC",
#  "Probability of recommending target combinations as ODC",
#  "Probability of recommending overly toxic combinations as ODC",
#  "Mean # of patients enrolled", "Proportion of patients allocated to true ODC(s)",
#  "Proportion stopped for safety", "Proportion stopped for futility",
#  "Observed DLT rate", "Observed response rate")

## ---- eval = FALSE------------------------------------------------------------
#  # empiric, normal prior
#  DLT_skeleton1 <- priorSkeletons(updelta=0.045, target=0.3, npos=5, ndose=9,
#                                  model = "empiric", prior = "normal", beta_mean=0)
#  DLT_skeleton2 <- priorSkeletons(updelta=0.06, target=0.3, npos=4, ndose=9,
#                                  model = "empiric", prior = "normal", beta_mean=0)
#  Efficacy_skeleton1 <- priorSkeletons(updelta=0.045, target=0.5, npos=5, ndose=9,
#                                       model = "empiric", prior = "normal", beta_mean=0)
#  Efficacy_skeleton2 <- priorSkeletons(updelta=0.06, target=0.5, npos=4, ndose=9,
#                                       model = "empiric", prior = "normal", beta_mean=0)
#  
#  conds <- orderings(DLT1=DLT_skeleton1, DLT2=DLT_skeleton2,
#                    ORR1=Efficacy_skeleton1, ORR2=Efficacy_skeleton2)

## ---- eval = FALSE------------------------------------------------------------
#  for (s in 1:length(SC)){
#    for (c in 1:length(conds)){
#      curr <- SIM_phase_I_II(nsim=1000, Nmax=conds[[c]]$N, DoseComb=SC[[s]],
#                             input_doseComb_forMat=c(3,3),
#                             input_type_forMat="matrix",
#                             input_Nphase=conds[[c]]$Nphase,
#                             input_DLT_skeleton=conds[[c]]$DLT,
#                             input_efficacy_skeleton=conds[[c]]$ORR,
#                             input_DLT_thresh=0.3, input_efficacy_thresh=0.3,
#                             input_cohortsize=1, input_corr=conds[[c]]$corr,
#                             input_early_stopping_safety_thresh=0.33,
#                             input_early_stopping_futility_thresh=0.2,
#                             input_model="empiric", input_para_prior="normal",
#                             input_beta_mean=0, input_beta_sd=sqrt(1.34),
#                             input_theta_mean=0, input_theta_sd=sqrt(1.34))
#      currTmp <- data.frame(currname, curr$prob_safe, curr$prob_target, curr$prob_toxic, curr$mean_SS, curr$mean_ODC,
#                            curr$prob_stop_safety, curr$prob_stop_futility, curr$mean_DLT, curr$mean_ORR)
#      currTmp$N <- conds[[c]]$N
#      currTmp$Skeleton <- ifelse(c<=27, 1, 2)
#      currTmp$Nphase <- conds[[c]]$Nphase
#      currTmp$corr <- conds[[c]]$corr
#      output <- rbind(output, currTmp)
#    }
#  }

## ---- eval = TRUE, fig.width=8, fig.height=6----------------------------------
sample_plot(dat, outcome = "curr.prob_target",
            outname = "Probability of ODC as target combinations",
            N = NULL, nR = 20, Skeleton = 1, corr = 0)

## ---- eval = TRUE, fig.width=8, fig.height=6----------------------------------
sample_plot(dat, outcome = "curr.mean_ODC",
            outname = "Proportion of patients allocated to true ODC(s)",
            N = 60, nR = 30, Skeleton = 1, corr = NULL)

