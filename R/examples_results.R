#' Output dataset for examples given list of inputs
#'
#' The dataset contains 1296 rows (6 scenarios, 4 different combinations of link functions and prior distributions, 
#'  2 sets of skeletons, 3 maximum number of patients, 3 toxicity and efficacy correlations, and 3 subset number 
#'  of patients) that represent each condition by 1000 simulations, along with 15 columns that contain 9 operating 
#'  characteristics, other 6 columns including Scenario, Model for link function and prior distribution, N for 
#'  maximum number of patients, Skeleton, Nphase for subset number of patients, and corr for toxicity and efficacy 
#'  correlation to separate each condition.
#'
#' @format A data frame with 1296 rows and 15 variables:
#' \itemize{
#'   \item Model: a string vector containing 4 models used in link functions. 
#'   \item{Scenario: a numeric vector containing 6 scenarios.
#'   \item Skeleton: a numeric vector containing 2 sets of skeletons used for both toxicity and efficacy.               
#'   \item N: a numeric vector containing 3 maximum number of patients.                     
#'   \item nR: a numeric vector containing 3 subset number of patients. 
#'   \item corr: a numeric vector containing 3 toxicity and efficacy correlations.                    
#'   \item prob_safe: a numeric vector for the probability of simulation trials with ODC identified as safe/ineffective combinations.      
#'   \item prob_target: a numeric vector for the probability of simulation trials with ODC identified as target combinations.        
#'   \item prob_toxic: a numeric vector for the probability of simulation trials with ODC identified as toxic combinations.         
#'   \item mean_SS: a numeric vector for the average number of patients enrolled.   
#'   \item mean_ODC: a numeric vector for the average proportion of patients allocated to target ODS(s).          
#'   \item prob_stop_safety: a numeric vector for the probability of simulation trials stopped early for safety.   
#'   \item prob_stop_futility: a numeric vector for the probability of simulation trials stopped early for futility. 
#'   \item mean_DLT: a numeric vector for the average observed DLT rate.           
#'   \item mean_ORR: a numeric vector for the average observed response rate. 
#' }
"examples_results"