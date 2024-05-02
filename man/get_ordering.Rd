\name{get_ordering}
\alias{get_ordering}
\title{
Complete orderings for combinations of two drugs
}
\description{
This function is to obtain complete orderings of both toxicity and efficacy for all drug combinations when considering the partial ordering issue for two combined drugs.
}
\usage{
get_ordering(doseComb_forMat, type_forMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{doseComb_forMat}{For 2 drugs, either a matrix with columns for all combinations and rows for two drugs when type_forMat is "comb", or a vector of length 2 indicating the number of levels of each drug when type_forMat is "matrix" should be input. For more than 2 drugs, a list of all possible orderings should be input when type_forMat is "self".}
  \item{type_forMat}{A character string to indicate the input type of dose combinations. The type_forMat is either "comb" for inputting all tested dose combinations for 2 drugs, "matrix" for entering number of levels of two drugs so that a matrix of combinations can be constructed, or "self" for more than 3 drugs with directly input of all possible orderings.}
}
\details{
Dose-toxicity and dose-efficacy curves for single drug are assumed as monotonically increasing as drug increases. After combined two drugs and traversing all combinations, the number of possible orderings exponentially expands when increasing dose levels for single drug or the dimensions of drug combinations, which will mask the information provided from orderings. Here, 6 typical complete orderings are utilized from practical design.

Take the \eqn{3 \times 3} matrix as an example so that there are total of \eqn{9} does combination.
  \deqn{B}
  \deqn{\quad \quad 1\quad    2\quad    3}
  \deqn{\quad 1\quad  d_1\quad  d_2\quad  d_3}
  \deqn{A\; 2\quad  d_4\quad  d_5\quad  d_6}
 \deqn{\quad 3\quad  d_7\quad  d_8\quad  d_9}

Across rows: \deqn{\pi_T(d_1)\le \pi_T(d_2)\le \pi_T(d_3)\le \pi_T(d_4)\le \pi_T(d_5)\le \pi_T(d_6) \le \pi_T(d_7)\le \pi_T(d_8) \le \pi_T(d_9)}
Up columns: \deqn{\pi_T(d_1)\le \pi_T(d_4)\le \pi_T(d_7)\le \pi_T(d_2)\le \pi_T(d_5)\le \pi_T(d_8) \le \pi_T(d_3)\le \pi_T(d_6) \le \pi_T(d_9)}
Up diagonals: \deqn{\pi_T(d_1)\le \pi_T(d_2)\le \pi_T(d_4)\le \pi_T(d_3)\le \pi_T(d_5)\le \pi_T(d_7) \le \pi_T(d_6)\le \pi_T(d_8) \le \pi_T(d_9)}
Down diagonals: \deqn{\pi_T(d_1)\le \pi_T(d_4)\le \pi_T(d_2)\le \pi_T(d_7)\le \pi_T(d_5)\le \pi_T(d_3) \le \pi_T(d_8)\le \pi_T(d_6) \le \pi_T(d_9)}
Alternating down-up diagonals: \deqn{\pi_T(d_1)\le \pi_T(d_2)\le \pi_T(d_4)\le \pi_T(d_7)\le \pi_T(d_5)\le \pi_T(d_3) \le \pi_T(d_6)\le \pi_T(d_8) \le \pi_T(d_9)}
Alternating up-down diagonals: \deqn{\pi_T(d_1)\le \pi_T(d_4)\le \pi_T(d_2)\le \pi_T(d_3)\le \pi_T(d_5)\le \pi_T(d_7) \le \pi_T(d_8)\le \pi_T(d_6) \le \pi_T(d_9)}
Finally, obtain the unique orderings if there are any duplicates.
}
\value{
A list of vectors are returned, with each vector representing a specific order of dose combinations in ascending order (dose combinations are denoted by natrual numbers).
}
\references{
Wages, N. A., & Conaway, M. R. (2013). Specifications of a continual reassessment method design for phase I trials of combined drugs. Pharmaceutical statistics, 12(4), 217-224. \doi{10.1002/pst.1575}
}

\examples{
# Input all dose combinations
doseComb <- matrix(c(10, 3, 20, 3, 30, 3, 40, 3, 20, 5, 20, 7), byrow=FALSE, nrow=2)
orderings <- get_ordering(doseComb_forMat = doseComb, type_forMat = "comb")

# Input number of dose levels
doseLevels <- c(2,3)
orderings <- get_ordering(doseComb_forMat = doseLevels, type_forMat = "matrix")

# Input customized orderings
orders <- list(c(1,2,3,4,5), c(2,1,4,3,5), c(1,2,4,3,5), c(2,1,3,4,5))
orderings <- get_ordering(doseComb_forMat = orders, type_forMat = "self")

# Extreme case: only one ordering will be input
orders <- list(c(1,2,3,4,5,6,7))
orderings <- get_ordering(doseComb_forMat = orders, type_forMat = "self")
}
