#' Title: Get 6 complete orderings for toxicity and efficacy for drug combinations
#'
#' Step 1: convert dose combinations into matrix
#' @param doseComb Either a numeric matrix or a numeric vector or a numeric list.
#' @param type A string.
#' condition 1: detailed dose combinations are directly given in matrix -> type = "comb"
#'                                        dose combinations
#'             dose A
#'             dose B
#' condition 2: directly given the number of levels for the two doses in vector -> type = "matrix"
#'
#' Step 2: get the 6 orderings
#' @return outMat -> Either a list or a matrix.
#' Note: each array refers to the index of drug combinations orderings)
#' @export
#'
#' @references
#' Wages NA, Conaway MR. Specifications of a continual reassessment method design for phase I trials of combined drugs. Pharmaceutical statistics. 2013 Jul;12(4):217-24.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3771354/}

doseComb_to_mat <- function(doseComb, type){
  if (type == "comb"){
    doseA <- dplyr::dense_rank(doseComb[1,])
    doseB <- dplyr::dense_rank(doseComb[2,])
    lA <- length(unique(doseA))
    lB <- length(unique(doseB))

    if (lA > lB){
      tempD <- doseA
      doseA <- doseB
      doseB <- tempD

      tempL <- lA
      lA <- lB
      lB <- tempL
    }

    outMat <- matrix(NA, nrow=lA, ncol=lB)
    for (i in 1:length(doseA)){
      outMat[doseA[i], doseB[i]] <- 1
    }

    n <- 1
    for (i in 1:lA){
      for (j in 1:lB){
        if (!is.na(outMat[i,j])){
          outMat[i,j] <- n
          n <- n+1
        }
      }
    }
    return(outMat)
  } else if (type == "matrix"){
    lA <- min(doseComb)
    lB <- max(doseComb)
    outMat <- matrix(1:(lA*lB), nrow=lA, ncol=lB, byrow=TRUE)
    return(outMat)
  } else if (type == "self"){
    return(doseComb)
  } else {
    stop("invalid name of type")
  }
}


#' Title
#'
#' @param doseComb_forMat Either a numeric matrix or a numeric vector or a numeric list.
#' @param type_forMat A string
#'
#' @return lst_out -> a list of vectors.
#' @export

get_ordering <- function(doseComb_forMat, type_forMat){
  orderMat <- doseComb_to_mat(doseComb = doseComb_forMat, type = type_forMat)

  if (type_forMat == "self"){
    return(orderMat)
  }
  
  #            B
  #       1    2    3
  #   1  d11  d12  d13
  # A 2  d21  d22  d23
  #   3  d31  d32  d33

  lst_out <- list()

  # 1. By rows
  out1 <- c(t(orderMat))
  lst_out[[1]] <- out1[!is.na(out1)]

  # 2. By columns
  out2 <- c(orderMat)
  lst_out[[2]] <- out2[!is.na(out2)]

  # 3. Up diagonals
  nr <- nrow(orderMat)
  nc <- ncol(orderMat)

  res_list <- vector(mode="list",length=nr+nc-1)
  for(r in 1:nr){
    for(c in 1:nc){
      sum_diag <- r+c-1
      res_list[[sum_diag]] <- c(res_list[[sum_diag]],list(c(r,c)))
    }
  }

  out3 <- NULL
  for(i in 1:(nr+nc-1)){
    ll <- length(res_list[[i]])
    if(ll>=2){
      res_list[[i]] <- res_list[[i]][order(sapply(res_list[[i]],'[[',2))]
    }
    for (j in 1:ll){
      out3 <- c(out3, orderMat[res_list[[i]][[j]][1], res_list[[i]][[j]][2]])
    }
  }
  lst_out[[3]] <- out3[!is.na(out3)]

  # 4. Down diagonals
  out4 <- NULL
  for(i in 1:(nr+nc-1)){
    ll <- length(res_list[[i]])
    if(ll>=2){
      res_list[[i]] <- res_list[[i]][order(sapply(res_list[[i]],'[[',1))]
    }
    for (j in 1:ll){
      out4 <- c(out4, orderMat[res_list[[i]][[j]][1], res_list[[i]][[j]][2]])
    }
  }
  lst_out[[4]] <- out4[!is.na(out4)]

  # 5. Alternating down-up diagonals
  out5 <- NULL
  for(i in 1:(nr+nc-1)){
    ind <- ifelse(i%%2==0, 2, 1)
    ll <- length(res_list[[i]])
    if(ll>=2){
      res_list[[i]] <- res_list[[i]][order(sapply(res_list[[i]],'[[',ind))]
    }
    for (j in 1:ll){
      out5 <- c(out5, orderMat[res_list[[i]][[j]][1], res_list[[i]][[j]][2]])
    }
  }
  lst_out[[5]] <- out5[!is.na(out5)]

  # 6. Alternating up-down diagonals
  out6 <- NULL
  for(i in 1:(nr+nc-1)){
    ind <- ifelse(i%%2==0, 1, 2)
    ll <- length(res_list[[i]])
    if(ll>=2){
      res_list[[i]] <- res_list[[i]][order(sapply(res_list[[i]],'[[',ind))]
    }
    for (j in 1:ll){
      out6 <- c(out6, orderMat[res_list[[i]][[j]][1], res_list[[i]][[j]][2]])
    }
  }
  lst_out[[6]] <- out6[!is.na(out6)]

  ### remove duplicates
  lst_out <- unique(lst_out)

  return(lst_out)
}

