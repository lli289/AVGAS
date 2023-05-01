#' Extracting specific columns from a data\cr
#'
#' @param X
#' @param varind
#' @param interaction.ind
#'
#' @return
#' @export
#'
#' @examples
Extract <- function(X, varind, interaction.ind){
  if (as.logical(any(duplicated(varind[which(varind!=0)])))){
    stop("There cannot be duplicated values in varind.")
  }
  ncoln <- length(varind)
  nrown <- nrow(X)
  nmain.p <- ncol(X)
  mainind <- varind[which(varind%in%1:nmain.p)]
  mainvars.matrix <- X[, mainind]
  interind <- varind[which(varind>nmain.p)]

  a <- interaction.ind[interind-nmain.p,]
  if (length(a) >1){
    intervars.matrix <- matrix(0, nrow = nrown,ncol = mydim(a)[1])
    for (i in 1:mydim(a)[1]) {
      if (mydim(a)[1]==1) intervars.matrix[,i] <- X[, a[1]]*X[,a[2]]
      else {
        intervars.matrix[,i] <- X[, a[i,1]]*X[,a[i,2]]
      }
    }
    colnames(intervars.matrix) <- paste0("X.", interind)
    data_extract <- cbind(mainvars.matrix, intervars.matrix)
  }
  else{
    data_extract <- mainvars.matrix
  }
  if (length(mainind) == 1){
    data_extract <- as.matrix(data_extract)
    colnames(data_extract)[which(colnames(data_extract)=="mainvars.matrix" )] <- paste0("X.", mainind)
  }
  if (length(mainind) == 1 && length(interind) ==0){
    colnames(data_extract) <- paste0("X.", mainind)
  }
  if (length(mainind) == 0){
    data_extract <- data_extract
  }
  if (length(mainind)>1){
    colnames(data_extract)[1:dim(mainvars.matrix)[2]] <- paste0("X.", mainind)
  }
  return(data_extract)
}

mydim <- function(x){
  if (is.matrix(x)) dim(x)
  else return(1)
}

mychoose <- function(k1){
  if (k1 == 1){
    return(1)
  }else{
    return(choose(k1,2))
  }
}

mysample <- function(x, size, replace = F, prob = NULL){
  if (length(x) == 1) return(x)
  if (length(x) > 1) return(sample(x, size, replace, prob))
}

sort_zeros <- function(vec){
  non_zeros <- vec[vec != 0]
  zeros <- vec[vec == 0]
  if (length(zeros) > 0){
    return(c(sort(non_zeros), zeros))
  }else{
    return(c(sort(non_zeros)))
  }
}

unique_rows <- function(matrix) {
  unique_matrix <- matrix[!duplicated(matrix[,ncol(matrix)]), ]
  return(unique_matrix)
}

predictor_match <- function(candidate.model, r1,r2, nmain.p, interaction.ind){
  model.match <- list()
  if (mydim(candidate.model)[1]==1){
    ee <- candidate.model[1:(r1+r2)]
    ee1 <- ee[which(ee%in% 1:nmain.p)]
    ee2 <- ee[which(ee > nmain.p)]
    if (length(ee2) == 0){
      model.match[[1]] <-c(paste0("X.", ee1))
    }
    if (length(ee1) == 0){
      model.match[[1]] <- c( paste0("X.", interaction.ind[ee2-nmain.p,1], "X.", interaction.ind[ee2-nmain.p,2]))
    }
    if (!length(ee2) ==0 & !length(ee1) ==0){
      model.match[[1]] <-c(paste0("X.", ee1), paste0("X.", interaction.ind[ee2-nmain.p,1], "X.", interaction.ind[ee2-nmain.p,2]))
    }
  }else{
    for (i in 1:nrow(candidate.model)) {
      ee <- candidate.model[i, 1:(r1+r2)]
      ee1 <- ee[which(ee%in% 1:nmain.p)]
      ee2 <- ee[which(ee > nmain.p)]
      if (length(ee2) == 0){
        model.match[[i]] <-c(paste0("X.", ee1))
      }
      if (length(ee1) == 0){
        model.match[[i]] <- c( paste0("X.", interaction.ind[ee2-nmain.p,1], "X.", interaction.ind[ee2-nmain.p,2]))
      }
      if (!length(ee2) ==0 & !length(ee1) ==0){
        model.match[[i]] <-c(paste0("X.", ee1), paste0("X.", interaction.ind[ee2-nmain.p,1], "X.", interaction.ind[ee2-nmain.p,2]))
      }
    }
  }
  return(model.match = model.match)
}
