#' Setting up initial candidate models\cr
#'
#' This function automatically ranks all candidate interaction effects under
#' Strong, Weak, or No heredity condition and obtains initial candidate models.
#'
#' @details
#' \itemize{
#' \item For inputs \code{pi1}, \code{pi2}, and \code{pi3}, the number needs to
#' satisfy the condition: \eqn{\pi_1+\pi_2+\pi_3=1-\pi_0} where \eqn{\pi_0}
#' is a numeric value between 0 and 1, the smaller the better.
#' \item For input \code{lambda}, the number needs to satisfy the condition:
#' \eqn{\lambda\geq 5.1/log(2)}.
#' \item For input \code{q}, this number is fixed throughout each generation.
#' }
#'
#' @param X Input data. An optional data frame, or numeric matrix of dimension
#' \code{n} by \code{nmain.p}. Note that the two-way interaction effects should not
#' be included in \code{X} because this function automatically generates the
#' corresponding two-way interaction effects if needed.
#' @param y Response variable. A \code{n}-dimensional vector, where \code{n} is the number
#' of observations in \code{X}.
#' @param heredity whether to enforce Strong, Weak, or No heredity.
#' Default is "Strong".
#' @param nmain.p A numeric value that represents the total number of main effects
#' in \code{X}.
#' @param r1 A numeric value indicating the maximum number of main effects.
#' @param r2 A numeric value indicating the maximum number of interaction effects.
#' @param interaction.ind A two-column numeric matrix containing all possible
#' two-way interaction effects. It must be generated outside of this function
#' using \code{t(utils::combn())}. See Example section for details.
#' @param pi1 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to the Details section.
#' @param pi2 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to the Details section.
#' @param pi3 A numeric value between 0 and 1, defined by users. Default is 0.32.
#' For guidance on selecting an appropriate value, please refer to the Details section.
#' @param lambda A numeric value defined by users. Default is 10.
#' For guidance on selecting an appropriate value, please refer to the Details section.
#' @param q A numeric value indicating the number of models in each generation (e.g.,
#' the population size). Default is 40.
#'
#' @return A \code{list} of output. The components are:
#' \item{initialize}{Initial candidate models. A numeric matrix of dimension \code{q} by
#' \code{r1+r2} where each row represents a fitted model. Duplicated models are allowed.}
#' \item{InterRank}{Rank of all candidate interaction effects. A two-column numeric
#' matrix. The first column contains indices of ranked two-way interaction effects, and the
#' second column contains its corresponding ABC score.}
#' \item{mainind.sel}{Selected main effects.  A \code{r1}-dimensional vector.}
#' \item{mainpool}{Ranked main effects in \code{X}.}
#' @export
#' @importFrom utils combn
#' @importFrom stats rnorm
#' @importFrom VariableScreening screenIID
#' @importFrom stats na.exclude
#'
#' @seealso \code{\link{ABC}}, \code{\link{Extract}}.
#' @examples # Under Strong heredity
#' set.seed(0)
#' nmain.p <- 4
#' interaction.ind <- t(combn(4,2))
#' X <- matrix(rnorm(50*4,1,0.1), 50, 4)
#' epl <- rnorm(50,0,0.01)
#' y<- 1+X[,1]+X[,2]+X[,1]*X[,2]+epl
#' p1 <- initial(X, y, nmain.p = 4, r1 = 3, r2 = 3,
#'     interaction.ind = interaction.ind, q = 5)
#' p1
initial <- function(X, y, heredity = "Strong",
                    nmain.p, r1, r2,
                    interaction.ind = NULL,
                    pi1 = 0.32, pi2 = 0.32, pi3 = 0.32,
                    lambda = 10, q = 40){
  if (is.null(interaction.ind)) stop("Interaction.ind is missing.
                                       Use t(combn()) to generate interaction matrix.")
  colnames(X) <- make.names(rep("","X",ncol(X)+1),unique=TRUE)[-1]
  max_model_size <- r1 + r2
  parents <- matrix(0, nrow = q, ncol = max_model_size)

  DCINFO <- VariableScreening::screenIID(X, y, method="DC-SIS")
  mainind <- as.numeric(gsub(".*?([0-9]+).*", "\\1",  colnames(X)[order(DCINFO$rank)]))
  Shattempind <- mainind[1:r1]

  # no heredity pool
  if (heredity == "No"){
    interpooltemp <- interaction.ind
  }

  # weak heredity pool
  if (heredity =="Weak"){
    for (i in 1:q) {
      df <- rbind(interaction.ind[interaction.ind[,1] %in% Shattempind,][order(stats::na.exclude(match(interaction.ind[,1], Shattempind))),],
                  interaction.ind[interaction.ind[,2] %in% Shattempind,][order(stats::na.exclude(match(interaction.ind[,2], Shattempind))),])
      interpooltemp <- df[!duplicated(df),]
    }
  }

  # strong heredity pool
  if (heredity =="Strong"){
    interpooltemp <- t(utils::combn(sort(Shattempind),2))
  }

  intercandidates.ind <- match(do.call(paste, as.data.frame(interpooltemp)), do.call(paste, as.data.frame(interaction.ind)))+nmain.p

  interscoreind <- list()
  for (i in 1:length(intercandidates.ind)) {
    interscoreind[[i]] <- ABC(X, y, heredity = heredity, nmain.p = nmain.p,
                              extract = "Yes", varind = c(Shattempind,intercandidates.ind[i]),
                              interaction.ind = interaction.ind,
                              pi1 = pi1, pi2 = pi2, pi3= pi3, lambda = lambda)
  }
  interscoreind <- interscoreind
  MA <- as.matrix(cbind(inter = intercandidates.ind, scores = interscoreind))
  if (dim(MA)[1] == 1){
    MB <- MA
  }else{
    MB <- MA[order(unlist(MA[,2]), na.last = TRUE),]
  }

  interind <- unlist(MB[,1])
  for (i in 1:dim(parents)[1]) {
    parents[i,c(1:length(Shattempind))] <- Shattempind
    if (length(unlist(MB[,1])) < dim(parents)[1]){
      parents[i,max(which(!parents[i,]==0))+1] <- rep_len(interind, length.out=dim(parents)[1])[i]
    }else{
      parents[i,max(which(!parents[i,]==0))+1] <- interind[i]
    }
  }
  parents <- list(
    initialize = parents,
    InterRank = MB,
    mainind.sel = Shattempind,
    mainpool = mainind
  )
  return(parents)
}
