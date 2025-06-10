#' @title Compute the Generalized Fréchet Distance Between Two Trajectories
#' 
#' @description Calculates the discrete Fréchet distance between two trajectories, which is 
#' used as the distance metric in clustering algorithms for longitudinal data.
#' 
#' @param traj1 A numeric matrix or data.frame representing the first trajectory.
#' The first column must be time points, and the remaining columns should be one
#' or more variables observed at each time point (e.g., Variable1, Variable2, ...).
#' Each row corresponds to a single time point.
#' @param traj2 A numeric matrix or data.frame representing the second trajectory.
#' The format should be the same as for \code{traj1}, where the first column is
#' time, and the subsequent columns are variables.
#' @param form A character string specifying the return format. 
#' Should be either `"scalar"` to return the scalar Fréchet distance,
#' or `"matrix"` to return the full dynamic programming matrix.
#' 
#' @return A numeric value or matrix.
#' If `form = "scalar"`, returns the Fréchet distance between the two trajectories as a single numeric value.
#' If `form = "matrix"`, returns the dynamic programming matrix used to compute the distance.
#' 
#' @details
#' This function is primarily used internally by clustering functions to evaluate 
#' the similarity between trajectories based on the Fréchet distance. 
#' It is used in the \code{mfkml} function and for generating the distance array
#' used in \code{SFclust} function.
#' 
#' @importFrom proxy dist
#' 
#' @examples
#' # Example trajectories with 3 variables
#' traj1 <- data.frame(
#'   Time = 1:4,
#'   Variable1 = c(1.2, 1.4, 1.6, 1.8),
#'   Variable2 = c(2.3, 2.1, 2.0, 1.9),
#'   Variable3 = c(3.1, 3.3, 3.5, 3.7)
#' ) 
#' traj2 <- data.frame(
#'   Time = 1:3,
#'   Variable1 = c(2.0, 2.2, 2.4),
#'   Variable2 = c(3.0, 2.9, 2.8),
#'   Variable3 = c(1.0, 1.1, 1.2)
#' )
#' 
#' # Compute Fréchet distance (scalar output)
#' fredist(traj1, traj2, form = "scalar")
#'
#' # Compute Fréchet distance matrix
#' fredist(traj1, traj2, form = "matrix")
#' 
#' @export
fredist <- function(traj1, traj2, form) {

  if (form != "scalar" & form != "matrix") {
    stop ("The 'form' argument must have 'scalar' or 'matrix' only\n")
  }
  
  m <- nrow(traj1)
  n <- nrow(traj2)
  
  Dmat <- proxy::dist(traj1, traj2)
  
  Fmat <- matrix(nrow = m, ncol = n)
  Fmat[1, 1] <- Dmat[1, 1]
  
  for (j in 2:n) {
    Fmat[1, j] <- max(Fmat[1, j - 1], Dmat[1, j])
  }
  
  for (i in 2:m) {
    Fmat[i, 1] <- max(Fmat[i - 1, 1], Dmat[i, 1])
    for (j in 2:n) {
      Fmat[i, j] <- max(min(Fmat[i - 1, j - 1],
                            Fmat[i - 1, j],
                            Fmat[i, j - 1]),
                        Dmat[i, j])
    }
  }
  
  return(if (form == "scalar") Fmat[m, n] else Fmat)
}