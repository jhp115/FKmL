# pairfmean: internal function, used only within mfkml() function
# Computes the pairwise weighted Fréchet mean of two trajectories.
# This function is used in the initial stage of clustering, where trajectory pairs are matched 
# (initial league matches) to compute weighted Fréchet means that serve as initial cluster centers.
# Different weights can be assigned to each trajectory by providing a data.frame 
# with ID and Weight columns to the mfkml() function. The pairfmean() function internally uses 
# this weighted data to compute pairwise Fréchet means.

# Inputs:
# - grp_dt: A long-format data.frame with columns:
#           ID (subject ID), Weight (observation weight), Time, and one or more Variable columns.
#           This data.frame is generated inside the mfkml() function.
#           If the input data.frame does not include a Weight column, 
#           all trajectories are treated as having equal weights. 
#           In that case, the mfkml() function internally adds a Weight column filled with 1s.
#           If different weights should be assigned to each trajectory,
#           the input must include a Weight column with the specified values.
#
# - nmat: An integer specifying the number of pairwise matches in the first round
#         (e.g., \code{floor(length(unique(ID)) / 2)}). 
#         If the number of trajectories is even, all are included in pairwise matches. 
#         If it is odd, the unmatched trajectory advances to the next round without pairing.
#
# Output:
# - A data.frame containing the Fréchet mean trajectories from each pair in the initial round.
#   If the number of IDs is odd, the remaining unmatched trajectory is retained and added at the end.

pairfmean <- function(grp_dt, nmat) {

  Mtraj <- grp_dt[F, ]
  idl <- unique(grp_dt[, 1])
  
  # Pairwise the weighted Fréchet means
  if (length(idl) >= 2) {
    for (idx in 1:nmat) {
      traj1 <- grp_dt[grp_dt[, 1] == idl[2 * idx - 1], -c(1:2)]
      traj2 <- grp_dt[grp_dt[, 1] == idl[2 * idx], -c(1:2)]
      weight <- c(grp_dt[grp_dt[, 1] == idl[2 * idx - 1], 2][1],
                  grp_dt[grp_dt[, 1] == idl[2 * idx], 2][1])
      m <- nrow(traj1)
      n <- nrow(traj2)
      
      # Get coupling sequence between two trajectories
      i <- m + 1
      j <- n + 1
      Fmat <- matrix(Inf, nrow = i, ncol = j)
      Fmat[-1, -1] <- fredist(traj1, traj2, form = 'matrix')
      csq <- matrix(c(m, n), ncol = 2, byrow = T)
      
      # Step down Fmat gradient
      while (i >= 2 | j >= 2) {
        if ((i - 3) * (n - 1) >= (j - 2) * (m - 1)) {
          switch (which.min(c(Fmat[i - 1, j], Fmat[i - 1, j - 1], Fmat[i, j - 1])),
                  {csq <- rbind(csq, c(i - 2, j - 1)); i <- i - 1},
                  {csq <- rbind(csq, c(i - 2, j - 2)); i <- i - 1; j <- j - 1},
                  {csq <- rbind(csq, c(i - 1, j - 2)); j <- j - 1})
        } else if ((j - 3) * (m - 1) >= (i - 2) * (n - 1)) {
          switch (which.min(c(Fmat[i, j - 1], Fmat[i - 1, j - 1], Fmat[i - 1, j])),
                  {csq <- rbind(csq, c(i - 1, j - 2)); j <- j - 1},
                  {csq <- rbind(csq, c(i - 2, j - 2)); i <- i - 1; j <- j - 1},
                  {csq <- rbind(csq, c(i - 2, j - 1)); i <- i - 1})
        } else {
          switch (which.min(c(Fmat[i - 1, j - 1], Fmat[i - 1, j], Fmat[i, j - 1])),
                  {csq <- rbind(csq, c(i - 2, j - 2)); i <- i - 1; j <- j - 1},
                  {csq <- rbind(csq, c(i - 2, j - 1)); i <- i - 1},
                  {csq <- rbind(csq, c(i - 1, j - 2)); j <- j - 1})
        }
      }
      
      csq <- csq[order(csq[, 1], csq[, 2]),][-1, ]
      
      # Calculate the weighted Fréchet mean with weights
      mtraj <- (weight[1] * traj1[csq[, 1], ] + weight[2] * traj2[csq[, 2], ]) / sum(weight)
      Mt <- data.frame(cbind(idl[2 * idx - 1], sum(weight), mtraj))
      names(Mt) <- names(Mtraj)
      Mtraj <- rbind(Mtraj, Mt)
    }
  } else {  
    Mtraj <- grp_dt # If there is only 1 trajectory in 'grp_dt', return the trajectory intactly
  }
  
  # If it was playoff then mean trajectories bind with original trajectories
  if (nmat * 2 < length(idl)) {
    Mtraj <- rbind(Mtraj, grp_dt[grp_dt[, 1] %in% idl[(nmat * 2 + 1):length(idl)], ])
  }
  
  return(Mtraj)
}

# fremean: internal function, used only within mfkml()
# Calculates the Fréchet mean for a set of trajectories.
# This function is used iteratively during the clustering process to update 
# cluster centers based on the current assignment of trajectories to clusters.

# Inputs:
# - grp_dt: A long-format data.frame with columns:
#           ID (subject ID), Weight (observation weight), Time, and one or more Variable columns.
#           This data.frame is generated inside the mfkml() function.
#           If the input data.frame does not include a Weight column, 
#           all trajectories are treated as having equal weights. 
#           In that case, the mfkml() function internally adds a Weight column filled with 1s.
#           If different weights should be assigned to each trajectory,
#           the input must include a Weight column with the specified values.
#
# Output:
# - A data.frame representing the Fréchet mean of all input trajectories.
#   The structure is the same as the input data.frame (excluding the ID and Weight columns),
#   where each row corresponds to a time point and contains the averaged variable values.
#   The final result is obtained through iterative pairwise Fréchet mean calculations
#   using a tournament-style process (playoffs and rounds) until only one representative trajectory remains.

fremean <- function(grp_dt) {
  
  nplayer <- length(unique(grp_dt[, 1]))
  
  if (nplayer >= 2) {
    
    # Set tournament to get the Fréchet mean
    # Shuffle samples maintaining long format
    ord <- sample(x = 1:nplayer, size = nplayer)
    grp_dt[, 1] <- factor(grp_dt[, 1])
    grp_dt[, 1] <- factor(grp_dt[, 1], levels = levels(grp_dt[, 1])[ord])
    reord_dt <- grp_dt[order(grp_dt[, 1]), ]
    reord_dt[, 1] <- as.integer(reord_dt[, 1])
    
    ## Matching trajectories on tournament ##
    idl <- unique(reord_dt[, 1])
    Mtraj <- reord_dt
    tourdeep <- floor(log2(nplayer))
    playoff <- nplayer - 2^tourdeep
    
    if (playoff > 0) {
      Mtraj <- pairfmean(grp_dt = reord_dt, nmat = playoff)
    }
    
    for (i in tourdeep:1) {
      nmat <- 2^(i - 1)
      Mtraj <- pairfmean(grp_dt = Mtraj, nmat = nmat)
    }
  } else {  # If there is only 1 trajectory in 'grp_dt', return the trajectory intactly
    Mtraj <- grp_dt
  }
  
  return(Mtraj[-c(1:2)])
}


#' @title Multidimensional Fréchet Distance-Based K-means for Longitudinal Data
#' 
#' @description Extends \code{kmlShape} to multidimensional (p ≥ 2) longitudinal data.
#' It performs scale adjustment and trajectory alignment across all variables prior to clustering
#' to reduce distortions caused by differences in time grids and amplitude scales. 
#' When variables exhibit substantially different ranges, standardization is required to prevent any single variable 
#' from disproportionately influencing the clustering outcome.
#' 
#' The clustering process follows an iterative K-means framework, where cluster assignments are updated based on Fréchet distances.
#' Cluster centers are computed using the weighted Fréchet mean, which accounts for variable weights assigned to individual trajectories.
#' This allows the mean to be adjusted according to the relative importance of each trajectory in the clustering process.
#' 
#' @param dt A long-format data.frame containing the following columns in the specified order:
#'   - \code{ID}: An identifier for each trajectory.
#'   - \code{Time}: The time points at which measurements were recorded (numeric or integer vector).
#'   - \code{Variable1}, \code{Variable2}, ... : The measured variables over time (numeric values).
#'   The data.frame should not include any missing values.
#'   See 'Details' for structure requirements.
#' @param clt_n An integer specifying the number of clusters. 
#' The number of unique trajectories must be greater than or equal to \code{clt_n}.
#' @param scales A numeric vector used for scaling the time and variable columns. The length of \code{scales} must be equal to \code{ncol(dt) - 1}, 
#'   where each value in \code{scales} corresponds to the scaling factor for the respective column (excluding the ID column).
#'   See 'Details' for structure requirements.
#' @param weight Specifies the weights used for calculating the weighted Fréchet mean. It can take one of the following forms:
#'   - A data.frame with two columns: \code{ID} and \code{Weight}, where each \code{Weight} value indicates the importance of the corresponding trajectory.
#'   - A numeric value of 1, indicating equal weights for all trajectories.
#'   See 'Details' for structure requirements.
#' @param maxIter The maximum number of iterations allowed before stopping if convergence is not reached. The default value is 50.
#' 
#' @details
#' The input dataset (\code{dt}) must contain only numeric values (except for the ID column) 
#' and must not include any missing values. 
#' Each variable should be measured at least three times per trajectory, 
#' since the method relies on trajectory shapes. 
#' Two observations per trajectory are insufficient to capture shape trends (e.g., increasing, decreasing, or stable).
#'
#' Because the Fréchet distance is sensitive to measurement units, proper scaling is essential when applying the \code{mfkml} function.
#' The \code{scales} vector contains scaling factors for time and each variable, 
#' which are used to rescale the corresponding columns.
#' This scaling prevents distortion due to differences in the units of time and variables,
#' allowing for more accurate shape-based comparisons.
#'
#' This function involves random sampling internally. 
#' For reproducible results, set the random seed before calling the function using \code{set.seed()}.
#'
#'@return A list with the following components:
#' \describe{
#'   \item{\code{Cluster}}{A data.frame containing the \code{ID} and \code{Cluster} columns, which indicate the final cluster assignment for each trajectory.}
#'   \item{\code{Center}}{A data.frame representing the final cluster centers, with columns for the cluster IDs, time points, and variable values.}
#'   \item{\code{Iteration}}{The number of iterations the algorithm performed before reaching convergence.}
#' }
#'
#' @export

mfkml <- function(dt, clt_n, scales, weight, maxIter = 50) {
  
  # Check input validity
  if (!is.data.frame(dt) | !all(sapply(dt[-1], is.numeric)) | any(is.na(dt))) {
    stop("The 'dt' argument must be a data.frame (ID, Time, Variable1, Variable2, ...)
         composed by all numeric columns except for ID
         and have complete cases but their measurement time no needs to be matched\n")
  }
  
  # Check number of measurements per trajectory
  tst_dt <- dplyr::group_by_at(dt, 1)
  tst_dt <- dplyr::mutate(tst_dt, fkml_Tms = seq(dplyr::n()))
  tst_dt <- dplyr::mutate(tst_dt, fkml_Max = max(fkml_Tms))
  tst_dt <- dplyr::select(dplyr::ungroup(tst_dt), fkml_Max)
  
  if (min(tst_dt) < 3) {
    stop("There are cases measured only once or twice in 'dt'.
         The measurement times must be more than three.\n")
  }
  
  # Check number of unique trajectories
  dt <- as.data.frame(dt)
  id <- unique(dt[, 1])
  
  if (length(id) < clt_n) {
    stop("Trajectories in 'dt' are less than clusters\n")
  }
  
  # Apply scaling to time and variable columns
  dt[-1] <- t(t(dt[-1]) * scales)
  
  # Merge weights with data
  if (is.data.frame(weight) == TRUE) {
    weight <- weight
  } else {
    if (weight == 1) {
      weight <- data.frame(ID = id, Weight = 1)
    } else if (!is.data.frame(weight) | ncol(weight) != 2 & !all.equal(weight, id)) {
      stop("The 'weight' argument must have '1' or data.frame of weight
         composed by the same sample as unique 'dt' (ID, Weight)\n")
    }
  }
  
  dt <- merge(weight, dt, by.x = 1, by.y = 1, all = TRUE)
  dt <- dt[order(dt[, 1], dt[, 3]), ]           
  
  # Initialization
  iter <- 0
  group <- 0
  exGroup <- 1
  
  # Initial random cluster assignment
  seeds <- sample(id, clt_n)
  centers <- dt[dt[, 1] %in% seeds, -2]
  centers[, 1] <- as.integer(factor(unlist(centers[, 1])))  # Convert initial IDs to cluster labels
  
  # Implement clustering
  options(warn = 2)
  while (!identical(group, exGroup) && iter < maxIter) {
    
    iter <- iter + 1
    exGroup <- group
    
    # Compute Fréchet distances from each trajectory to each cluster center
    matrixDist <- numeric()
    for (cent_i in 1:clt_n) {
      
      # Compute Fréchet distance between an individual trajectory and the cluster center
      fDist <- function(i) {
        fredist(traj1 = dt[dt[, 1] == i, -c(1:2)],
                traj2 = centers[centers[, 1] == cent_i, -1],
                form = 'scalar')
      }
      
      distToMi <- sapply(id, fDist)
      matrixDist <- cbind(matrixDist, distToMi)
      
    }
    
    # Assign each trajectory to the closest cluster center
    group <- as.integer(apply(matrixDist, 1, which.min))  
    names(group) <- id
    groupl <- group[match(dt[, 1], names(group))]
    
    # Compute new cluster centers using the weighted Fréchet mean
    fMean <- function(cent_i) {
      data.frame(Cluster = cent_i,
                 fremean(dt[groupl == cent_i, ]))
    }
    
    centers <- do.call(rbind, lapply(1:clt_n, fMean))
  }
  
  options(warn = 0)
  if (iter == maxIter) {
    warning("The maximum number of iteration has been reached, the algorithm did not converge\n")
  }
  
  # Restore scaled value
  centers[-1] <- t(t(centers[-1]) / scales)  
  
  # Reorder cluster labels by cluster size (largest first)
  reOrder <- rank(-table(group), ties.method = 'first')
  group_re <- reOrder[group]
  group_re <- as.factor(group_re)
  
  centers[, 1] <- reOrder[centers[, 1]]
  centers <- centers[order(centers[, 1], centers[, 2]), ]
  centers$Cluster <- as.factor(centers$Cluster)
  
  return(list(Cluster = data.frame(ID = names(group), Cluster = group_re),
              Center = centers,
              Iteration = iter))
}

