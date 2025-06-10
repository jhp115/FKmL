#' @title Compute Distance Array for Multidimensional Functional Data
#' 
#' @description This function standardizes multidimensional functional data using provided scaling factors,
#' computes pairwise Fréchet distances between trajectories for each variable, and returns a
#' distance array (3-dimensional array of distance matrices).
#' 
#' @param dt A long-format data.frame containing the following columns in the specified order:
#'   - \code{ID}: An identifier for each trajectory.
#'   - \code{Time}: The time points at which measurements were recorded (numeric or integer vector).
#'   - \code{Variable1}, \code{Variable2}, ... : The measured variables over time (numeric values).
#'   The data.frame should not include any missing values.
#'   See 'Details' for structure requirements.
#' @param time_scale A single numeric value used to scale the \code{Time} column. 
#' This ensures that time is appropriately weighted relative to the variables.
#' @param var_scales A numeric vector of scaling factors for the measured variables. 
#' Its length must be equal to \code{ncol(dt) - 2}.
#'
#' @return A numeric value or matrix.
#' If `form = "scalar"`, returns the Fréchet distance between the two trajectories as a single numeric value.
#' If `form = "matrix"`, returns the dynamic programming matrix used to compute the distance.
#' 
#' @details
#' The \code{dist.array} function first applies scaling to the \code{Time} and each measured variable.
#' Then, it computes pairwise Fréchet distances between trajectories for each variable separately.
#' The output is a 3-dimensional array in which each slice corresponds to a variable-specific distance matrix.
#'
#' Unlike the \code{mfkml} function, which requires at least three measurements across time for each trajectory,
#' the SFKmL ((Sparse multi-dimensional Fréchet distance-based K-medoids for Longitudinal data), which uses \code{dist.array}, 
#' allows for trajectories with missing values, as long as each variable has at least three time points for each trajectory. 
#' Therefore, \code{dt} may include missing values.
#'
#' @return A 3-dimensional array of pairwise distances with dimensions \code{[n, n, p]}, where:
#'   \describe{
#'     \item{n}{Number of unique trajectories.}
#'     \item{p}{Number of variables.}
#'   }
#' Each slice \code{[, , k]} is a distance matrix for variable \code{k}.
#'
#' @importFrom abind abind
#' 
#' @export

dist.array <- function(dt, time_scale, var_scales){
  
  # Check input validity
  if (!is.data.frame(dt) | !all(sapply(dt[-1], is.numeric)) | any(is.na(dt[, c(1:2)]))) {
    stop (
      "The 'dt' argument must be a data.frame containing the columns ID, Time, Variable1, Variable2, etc., 
    with all columns being numeric except for ID. 
    Additionally, the ID and Time columns must not have any missing values.\n")
  }
  
  # Check for trajectories where any variable is measured fewer than three times
  uniq_id <- unique(dt[, 1])
  sumnonisna <- function(x) {sum(!is.na(x))}
  
  naid <- c()
  for (i in uniq_id) {
    dd1 <- dt[which(dt[, 1] %in% i), ]
    dd1 <- dd1[, -c(1:2)]
    vv <- apply(dd1, 2, sumnonisna)
    if (sum(vv < 3) > 0) naid <- c(naid, i)
  }
  
  if (length(naid) > 0) 
    stop ("All variables must be measured at least three times.")
  
  # Check number of measurements per trajectory
  tst_dt <- dplyr::group_by_at(dt, 1)
  tst_dt <- dplyr::mutate(tst_dt, fkml_Tms = seq(dplyr::n()))
  tst_dt <- dplyr::mutate(tst_dt, fkml_Max = max(fkml_Tms))
  tst_dt <- dplyr::select(dplyr::ungroup(tst_dt), fkml_Max)
  
  if (min(tst_dt) < 3) {
    stop("There are cases measured only once or twice in 'dt'.
         The measurement times must be more than three\n")
  }
  
  # Check scaling parameter lengths
  num.var <- dim(dt)[2] - 2     # the number of variables
  
  scales <- c(time_scale, var_scales)
  
  if (num.var + 1 != length(scales)) {
    stop("the number of scale parameters does not match with the number of variables")
  }
  
  # Apply scaling to time and variable columns
  dt[-1] <- t(t(dt[-1]) * scales)
  
  # Pairwise Fréchet distance computation for variable k
  pairwise <- function (i, j, k) {
    fredist(traj1 = dt[dt[, 1] == i, c(2, k+2)], traj2 = dt[dt[, 1] == j, c(2, k+2)], form = "scalar")
  }
  
  id <- unique(dt[, 1])
  dist.mat <- function(p) {outer(id, id, Vectorize(pairwise), k = p)}
  
  # Construct 3-dimensional distance array
  dist.ary <- NULL
  for (k in 1:num.var) {
    dist.ary <- abind(dist.ary, dist.mat(k), along = 3)
  }
  
  return(dist.ary)
}

# UpdateCs.SFclust: internal function, used within SFclust() function
# Performs k-medoids clustering given a (possibly weighted) distance array.
# The function supports both 3-dimensional arrays (multiple variables) and 2D matrices (single variable).
# When weights (wss) are provided, a weighted sum of distances is used.
# Returns a vector of cluster assignments.

UpdateCs.SFclust <- function(k, dist.ary, wss = NULL) {
  
  # Check if all weights are zero
  
  if (is.null(wss) == FALSE & all(wss == 0)) {
    stop("All weights are zero, so no further iterations can be performed.\n",
         "Please consider a l1bound parameter.")
  }
  
  if (length(dim(dist.ary)) == 3) {
    
    num.var <- dim(dist.ary)[3]
    
    # If no weights are given, assign equal weight scaled by 1/sqrt(num.var)
    if (is.null(wss)) wss <- rep(1/sqrt(num.var),  num.var)
    
    # Weighted distance array
    w.dist.ary <- dist.ary
    
    if (sum(wss == 0) > 0) {
      # When only a subset of weights are non-zero, subset w.dist.ary without collapsing the third dimension
      w.dist.ary <- w.dist.ary[,, which(wss!=0), drop = FALSE]
      wss <- wss[wss != 0]
    }
    
    # Apply weights to the distance arrays
    for (p in 1:length(wss)) w.dist.ary[,, p] <- w.dist.ary[,, p] * wss[p]
    
    id <- 1:dim(dist.ary)[1]
    
    # Initial random cluster assignment
    clust.assign <- ceiling(runif(id, 0, k))
    clust.assign.old <- ceiling(runif(id, 0, k))
    
    iter = 0
    while (!identical(clust.assign.old,  clust.assign)) {
      iter = iter + 1
      clust.assign.old <- clust.assign
      ctrd.idx <- NULL
      
      for (kk in 1:k) {
        clust.idx <- which(clust.assign == kk)
        
        
        if (length(clust.idx) == 1) {
          # If only one observation in cluster, assign it as the center
          ctrd.idx[kk] <- clust.idx
        } else if (length(clust.idx) == 2) {
          # Arbitrarily pick the first if two observations
          ctrd.idx[kk] <- clust.idx[1]
        } else {
          # Compute total distances within cluster and choose the medoid
          dist.mat.clust <- apply(w.dist.ary[clust.idx, clust.idx, ], c(1, 2), sum)
          min.idx <- which.min(apply(dist.mat.clust, 2, sum))
          ctrd.idx[kk] <- clust.idx[min.idx]
        }
      }
      
      # Assign trajectories to the closest cluster medoid
      dist.mat.ctrd.vs.all <- apply(w.dist.ary[, ctrd.idx, ], c(1, 2), sum)
      clust.assign <- apply(dist.mat.ctrd.vs.all, 1, which.min)
      
    }
  }
  
  if (length(dim(dist.ary)) == 2) {
    
    id <- 1:dim(dist.ary)[1]
    
    # Initial random cluster assignment
    clust.assign <- ceiling(runif(id, 0, k))
    clust.assign.old <- ceiling(runif(id, 0, k))
    
    iter = 0
    while (!identical(clust.assign.old, clust.assign)) {
      iter = iter + 1
      clust.assign.old <- clust.assign
      ctrd.idx <- NULL
      
      for (kk in 1:k) {
        clust.idx <- which(clust.assign == kk)
        
        if (length(clust.idx) == 1) {
          ctrd.idx[kk] <- clust.idx
        } else if (length(clust.idx) == 2) {
          ctrd.idx[kk] <- clust.idx[1]
        } else {
          dist.mat.clust <- dist.ary[clust.idx, clust.idx]
          min.idx <- which.min(apply(dist.mat.clust, 2, sum))
          ctrd.idx[kk] <- clust.idx[min.idx]
        }
      }
      
      dist.mat.ctrd.vs.all <- dist.ary[, ctrd.idx]
      clust.assign <- apply(dist.mat.ctrd.vs.all, 1, which.min)
      
    }
  }
  
  return(clust.assign)
}


# bcss.feature: internal function, used within SFclust() function
# Computes the Between-Cluster Sum of Squares (BCSS) for a given feature (variable).
# This is used to assess the contribution of each feature to the clustering structure.
# Returns the BCSS value for the p-th feature.

bcss.feature <- function(p, tot.num.obs, clust.num.obs, dist.ary, clust.assign) {
  
  # Compute total inertia (sum of all pairwise distances for the p-th feature)
  tot.inertia <- sum(dist.ary[,, p])/(2*tot.num.obs)
  clust.inertia <- 0
  
  # Compute within-cluster inertia for the p-th feature
  for (kk in 1:length(unique(clust.assign))) {
    clust.idx <- which(clust.assign == kk)
    clust.inertia <- clust.inertia + sum(dist.ary[clust.idx, clust.idx, p])/(2*clust.num.obs[kk])
  }
  
  # BCSS = Total inertia - Within-cluster inertia
  bcss <- tot.inertia - clust.inertia
  return(bcss)
}

# optim.w: internal function, used within SFclust() function
# Computes the optimal weight vector `w` for a given input vector `a` under an L1 norm constraint.
# The resulting weight vector is L2-normalized and satisfies a target L1 norm (approximately equal to `S`).
# Returns a numeric vector representing the optimized and normalized weights.

optim.w <- function(a, S) {
  
  # Sort input vector in decreasing order
  sa <- sort(a, decreasing = TRUE)
  p <- length(a)
  
  # Initialization
  lam1 <- rep(Inf, p)
  l1norm <- rep(Inf, p)
  l2norm <- rep(Inf, p)
  
  cnt  <-  0
  tot <- c()
  
  # Loop over possible truncation levels
  for (l in p:1) {
    cnt <- cnt + 1
    
    # t1: sum of top-l values, t2: sum of squared top-l values
    t1 <- sum(sa[1:l])
    t2 <- sum(sa[1:l]^2)
    
    # Solve quadratic to find lambda that controls L1 norm
    bb <- t1^2/l - (t1^2 - t2*S^2)/(l - S^2)
    b <- t1 * (l - S^2)
    aval <- l * (l - S^2)
    cval <- t1^2 - t2*S^2
    dett <- b^2 - aval*cval  # discriminant
    
    # If lambda is not real, skip
    if (bb < 0) {
      next
    }
    
    # Compute lambda and corresponding weight vector
    lam1[cnt] <- t1/l - 1/sqrt(l)*sqrt(bb)
    w_temp  =  pmax(a-lam1[cnt], 0)
    w <- w_temp/sqrt(sum(w_temp^2))
    
    # Store norms and weight vector
    l1norm[cnt] <- sum(w)
    l2norm[cnt] <- sum(w^2)
    tot <- rbind(tot, c(lam1[cnt], l2norm[cnt], l1norm[cnt], w))
  }
  
  tot <- data.frame(tot)
  colnames(tot)[1:3] <- c("lambda", "l2", "l1")
  colnames(tot)[4:(3+p)] <- paste0("x", seq(1, p))
  
  # Find solution closest to target L1 norm S
  l1.bound <- tot$l1[!is.nan(tot$l1)]
  soln.row <- which.min(abs(l1.bound - S))
  
  return(as.numeric(tot[soln.row, -(1:3)]))
}

# UpdateWs.SFclust: internal function, used within SFclust() function
# Computes updated feature weight vector (`wss`) based on BCSS (between-cluster sum of squares)
# for each feature and enforces sparsity via L1 norm constraint using `optim.w`.
#
# Output: A list with:
# - wss: Normalized feature weights satisfying the L1 norm constraint.
# - bcss: Raw BCSS values per feature before optimization.

UpdateWs.SFclust <- function(clust.assign, dist.ary, l1bound, num.var) {
  
  # Total number of observations
  tot.num.obs <- length(clust.assign)
  
  # Number of observations in each cluster
  clust.num.obs <- as.numeric(table(clust.assign))
  
  # Compute BCSS for each feature using `bcss.feature`
  bcss.perfeature <- sapply(1:num.var, bcss.feature, 
                            tot.num.obs =  tot.num.obs, 
                            clust.num.obs =  clust.num.obs, 
                            dist.ary = dist.ary, 
                            clust.assign = clust.assign)
  
  # Compute optimal weights using BCSS and target L1 norm
  wss <- optim.w(bcss.perfeature, l1bound)
  
  return(list(wss = wss, bcss = bcss.perfeature))
}


#' @title Sparse Fréchet Distance-Based K-medoids for Longitudinal Data
#'
#' @description Performs clustering on longitudinal trajectories using a sparse feature weighting
#' scheme and Fréchet distance. The method iteratively updates cluster assignments 
#' and feature weights subject to an \eqn{\ell_1} norm constraint.
#' 
#' @param k The number of clusters.
#' @param l1bound A bound on the \eqn{\ell_1} norm for the weight updates. It must lie between 1 and the square root of the number of variables.
#' @param dist.ary A 3-dimensional array of pairwise Fréchet distances. The array should be of shape (n, n, p), where \code{n} is the number of trajectories and 
#' \code{p} is the number of variables. Each \code{dist.ary[,,j]} stores the pairwise distances for the \code{j}-th variable.
#' @param maxIter The maximum number of iterations before stopping if convergence is not reached. Default is 20.
#' @param eps A small positive threshold for convergence. The algorithm stops when the change in weights becomes smaller than this threshold. Default is 1e-4.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{clust}{A vector of cluster assignments for each trajectory.}
#'   \item{final.weight}{The final weight vector after the last iteration, reflecting the contribution of each variable to the clustering process.}
#'   \item{weight.history}{A matrix of weight values at each iteration, showing how the feature weights evolved.}
#'   \item{criteria}{A vector of convergence criteria values for each iteration, quantifying the change in weights.}
#'   \item{iteration}{The number of iterations performed before convergence or reaching \code{maxIter}.}
#' }
#'
#' @details
#' The function assumes that the input \code{dist.array} contains pairwise distances between trajectories for each variable, 
#' using the generalized Fréchet distance. Clustering is performed via a k-medoids algorithm, 
#' and feature weights are updated using between-cluster sum of squares (BCSS) with sparsity control. 
#' If the number of variables is one, only clustering is performed, and no variable weighting is applied.
#' This function involves random sampling internally. For reproducible results, set the random seed before calling the function using \code{set.seed()}.
#' 
#' @export

SFclust <- function(k, l1bound, dist.ary, maxIter = 20, eps = 1e-4) {
  
  # Check if the number of trajectories is less than the number of clusters
  if (dim(dist.ary)[1] < k) {
    stop("Trajectories are less than clusters.\n")
  }
  
  # Check if l1bound is within the valid range
  if (l1bound < 1 || l1bound > sqrt(dim(dist.ary)[3])) {
    stop(sprintf("Error: l1bound must be between 1 and sqrt(number of variables) (%.4f).\n", sqrt(dim(dist.ary)[3])))
  }  
  
  # Determine the number of variables from the distance array
  if (length(dim(dist.ary)) == 3) {
    num.var <- dim(dist.ary)[3]
  } else if (length(dim(dist.ary)) == 2) {  # Handle the case when only one variable is present
    num.var <- 1
  }
  
  # Initialize cluster assignments
  css <- UpdateCs.SFclust(k, dist.ary)
  
  # If only one variable is provided, return cluster assignment only
  if (num.var == 1) {
    message("There is only one variable; weight update is skipped.")
    return(list(clust = css))
  }
  
  # Initialize variable weights
  wss <- UpdateWs.SFclust(css, dist.ary, l1bound, num.var)$wss
  
  # Check for invalid weights (NA)
  if (any(is.na(wss))) {
    stop(paste0("Error: Weight contains NA values after initialization, so no further iterations can be performed.\n",
                "Please consider a different l1bound instead of l1bound = ", l1bound, "."))
  }
  
  # Update clusters and weights until convergence
  count <- 0
  criteria <- NULL
  wss.old <- rep(1, length(wss))
  Wss <- NULL
  
  while (count < maxIter & (sum(abs(wss-wss.old)^2)/(num.var*sum(abs(wss.old)))) > eps) {
    
    count <- count + 1
    wss.old <- wss
    
    # Update cluster assignments with current weights
    css <- UpdateCs.SFclust(k, dist.ary, wss)
    
    # Update weights based on new assignments
    wss <- UpdateWs.SFclust(css, dist.ary, l1bound, num.var)$wss
    
    if (any(is.na(wss))) {
      stop(paste("Error: wss contains NA values at iteration", count))
    }
    
    # Track weight history and convergence criteria
    Wss <- rbind(Wss, wss)
    criteria <- c(criteria, sum(abs(wss-wss.old)^2)/(num.var*sum(abs(wss.old))))
  }
  
  return(list(clust = css, final.weight = wss, weight.history = Wss, criteria = criteria, iteration = count))
}


# permute.distance.array: internal function, used only within SFclust.permute() function
# This function performs permutation on the lower triangle of the distance matrix for each variable.
# It shuffles the pairwise distances between trajectories within each variable to create a permuted distance array.
# The function ensures that the permuted distance matrix retains the original structure while randomizing the lower triangle.
# Returns a 3-dimensional array of permuted pairwise distances with the same structure as the input distance array.

permute.distance.array <- function(dist.ary) {
  
  dist.ary.perm <- array(0, dim = dim(dist.ary))
  
  for (k in 1:dim(dist.ary.perm)[3]) {
    dist.ary.perm.p <- dist.ary[,, k]
    dist.ary.perm.p.low.elem <- dist.ary.perm.p[lower.tri(dist.ary.perm.p)]
    
    dist.ary.perm.p.low.samp <- sample(dist.ary.perm.p.low.elem)
    
    dist.ary.perm.p[lower.tri(dist.ary.perm.p)] <- dist.ary.perm.p.low.samp
    t.dist.ary.perm.p <- t(dist.ary.perm.p)
    t.dist.ary.perm.p[lower.tri(t.dist.ary.perm.p)] <- dist.ary.perm.p.low.samp
    dist.ary.perm[,,k] <- t.dist.ary.perm.p
    
  }
  
  return(dist.ary.perm)
}


#' @title Perform Permutation-Based Clustering Evaluation for SFclust
#'
#' @description Performs a permutation-based analysis to evaluate clustering results across different 
#' values of the \eqn{\ell_1} norm constraint (`s`). This function is designed to help determine the 
#' most appropriate \eqn{\ell_1} norm value by comparing the observed clustering outcome with those 
#' obtained under random permutations.
#' 
#' The function computes gap statistics for each \eqn{\ell_1} norm constraint value based on permuted 
#' versions of the input distance array, and identifies the optimal `s` as the one 
#' maximizing the gap statistic. Two ggplot objects are returned to visualize the gap patterns.
#'  
#' @param dist.ary A 3-dimensional distance array representing pairwise distances 
#' between trajectories across multiple variables. Follows the same format used in `SFclust`.
#' @param k An integer specifying the number of clusters.
#' @param nperms An integer specifying the number of permutations to perform.
#' @param l1b A numeric vector of \eqn{\ell_1} norm constraint values to test during clustering. These values 
#' control the sparsity of the weights during clustering.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{totss}{A numeric vector of total within-cluster sum of squared distances for each \eqn{\ell_1} norm value.}
#'   \item{permtotss}{A matrix of total sum of squared distances for each permutation and each \eqn{\ell_1} norm value.}
#'   \item{nnonzerowss}{A numeric vector of the number of nonzero weights for each \eqn{\ell_1} norm value.}
#'   \item{gaps}{A numeric vector of gap statistics: the difference between observed and permuted clustering results.}
#'   \item{sdgaps}{A numeric vector of standard deviations of the gaps across permutations.}
#'   \item{l1bounds}{A vector of \eqn{\ell_1} norm constraint values that were successfully processed without error.}
#'   \item{bestl1b}{The \eqn{\ell_1} norm constraint value that yielded the largest gap.}
#'   \item{failed_j}{Indices of `l1b` values that caused errors during the clustering process.}
#'   \item{failed_l1b}{The actual \eqn{\ell_1} norm values that caused errors.}
#'   \item{gapplot.l1b}{A ggplot object showing the gap statistics plotted against \eqn{\ell_1} norm constraint values.}
#'   \item{gapplot.nnz}{A ggplot object showing the gap statistics plotted against the number of nonzero weights.}
#' }
#'
#' @details
#' This function helps assess the robustness of clustering structure and select an optimal level of sparsity.
#' If any clustering attempt fails (e.g., due to convergence issues or weight update errors), the corresponding 
#' `l1b` values are reported in `failed_l1b` and `failed_j`.
#' This function returns two ggplot objects (`gapplot.l1b` and `gapplot.nnz`) that can be used to visualize the
#' gap statistics. These are not automatically printed, allowing users to decide when and how to display them.
#' This function involves random sampling internally. For reproducible results, set the random seed before calling the function using \code{set.seed()}.
#' 
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal
#' @importFrom stats runif na.omit sd
#' 
#' @export  

SFclust.permute <- function(dist.ary, k, nperms, l1b) {
  
  # Check if distance array  has three dimensions  
  if (length(dim(dist.ary)) != 3) {
    stop("There is only one variable: no need to tune parameter w")
  }
  
  # Clustering with original distance array for each l1b
  nnonzerowss <- NULL
  totss <- NULL
  failed_j <- c()
  
  for (j in 1:length(l1b)) {
    out <- tryCatch(
      SFclust(k, l1b[j], dist.ary = dist.ary),
      error = function(e) {
        message("Error at l1b[", j, "] :", e$message)
        return(NULL)
      }
    )
    if (is.null(out)) {
      failed_j <- c(failed_j, j)
      next
    }
    
    nnonzerowss <- c(nnonzerowss, sum(out$final.weight != 0))
    bcss <- UpdateWs.SFclust(out$clust, dist.ary, l1b[j], dim(dist.ary)[3])$bcss
    totss <- c(totss, sum(out$final.weight * bcss))
  }
  
  # Create permuted distance arrays
  permx <- list()
  seed.vec <- round(runif(nperms, 1, 1000000))
  
  for (kk in 1:nperms) {
    set.seed(seed.vec[kk])
    permx[[kk]] <- permute.distance.array(dist.ary)
  }
  
  permtotss <- matrix(NA, nrow = length(setdiff(1:length(l1b), failed_j)), ncol = nperms)
  
  for (kk in 1:nperms) {
    dist.ary.pm <- permx[[kk]]
    perm.totss <- numeric(length(l1b[setdiff(1:length(l1b), failed_j)]))
    perm.totss[] <- NA
    
    for (j in setdiff(1:length(l1b), failed_j)) {
      perm.out <- tryCatch(
        SFclust(k, l1b[j], dist.ary = dist.ary.pm),
        error = function(e) {
          message("Error at permuted dist.ary for l1b[", j, "] (l1b = ", l1b[j], ") :", e$message)
          return(NULL) 
        }
      )
      if (is.null(perm.out)) {
        failed_j <- c(failed_j, j)
        perm.totss[j] <- NA
        next
      }
      
      perm.bcss <- UpdateWs.SFclust(perm.out$clust, dist.ary.pm, l1b[j], dim(dist.ary.pm)[3])$bcss
      perm.totss[j] <- sum(perm.out$final.weight * perm.bcss)
    }
    
    permtotss[, kk] <- perm.totss
  }
  
  permtotss <- na.omit(permtotss)
  
  # Compute gap statistic
  gaps <- log(totss[setdiff(1:length(l1b), failed_j)]) - apply(log(permtotss), 1, mean)
  
  # Create gap plots (no printing here)
  gapplot.l1b <- ggplot(data = data.frame(x = l1b[setdiff(1:length(l1b), failed_j)], y = gaps), aes(x = x, y = y)) + 
    geom_line() + 
    geom_point() +
    labs(x = "l1bound (s)", y = "Gap", title = paste("k = ", k)) +
    theme_minimal()
  
  gapplot.nnz <- ggplot(data = data.frame(x = nnonzerowss, y = gaps), aes(x = x, y = y)) + 
    geom_point() + 
    labs(x = "Number of nonzero weights", y = "Gap", title = paste("k = ", k)) +
    theme_minimal()
  
  return(invisible(list(totss = totss, permtotss = permtotss, nnonzerowss = nnonzerowss, 
                        gaps = gaps, sdgaps = apply(log(permtotss), 1, sd), 
                        l1bounds = l1b[setdiff(1:length(l1b), failed_j)], bestl1b = l1b[which.max(gaps)],
                        failed_j = failed_j, failed_l1b = l1b[failed_j],
                        gapplot.l1b = gapplot.l1b, gapplot.nnz = gapplot.nnz)))
}
