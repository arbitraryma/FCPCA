

library(base)
library(fungible)


# Function to compute the weighted cross-covariance matrices for different clusters, and then is used to construct the common projection axes 

compute_weighted_cross_cov <- function(matrixU, m, cross_cov_lag0_lag1, cross_cov_lag0_lag2) {
  n <- nrow(matrixU)  # Number of MTS objects
  k <- ncol(matrixU)  # Number of clusters
  
  # Initialize lists to store weighted cross-covariance matrices for each cluster
  weighted_cross_cov_lag0_lag1 <- vector("list", k)
  weighted_cross_cov_lag0_lag2 <- vector("list", k)
  
  for (cluster in 1:k) {  # Iterate over each cluster
    # Get the weights for the current cluster, using the fuzziness parameter m
    weights <- matrixU[, cluster]^m
    
    # Normalize weights to prevent division by zero
    total_weight <- sum(weights)
    if (total_weight == 0) {
      stop("Total weight for cluster ", cluster, " is zero. Check your input matrixU or fuzziness parameter.")
    }
    
    # Compute the weighted cross-covariance for lag 0-lag 1 for this cluster
    weighted_cross_cov_lag0_lag1[[cluster]] <- Reduce(`+`, 
                                                      Map(function(cov_matrix, weight) weight * cov_matrix, cross_cov_lag0_lag1, weights)) / total_weight
    
    # Compute the weighted cross-covariance for lag 0-lag 2 for this cluster
    weighted_cross_cov_lag0_lag2[[cluster]] <- Reduce(`+`, 
                                                      Map(function(cov_matrix, weight) weight * cov_matrix, cross_cov_lag0_lag2, weights)) / total_weight
  }
  
  # Return a list containing the weighted cross-covariance matrices for each cluster
  return(list(
    weighted_cross_cov_lag0_lag1 = weighted_cross_cov_lag0_lag1, 
    weighted_cross_cov_lag0_lag2 = weighted_cross_cov_lag0_lag2
  ))
}



# This function standardizes each multivariate time series (MTS) object by removing the mean of each variable (column-wise) 
# and possibly scaling by some measure of spread (e.g., standard deviation). 
standardize_mts <- function(mts_data) {
  # mts_data is a list where each element is an MTS object (a matrix)
  lapply(mts_data, function(ts) {
    scale(ts, center = TRUE, scale = TRUE) # Subtracts the column means
  })
}


# This function calculates the cross-covariance matrices for each MTS object in a list. 
# It constructs a block matrix including both the lag-0 covariance matrix and the lag-1 cross-covariance matrix:

cross_covariance_lag1 <- function(mts_data) {
  # mts_data is a multivariate time series (matrix or data frame)
  
  n <- nrow(mts_data) # number of time points
  p <- ncol(mts_data) # number of variables
  
  # Check if there are enough time points
  if (n < 2) {
    stop("The multivariate time series should have at least 2 time points.")
  }
  
  # Compute Gamma_0 (lag-0 covariance matrix), biased
  Gamma_0 <- cov(mts_data)
  
  # Compute Gamma_1 (lag-1 cross-covariance matrix), biased
  X_t <- mts_data[1:(n-1), ]    # Data from time t
  X_t_plus_1 <- mts_data[2:n, ]  # Data from time t+1
  
  # Cross-covariance between X_t and X_{t+1}, biased
  Gamma_1 <- cov(X_t, X_t_plus_1)
  
  # Construct the block matrix
  block_matrix <- rbind(
    cbind(Gamma_0, Gamma_1),   # First row: [Gamma(0), Gamma(1)]
    cbind(t(Gamma_1), Gamma_0) # Second row: [Gamma(-1), Gamma(0)]
  )
  
  output <- list(block_matrix, X_t, X_t_plus_1)
  return(output)
}



# this is the function when we consider lag-0 and lag-2 

cross_covariance_lag2 <- function(mts_data) {
  # mts_data is a multivariate time series (matrix or data frame)
  
  n <- nrow(mts_data) # number of time points
  p <- ncol(mts_data) # number of variables
  
  # Check if there are enough time points
  if (n < 3) { # Need at least 3 time points for lag 2
    stop("The multivariate time series should have at least 3 time points.")
  }
  
  # Compute Gamma_0 (lag-0 covariance matrix), biased
  Gamma_0 <- cov(mts_data)
  
  # Compute Gamma_2 (lag-2 cross-covariance matrix), biased
  X_t <- mts_data[1:(n-2), ]    # Data from time t
  X_t_plus_2 <- mts_data[3:n, ]  # Data from time t+2
  
  # Cross-covariance between X_t and X_{t+2}, biased
  Gamma_2 <- cov(X_t, X_t_plus_2)
  
  # Construct the block matrix with lag 2
  block_matrix <- rbind(
    cbind(Gamma_0, Gamma_2),   # First row: [Gamma(0), Gamma(2)]
    cbind(t(Gamma_2), Gamma_0) # Second row: [Gamma(-2), Gamma(0)]
  )
  
  output <- list(block_matrix, X_t, X_t_plus_2)
  return(output)
}


combine_xt_xt1 <- function(X_t, X_t_plus_1) {
  # X_t: (n-1) x p matrix for X_t
  # X_t_plus_1: (n-1) x p matrix for X_t_plus_1
  
  # Combine X_t and X_t_plus_1 into a single matrix of dimension (n-1) x 2p
  combined_matrix <- cbind(X_t, X_t_plus_1)  # Combine along columns
  
  return(combined_matrix)  # Return the combined matrix
}




# this function helps to determine the selection of the hyper-paramter m, where the smallest S is preferred 
compute_S <- function(reconstruction_error, projection_axes,n, k) {
  # 1) Numerator: total reconstruction error
  numerator <- sum(reconstruction_error)
  
  # 2) Build the projection matrices P_c(lag) for each cluster c and each lag
  
  P_list <- vector("list", k)
  
  for (c in seq_len(k)) {
    # Extract the principal axes for lag 1 and lag 2
    Clag1 <- projection_axes[[c]]$lag01  # C^(c,1), dimension p x k1
    Clag2 <- projection_axes[[c]]$lag02  # C^(c,2), dimension p x k2
    
    # Form the projection matrices P_c(1) = C^(c,1) * C^(c,1)^T
    P_lag1 <- Clag1 %*% t(Clag1)
    P_lag2 <- Clag2 %*% t(Clag2)
    
    # Store them
    P_list[[c]] <- list(P_lag1 = P_lag1,
                        P_lag2 = P_lag2)
  }
  
  # 3) Define a function to measure subspace distance across both lags
  #    Here, we sum the Frobenius norms for lag1 and lag2
  subspace_dist <- function(Pc, Pm) {
    d_lag1 <- norm(Pc$P_lag1 - Pm$P_lag1, type = "F")
    d_lag2 <- norm(Pc$P_lag2 - Pm$P_lag2, type = "F")
    return(d_lag1 + d_lag2)
  }
  
  # 4) Find the minimum pairwise distance among all clusters
  min_dist <- Inf
  for (c in seq_len(k)) {
    for (m in seq_len(k)) {
      if (m > c) {
        dist_val <- subspace_dist(P_list[[c]], P_list[[m]])
        if (dist_val < min_dist) {
          min_dist <- dist_val
        }
      }
    }
  }
  
  # 5) The denominator is the smallest distance among cluster subspaces
  denominator <- min_dist
  
  # 6) Final S value
  S_value <- numerator / (n*denominator)
  list(S_value = S_value, P_list = P_list)
}





# this is the main function for the fuzzy CPCA method  with 3 times replications
# and pre-defined number of principal components, the final result is with the lowest weighted sum reconstruction error 






###############################################################################
## FCPCA with Auto/Grid Search over k and m (select by lowest S)
###############################################################################

# Core single-run FCPCA (standard L2, no exponential), returns S_value via compute_S 
.fcpca_run_once <- function(ts,
                            k,
                            m = 1.5,
                            startU = NULL,
                            conver = 1e-3,
                            maxit = 1000,
                            verbose = FALSE) {
  
  ts <- standardize_mts(ts)
  
  n <- length(ts)
  # lag-1 and lag-2 blocks + combined [X_t, X_{t+1}] and [X_t, X_{t+2}]
  res1  <- lapply(ts, cross_covariance_lag1)
  res2  <- lapply(ts, cross_covariance_lag2)
  
  sigma  <- lapply(res1, function(x) x[[1]])
  sigma2 <- lapply(res2, function(x) x[[1]])
  
  comb1  <- mapply(combine_xt_xt1, lapply(res1, `[[`, 2), lapply(res1, `[[`, 3), SIMPLIFY = FALSE)
  comb2  <- mapply(combine_xt_xt1, lapply(res2, `[[`, 2), lapply(res2, `[[`, 3), SIMPLIFY = FALSE)
  
  # init memberships
  if (is.null(startU)) {
    U <- matrix(runif(n * k), nrow = n); U <- U / rowSums(U)
  } else U <- startU
  
  # first weighted covs to size the subspaces (95% energy per lag)
  wcov <- compute_weighted_cross_cov(U, m, sigma, sigma2)
  
  svd01_1 <- svd(wcov$weighted_cross_cov_lag0_lag1[[1]])
  k1 <- which(cumsum(svd01_1$d)/sum(svd01_1$d) >= 0.95)[1]; if (is.na(k1)) k1 <- length(svd01_1$d)
  
  svd02_1 <- svd(wcov$weighted_cross_cov_lag0_lag2[[1]])
  k2 <- which(cumsum(svd02_1$d)/sum(svd02_1$d) >= 0.95)[1]; if (is.na(k2)) k2 <- length(svd02_1$d)
  
  iteration <- 0
  diffU <- Inf
  prev_obj <- Inf
  cur_obj  <- Inf
  
  repeat {
    iteration <- iteration + 1
    
    # update weighted covariances
    wcov <- compute_weighted_cross_cov(U, m, sigma, sigma2)
    
    # update projection axes (fixed dims k1, k2 chosen above)
    projAx <- vector("list", k)
    for (c in 1:k) {
      s1 <- svd(wcov$weighted_cross_cov_lag0_lag1[[c]])
      s2 <- svd(wcov$weighted_cross_cov_lag0_lag2[[c]])
      U1 <- s1$u[, 1:k1, drop = FALSE]
      U2 <- s2$u[, 1:k2, drop = FALSE]
      projAx[[c]] <- list(lag01 = U1, lag02 = U2)
    }
    
    # reconstruction error matrix (n x k)
    R <- matrix(0, nrow = n, ncol = k)
    for (i in 1:n) for (c in 1:k) {
      rec1 <- comb1[[i]] %*% projAx[[c]]$lag01 %*% t(projAx[[c]]$lag01)
      rec2 <- comb2[[i]] %*% projAx[[c]]$lag02 %*% t(projAx[[c]]$lag02)
      R[i, c] <- sum((comb1[[i]] - rec1)^2) + sum((comb2[[i]] - rec2)^2)
    }
    
    # standard fuzzy update
    U_new <- matrix(0, nrow = n, ncol = k)
    for (i in 1:n) for (c in 1:k) {
      denom <- sum((R[i, c] / (R[i, ] + 1e-12))^(1/(m - 1)))
      U_new[i, c] <- 1 / denom
    }
    
    diffU   <- max(abs(U_new - U))
    U       <- U_new
    prev_obj <- cur_obj
    cur_obj  <- sum((U^m) * R)
    
    if (verbose) cat(sprintf("Iter %d: Obj=%.6f, dU=%.6f\n", iteration, cur_obj, diffU))
    if ((abs(prev_obj - cur_obj) < conver && diffU < conver) || iteration >= maxit) break
  }
  
  final_wcov <- compute_weighted_cross_cov(U, m, sigma, sigma2)
  S_val <- compute_S(sum(U^m * R), final_wcov, n, k)
  
  list(
    membership_matrix           = U,
    projection_axes             = projAx,
    iterations                  = iteration,
    converged                   = (iteration < maxit),
    reconstruction_error        = cur_obj,
    reconstruction_error_matrix = R,
    hard_cluster                = max.col(U),
    weighted_covariances        = final_wcov,
    n_components_lag01          = k1,
    n_components_lag02          = k2,
    S_value                     = S_val,
    k                           = k,
    m                           = m,
    sigma                       = sigma,
    sigma2                      = sigma2,
    combined_list               = comb1,
    combined_list2              = comb2
  )
}

# Public fcpca(): single (k,m) OR grid over k,m (auto), parallel optional; selects min S_value
fcpca <- function(ts,
                  k = NULL,
                  m = NULL,
                  startU = NULL,
                  conver = 1e-3,
                  maxit = 1000,
                  verbose = TRUE,
                  parallel = TRUE) {
  
  # auto/grid ranges
  auto_mode <- is.null(k) || is.null(m) || length(k) > 1 || length(m) > 1
  if (auto_mode) {
    k_range <- if (is.null(k)) 2:6 else as.integer(k)
    m_range <- if (is.null(m)) c(1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5) else as.numeric(m)
  } else {
    k_range <- as.integer(k)
    m_range <- as.numeric(m)
  }
  grid <- expand.grid(k = k_range, m = m_range)
  
  run_once <- function(idx) {
    kk <- grid$k[idx]; mm <- grid$m[idx]
    if (verbose) cat(sprintf("FCPCA run: k=%d, m=%.2f\n", kk, mm))
    .fcpca_run_once(ts, k = kk, m = mm, startU = startU, conver = conver, maxit = maxit, verbose = FALSE)
  }
  
  if (auto_mode) {
    # parallel map (optional)
    results_list <- if (parallel) {
      if (!requireNamespace("furrr", quietly = TRUE)) stop("Package 'furrr' required for parallel=TRUE. Install it or set parallel=FALSE.")
      furrr::future_map(seq_len(nrow(grid)), run_once, .options = furrr::furrr_options(seed = TRUE))
    } else {
      # sequential with simple early stopping on exploding S (same spirit as RFCPCA_E)
      out <- list(); best_S <- Inf
      for (i in seq_len(nrow(grid))) {
        res <- run_once(i)
        out[[i]] <- res
        if (!is.na(res$S_value) && res$S_value < best_S) best_S <- res$S_value
        if (!is.infinite(best_S) && res$S_value > best_S * 100) {
          if (verbose) cat("Early stopping triggered (S spiked).\n")
          break
        }
      }
      out
    }
    
    S_vals <- sapply(results_list, function(x) x$S_value)
    best_idx <- which.min(S_vals)
    best <- results_list[[best_idx]]
    
    best$optimal_k   <- grid$k[best_idx]
    best$optimal_m   <- grid$m[best_idx]
    best$all_results <- results_list
    best
    
  } else {
    # single setting
    res <- .fcpca_run_once(ts, k = k_range[1], m = m_range[1], startU = startU,
                           conver = conver, maxit = maxit, verbose = verbose)
    res$optimal_k   <- k_range[1]
    res$optimal_m   <- m_range[1]
    res$all_results <- list(res)
    res
  }
}
