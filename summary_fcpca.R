# This function is the fuzzy CPCA clustering algorithm, all functions are wrote in this file

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
  return(S_value)
}


# this is the main function for the fuzzy CPCA method  with 3 times replications
# and pre-defined number of principal components, the final result is with the lowest weighted sum reconstruction error 


fcpca <- function(ts, k, m = 1.5, startU = NULL, conver = 1e-3, maxit = 1000, verbose = TRUE, replicates = 1) {
  # ts: List of N MTS objects
  # k: Number of clusters
  # m: Fuzziness parameter
  # startU: Initial membership matrix
  # conver: Convergence criterion
  # maxit: Maximum number of iterations
  # verbose: Whether to display iteration details
  # replicates: Number of times to repeat clustering
  
  # Standardize MTS objects
  ts <- standardize_mts(ts)
  
  n <- length(ts)  # Number of MTS objects
  p <- ncol(ts[[1]])  # Dimensionality of MTS
  
  # Apply the cross_covariance function to each MTS object
  results <- lapply(ts, cross_covariance_lag1) # this helps to obtain X_t, X_t+1, and lag-1 block matrix
  results2 <-  lapply(ts, cross_covariance_lag2) # this helps to obtain X_t, X_t+2, and lag-2 block matrix
  
  # Extract the block matrices (sigma), X_t (tsxt), and X_t_plus_1 (tsxt1) into separate lists
  
  sigma <- lapply(results, function(x) x[[1]]) # with dimension 2p * 2p, lag1 cross-covariance 
  
  sigma2 <-  lapply(results2, function(x) x[[1]]) # with dimension 2p * 2p, lag2 cross-covariance 
  
  tsxt <- lapply(results, function(x) x[[2]])  #  X_t, with length n-1
  tsxt1 <- lapply(results, function(x) x[[3]]) #  X_t+1,  with length n-1
  
  combined_list <- mapply(combine_xt_xt1, tsxt, tsxt1, SIMPLIFY = FALSE) # with (X_t,X_t+1) this gives dimension (n-1)*2p 
  
  tsxt2 <- lapply(results2, function(x) x[[2]])  # X_t, with length n-2
  tsxt12 <- lapply(results2, function(x) x[[3]]) # X_t+2, with length n-2
  combined_list2 <- mapply(combine_xt_xt1, tsxt2, tsxt12, SIMPLIFY = FALSE) # this is (X_t, X_t+2), with dimension (n-2) * 2p
  
  
  # Compute weighted covariance matrices for initial PCs
  if (is.null(startU)) {
    U <- matrix(runif(n * k), nrow = n)
    U <- U / rowSums(U) # normalize rows to sum to 1
  } else {
    U <- startU
  }
  weighted_covariances <- compute_weighted_cross_cov(U, m, sigma, sigma2)
  
  # Compute SVD for lag-1 and lag-2 cross-covariances to determine fixed PCs
  svd_result_lag01 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag1[[1]])
  n_components_lag01 <- which(cumsum(svd_result_lag01$d) / sum(svd_result_lag01$d) >= 0.95)[1]
  
  
  
  svd_result_lag02 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag2[[1]])
  n_components_lag02 <- which(cumsum(svd_result_lag02$d) / sum(svd_result_lag02$d) >= 0.95)[1]
    
  
  # Repeat clustering multiple times
  best_result <- NULL
  best_error <- Inf
  
  for (replicate in 1:replicates) {
    if (verbose) cat(sprintf("\nReplicate %d:\n", replicate))
    
    # Initialize membership matrix U, for each replicate, a different U is generated in case a very poor initialization
    if (is.null(startU)) {
      U <- matrix(runif(n * k), nrow = n)
      U <- U / rowSums(U) # normalize rows to sum to 1
    } else {
      U <- startU
    }
    
    # Initialize variables for convergence
    iteration <- 0
    diff <- Inf
    prev_error <- Inf
    current_error <- Inf
    
    # Start iteration
    while (iteration < maxit && (abs(prev_error - current_error) > conver || diff > conver)) {
      iteration <- iteration + 1
      
      # Compute weighted covariance matrices
      weighted_covariances <- compute_weighted_cross_cov(U, m, sigma, sigma2)
      
      # Update projection axes using fixed number of PCs
      projection_axes <- vector("list", k)
      for (c in 1:k) {
        # SVD for lag-1
        svd_result_lag01 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag1[[c]])
        proj_axes_lag01 <- svd_result_lag01$u[, 1:n_components_lag01]
        
        # SVD for lag-2
        svd_result_lag02 <- svd(weighted_covariances$weighted_cross_cov_lag0_lag2[[c]])
        proj_axes_lag02 <- svd_result_lag02$u[, 1:n_components_lag02]
        
        projection_axes[[c]] <- list(lag01 = proj_axes_lag01, lag02 = proj_axes_lag02)
      }
      
      # Compute reconstruction errors
      total_reconstruction_error <- matrix(0, nrow = n, ncol = k)
      for (i in 1:n) {
        for (c in 1:k) {
          reconstructed_tt1 <- combined_list[[i]] %*% projection_axes[[c]]$lag01 %*% t(projection_axes[[c]]$lag01)
          reconstruction_error_tt1 <- norm(combined_list[[i]] - reconstructed_tt1, "F")^2
          
          reconstructed_tt2 <- combined_list2[[i]] %*% projection_axes[[c]]$lag02 %*% t(projection_axes[[c]]$lag02)
          reconstruction_error_tt2 <- norm(combined_list2[[i]] - reconstructed_tt2, "F")^2
          
          total_reconstruction_error[i, c] <- reconstruction_error_tt1 + reconstruction_error_tt2
        }
      }
      
      # Update membership matrix
      U_new <- matrix(0, nrow = n, ncol = k)
      for (i in 1:n) {
        for (c in 1:k) {
          sum_term <- sum((total_reconstruction_error[i, c] / (total_reconstruction_error[i, ] + 1e-10))^(1 / (m - 1)))
          U_new[i, c] <- 1 / sum_term
        }
      }
      
      # Update convergence metrics
      diff <- max(abs(U - U_new))
      prev_error <- current_error
      U <- U_new
      current_error <- sum(U^m * total_reconstruction_error)
      
      if (verbose) {
        cat(sprintf("Iteration %d: Reconstruction Error = %.6f, Max Change in U = %.6f\n", 
                    iteration, current_error, diff))
      }
    }
    
    # Check if this replicate is the best
    if (current_error < best_error) {
      best_error <- current_error
      best_result <- list(
        membership_matrix = U,
        projection_axes = projection_axes,
        iterations = iteration,
        converged = iteration < maxit,
        reconstruction_error = current_error,
        reconstruction_error_matrix = total_reconstruction_error,
        hard_cluster = max.col(U),
        weighted_covariances = weighted_covariances,
        S = compute_S(current_error,projection_axes,n,k)
      )
    }
  }
  
  # Return the best result
  return(best_result)
}
















