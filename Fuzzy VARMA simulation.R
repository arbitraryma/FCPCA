library(gtools)
library(base)
library(MTS)
library(dtwclust)
library(mlmts)
# function to generate high dimensional coefficient matrix for VAR process 
generate_stationary_matrix_optimized <- function(dim, lags = 1, tol = 1e-6) {
  phi <- array(0, dim = c(dim, dim, lags))
  for (lag in 1:lags) {
    mat <- matrix(runif(dim^2, -1, 1), nrow = dim)
    spectral_radius <- max(abs(eigen(mat, only.values = TRUE)$values))
    if (spectral_radius >= 1 - tol) {
      mat <- mat / (spectral_radius + tol)  # Scale to make spectral radius < 1
    }
    phi[, , lag] <- mat
  }
  return(phi)
}

# function to generate high dimensional coefficient matrix for VMA process 

generate_invertible_matrix_optimized <- function(dim, lags = 1, tol = 1e-6) {
  theta <- array(0, dim = c(dim, dim, lags))
  for (lag in 1:lags) {
    mat <- matrix(runif(dim^2, -1, 1), nrow = dim)
    spectral_radius <- max(abs(eigen(mat, only.values = TRUE)$values))
    if (spectral_radius >= 1 - tol) {
      mat <- mat / (spectral_radius + tol)  # Scale to make spectral radius < 1
    }
    theta[, , lag] <- mat
  }
  return(theta)
}




# Assign objects to a group based on 0.65 threshold
assign_groups <- function(membership_matrix, threshold = 0.7) {
  # Args:
  #   membership_matrix: A matrix of membership values (rows = objects, cols = clusters)
  #   threshold: The minimum membership value required to assign an object to a group
  #
  # Returns:
  #   A vector of group assignments (1 to k for assigned objects, NA for unassigned objects)
  
  # Apply the threshold to each row of the membership matrix
  group_assignment <- apply(membership_matrix, 1, function(row) {
    if (max(row) >= threshold) {
      which.max(row)  # Assign to the group with the highest membership
    } else {
      NA  # Uncertain or not assigned to any group
    }
  })
  
  return(group_assignment)
}





calculate_accuracy <- function(predicted_labels, true_labels) {
  # Args:
  #   predicted_labels: A vector of cluster assignments (from clustering algorithm)
  #   true_labels: A vector of ground truth cluster labels
  #
  # Returns:
  #   Accuracy as a percentage (between 0 and 100)
  
  # Ensure both vectors are the same length
  if (length(predicted_labels) != length(true_labels)) {
    stop("Predicted and true labels must have the same length.")
  }
  
  # Identify indices where the ground truth expects NA
  expected_na_indices <- which(is.na(true_labels))
  
  # Calculate the number of correctly classified expected NA
  correct_na_classifications <- sum(is.na(predicted_labels[expected_na_indices]))
  
  # Valid indices for permutation-based comparison (exclude expected NA)
  valid_indices <- setdiff(seq_along(true_labels), expected_na_indices)
  valid_predicted_labels <- predicted_labels[valid_indices]
  valid_true_labels <- true_labels[valid_indices]
  
  # If there are no valid labels to compare, return NA
  if (length(valid_true_labels) == 0 || length(valid_predicted_labels) == 0) {
    return(NA)
  }
  
  # Unique cluster labels (excluding NAs)
  unique_true <- unique(na.omit(valid_true_labels))
  unique_predicted <- unique(na.omit(valid_predicted_labels))
  
  # Generate all permutations of the true labels
  permutations <- gtools::permutations(length(unique_true), length(unique_true), unique_true)
  
  # Initialize best accuracy
  best_accuracy <- 0
  
  # Loop through each permutation to calculate accuracy
  for (perm in 1:nrow(permutations)) {
    # Create a mapping between predicted labels and the permuted true labels
    mapping <- setNames(permutations[perm, ], unique_predicted)
    
    # Map predicted labels to the permuted true labels
    remapped_labels <- sapply(valid_predicted_labels, function(x) {
      if (!is.na(x) && !is.null(mapping[x])) {
        mapping[x]
      } else {
        NA
      }
    })
    
    # Count correct classifications for valid indices
    correct_predictions <- sum(remapped_labels == valid_true_labels, na.rm = TRUE)
    
    # Include correct NA classifications
    total_correct <- correct_predictions + correct_na_classifications
    
    # Total labels include all true labels (both valid and expected NA)
    total_labels <- length(true_labels)
    
    # Calculate accuracy
    accuracy <- total_correct / total_labels 
    
    # Update best accuracy
    best_accuracy <- max(best_accuracy, accuracy)
  }
  
  return(best_accuracy)
}




# The simulation part, the complete code is not provided, we just show the methods, parameters we use in our simulation, if you would like to test, feel free to add the loops 
# and replicates and check the results. 
labels <- c(rep(1,10),rep(2,10),NA,NA)

rep <- 100
 
f <- c(20, 60, 100)

for (j in 1:length(f)) {
  e  <- f[j]
  # for this simulation, we fix the coefficient matrix to generate the var and ma, VARMA process
  m <- drop(generate_stationary_matrix_optimized(e))
  n <- drop(generate_invertible_matrix_optimized(e))
  w <- generate_stationary_matrix(e)
  Sigma <- diag(e)
  Sigma_VAR <- diag(e) * 0.1  # Smaller variance for VAR
  Sigma_VMA <- diag(e) * 0.5  # Larger variance for VMA
  
  Sigma_VARMA <- diag(e) * 0.3 
  
 #   lengths1 <- sample(200:600, 10)
 #   lengths2 <- sample(200:600, 10)
 #    lengths3 <- sample(200:600, 2)
      
     lengths1 <- rep(200,10)
       lengths2 <- rep(200,10)
     lengths3 <- rep(200,2)

   #   lengths1 <- sample(400, 10)
 #   lengths2 <- sample(400, 10)
 #    lengths3 <- sample(400, 2)
      
    
    data_list_var1 <- lapply(lengths1, function(len) {
      MTS::VARMAsim(len, arlags = 1, phi = m * 0.8 , sigma = Sigma_VAR )$series
    })
    
    
    data_list_var2 <- lapply(lengths2, function(len) {
      MTS::VARMAsim(len, malags = 1, theta = n * 2  , sigma = Sigma_VMA )$series
    })
    
 
    
    data_list_var3 <- lapply(lengths3, function(len) {
      phi_VARMA <- m * 0.8   # Add interaction and scale
      theta_VARMA <- n * 2
      Sigma_VARMA <- diag(e) * 0.3  # Adjust noise
      ts <- MTS::VARMAsim(len, arlags = 1, malags = 1, phi = phi_VARMA, theta = theta_VARMA, sigma = Sigma_VARMA)$series
      ts + matrix(rnorm(length(ts), 0, 0.5), ncol = e)  # Add noise
    })
    
    ts <- c(data_list_var1, data_list_var2, data_list_var3)    
    
 
    roblag1 <- fcpca(ts,2,1.2,replicates = 1)
    roblag2 <- fcpca(ts,2,1.4,replicates = 1)
    roblag3 <- fcpca(ts,2,1.6,replicates = 1)
     roblag4 <- fcpca(ts,2,1.8,replicates = 1)
     roblag5 <- fcpca(ts,2,2,replicates = 1)
    roblag6 <- fcpca(ts,2,2.2,replicates = 1)
 
    
    roblag10 <- assign_groups(roblag1$membership_matrix)
    roblag20 <- assign_groups(roblag2$membership_matrix)
    roblag30 <- assign_groups(roblag3$membership_matrix)
    roblag40 <- assign_groups(roblag4$membership_matrix)
    roblag50 <- assign_groups(roblag5$membership_matrix)
    roblag60 <- assign_groups(roblag6$membership_matrix)
    
 
    calculate_accuracy(roblag10, labels )
    calculate_accuracy(roblag20, labels )
    calculate_accuracy(roblag30, labels )
    calculate_accuracy(roblag40, labels )
    calculate_accuracy(roblag50, labels )
    calculate_accuracy(roblag60, labels )
    
    
    roblag11 <- assign_groups(roblag1$membership_matrix,threshold = 0.6)
    roblag21 <- assign_groups(roblag2$membership_matrix,threshold = 0.6)
    roblag31 <- assign_groups(roblag3$membership_matrix,threshold = 0.6)
    roblag41 <- assign_groups(roblag4$membership_matrix,threshold = 0.6)
    roblag51 <- assign_groups(roblag5$membership_matrix,threshold = 0.6)
    roblag61 <- assign_groups(roblag6$membership_matrix,threshold = 0.6)
    
     calculate_accuracy(roblag11 , labels )
     calculate_accuracy(roblag21, labels )
     calculate_accuracy(roblag31, labels )
     calculate_accuracy(roblag41, labels )
     calculate_accuracy(roblag51, labels )
     calculate_accuracy(roblag61, labels )
    
    
    
    
    # vpca method as a comparison
     
    vp1 <- vpca_clustering(ts,2,m = 1.2,crisp = FALSE)
    vp2 <- vpca_clustering(ts,2,m = 1.4,crisp = FALSE)
    vp3 <- vpca_clustering(ts,2,m = 1.6,crisp = FALSE)
    vp4 <- vpca_clustering(ts,2,m = 1.8,crisp = FALSE)
    vp5 <- vpca_clustering(ts,2,m = 2,crisp = FALSE)
    vp6 <- vpca_clustering(ts,2,m = 2.2,crisp = FALSE)
    
    vp11 <- assign_groups(t(vp1$U))
    vp21 <- assign_groups(t(vp2$U))
    vp31 <- assign_groups(t(vp3$U))
    vp41 <- assign_groups(t(vp4$U))
    vp51 <- assign_groups(t(vp5$U))
    vp61 <- assign_groups(t(vp6$U))
    
    
    calculate_accuracy(vp11, labels )
    calculate_accuracy(vp21, labels )
     calculate_accuracy(vp31, labels )
     calculate_accuracy(vp41, labels )
     calculate_accuracy(vp51, labels )
     calculate_accuracy(vp61, labels )
    
    
   
    vp12 <- assign_groups(t(vp1$U),threshold = 0.6)
    vp22 <- assign_groups(t(vp2$U),threshold = 0.6)
    vp32 <- assign_groups(t(vp3$U),threshold = 0.6)
    vp42 <- assign_groups(t(vp4$U),threshold = 0.6)
    vp52 <- assign_groups(t(vp5$U),threshold = 0.6)
    vp62 <- assign_groups(t(vp6$U),threshold = 0.6)
    
    
     calculate_accuracy(vp12, labels )
     calculate_accuracy(vp22, labels )
     calculate_accuracy(vp32, labels )
     calculate_accuracy(vp42, labels )
    calculate_accuracy(vp52, labels )
    calculate_accuracy(vp62, labels )
    
    
    
    
    
    
    
    # here we perform the FCMD clustering method, using the 'dtw_basic" to make it faster 
    fcmd1 <- tsclust(series  = ts, type = "fuzzy",  k  = 2,  distance  = "dtw_basic", centroid  = "fcmd",              
      trace     = FALSE , fuzzy_control( fuzziness = 1.2))
    fcmd2 <- tsclust(series  = ts, type = "fuzzy",  k  = 2,  distance  = "dtw_basic", centroid  = "fcmd",              
                     trace     = FALSE , fuzzy_control( fuzziness = 1.4))
    fcmd3 <- tsclust(series  = ts, type = "fuzzy",  k  = 2,  distance  = "dtw_basic", centroid  = "fcmd",              
                     trace     = FALSE , fuzzy_control( fuzziness = 1.6))
    fcmd4 <- tsclust(series  = ts, type = "fuzzy",  k  = 2,  distance  = "dtw_basic", centroid  = "fcmd",              
                     trace     = FALSE , fuzzy_control( fuzziness = 1.8))
    fcmd5 <- tsclust(series  = ts, type = "fuzzy",  k  = 2,  distance  = "dtw_basic", centroid  = "fcmd",              
                     trace     = FALSE , fuzzy_control( fuzziness = 2))
    fcmd6 <- tsclust(series  = ts, type = "fuzzy",  k  = 2,  distance  = "dtw_basic", centroid  = "fcmd",              
                     trace     = FALSE , fuzzy_control( fuzziness = 2.2))
    
    
    
    fcmd11 <- assign_groups(fcmd1@fcluster)
    fcmd21 <- assign_groups(fcmd2@fcluster )
    fcmd31 <- assign_groups(fcmd3@fcluster)
    fcmd41 <- assign_groups(fcmd4@fcluster)
    fcmd51 <- assign_groups(fcmd5@fcluster)
    fcmd61 <- assign_groups(fcmd6@fcluster)
    
    
     calculate_accuracy(fcmd11, labels )
     calculate_accuracy(fcmd21, labels )
     calculate_accuracy(fcmd31, labels )
     calculate_accuracy(fcmd41, labels )
     calculate_accuracy(fcmd51, labels )
     calculate_accuracy(fcmd61, labels )
    
    
    
    fcmd12 <- assign_groups(fcmd1@fcluster,0.6)
    fcmd22 <- assign_groups(fcmd2@fcluster,0.6)
    fcmd32 <- assign_groups(fcmd3@fcluster,0.6)
    fcmd42 <- assign_groups(fcmd4@fcluster,0.6)
    fcmd52 <- assign_groups(fcmd5@fcluster,0.6)
    fcmd62 <- assign_groups(fcmd6@fcluster,0.6)
    
    
     calculate_accuracy(fcmd12, labels )
     calculate_accuracy(fcmd22, labels )
     calculate_accuracy(fcmd32, labels )
     calculate_accuracy(fcmd42, labels )
     calculate_accuracy(fcmd52, labels )
     calculate_accuracy(fcmd62, labels )
    
    

  
}

















