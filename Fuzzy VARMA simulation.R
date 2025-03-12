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











# The simulation part 
labels <- c(rep(1,10),rep(2,10),NA,NA)

rep <- 100
 
ac <- matrix(0,nrow = 3, ncol = rep)
arob1 <- matrix(0,nrow = 3, ncol = rep)
ac1 <- matrix(0,nrow = 3, ncol = rep)
autoselectac <- matrix(0,nrow = 3, ncol = rep)
auto_iniac <- matrix(0,nrow = 3, ncol = rep)
rob6 <- matrix(0,nrow = 3, ncol = rep)

ac1 <- matrix(0,nrow = 3, ncol = rep)
arob11 <- matrix(0,nrow = 3, ncol = rep)
ac11 <- matrix(0,nrow = 3, ncol = rep)
autoselectac1 <- matrix(0,nrow = 3, ncol = rep)
auto_iniac1 <- matrix(0,nrow = 3, ncol = rep)
rob61 <- matrix(0,nrow = 3, ncol = rep)



vp11_ta <- matrix(0,nrow = 3, ncol = rep)
vp21_ta <- matrix(0,nrow = 3, ncol = rep)
vp31_ta <- matrix(0,nrow = 3, ncol = rep)
vp41_ta <- matrix(0,nrow = 3, ncol = rep)
vp51_ta <- matrix(0,nrow = 3, ncol = rep)
vp61_ta <- matrix(0,nrow = 3, ncol = rep)



vp12_ta <- matrix(0,nrow = 3, ncol = rep)
vp22_ta <- matrix(0,nrow = 3, ncol = rep)
vp32_ta <- matrix(0,nrow = 3, ncol = rep)
vp42_ta <- matrix(0,nrow = 3, ncol = rep)
vp52_ta <- matrix(0,nrow = 3, ncol = rep)
vp62_ta <- matrix(0,nrow = 3, ncol = rep)


fcmd11_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd21_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd31_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd41_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd51_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd61_ta <- matrix(0,nrow = 3, ncol = rep)


fcmd12_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd22_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd32_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd42_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd52_ta <- matrix(0,nrow = 3, ncol = rep)
fcmd62_ta <- matrix(0,nrow = 3, ncol = rep)


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
  
  for (l in 61:65) {
 #   lengths1 <- sample(200:600, 10)
 #   lengths2 <- sample(200:600, 10)
 #    lengths3 <- sample(200:600, 2)
      
     lengths1 <- rep(200,10)
       lengths2 <- rep(200,10)
     lengths3 <- rep(200,2)
    
    
    data_list_var1 <- lapply(lengths1, function(len) {
      MTS::VARMAsim(len, arlags = 1, phi = m * 0.8 , sigma = Sigma_VAR )$series
    })
    
    
    data_list_var2 <- lapply(lengths2, function(len) {
      MTS::VARMAsim(len, malags = 1, theta = n * 2  , sigma = Sigma_VMA )$series
    })
    
 
    
    data_list_var3 <- lapply(lengths3, function(len) {
     # phi_VARMA <- m * 0.5 + 0.2 * n  # Add interaction and scale
     # theta_VARMA <- n * 1.5 + 0.2 * m
     # Sigma_VARMA <- diag(e) * 0.3  # Adjust noise
      phi_VARMA <- m * 0.8   # Add interaction and scale
      theta_VARMA <- n * 2
      Sigma_VARMA <- diag(e) * 0.3  # Adjust noise
      ts <- MTS::VARMAsim(len, arlags = 1, malags = 1, phi = phi_VARMA, theta = theta_VARMA, sigma = Sigma_VARMA)$series
      ts + matrix(rnorm(length(ts), 0, 0.5), ncol = e)  # Add noise
    })
    
    ts <- c(data_list_var1, data_list_var2, data_list_var3)
    # ts <- c(data_list_var1, data_list_var2, data_list_var3,data_list_var4)
    
    
   # auto <- cpca(ts,2)
    roblag1 <- fcpca(ts,2,1.2,replicates = 1)
    roblag2 <- fcpca(ts,2,1.4,replicates = 1)
    roblag3 <- fcpca(ts,2,1.6,replicates = 1)
     roblag4 <- fcpca(ts,2,1.8,replicates = 1)
     roblag5 <- fcpca(ts,2,2,replicates = 1)
    roblag6 <- fcpca(ts,2,2.2,replicates = 1)
    # rob <-  robcpca(ts,2)
    
    roblag10 <- assign_groups(roblag1$membership_matrix)
    roblag20 <- assign_groups(roblag2$membership_matrix)
    roblag30 <- assign_groups(roblag3$membership_matrix)
    roblag40 <- assign_groups(roblag4$membership_matrix)
    roblag50 <- assign_groups(roblag5$membership_matrix)
    roblag60 <- assign_groups(roblag6$membership_matrix)
    
 
    ac[j,l] <-  calculate_accuracy(roblag10, labels )
    arob1[j,l] <- calculate_accuracy(roblag20, labels )
    ac1[j,l] <-  calculate_accuracy(roblag30, labels )
    autoselectac[j,l] <-  calculate_accuracy(roblag40, labels )
    auto_iniac[j,l] <-  calculate_accuracy(roblag50, labels )
    rob6[j,l] <-  calculate_accuracy(roblag60, labels )
    
    
    roblag11 <- assign_groups(roblag1$membership_matrix,threshold = 0.6)
    roblag21 <- assign_groups(roblag2$membership_matrix,threshold = 0.6)
    roblag31 <- assign_groups(roblag3$membership_matrix,threshold = 0.6)
    roblag41 <- assign_groups(roblag4$membership_matrix,threshold = 0.6)
    roblag51 <- assign_groups(roblag5$membership_matrix,threshold = 0.6)
    roblag61 <- assign_groups(roblag6$membership_matrix,threshold = 0.6)
    
    ac1[j,l] <-  calculate_accuracy(roblag11 , labels )
    arob11[j,l] <- calculate_accuracy(roblag21, labels )
    ac11[j,l] <-  calculate_accuracy(roblag31, labels )
    autoselectac1[j,l] <-  calculate_accuracy(roblag41, labels )
    auto_iniac1[j,l] <-  calculate_accuracy(roblag51, labels )
    rob61[j,l] <-  calculate_accuracy(roblag61, labels )
    
    
    
    
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
    
    
    vp11_ta[j,l] <- calculate_accuracy(vp11, labels )
    vp21_ta[j,l] <- calculate_accuracy(vp21, labels )
    vp31_ta[j,l] <- calculate_accuracy(vp31, labels )
    vp41_ta[j,l] <- calculate_accuracy(vp41, labels )
    vp51_ta[j,l] <- calculate_accuracy(vp51, labels )
    vp61_ta[j,l] <- calculate_accuracy(vp61, labels )
    
    
   
    vp12 <- assign_groups(t(vp1$U),threshold = 0.6)
    vp22 <- assign_groups(t(vp2$U),threshold = 0.6)
    vp32 <- assign_groups(t(vp3$U),threshold = 0.6)
    vp42 <- assign_groups(t(vp4$U),threshold = 0.6)
    vp52 <- assign_groups(t(vp5$U),threshold = 0.6)
    vp62 <- assign_groups(t(vp6$U),threshold = 0.6)
    
    
    vp12_ta[j,l] <- calculate_accuracy(vp12, labels )
    vp22_ta[j,l] <- calculate_accuracy(vp22, labels )
    vp32_ta[j,l] <- calculate_accuracy(vp32, labels )
    vp42_ta[j,l] <- calculate_accuracy(vp42, labels )
    vp52_ta[j,l] <- calculate_accuracy(vp52, labels )
    vp62_ta[j,l] <- calculate_accuracy(vp62, labels )
    
    
    
    
    
    
    
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
    
    
    fcmd11_ta[j,l] <- calculate_accuracy(fcmd11, labels )
    fcmd21_ta[j,l] <- calculate_accuracy(fcmd21, labels )
    fcmd31_ta[j,l] <- calculate_accuracy(fcmd31, labels )
    fcmd41_ta[j,l] <- calculate_accuracy(fcmd41, labels )
    fcmd51_ta[j,l] <- calculate_accuracy(fcmd51, labels )
    fcmd61_ta[j,l] <- calculate_accuracy(fcmd61, labels )
    
    
    
    fcmd12 <- assign_groups(fcmd1@fcluster,0.6)
    fcmd22 <- assign_groups(fcmd2@fcluster,0.6)
    fcmd32 <- assign_groups(fcmd3@fcluster,0.6)
    fcmd42 <- assign_groups(fcmd4@fcluster,0.6)
    fcmd52 <- assign_groups(fcmd5@fcluster,0.6)
    fcmd62 <- assign_groups(fcmd6@fcluster,0.6)
    
    
    fcmd12_ta[j,l] <- calculate_accuracy(fcmd12, labels )
    fcmd22_ta[j,l] <- calculate_accuracy(fcmd22, labels )
    fcmd32_ta[j,l] <- calculate_accuracy(fcmd32, labels )
    fcmd42_ta[j,l] <- calculate_accuracy(fcmd42, labels )
    fcmd52_ta[j,l] <- calculate_accuracy(fcmd52, labels )
    fcmd62_ta[j,l] <- calculate_accuracy(fcmd62, labels )
    
    
    
    print(l)
    
  }
  
}


###############################################################################
# 1) Create a named list of all accuracy matrices
###############################################################################
accuracy_matrices <- list(
  # Newly added ones
  ac               = ac,
  arob1            = arob1,
  ac1              = ac1,
  autoselectac     = autoselectac,
  auto_iniac       = auto_iniac,
  rob6             = rob6,
  
  ac1_roblag       = ac1,         # NOTE: you already have ac1 above
  arob11           = arob11,
  ac11             = ac11,
  autoselectac1    = autoselectac1,
  auto_iniac1      = auto_iniac1,
  rob61            = rob61,
  
  # VPCA-based ones
  vp11_ta          = vp31_ta,
  vp21_ta          = vp31_ta,
  vp31_ta          = vp31_ta,
  vp41_ta          = vp41_ta,
  vp51_ta          = vp51_ta,
  vp61_ta          = vp61_ta,
  
  
  vp12_ta          = vp12_ta,
  vp22_ta          = vp22_ta,
  vp32_ta          = vp32_ta,
  vp42_ta          = vp42_ta,
  vp52_ta          = vp52_ta,
  vp62_ta          = vp62_ta,
  
  # FCMD-based ones
  fcmd11_ta        = fcmd11_ta,
  fcmd21_ta        = fcmd21_ta,
  fcmd31_ta        = fcmd31_ta,
  fcmd41_ta        = fcmd41_ta,
  fcmd51_ta        = fcmd51_ta,
  fcmd61_ta        = fcmd61_ta,
  
  fcmd12_ta        = fcmd12_ta,
  fcmd22_ta        = fcmd22_ta,
  fcmd32_ta        = fcmd32_ta,
  fcmd42_ta        = fcmd42_ta,
  fcmd52_ta        = fcmd52_ta,
  fcmd62_ta        = fcmd62_ta
)

###############################################################################
# 2) Compute row means for the first 5 columns of each matrix
###############################################################################
row_means_list <- lapply(accuracy_matrices, function(mat) {
  # Some matrices may have fewer than 5 columns, so take min(ncol(mat), 5)
  cols_to_take <- min(ncol(mat), 5)
  rowMeans(mat[, 51:55, drop = FALSE], na.rm = TRUE)
})




sd_list <-  lapply(accuracy_matrices, function(mat) {
    # Some matrices may have fewer than 5 columns, so take min(ncol(mat), 5)
    
    sd(mat[, 1:50, drop = FALSE], na.rm = TRUE)
  })


###############################################################################
# 3) (Optional) Convert the list of vectors to a data frame for easier viewing
###############################################################################
row_means_df <- as.data.frame(row_means_list)

# Preview the first few rows
head(row_means_df)

# row_means_df now has one column per matrix, each containing the row means 
# of the first 5 columns for that matrix.


sd_df_1 <- as.data.frame(sd_list)

head(sd_df_1)



fcpca_times <- numeric(6)  # to store execution times
fcpca_m_values <- c(1.2, 1.4, 1.6, 1.8, 2, 2.2)

for (i in seq_along(fcpca_m_values)) {
  start_time <- Sys.time()
  # roblagX <- fcpca(ts, 2, m_value, replicates = 1)
  fcpca(ts, 2, fcpca_m_values[i], replicates = 1)
  end_time <- Sys.time()
  
  fcpca_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
}

# 2) vpca_clustering calls (vp1 to vp6)

vpca_times <- numeric(6)
vpca_m_values <- c(1.2, 1.4, 1.6, 1.8, 2, 2.2)

for (i in seq_along(vpca_m_values)) {
  start_time <- Sys.time()
  # vpX <- vpca_clustering(ts, 2, m = m_value, crisp = FALSE)
  vpca_clustering(ts, 2, m = vpca_m_values[i], crisp = FALSE)
  end_time <- Sys.time()
  
  vpca_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
}

# 3) tsclust calls (fcmd1 to fcmd6)
tsclust_times <- numeric(6)
tsclust_m_values <- c(1.2, 1.4, 1.6, 1.8, 2, 2.2)

for (i in seq_along(tsclust_m_values)) {
  start_time <- Sys.time()
  # fcmdX <- tsclust(series = ts, type = "fuzzy", k = 2, 
  #                  distance = "dtw_basic", centroid = "fcmd",
  #                  trace = FALSE, fuzzy_control(fuzziness = m_value))
  tsclust(
    series = ts,
    type = "fuzzy",
    k = 2,
    distance = "dtw",
    centroid = "fcmd",
    trace = FALSE,
    fuzzy_control(fuzziness = tsclust_m_values[i])
  )
  end_time <- Sys.time()
  
  tsclust_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
}


###############################################################################
# compare execution times
###############################################################################
cat("fcpca times (in seconds):       ", fcpca_times, "\n")
cat("vpca_clustering times (seconds):", vpca_times,  "\n")
cat("tsclust times (seconds):        ", tsclust_times,"\n")










# realization of 2 dimension 




library(MTS)

# --- 1) Generate one realization of each process (example) ---
e <- 2  # dimension
set.seed(1)
m <- matrix(c(0.4, 0.1, -0.1, 0.3), nrow=e, byrow=TRUE)
n <- matrix(c(0.2, 0.05, -0.05, 0.2), nrow=e, byrow=TRUE)

ts_VAR <- VARMAsim(n=100, arlags=1, phi=m*0.8, sigma=diag(e)*0.1)$series
ts_VMA <- VARMAsim(n=100, malags=1, theta=n*2, sigma=diag(e)*0.5)$series
ts_VARMA <- VARMAsim(n=100, arlags=1, malags=1, phi=m*0.8, theta=n*2, sigma=diag(e)*0.3)$series

par(mfrow = c(2,1),          # 2 rows, 1 column
    mar = c(5,5,2,1),        # bigger left margin
    mgp = c(3.2, 1, 0))      # move labels outward slightly

#    We'll do dimension 1 first:
all_vals_dim1 <- c(ts_VAR[,1], ts_VMA[,1], ts_VARMA[,1])
ylim_dim1 <-  c(min(all_vals_dim1),max(all_vals_dim1)+0.8)

# 2) Do the same for dimension 2 (if e=2):
all_vals_dim2 <- c(ts_VAR[,2], ts_VMA[,2], ts_VARMA[,2])
ylim_dim2 <-  c(min(all_vals_dim2),max(all_vals_dim1)+0.8)

par(mfrow = c(2,1))

# Dimension 1
plot(ts_VAR[,1], type="l", col="seagreen4", lty=1,
     xlab="Time", ylab="Value",
     main="Realization 1",
     ylim = ylim_dim1)  # Use the combined range
lines(ts_VMA[,1],   col="orange4",  lty=1)
lines(ts_VARMA[,1], col="palevioletred3", lty=2)
legend("topleft", legend=c("VAR(1)","VMA(1)","VARMA(1,1)"),
       col=c("seagreen4","orange4","palevioletred3"), lty=c(1,1,2), bty="n")

# Dimension 2
plot(ts_VAR[,2], type="l", col="seagreen4", lty=1,
     xlab="Time", ylab="Value",
     main="Realization 2",
     ylim = ylim_dim2)  # Use the combined range
lines(ts_VMA[,2],   col="orange4",  lty=1)
lines(ts_VARMA[,2], col="palevioletred3", lty=2)
legend("topleft", legend=c("VAR(1)","VMA(1)","VARMA(1,1)"),
       col=c("seagreen4","orange4","palevioletred3"), lty=c(1,1,2), bty="n")



















par(mfrow = c(2,2), mar = c(4,4,2,1))

acf(ts_VAR[,1], main = "ACF: VAR(1), dim1")
pacf(ts_VAR[,1], main = "PACF: VAR(1), dim1")

acf(ts_VAR[,2], main = "ACF: VAR(1), dim2")
pacf(ts_VAR[,2], main = "PACF: VAR(1), dim2")


par(mfrow = c(2,2))
acf(ts_VMA[,1], main = "ACF: VMA(1), dim1")
pacf(ts_VMA[,1], main = "PACF: VMA(1), dim1")
acf(ts_VMA[,2], main = "ACF: VMA(1), dim2")
pacf(ts_VMA[,2], main = "PACF: VMA(1), dim2")

# ACF/PACF for VARMA(1,1)
par(mfrow = c(2,2))
acf(ts_VARMA[,1], main = "ACF: VARMA(1,1), dim1")
pacf(ts_VARMA[,1], main = "PACF: VARMA(1,1), dim1")
acf(ts_VARMA[,2], main = "ACF: VARMA(1,1), dim2")
pacf(ts_VARMA[,2], main = "PACF: VARMA(1,1), dim2")


















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
  
  for (l in 68:63) {
    lengths1 <- sample(200:600, 10)
    lengths2 <- sample(200:600, 10)
    lengths3 <- sample(200:600, 2)
    
  
    
    data_list_var1 <- lapply(lengths1, function(len) {
      MTS::VARMAsim(len, arlags = 1, phi = m * 0.8 , sigma = Sigma_VAR )$series
    })
    
    
    data_list_var2 <- lapply(lengths2, function(len) {
      MTS::VARMAsim(len, malags = 1, theta = n * 2  , sigma = Sigma_VMA )$series
    })
    
    
    
    data_list_var3 <- lapply(lengths3, function(len) {
      # phi_VARMA <- m * 0.5 + 0.2 * n  # Add interaction and scale
      # theta_VARMA <- n * 1.5 + 0.2 * m
      # Sigma_VARMA <- diag(e) * 0.3  # Adjust noise
      phi_VARMA <- m * 0.8   # Add interaction and scale
      theta_VARMA <- n * 2
      Sigma_VARMA <- diag(e) * 0.3  # Adjust noise
      ts <- MTS::VARMAsim(len, arlags = 1, malags = 1, phi = phi_VARMA, theta = theta_VARMA, sigma = Sigma_VARMA)$series
      ts + matrix(rnorm(length(ts), 0, 0.5), ncol = e)  # Add noise
    })
    
    ts <- c(data_list_var1, data_list_var2, data_list_var3)
    # ts <- c(data_list_var1, data_list_var2, data_list_var3,data_list_var4)
    
    
    # auto <- cpca(ts,2)
    roblag1 <- fcpca(ts,2,1.2,replicates = 1)
    roblag2 <- fcpca(ts,2,1.4,replicates = 1)
    roblag3 <- fcpca(ts,2,1.6,replicates = 1)
    roblag4 <- fcpca(ts,2,1.8,replicates = 1)
    roblag5 <- fcpca(ts,2,2,replicates = 1)
    roblag6 <- fcpca(ts,2,2.2,replicates = 1)
    # rob <-  robcpca(ts,2)
    
    roblag10 <- assign_groups(roblag1$membership_matrix)
    roblag20 <- assign_groups(roblag2$membership_matrix)
    roblag30 <- assign_groups(roblag3$membership_matrix)
    roblag40 <- assign_groups(roblag4$membership_matrix)
    roblag50 <- assign_groups(roblag5$membership_matrix)
    roblag60 <- assign_groups(roblag6$membership_matrix)
    
    
    ac[j,l] <-  calculate_accuracy(roblag10, labels )
    arob1[j,l] <- calculate_accuracy(roblag20, labels )
    ac1[j,l] <-  calculate_accuracy(roblag30, labels )
    autoselectac[j,l] <-  calculate_accuracy(roblag40, labels )
    auto_iniac[j,l] <-  calculate_accuracy(roblag50, labels )
    rob6[j,l] <-  calculate_accuracy(roblag60, labels )
    
    
    roblag11 <- assign_groups(roblag1$membership_matrix,threshold = 0.6)
    roblag21 <- assign_groups(roblag2$membership_matrix,threshold = 0.6)
    roblag31 <- assign_groups(roblag3$membership_matrix,threshold = 0.6)
    roblag41 <- assign_groups(roblag4$membership_matrix,threshold = 0.6)
    roblag51 <- assign_groups(roblag5$membership_matrix,threshold = 0.6)
    roblag61 <- assign_groups(roblag6$membership_matrix,threshold = 0.6)
    
    ac1[j,l] <-  calculate_accuracy(roblag11 , labels )
    arob11[j,l] <- calculate_accuracy(roblag21, labels )
    ac11[j,l] <-  calculate_accuracy(roblag31, labels )
    autoselectac1[j,l] <-  calculate_accuracy(roblag41, labels )
    auto_iniac1[j,l] <-  calculate_accuracy(roblag51, labels )
    rob61[j,l] <-  calculate_accuracy(roblag61, labels )
    
    
    
    
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
    
    
    fcmd11_ta[j,l] <- calculate_accuracy(fcmd11, labels )
    fcmd21_ta[j,l] <- calculate_accuracy(fcmd21, labels )
    fcmd31_ta[j,l] <- calculate_accuracy(fcmd31, labels )
    fcmd41_ta[j,l] <- calculate_accuracy(fcmd41, labels )
    fcmd51_ta[j,l] <- calculate_accuracy(fcmd51, labels )
    fcmd61_ta[j,l] <- calculate_accuracy(fcmd61, labels )
    
    
    
    fcmd12 <- assign_groups(fcmd1@fcluster,0.6)
    fcmd22 <- assign_groups(fcmd2@fcluster,0.6)
    fcmd32 <- assign_groups(fcmd3@fcluster,0.6)
    fcmd42 <- assign_groups(fcmd4@fcluster,0.6)
    fcmd52 <- assign_groups(fcmd5@fcluster,0.6)
    fcmd62 <- assign_groups(fcmd6@fcluster,0.6)
    
    
    fcmd12_ta[j,l] <- calculate_accuracy(fcmd12, labels )
    fcmd22_ta[j,l] <- calculate_accuracy(fcmd22, labels )
    fcmd32_ta[j,l] <- calculate_accuracy(fcmd32, labels )
    fcmd42_ta[j,l] <- calculate_accuracy(fcmd42, labels )
    fcmd52_ta[j,l] <- calculate_accuracy(fcmd52, labels )
    fcmd62_ta[j,l] <- calculate_accuracy(fcmd62, labels )
    
    
    
    print(l)
    
  }
  
}





accuracy_matrices <- list(
  # Newly added ones
  ac               = ac,
  arob1            = arob1,
  ac1              = ac1,
  autoselectac     = autoselectac,
  auto_iniac       = auto_iniac,
  rob6             = rob6,
  
  ac1_roblag       = ac1,         # NOTE: you already have ac1 above
  arob11           = arob11,
  ac11             = ac11,
  autoselectac1    = autoselectac1,
  auto_iniac1      = auto_iniac1,
  rob61            = rob61,
  
  # VPCA-based ones
  vp11_ta          = vp31_ta,
  vp21_ta          = vp31_ta,
  vp31_ta          = vp31_ta,
  vp41_ta          = vp41_ta,
  vp51_ta          = vp51_ta,
  vp61_ta          = vp61_ta,
  
  
  vp12_ta          = vp12_ta,
  vp22_ta          = vp22_ta,
  vp32_ta          = vp32_ta,
  vp42_ta          = vp42_ta,
  vp52_ta          = vp52_ta,
  vp62_ta          = vp62_ta,
  
  # FCMD-based ones
  fcmd11_ta        = fcmd11_ta,
  fcmd21_ta        = fcmd21_ta,
  fcmd31_ta        = fcmd31_ta,
  fcmd41_ta        = fcmd41_ta,
  fcmd51_ta        = fcmd51_ta,
  fcmd61_ta        = fcmd61_ta,
  
  fcmd12_ta        = fcmd12_ta,
  fcmd22_ta        = fcmd22_ta,
  fcmd32_ta        = fcmd32_ta,
  fcmd42_ta        = fcmd42_ta,
  fcmd52_ta        = fcmd52_ta,
  fcmd62_ta        = fcmd62_ta
)

###############################################################################
# 2) Compute row means for the first 5 columns of each matrix
###############################################################################
row_means_list <- lapply(accuracy_matrices, function(mat) {
  # Some matrices may have fewer than 5 columns, so take min(ncol(mat), 5)
  cols_to_take <- min(ncol(mat), 5)
  rowMeans(mat[, 61:70, drop = FALSE], na.rm = TRUE)
})

as.data.frame(row_means_list)
