library(R.matlab)
data <- readMat("/Users/maz0b/Desktop/2425FALL/fuzzy clustering FCPCA/real data/dataset.mat")
EEG_data <- data$EEGsample
# Suppose 'EEG_data' is a 3D array with dimensions [2022, 30, 384]
# Initialize an empty list of the correct length
EEGsample <- vector("list", length = 2022)

# Loop through each sample (1 to 2022)
for (i in seq_len(2022)) {
  # Extract the i-th sample from EEG_data
  # EEG_data[i, , ] is a 30 x 384 matrix
  # We transpose it to get 384 x 30
  EEGsample[[i]] <- t(EEG_data[i, , ])
}

# Now 'EEGsample' is a list of length 2022
# Each element is a 384 x 30 matrix
str(EEGsample)

label <- data$substate


one_sub <- data$subindex

# Here, we only provide the fuzzy clustering for subject 11, but you can just change the subject_of_interest below from 1 to 11 to scrutinize other subjects' results.

subject_of_interest <- 11

# 1. Identify which samples belong to this subject
idx <- which(data$subindex == subject_of_interest)

# 2. Extract those samples from the EEGsample list
EEGsample_subject <- EEGsample[idx]

# 3. Extract the corresponding labels (if needed)
label_subject <- data$substate[idx]

# (Optional) Check the results
length(EEGsample_subject)
length(label_subject)
unique(data$subindex[idx])  # Should only show the subject_of_interest

a <- fcpca(EEGsample_subject,2,1.1,replicates = 1)
true_accuracy(a$hard_cluster,label_subject)  # 0.9336283
a$S 



e <- tsclust(EEGsample_subject, type      = "fuzzy",               # fuzzy clustering
             k         = 2,                     # number of clusters (choose based on your goal)
             distance  = "gak",           # distance measure; can also try "sdtw", "gak", etc.
             centroid  = "fcmd",               # centroid type (often good for DTW-based methods)
             trace     = FALSE ,                 # suppress progress messages
             fuzzy_control(
               # notice the vector
               fuzziness = 1.2
             )
)


e@fcluster


true_accuracy(max.col((e@fcluster)),label_subject)


e <- vpca_clustering(EEGsample_subject,2,1.1,crisp = FALSE)
true_accuracy(e$U,label_subject)





# The way to generate the plots m versus CVIs and RIs

# Define the range of m values you want to test
m_values <- seq(1.1, 2.2, by = 0.1)

# Initialize vectors to hold RI and CVI
RI_vec <- numeric(length(m_values))
CVI_vec <- numeric(length(m_values))

# Loop over each m value
for (i in seq_along(m_values)) {
  m_val <- m_values[i]
  
  # Run fcpca() with fuzziness parameter m_val
  # 'a' will typically contain clustering results, membership matrix, etc.
  a <- fcpca(EEGsample_subject, 2, m_val)
  
  # Compute RI (hard cluster accuracy)
  RI_vec[i] <- true_accuracy(a$hard_cluster, label_subject)
  
  # Extract CVI (S) from the result
  CVI_vec[i] <- a$S
}

# 1. Combine results into a data frame
results_df <- data.frame(
  m   = m_values,
  CVI = CVI_vec,
  RI  = RI_vec
)

# 2. Log-transform CVI (for plotting on a similar scale)
results_df$logCVI <- log(results_df$CVI)

# 3. Normalize logCVI and RI to [0, 1] for plotting on the same y-axis
logCVI_min <- min(results_df$logCVI)
logCVI_max <- max(results_df$logCVI)
results_df$norm_logCVI <- (results_df$logCVI - logCVI_min) / (logCVI_max - logCVI_min)

RI_min <- min(results_df$RI)
RI_max <- max(results_df$RI)
results_df$norm_RI <- (results_df$RI - RI_min) / (RI_max - RI_min)

# 4. Plot both metrics in the same figure
library(ggplot2)

ggplot(results_df, aes(x = m)) +
  
  # --- log(CVI) trend ---
  geom_line(aes(y = norm_logCVI, color = "log(CVI)"), size = 1) +
  geom_point(aes(y = norm_logCVI, color = "log(CVI)"), size = 3) +
  # Annotate actual CVI values above each point
  geom_text(
    aes(y = norm_logCVI, label = round(CVI, 2), color = "log(CVI)"),
    vjust = -1, size = 4
  ) +
  
  # --- RI trend ---
  geom_line(aes(y = norm_RI, color = "RI"), size = 1) +
  geom_point(aes(y = norm_RI, color = "RI"), size = 3) +
  # Annotate actual RI values below each point
  geom_text(
    aes(y = norm_RI, label = round(RI, 2), color = "RI"),
    vjust = 1.5, size = 4
  ) +
  
  # --- Labels, theme, and legend ---
  labs(
    x = "Fuzziness Parameter (m)",
    y = NULL,
    color = "Metric"  # Legend title
  ) +
  theme_minimal(base_size = 14) +
  
  # Remove y-axis text, ticks, and grid lines
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  # Ensure x-axis has ticks at each m-value
  scale_x_continuous(breaks = results_df$m) +
  
  # Use the same color palette as before
  scale_color_manual(
    values = c(
      "log(CVI)" = "#0072B2",  # Blue
      "RI"       = "#D55E00"   # Orange
    )
  )





# plots for the membership matrix 
library(tidyr)
library(ggplot2)

# 1. Convert your membership matrix to a data frame
U_df <- as.data.frame(a$membership_matrix)

# 2. Create a 'sample' identifier
U_df$sample <- factor(seq_len(nrow(U_df)))  # 1 to 226 as factors

# 3. Reshape data to "long" format
U_long <- gather(U_df, key = "group", value = "membership", -sample)

# 4. Rename factor levels (and set order)
#    'drowsy state' is first, 'alert state' is second
U_long$group <- factor(U_long$group,
                       levels = c("V1", "V2"), 
                       labels = c("drowsy state", "alert state")
)

# 5. Plot as stacked bars
ggplot(U_long, aes(x = sample, y = membership, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Membership Degree") +
  
  # Use a custom fill scale: pink for drowsy, cyan for alert
  scale_fill_manual(
    name   = "group",  # Legend title (optional)
    values = c("drowsy state" = "lightpink",
               "alert state"  = "lightcyan3"),
    # The breaks below ensure the legend order is exactly as specified
    breaks = c("drowsy state", "alert state")
  ) +
  
  # Minimal theme + remove additional grid lines / background
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.background  = element_blank(),
    legend.background = element_blank(),
    legend.key        = element_blank()
  )

