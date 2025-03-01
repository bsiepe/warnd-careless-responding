# =====================================================================
# Title:        Functions for careless responding indicators
# Author:       Joran Geeraerts
# Date:         2025-02-24
# Description:  This code contains functions to calculate:
#               Mahalanobis distance, robust PCA orthogonal distance
#               Psychometric synonym and antonym violations
#               SD of response times
# Dependencies: {dplyr}, {tidyverse}, {rrcov}, {Hmisc}
# =====================================================================

# ---- SETUP ----
# Clear workspace
rm(list = ls())

# Set working directory
setwd("C:\\Users\\Joran\\OneDrive - Vrije Universiteit Brussel\\SIDE PROJECTS\\Careless responding - WARN-D\\Data analyses")

# Libraries
library(dplyr)
library(tidyr)
library(rrcov) # pcaHubert()
library(Hmisc) # rcorr()

# Read data
filepath <- 'Data\\WARND_stage2_c1234_trainingset_before_15min_v2024_07_16.csv'
data <- read.csv(filepath)

# ---- VARIABLES USED TO CREATE AND TEST THE FUNCTIONS ----
## Items to calculate indices on
itemsToCalc <- c("sad_d", "stressed_d", "overwhelm_d", "nervous_d", "ruminate_d", "irritable_d", "tired_d", "cheerful_d", "motivated_d", "relaxed_d")

## Response time variables for selection of items
rt_itemsToCalc <- paste0("rt_", itemsToCalc)

## Cut-offs
cor_threshold <- .45

synonym_difference_threshold <- 5
antonym_maxvalue_threshold <- 5

# ---- MAHALANOBIS DISTANCE ----
# Requires dataset and vector with names of items to include
calc_indicator_mahalanobis <- function(data, itemsToCalc) {
  ## Note: Currently we provide missing values for distances that are not able to compute due to, for instance, singular matrices
  ## Calculate Mahalanobis distances for each participant (splits data for each id, applies function to each subset, then binds rows of resulting distances)
  mahalanobis.results <- do.call(rbind, lapply(split(data, data$external_id), function(df) {
    # Select the relevant question columns
    question_data <- df |>
      select(all_of(itemsToCalc)) |>
      mutate(across(everything(), as.numeric)) |> # Ensure numeric data
      na.omit() # Remove rows with missing data
    
    # Calculate Mahalanobis distances manually
    if (nrow(question_data) < 2) {
      # If there are fewer than 2 rows, covariance cannot be computed properly
      distances <- rep(NA, nrow(df))
    } else {
      # Calculate mean vector and covariance matrix
      mu <- colMeans(question_data)
      cov_matrix <- cov(question_data)
      
      # Calculate Mahalanobis distance with error handling
      distances <- tryCatch({
        apply(question_data, 1, function(x) {
          sqrt(t(x - mu) %*% solve(cov_matrix) %*% (x - mu)) # Mahalanobis distance
        })
      }, error = function(e) {
        rep(NA, nrow(question_data))
      })
    }
    
    # Store results
    data.frame(
      external_id = df$external_id, 
      counter = df$counter,
      mahalanobis_dist = c(distances, rep(NA, nrow(df) - length(distances)))
    )
    
  }))
  
  return(mahalanobis.results)
}


# ---- ROBUST PCA ORTHOGONAL DISTANCE ----
# Requires dataset and vector with names of items to include
calc_indicator_rob_PCA_orthogonal_distance <- function(data, itemsToCalc) {
  ## Note: The robust PCA calculates the optimal number of principal components to compute for every participants
  ## I'm currently not sure if the values can be compared between participants, as it could impact the size of the orthogonal distance
  robustpca.results <- do.call(rbind, lapply(split(data, data$external_id), function(df) {
    # Select the relevant question columns and convert to matrix
    question_data <- df |>
      select(all_of(itemsToCalc)) |>
      mutate(across(everything(), as.numeric)) 
    
    # Identify complete rows (non-NA rows)
    complete_rows <- complete.cases(question_data)
    
    # Apply PcaHubert only to complete rows
    question_data_complete <- question_data[complete_rows, ] |>
      data.matrix()
    
    # Initialize distances and kValues with NA, matching original row count
    distances <- rep(NA, nrow(df))
    kValues <- rep(NA, nrow(df))
    
    if (nrow(question_data_complete) > 1) { # Check if enough data remains
      # Use tryCatch to handle potential errors during PcaHubert
      tryCatch({
        output <- PcaHubert(question_data_complete, kmax = 15)
        
        # Fill in distances and kValues only for complete rows
        distances[complete_rows] <- output$od
        kValues[complete_rows] <- output$k
        
      }, error = function(e) {
        message("Error in PcaHubert calculation for external_id: ", unique(df$external_id), " - Replacing with NA")
        # NAs are already in place, so no need to modify distances/kValues
      })
    }
    
    # Store results
    data.frame(
      external_id = df$external_id, 
      counter = df$counter,
      robustpca_dist = distances,
      kValues = kValues
    )
  }))
  
  return(robustpca.results)
}

# ---- PSYCHOMETRIC SYNONYM VIOLATIONS ----
# Requires dataset, vector with names of items to include, threshold for correlation to be considered synonym, 
# and threshold for difference to be considered violation
calc_psychometric_synonym_violations <- function(data, itemsToCalc, cor_threshold, synonym_difference_threshold) {
  ## Select relevant data to calculate correlations on
  question_data <- data |>
    select(external_id, counter, all_of(itemsToCalc)) |>
    mutate(across(all_of(itemsToCalc), as.numeric)) |> # Ensure numeric data
    na.omit() # Remove rows with missing data
  
  ## Calculate correlations
  corr.temp <- rcorr(question_data |> select(3:ncol(.)) |> data.matrix(), type = 'spearman') # Correlation matrix
  
  flattenCorrMatrix <- function(corrmat, pmat) { # Flatten the correlation matrix
    ut <- upper.tri(corrmat)
    data.frame(
      row = rownames(corrmat)[row(corrmat)[ut]],
      column = rownames(corrmat)[col(corrmat)[ut]],
      cor  =(corrmat)[ut],
      p = pmat[ut]
    )
  }
  corr.temp <- flattenCorrMatrix(corr.temp$r, corr.temp$P) |> arrange(cor)
  
  ## Filter correlations that reach cut-off to be considered synonyms
  synonym_pairs <- corr.temp |> 
    filter(cor >= cor_threshold) |>
    rename('word1' = row ,
           'word2' = column) |>
    select(1:2)
  
  ## Calculate violations per row
  psychometric.synonym.count <- matrix(nrow = nrow(question_data), ncol = nrow(synonym_pairs)) # Create empty matrix to fill
  
  for (row in 1:nrow(synonym_pairs)){ # Loop over pairs
    # Select two synonym items from the relevant response data
    df.temp <- question_data |> 
      select(as.character(synonym_pairs[row, 1]), as.character(synonym_pairs[row, 2]))
    
    # Check if the difference between scores on those items reach the threshold
    psychometric.synonym.count[, row] <- abs(df.temp[, 1] - df.temp[, 2]) >= synonym_difference_threshold
  }
  psychometric.synonym.count <- cbind(question_data$external_id, question_data$counter, # Bind id, counter, and synonym violation counts
                                      rowSums(psychometric.synonym.count)) |> 
    as.data.frame() |> 
    rename('external_id' = V1, 'counter' = V2, 'synoViolations' = V3)
  
  return(psychometric.synonym.count)
}

# ---- PSYCHOMETRIC ANTONYM VIOLATIONS ----
# Requires dataset, vector with names of items to include, negative threshold for correlation to be considered antonym, 
# and threshold that both items need to reach to be considered a violation
calc_psychometric_antonym_violations <- function(data, itemsToCalc, cor_threshold, antonym_maxvalue_threshold) {
  ## Select relevant data to calculate correlations on
  question_data <- data |>
    select(external_id, counter, all_of(itemsToCalc)) |>
    mutate(across(all_of(itemsToCalc), as.numeric)) |> # Ensure numeric data
    na.omit() # Remove rows with missing data
  
  ## Calculate correlations
  corr.temp <- rcorr(question_data |> select(3:ncol(.)) |> data.matrix(), type = 'spearman') # Correlation matrix
  
  flattenCorrMatrix <- function(corrmat, pmat) { # Flatten the correlation matrix
    ut <- upper.tri(corrmat)
    data.frame(
      row = rownames(corrmat)[row(corrmat)[ut]],
      column = rownames(corrmat)[col(corrmat)[ut]],
      cor  =(corrmat)[ut],
      p = pmat[ut]
    )
  }
  corr.temp <- flattenCorrMatrix(corr.temp$r, corr.temp$P) |> arrange(cor)
  
  ## Filter correlations that reach cut-off to be considered antonyms
  antonym_pairs <- corr.temp |> 
    filter(cor <= cor_threshold) |>
    rename('word1' = row ,
           'word2' = column) |>
    select(1:2)
  
  ## Calculate violations per row
  psychometric.antonym.count <- matrix(nrow = nrow(question_data), ncol = nrow(antonym_pairs)) # Create empty matrix to fill
  
  for (row in 1:nrow(antonym_pairs)){ # Loop over pairs
    # Select two synonym items from the relevant response data
    df.temp <- question_data |> 
      select(as.character(antonym_pairs[row, 1]), as.character(antonym_pairs[row, 2]))
    
    # Check if both items score above the threshold
    psychometric.antonym.count[, row] <- df.temp[, 1] > antonym_maxvalue_threshold & df.temp[, 2] > antonym_maxvalue_threshold
  }
  
  psychometric.antonym.count <- cbind(question_data$external_id, question_data$counter, # Bind id, counter, and antonym violation counts
                                      rowSums(psychometric.antonym.count)) |> 
    as.data.frame() |> 
    rename('external_id' = V1, 'counter' = V2, 'antoViolations' = V3)
  return(psychometric.antonym.count)
}

# ---- SD OF RESPONSE TIMES ----
# Requires dataset, vector with names of response times of included items
calc_sd_response_times <- function(data, rt_itemsToCalc) {
  # Apply as.POSIXct to the specified columns in the data frame, then make it numeric to interpret it as seconds (needed for differences later)
  data[rt_itemsToCalc] <- lapply(data[rt_itemsToCalc], as.POSIXlt, format = "%Y-%m-%d %H:%M:%OS")
  data[rt_itemsToCalc] <- lapply(data[rt_itemsToCalc], as.numeric)
  
  # Reorder each row based on the response times in rt_itemsToCalc
  # For each row, order the response times
  ordered_rows <- t(apply(data[rt_itemsToCalc], 1, function(x) x[order(x)]))
  
  # Convert to data frame and maintain column names
  ordered_rows_df <- as.data.frame(ordered_rows)
  
  # Calculate time differences between consecutive items for each row
  time_diffs <- t(apply(ordered_rows_df, 1, function(x) diff(x)))
  
  # Convert time_diffs into a data frame
  time_diffs_df <- as.data.frame(time_diffs)
  
  # Calculate the standard deviation of time differences for each respondent
  sd_time_diff <- apply(time_diffs_df, 1, sd)
  
  # Return the final result
  return(sd_time_diff)
}