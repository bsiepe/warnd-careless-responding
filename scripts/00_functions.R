# -------------------------------------------------------------------------
# Indices functions -------------------------------------------------------
# -------------------------------------------------------------------------


# Standard deviations within-assessment -----------------------------------
# calc_within_assessment_sd <- function(data,
#                                  items = ema_items,
#                                  id_col = user_id,
#                                  answer_col = answer_id) {
#   data |>
#     dplyr::filter(item %in% items) |>
#     dplyr::group_by({{id_col}}, counter) |>
#     dplyr::summarise(sd = sd({{answer_col}}, na.rm = TRUE)) |>
#     dplyr::ungroup()
# }

calc_within_assessment_sd <- function(data,
                                      items = ema_items,
                                      id_col = user_id) {
  data |>
    dplyr::rowwise() |>
    dplyr::mutate(assessment_sd = sd(dplyr::c_across(dplyr::all_of(items)), na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::select({{ id_col }}, counter, assessment_sd)
}



# Percentage of items at the mode -----------------------------------------
calc_mode <- function(x) {
  unique_x <- unique(x)
  unique_x[which.max(tabulate(match(x, unique_x)))]
}


# calc_mode_percentage <- function(data,
#                             items = ema_items,
#                             id_col = user_id,
#                             answer_col = answer_id) {
#   data |>
#     dplyr::filter(item %in% items) |>
#     dplyr::group_by({{id_col}}, counter) |>
#     dplyr::mutate(mode = calc_mode({{answer_col}})) |>
#     dplyr::mutate(is_mode = ifelse({{answer_col}} == mode, 1, 0)) |>
#     dplyr::summarise(mode_pct = mean(is_mode, na.rm = TRUE)) |>
#     dplyr::ungroup()
# }

calc_mode_percentage <- function(data,
                                 items = ema_items,
                                 id_col = user_id) {
  data |>
    dplyr::rowwise() |>
    dplyr::mutate(mode = calc_mode(dplyr::c_across(dplyr::all_of(items))),
           is_mode = mean(dplyr::c_across(dplyr::all_of(items)) == mode, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::select({{ id_col }}, counter, mode_pct = is_mode) |> 
    dplyr::mutate(mode_pct = if_else(is.nan(mode_pct), NA, mode_pct))
}



# Time per item -----------------------------------------------------------





# Mahalanobis Distance ----------------------------------------------------
calc_indicator_mahalanobis <- function(data, items) {
  # Note: Currently we provide missing values for distances that are not able to compute due to, for instance, singular matrices
  # Calculate Mahalanobis distances for each participant
  # (splits data for each id, applies function to each subset, then binds rows of resulting distances)
  
  mahalanobis.results <- do.call(rbind, lapply(split(data, data$external_id), function(df) {
    # Select the relevant question columns
    question_data <- df |>
      dplyr::select(external_id, counter, all_of(items)) |>
      dplyr::mutate(dplyr::across(all_of(items), as.numeric)) |>
      stats::na.omit()
    
    # Calculate Mahalanobis distances manually
    if (nrow(question_data) < 2) {
      # If there are fewer than 2 rows, covariance cannot be computed properly
      distances <- rep(NA, nrow(question_data))
    } else {
      # Calculate mean vector and covariance matrix
      mu <- colMeans(question_data[,items])
      cov_matrix <- cov(question_data[,items])
      
      # Calculate Mahalanobis distance with error handling
      distances <- tryCatch({
        apply(question_data[,items], 1, function(x) {
          sqrt(t(x - mu) %*% solve(cov_matrix) %*% (x - mu)) # Mahalanobis distance
        })
      }, error = function(e) {
        rep(NA, nrow(question_data))
      })
    }
    
    # Store results
    data.frame(
      external_id = question_data$external_id, 
      counter = as.integer(question_data$counter),
      mahalanobis_dist = distances
    )
    
  }))
  
  return(mahalanobis.results)
}



# Robust PCA --------------------------------------------------------------
calc_indicator_rob_PCA_orthogonal_distance <- function(data, items) {
  ## Note: The robust PCA calculates the optimal number of principal components to compute for every participants
  ## I'm currently not sure if the values can be compared between participants, as it could impact the size of the orthogonal distance
  robustpca.results <- do.call(rbind, lapply(split(data, data$external_id), function(df) {
    # Select the relevant question columns and convert to matrix
    question_data <- df |>
      dplyr::select(external_id, counter, dplyr::all_of(items)) |>
      dplyr::mutate(dplyr::across(dplyr::all_of(items), as.numeric)) |> 
      na.omit() |> 
      data.matrix()
    
    # Initialize distances and kValues with NA, matching row count of data with NA's removed
    distances <- rep(NA, nrow(question_data))
    kValues <- rep(NA, nrow(question_data))
    

    if (nrow(question_data) > 1) { # Check if enough data remains
      # Use tryCatch to handle potential errors during PcaHubert
      tryCatch({
        output <- PcaHubert(question_data, kmax = 15)
        
        # Fill in distances and kValues
        distances <- output$od
        kValues<- output$k
        
      }, error = function(e) {
        message("Error in PcaHubert calculation for external_id: ", unique(df$external_id), " - Replacing with NA")
        # NAs are already in place, so no need to modify distances/kValues
      })
    }
    
    # Store results
    data.frame(
      external_id = question_data[,"external_id"], 
      counter = question_data[,"counter"],
      robustpca_dist = distances,
      kValues = kValues
    )
  }))
  # fix issue with external_id column
  robustpca.results <- robustpca.results |>
    dplyr::select(!external_id) |> 
    tibble::rownames_to_column(var = "external_id") |> 
    dplyr::mutate(external_id = gsub("\\..*", "", external_id)) |> 
    dplyr::mutate(counter = as.integer(counter))
  
  return(robustpca.results)
}


# Psychometrics synonym violation -----------------------------------------
calc_psychometric_synonym_violations <- function(data, 
                                                 items, 
                                                 cor_threshold, 
                                                 synonym_difference_threshold) {
  
  ## Select relevant data to calculate correlations on
  question_data <- data |>
    dplyr::select(external_id, counter, all_of(items)) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(items), as.numeric)) |> # Ensure numeric data
    stats::na.omit() 
  
  n_question_cols <- ncol(question_data)
  
  ## Calculate correlations
  corr.temp <- rcorr(question_data |> 
                       dplyr::select(3:n_question_cols) |>  
                       data.matrix(), type = 'spearman') # Correlation matrix
  
  flattenCorrMatrix <- function(corrmat, pmat) { # Flatten the correlation matrix
    ut <- upper.tri(corrmat)
    data.frame(
      row = rownames(corrmat)[row(corrmat)[ut]],
      column = rownames(corrmat)[col(corrmat)[ut]],
      cor = (corrmat)[ut],
      p = pmat[ut]
    )
  }
  corr.temp <- flattenCorrMatrix(corr.temp$r, corr.temp$P) |>  
    dplyr::arrange(cor)
  
  ## Filter correlations that reach cut-off to be considered synonyms
  synonym_pairs <- corr.temp |>  
    dplyr::filter(cor >= cor_threshold)  |> 
    dplyr::rename('word1' = row ,
           'word2' = column) |> 
    dplyr::select(1:2)
  
  ## Calculate violations per row
  psychometric.synonym.count <- matrix(nrow = nrow(question_data), ncol = nrow(synonym_pairs)) # Create empty matrix to fill
  
  for (row in 1:nrow(synonym_pairs)){ # Loop over pairs
    # Select two synonym items from the relevant response data
    df.temp <- question_data |>  
      dplyr::select(as.character(synonym_pairs[row, 1]), as.character(synonym_pairs[row, 2]))
    
    # Check if the difference between scores on those items reach the threshold
    psychometric.synonym.count[, row] <- abs(df.temp[, 1] - df.temp[, 2]) >= synonym_difference_threshold
  }
  psychometric.synonym.count <- cbind(question_data$external_id, question_data$counter, # Bind id, counter, and synonym violation counts
                                      rowSums(psychometric.synonym.count)) |>  
    as.data.frame() |>  
    dplyr::rename('external_id' = V1, 'counter' = V2, 'syno_violations' = V3) |> 
    dplyr::mutate(counter = as.numeric(counter),
                  syno_violations = as.numeric(syno_violations))
  
  return(psychometric.synonym.count)
}



# Psychometric antonym violation ------------------------------------------
calc_psychometric_antonym_violations <- function(data, 
                                                 items, 
                                                 cor_threshold, 
                                                 antonym_maxvalue_threshold) {
  ## Select relevant data to calculate correlations on
  question_data <- data |>
    dplyr::select(external_id, counter, dplyr::all_of(items)) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(items), as.numeric)) |> # Ensure numeric data
    stats::na.omit() 
  
  n_question_cols <- ncol(question_data)
  
  ## Calculate correlations
  corr.temp <- rcorr(question_data |> 
                       dplyr::select(3:n_question_cols) |> 
                       data.matrix(), type = 'spearman') # Correlation matrix
  
  flattenCorrMatrix <- function(corrmat, pmat) { # Flatten the correlation matrix
    ut <- upper.tri(corrmat)
    data.frame(
      row = rownames(corrmat)[row(corrmat)[ut]],
      column = rownames(corrmat)[col(corrmat)[ut]],
      cor  =(corrmat)[ut],
      p = pmat[ut]
    )
  }
  corr.temp <- flattenCorrMatrix(corr.temp$r, corr.temp$P) |> 
    dplyr::arrange(cor)
  
  ## Filter correlations that reach cut-off to be considered antonyms
  antonym_pairs <- corr.temp |> 
    dplyr::filter(cor <= cor_threshold) |>
    dplyr::rename('word1' = row ,
           'word2' = column) |>
    dplyr::select(1:2)
  
  ## Calculate violations per row
  psychometric.antonym.count <- matrix(nrow = nrow(question_data), ncol = nrow(antonym_pairs)) # Create empty matrix to fill
  
  for (row in 1:nrow(antonym_pairs)){ # Loop over pairs
    # Select two synonym items from the relevant response data
    df.temp <- question_data |> 
      dplyr::select(as.character(antonym_pairs[row, 1]), as.character(antonym_pairs[row, 2]))
    
    # Check if both items score above the threshold
    psychometric.antonym.count[, row] <- df.temp[, 1] > antonym_maxvalue_threshold & df.temp[, 2] > antonym_maxvalue_threshold
  }
  
  psychometric.antonym.count <- cbind(question_data$external_id, question_data$counter, # Bind id, counter, and antonym violation counts
                                      rowSums(psychometric.antonym.count)) |> 
    as.data.frame() |> 
    dplyr::rename('external_id' = V1, 'counter' = V2, 'anto_violations' = V3) |> 
    dplyr::mutate(counter = as.numeric(counter),
                  anto_violations = as.numeric(anto_violations))
  return(psychometric.antonym.count)
}


# SD Response Times -------------------------------------------------------
calc_sd_response_times <- function(data, 
                                   items) {
  
  # preserve identifying columns
  id_cols <- data[, c("external_id", "counter")]
  
  #  make numeric to interpret it as seconds (needed for differences later)
  data[items] <- lapply(data[items], as.POSIXlt, format = "%Y-%m-%d %H:%M:%OS")
  data[items] <- lapply(data[items], as.numeric)
  
  # reorder each row based on the response times
  ordered_rows <- t(apply(data[items], 1, function(x) x[order(x)]))
  ordered_rows_df <- as.data.frame(ordered_rows)
  
  # Calculate time differences between consecutive items for each row
  time_diffs <- t(apply(ordered_rows_df, 1, function(x) diff(x)))
  
  time_diffs_df <- as.data.frame(time_diffs)
  
  # SDs of time differences for each respondent
  sd_time_diff <- apply(time_diffs_df, 1, sd)
  
  df_sd_time_diff <- data.frame(id_cols, sd_time_diff)
  
  
  # Return the final result
  return(df_sd_time_diff)
}



# Longstring Proportions --------------------------------------------------
calc_longstring_proportions <- function(data,
                                        items,
                                        rt_items){
  
  # preserve identifying columns
  id_cols <- data[, c("external_id", "counter")]
  
  # ensure items and rt_items are of the same length and match
  stopifnot(length(items) == length(rt_items))
    rt_truncated <- gsub("rt_", "", rt_items)
  stopifnot(identical(rt_truncated, items))
  
  #  make numeric to interpret it as seconds (needed for differences later)
  data[rt_items] <- lapply(data[rt_items], function(x) {
    as.numeric(as.POSIXct(x, format = "%Y-%m-%d %H:%M:%OS"))
  })
  
  
  # separate items
  response_times <- data[, rt_items]
  actual_responses <- data[, items]
  
  # order each row based on response times
  ordered_responses <- t(sapply(seq_len(nrow(data)), function(i) {
    # order indices for the current row based on response times
    ord_idx <- order(as.numeric(response_times[i, ]))
    # Return responses ordered by these indices
    actual_responses[i, ord_idx]
  }))
  # convert ordered_responses to a data frame
  ordered_responses_df <- as.data.frame(ordered_responses)
  
  # unlist 
  ordered_responses_df <- as.data.frame(do.call(cbind, lapply(ordered_responses_df, unlist)))
 
  # calculate maximum consecutive identical responses, handling NA
  max_consecutive <- function(row) {
    if (all(is.na(row))) {
      return(NA)
    }
    
    max_count <- 1
    current_count <- 1
    for (i in 2:length(row)) {
      # Skip comparison if either value is NA
      if (is.na(row[i]) | is.na(row[i - 1])) {
        current_count <- 1
        next
      }
      if (row[i] == row[i - 1]) {
        current_count <- current_count + 1
        if (current_count > max_count) {
          max_count <- current_count
        }
      } else {
        current_count <- 1
      }
    }
    return(max_count)
  }
  
  # Calculate maximum consecutive identical responses for each row
  max_consecutives <- apply(ordered_responses_df, 1, max_consecutive)
  
  # Calculate proportions
  proportions <- max_consecutives / ncol(ordered_responses_df)
  
  # Combine results with identifying columns
  result <- data.frame(id_cols, longstring_proportion = proportions)
  
  return(result)
   
  
  
}



# Visualization helpers ---------------------------------------------------
theme_bs <- function(){
  # add google font
  sysfonts::font_add_google("News Cycle", "news")
  # use showtext
  showtext::showtext_auto()
  
  # theme
  ggplot2::theme_minimal(base_family = "news") +
    ggplot2::theme(
      # remove minor grid
      panel.grid.minor = ggplot2::element_blank(),
      # Title and Axis Texts
      plot.title = ggplot2::element_text(face = "plain",
                                         size = 22,
                                         hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 16,
                                            hjust = 0.5),
      axis.text.x = ggplot2::element_text(face = "plain", size = 18),
      axis.title.x = ggplot2::element_text(face = "plain", size = 18),
      axis.text.y = ggplot2::element_text(face = "plain", size = 18),
      axis.title.y = ggplot2::element_text(face = "plain", size = 18),
      axis.line = element_line(colour = "#6d6d6e"),
      
      # Faceting
      strip.text = ggplot2::element_text(face = "plain",
                                         size = 20,
                                         hjust = 0.5),
      strip.text.x.top = ggplot2::element_text(face = "plain", 
                                               size = 20,
                                               hjust = 0.5),
      # strip.text.y = element_blank(),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      # Grid
      panel.grid = ggplot2::element_line(colour = "#F3F4F5"),
      # Legend
      legend.title = ggplot2::element_text(face = "plain"),
      legend.position = "top",
      legend.justification = 1,
      # Panel/Facets
      panel.spacing.x = ggplot2::unit(1.25, "lines"),
      panel.spacing.y = ggplot2::unit(1.25, "lines"),
      # Remove vertical grid lines
      panel.grid.major.x = ggplot2::element_blank()
      
    )
}



