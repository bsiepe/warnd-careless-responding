# <!-- # OLD -->
# 
# <!-- ## Mahalanobis Distance -->
# 
# <!-- Compute the Mahalanobis distance:  -->
# 
# <!-- ```{r mahalanobis-distance} -->
# <!-- indices_mahalanobis <- calc_indicator_mahalanobis(subset_data,  -->
# <!--                                                   items = ema_items) -->
# <!-- ``` -->
# 
# 
# <!-- ## Robust PCA Orthogonal Distance -->
# 
# <!-- ```{r robust-pca} -->
# <!-- indices_robust_pca <- calc_indicator_rob_PCA_orthogonal_distance(subset_data,  -->
# <!--                                                                  items = ema_items) -->
# <!-- ``` -->
# 
# 
# <!-- ## Psychometric syntonym violations -->
# 
# <!-- ```{r psychometrics-syntonym} -->
# <!-- cor_threshold <- .45 -->
# <!-- synonym_difference_threshold <- 5 -->
# 
# <!-- indices_psychometric_syntonym <- calc_psychometric_synonym_violations(subset_data, -->
# <!--                                                                       items = ema_items,  -->
# <!--                                                                       cor_threshold = cor_threshold,  -->
# <!--                                                                       synonym_difference_threshold = synonym_difference_threshold) -->
# 
# <!-- ``` -->
# <!-- Do we already need a threshold here? -->
# 
# 
# <!-- ## Psychometric antonym violations -->
# 
# <!-- ```{r psychometrics-antonym} -->
# <!-- cor_threshold <- .45 -->
# <!-- antonym_maxvalue_threshold <- 5 -->
# 
# <!-- indices_psychometric_antonym <- calc_psychometric_antonym_violations(subset_data,  -->
# <!--                                                                      items = ema_items,  -->
# <!--                                                                      cor_threshold = cor_threshold,  -->
# <!--                                                                      antonym_maxvalue_threshold = antonym_maxvalue_threshold) -->
# 
# <!-- ``` -->
# 
# 
# 
# <!-- ## SD Response Times -->
# 
# <!-- ```{r sd-response-times} -->
# <!-- indices_sd_response_times <- calc_sd_response_times(subset_data,  -->
# <!--                                                     items = rt_ema_items) -->
# <!-- ``` -->
# 
# 
# 
# <!-- ## SD Within-Assessment -->
# 
# <!-- ```{r sd-within-assessment} -->
# <!-- indices_sd_within_assessment <- calc_within_assessment_sd(subset_data, -->
# <!--                                                           items = ema_items,  -->
# <!--                                                           id_col = "external_id") -->
# 
# <!-- ``` -->
# 
# 
# 
# <!-- ## Percentage of Items at the Mode -->
# 
# <!-- ```{r items-at-mode} -->
# <!-- indices_items_at_mode <- calc_mode_percentage(subset_data,  -->
# <!--                                               items = ema_items, -->
# <!--                                               id_col = "external_id") -->
# 
# 
# <!-- ``` -->
# # Combine all indices

# ```{r combine-indices}
# df_indices <- indices_items_at_mode |> 
#   left_join(indices_sd_within_assessment, by = c("external_id", "counter")) |> 
#   left_join(indices_sd_response_times, by = c("external_id", "counter")) |> 
#   left_join(indices_psychometric_antonym, by = c("external_id", "counter")) |> 
#   left_join(indices_psychometric_syntonym, by = c("external_id", "counter")) |> 
#   left_join(indices_robust_pca, by = c("external_id", "counter")) |> 
#   left_join(indices_mahalanobis, by = c("external_id", "counter"))
# 
# 
# ```