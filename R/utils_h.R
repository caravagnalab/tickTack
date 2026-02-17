#' get_clonal_peaks Function
#'
#' Obtain the theoretical peaks that would be observed in a VAF spectrum.
#' @param k karyotype
#' @param purity peaks
#' @keywords peaks
#' @export

get_clonal_peaks = function(k, purity) {
  multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
  major <- multiplicities[1]
  n_tot <- sum(multiplicities)
  # get only Major and 1
  multiplicities <- c(1, major)
  peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
  return(sort(peaks))
}


#' karyo_to_int Function
#'
#' This function allows you to change the string representation of the CN events into numeric values 0 - 1.
#' @param k karyotype string
#' @keywords karyotype
#' @export

karyo_to_int <- function(k) {
  if (k == "2:1") return(1)
  return(0)
}



#' Extract mutations.
#'
#' @description Getter to obtain mutation calls from an object.
#'
#' @param x A CNAqc object.
#' @param cna \code{"clonal"} for clonal CNAs, \code{"subclonal"} for subclonal CNAs.
#' @param type \code{"SNV"} for single-nucleotide variants, \code{"indel"} for insertion-deletions.
#'
#' @return A tibble with the data.
#' @export
#'

Mutations = function(x, cna = c("clonal", "subclonal"), type = c("SNV", "indel"))
{
  # stopifnot(inherits(x, 'cnaqc'))
  
  clonal = NULL
  # if("clonal" %in% cna)
  clonal = x$mutations %>% dplyr::mutate(cna = 'clonal')
  
  subclonal = NULL
  # if(("subclonal" %in% cna) & x$has_subclonal_CNA)
  # subclonal = x$cna_subclonal$mutations %>% Reduce(f = dplyr::bind_rows) %>% dplyr::mutate(cna = 'subclonal')
  
  mutations = dplyr::bind_rows(clonal, subclonal) %>%
    dplyr::select(.data$chr, .data$from, .data$to, .data$ref, .data$alt, .data$NV, .data$DP, .data$VAF, dplyr::everything())
  # %>%
  #   dplyr::filter(type %in% !!type)
  
  if((mutations %>% nrow())== 0) cli::cli_alert_danger("No mutations with these parameters: CNA {.field {cna}}, type {.field {type}}.")
  
  return(mutations)
}

#' Extract CNAs.
#'
#' @description Getter to obtain copy number calls from an object.
#'
#' @param x A CNAqc object.
#' @param type \code{"clonal"} for clonal CNAs, \code{"subclonal"} for subclonal CNAs.
#'
#' @return A tibble with the data.
#' @export
#'

CNA = function(x, type = c("clonal", "subclonal"))
{
  # stopifnot(inherits(x, 'cnaqc'))
  
  clonal = NULL
  # if("clonal" %in% type)
  clonal = x$cna
  
  subclonal = NULL
  # if(("subclonal" %in% type) & x$has_subclonal_CNA)
  #   subclonal = x$cna_subclonal
  
  cna = dplyr::bind_rows(clonal, subclonal) %>%
    dplyr::select(.data$chr, .data$from, .data$to, dplyr::starts_with('Major'), dplyr::starts_with('minor'), .data$CCF, dplyr::everything())
  
  if((cna %>% nrow())== 0) cli::cli_alert_danger("No CNAs with these parameters: {.field {cna}}.")
  
  return(cna)
}




#' Model Saelection.
#'
#' @description Perform model selection data preparation.
#'
#' @param fit results from tickTack_h inference.
#' @param purity sample purity
#' @param coverage sample coverage
#' 
#' @return A tibble with the data.
#' @export
#'

predict_best_number_of_clusters <- function(fit, purity, coverage){
  
  model_path <- system.file("xgboost", "model.bin", 
                            package = "tickTack", mustWork = TRUE)
  model <- xgboost::xgb.load(model_path)
  features <- c("BIC_w", "ICL_w", "rank_AIC", "rank_ICL", "rank_LOO", 
                "rank_BIC", "ICL", "LOO", "AIC", "BIC", "mean_cluster_separation", 
                "mean_cluster_overlap", "min_cluster_separation", "max_cluster_overlap", 
                "mean_density_mut", "mutations_density", "coverage", 
                "purity", "n_cna", "clock_sd", "clock_range", "clock_ci_width")
  tbl_info <- process_file_prediction(fit)
  df <- tbl_info %>% dplyr::mutate(mutations_density = .data$mean_density_mut) %>% 
    dplyr::mutate(min_AIC = min(.data$AIC), min_BIC = min(.data$BIC), 
                  min_ICL = min(.data$ICL), min_LOO = max(.data$LOO), 
                  rank_AIC = dplyr::dense_rank(.data$AIC), rank_BIC = dplyr::dense_rank(.data$BIC), 
                  rank_ICL = dplyr::dense_rank((.data$ICL)), rank_LOO = dplyr::dense_rank(dplyr::desc(.data$LOO))) %>% 
    dplyr::ungroup()
  df <- df %>% dplyr::mutate(score_geom = (.data$min_cluster_separation) - 
                               (.data$max_cluster_overlap)) %>% dplyr::mutate(max_score = max(.data$score_geom))
  df <- df %>% dplyr::rowwise() %>% dplyr::mutate(clock_vals = list(as.numeric(strsplit(.data$inferred_clock, 
                                                                                        "_")[[1]])), clock_low = list(as.numeric(strsplit(.data$inferred_low, 
                                                                                                                                          "_")[[1]])), clock_high = list(as.numeric(strsplit(.data$inferred_high, 
                                                                                                                                                                                             "_")[[1]])), clock_sd = stats::sd(.data$clock_vals), 
                                                  clock_range = max(.data$clock_vals) - min(.data$clock_vals), 
                                                  clock_ci_width = mean(.data$clock_high - .data$clock_low))
  df_single_sample <- df %>% dplyr::mutate(BIC_w = .data$BIC * 
                                             1, ICL_w = .data$ICL * 50) %>% dplyr::mutate(coverage = coverage, 
                                                                                          purity = purity, n_cna = fit$results$data$accepted_cna %>% 
                                                                                            nrow())
  data_pred <- as.matrix(df_single_sample[, features])
  df_single_sample$prob_correct <- stats::predict(model, data_pred)
  best_K_pred <- df_single_sample

  return(best_K_pred)
  
}



#' Model Selection function.
#'
#' @description Perform model selection data preparation.
#'
#' @param fit results from tickTack_h inference.
#' 
#' @return A tibble with the data.
#' @export
#'
process_file_prediction <- function(fit) {
  tryCatch({
    results <- fit$results_timing
    results_model_selection <- model_selection_h(fit$results_timing) #fit$results_model_selection
    summarized_results <- results_model_selection$best_fit$summarized_results
    best_K <- results_model_selection$best_K
    model_selection_tibble <- results_model_selection$model_selection_tibble
    entropy <- results_model_selection$entropy_list
    K = nrow(results_model_selection$model_selection_tibble)
    cna <- results$data$accepted_cna %>% dplyr::rowwise() %>% 
      dplyr::mutate(from = strsplit(.data$segment_name, 
                                    split = "_")[[1]][2], to = strsplit(.data$segment_name, 
                                                                        split = "_")[[1]][3]) %>% dplyr::mutate(len = as.numeric(.data$to) - 
                                                                                                                  as.numeric(.data$from))
    n_mut <- length(results$data$input_data$NV)
    mean_density_mut <- n_mut/sum(cna$len)
    results_model_selection$model_selection_tibble$mean_density_mut = mean_density_mut
    Kmax <- nrow(results_model_selection$model_selection_tibble)
    mean_sep <- numeric(Kmax)
    mean_ov <- numeric(Kmax)
    min_sep <- numeric(Kmax)
    max_ov <- numeric(Kmax)
    tau_inferred <- numeric(Kmax)
    tau_inferred_low <- numeric(Kmax)
    tau_inferred_high <- numeric(Kmax)
    entropy_per_segment <- numeric(Kmax)
    entropy_per_segment_norm <- numeric(Kmax)
    for (i in seq_len(Kmax)) {
      best_k = i
      best_fit = fit$results_timing$draws_and_summary[[i]]$summarized_results
      row = results_model_selection$model_selection_tibble[i, 
      ]
      unique_intervals <- unique(best_fit %>% dplyr::select(.data$clock_low, 
                                                            .data$clock_high, .data$clock_mean))
      K <- nrow(unique_intervals)
      inferred_clusters = best_fit$clock_mean
      tau_inferred[i] <- paste(inferred_clusters, collapse = "_")
      tau_inferred_low[i] <- paste(best_fit$clock_low, 
                                   collapse = "_")
      tau_inferred_high[i] <- paste(best_fit$clock_high, 
                                    collapse = "_")
      if (is.null(K) || K < 2) {
        mean_sep[i] <- NA
        mean_ov[i] <- NA
        min_sep[i] <- NA
        max_ov[i] <- NA
        next
      }
      intervals <- unique_intervals %>% dplyr::mutate(mean = .data$clock_mean, 
                                                      low = .data$clock_low, high = .data$clock_high)
      dist_mean_mat <- matrix(NA_real_, K, K)
      dist_overlap_mat <- matrix(NA_real_, K, K)
      mean_vec <- intervals$mean
      low_vec <- intervals$low
      high_vec <- intervals$high
      dist_mean_mat <- abs(outer(mean_vec, mean_vec, "-"))
      low_max <- pmax(outer(low_vec, low_vec, pmax))
      high_min <- pmin(outer(high_vec, high_vec, pmin))
      overlap <- pmax(0, high_min - low_max)
      union_len <- pmax(outer(high_vec, high_vec, pmax)) - 
        pmin(outer(low_vec, low_vec, pmin))
      dist_overlap_mat <- ifelse(union_len == 0, 0, overlap/union_len)
      diag(dist_mean_mat) <- 0
      diag(dist_overlap_mat) <- 1
      closest_mean_distance <- apply(dist_mean_mat, 1, 
                                     function(x) sort(x)[2])
      closest_overlap <- apply(dist_overlap_mat, 1, function(x) sort(x, 
                                                                     decreasing = TRUE)[2])
      mean_sep[i] <- mean(closest_mean_distance)
      mean_ov[i] <- mean(closest_overlap)
      min_sep[i] <- min(closest_mean_distance)
      max_ov[i] <- max(closest_overlap)
      entropy_per_segment[i] <- entropy$entropy_per_segment_matrix[as.character(i)][[1]][[1]][1]
      entropy_per_segment_norm[i] <- entropy$entropy_per_segment_matrix_norm[as.character(i)][[1]][[1]][1]
    }
    tbl <- results_model_selection$model_selection_tibble
    tbl$mean_cluster_separation <- mean_sep
    tbl$mean_cluster_overlap <- mean_ov
    tbl$min_cluster_separation <- min_sep
    tbl$max_cluster_overlap <- max_ov
    tbl$inferred_clock <- tau_inferred
    tbl$inferred_low <- tau_inferred_low
    tbl$inferred_high <- tau_inferred_high
    tbl$entropy_per_segment <- entropy_per_segment
    tbl$entropy_per_segment_norm <- entropy_per_segment_norm
    return(tbl)
    
  }, error = function(e) {
    message(sprintf("Error", e$message))
    NULL
  })
  
}
