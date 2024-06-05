#' Fit timing of clonal peaks in cancer genome sequencing data
#'
#' This function fits the timing of clonal peaks in cancer genome sequencing data using either a beta-binomial or binomial model.
#'
#' @param segments A data frame containing segment information with columns `chr`, `from`, `to`, `Major`, and `minor`.
#' @param mutations A data frame containing mutation information with columns `chr`, `from`, `to`, `DP`, and `NV`.
#' @param purity A numeric value representing the tumor purity.
#' @param possible_k A character vector of possible karyotypes in the format "Major:minor". Default is c("2:1", "2:2", "2:0").
#' @param alpha A numeric value for the significance level. Default is 0.05.
#' @param min_mutations_number An integer specifying the minimum number of mutations required for analysis. Default is 2.
#' @param beta_binomial A logical value indicating whether to use the beta-binomial model. Default is FALSE.
#' @param beta_binomial_disp A numeric value for the beta-binomial dispersion parameter. Default is 0.01.
#'
#' @return A list containing two tibbles: `inference_results` and `summarized_results`. Returns NULL if no results are obtained.
#' @export
fit <- function(segments, mutations, purity,
                       possible_k = c("2:1", "2:2", "2:0"),
                       alpha = .05,
                       min_mutations_number = 2,
                       beta_binomial = FALSE,
                       beta_binomial_disp = 0.01) {

  # Load the appropriate model based on the beta_binomial parameter
  if (beta_binomial) {
    model <- get_model("timing_betabinomial")
  } else {
    model <- get_model("timing_binomial")
  }

  # Drop segments with NA values in Major or minor columns
  segments <- segments %>%
    tidyr::drop_na(Major, minor)

  n_segments <- nrow(segments)
  inference_results <- dplyr::tibble()
  summarized_results <- dplyr::tibble()

  # Loop through each segment
  for (segment_idx in 1:n_segments) {
    segment <- segments[segment_idx, ]
    chr <- segment$chr
    segment_id <- paste(chr, segment$from, segment$to, sep = "_")

    # Get karyotype
    Major <- segment$Major
    minor <- segment$minor
    k <- paste(Major, minor, sep=':')

    # Get clonal peaks for the karyotype
    peaks <- get_clonal_peaks(k, purity)

    if (k %in% possible_k) {
      # Filter mutations within the segment and drop rows with NA in DP column
      segment_mutations <- mutations %>%
        dplyr::filter(.data$chr == segment$chr, .data$from > segment$from, .data$to < segment$to) %>%
        tidyr::drop_na(DP)

      accepted_mutations <- data.frame()

      if (nrow(segment_mutations) > 0) {
        # Define the confidence interval
        probs <- c(alpha / 2, 1 - alpha / 2)
        DP <- segment_mutations$DP
        NV <- segment_mutations$NV

        # Check if mutations fall within the confidence interval
        accepted_idx <- lapply(1:length(DP), function(i) {
          for (p in peaks) {
            if (beta_binomial) {
              quantiles <- TailRank::qbb(probs, DP[i], p * (1 - beta_binomial_disp) / beta_binomial_disp, (1 - p) * (1 - beta_binomial_disp) / beta_binomial_disp)
            } else {
              quantiles <- stats::qbinom(probs, DP[i], p)
            }
            if ((NV[i] >= quantiles[1]) && (NV[i] <= quantiles[2])) {
              return(i)
            }
          }
        }) %>% unlist()

        # Get the accepted mutations
        accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx])
      }

      # Fit the model if the number of accepted mutations is greater than or equal to the minimum required
      if (nrow(accepted_mutations) >= min_mutations_number) {
        cli::cli_alert_info("Fitting segment with index {.val {segment_idx}}")

        # Prepare input data for the model
        input_data <- list(
          N = nrow(accepted_mutations),
          NV = accepted_mutations$NV,
          DP = accepted_mutations$DP,
          peaks = peaks,
          beta_dispersion = beta_binomial_disp
        )

        fit <- model$sample(data = input_data, iter_warmup = 2000, iter_sampling = 2000, chains = 8, parallel_chains = 8)

        # Compute tau posteriors
        tau_posteriors <- get_tau_posteriors(fit, k)$tau %>% unname() %>% as.numeric()
        tau_low <- stats::quantile(tau_posteriors, alpha / 2) %>% unname()
        tau_high <- stats::quantile(tau_posteriors, 1 - alpha / 2) %>% unname()
        tau_mean <- mean(tau_posteriors)

        # Store inference results
        inference_results <- dplyr::bind_rows(inference_results, dplyr::tibble(tau = tau_posteriors, segment = segment_idx, karyotype = k, chr = chr, segment_id = segment_id))
        summarized_results <- dplyr::bind_rows(summarized_results, dplyr::tibble(tau_low = tau_low, tau_mean = tau_mean, tau_high = tau_high, segment = segment_idx, karyotype = k, chr = chr, segment_id = segment_id))
      }
    }
  }

  # Return NULL if no results are obtained
  if (nrow(inference_results) == 0) {
    cli::cli_alert_danger("Inference concluded without errors but with no results.")
    return(NULL)
  }

  # Return the results
  return(list(inference_results = inference_results, summarized_results = summarized_results))
}
