
#' Fit the tickTack Model for Timing CNA Events
#'
#' @description
#' Fits the tickTack variational inference model to a CNAqc object to estimate
#' the timing of copy number alteration (CNA) events. The function smooths
#' nearby segments, filters mutations and segments based on quality criteria,
#' and runs variational inference across a range of cluster numbers \code{K},
#' selecting the best model via BIC.
#'
#' @param x A CNAqc object containing somatic mutations and copy number segments,
#'   along with purity information. The object must have mutations accessible via
#'   \code{Mutations(x)} and CNAs via \code{CNA(x)}.
#' @param tolerance Numeric. Relative tolerance for the variational inference
#'   convergence criterion. Default is \code{1e-4}.
#' @param possible_k Character vector. Karyotypes to consider during filtering.
#'   Each entry should be a string of the form \code{"major:minor"}
#'   (e.g., \code{"2:1"}, \code{"2:2"}, \code{"2:0"}). Default is
#'   \code{c("2:1", "2:2", "2:0")}.
#' @param alpha Numeric. Significance level used when filtering mutations within
#'   expected VAF peaks. Default is \code{0.05}.
#' @param min_mutations_number Integer. Minimum number of mutations required for
#'   a segment to be included in the inference. Default is \code{10}.
#' @param max_distance_smooth Numeric. Maximum genomic distance (in base pairs)
#'   between adjacent segments for them to be merged during smoothing.
#'   Default is \code{5e6}.
#' @param min_segment_length Numeric. Minimum segment length (in base pairs)
#'   for a CNA segment to be retained. Default is \code{1e6}.
#' @param n_components Integer or \code{NULL}. Maximum number of timing clusters
#'   \eqn{K} to evaluate. If \code{NULL} (default), the range is automatically
#'   set to \code{1:floor(sqrt(S))}, where \code{S} is the number of
#'   accepted segments.
#'
#' @return The input CNAqc object \code{x}, augmented with a \code{results_timing}
#'   field. This field is a named list containing:
#'   \describe{
#'     \item{\code{data}}{Filtered input data, including accepted CNA segments
#'       and the prepared Stan input list.}
#'     \item{\code{draws_and_summary}}{Named list (by \code{K}) of variational
#'       inference draws and a per-segment summary tibble with columns
#'       \code{clock_mean}, \code{clock_median}, \code{clock_low} (5th percentile),
#'       and \code{clock_high} (95th percentile) of the assigned timing clock \eqn{\tau}.}
#'     \item{\code{log_lik_matrix_list}}{Named list (by \code{K}) of log-likelihood
#'       draw matrices used for BIC computation.}
#'     \item{\code{elbo_iterations}}{Named list (by \code{K}) of data frames
#'       tracking ELBO values across variational inference iterations.}
#'     \item{\code{bic_values}}{Numeric vector of BIC values for each \code{K}
#'       in \code{range_k}.}
#'     \item{\code{best_K}}{Integer. The value of \code{K} that minimises the BIC.}
#'   }
#'
#' @details
#' The function proceeds in the following steps:
#' \enumerate{
#'   \item Nearby CNA segments are merged using \code{\link[CNAqc]{smooth_segments}}.
#'   \item Mutations and segments are filtered by karyotype, minimum mutation
#'     count, VAF peak membership, and minimum segment length.
#'   \item For each \code{K} in \code{range_k}, a tickTack Stan model is fit
#'     via ADVI (automatic differentiation variational inference).
#'   \item Each segment is hard-assigned to the timing clock \eqn{\tau_k} with
#'     the highest posterior probability.
#'   \item Model selection is performed by computing BIC as
#'     \eqn{-2 \cdot \bar{\ell} + (2K - 1) \cdot \log(N)}, where
#'     \eqn{\bar{\ell}} is the mean total log-likelihood across draws and
#'     \eqn{N} is the total number of mutations.
#'   \item The \code{K} minimising BIC is stored as \code{best_K}.
#' }
#'
#' If inference fails to converge for a given \code{K}, a warning is emitted
#' and results for that \code{K} are omitted.
#'
#' @seealso
#' \code{\link[CNAqc]{smooth_segments}}, \code{\link{prepare_input_data}},
#' \code{\link{get_model}}
#'
#' @examples
#' \dontrun{
#' # x is a CNAqc object with purity, mutations, and CNA segments
#' x_timed <- fit_tickTack(
#'   x,
#'   tolerance          = 1e-4,
#'   possible_k         = c("2:1", "2:2", "2:0"),
#'   alpha              = 0.05,
#'   min_mutations_number = 10,
#'   max_distance_smooth  = 5e6,
#'   min_segment_length   = 1e6,
#'   n_components       = NULL
#' )
#'
#' # Inspect the best number of clocks
#' x_timed$results_timing$best_K
#'
#' # Per-segment timing summary for the best K
#' best <- x_timed$results_timing$best_K
#' x_timed$results_timing$draws_and_summary[[as.character(best)]]$summarized_results
#' }
#'
#' @importFrom dplyr filter left_join tibble
#' @importFrom utils read.csv
#'
#' @export
fit_tickTack = function(x,
                    tolerance = 1e-4,
                    possible_k = c("2:1", "2:2", "2:0"),
                    alpha = .05,
                    min_mutations_number = 10,
                    max_distance_smooth = 5e6, #1e7 in GEL
                    min_segment_length = 1e6,#,
                    n_components = NULL
) {
  # stopifnot(inherits(x, 'cnaqc'))
  if ("snvs" %in% names(x) & !"mutations" %in% names(x)) {x$mutations = x$snvs; x$snvs = NULL}

  # 1. Smooth segment close to each other
  x <- CNAqc::smooth_segments(x, maximum_distance = max_distance_smooth)

  # Extract Mutations
  mutations <- Mutations(x)

  if(!nrow(mutations)*ncol(mutations)){
    stop("No mutations have been called on this CNAqc object.")
  }

  # Extract CNAs
  segments <- CNA(x)

  # Filter out segments too short
  segments = segments %>% dplyr::filter((to - from) >= min_segment_length)


  if(!nrow(segments)*ncol(segments)){
    stop("No CNA events have been called on this CNAqc object.")
  }

  purity = x$purity

  # Filter data by minimum numver of mutations, karyotypes and filter mutations within expected peaks
  accepted_data <- prepare_input_data(mutations, segments, purity, possible_k = possible_k, alpha = alpha, min_mutations_number = min_mutations_number)

  input_data = accepted_data$input_data
  accepted_cna = accepted_data$accepted_cna

  message(" ", input_data$S, " segments will be included in the inference")

  if (is.null(n_components)) {
    range_k = 1:floor(sqrt(input_data$S))
  } else {
    range_k = 1:min(n_components,input_data$S)
  }

  message("Performing inference with maximum ", max(range_k), " clusters.\nInsert a specific set of values in the <n_components> parameter if a different value is desired! ")

  # if (input_data$S==1 & (tolerance>=0.0001|tolerance<=0.01)){
  #   message("Performing inference with ", range_k, " component. Decreasing tolerance to 0.01")
  #   tolerance = 0.001
  # }

  # before inference add K to the list obtained as input_data

  fit_single_k = function(input_data, K, tolerance) {
    m = get_model("tickTack")
    input_data$K = K
    fit = m$variational(input_data, tol_rel_obj = tolerance, save_latent_dynamics = TRUE)
    fit
  }

  assign_tau_to_segments = function(fit, S, K) {
    tau_draws    <- fit$draws("tau", format = "matrix")          # (draws, K)
    probs_draws  <- fit$draws("seg_probs", format = "matrix")    # (draws, S*K) — needs reshape
    n_draws <- nrow(tau_draws)

    # Reshape seg_probs to (draws, S, K)
    probs_arr <- array(probs_draws, dim = c(n_draws, S, K))

    # Average probabilities over draws → (S, K)
    mean_probs  <- apply(probs_arr, c(2, 3), mean)

    # Hard assignment per segment
    hard_assign <- apply(mean_probs, 1, which.max)   # length S

    # Taus per clock
    mean_tau <- apply(tau_draws, 2, mean)             # length K
    median_tau <- apply(tau_draws, 2, median)
    q5 = function(x) {quantile(x, .05)}
    q95 = function(x) {quantile(x, .95)}
    low_tau = apply(tau_draws, 2, q5)
    high_tau = apply(tau_draws, 2, q95)


    # Tau assigned to each segment
    dplyr::tibble(
      segment_id = 1:S,
      clock_mean = mean_tau[hard_assign],
      clock_median = median_tau[hard_assign],
      clock_low = low_tau[hard_assign],
      clock_high = high_tau[hard_assign]
    )
  }

  draws_and_summary = c()
  elbo_iterations = list()
  log_lik_matrix_list = list()

  for (K in range_k){

    res <- tryCatch({

      input_data$K = K
      fit = fit_single_k(input_data, K, tolerance)
      clock_assignment = assign_tau_to_segments(fit, input_data$S, input_data$K)

      summarized_results <- accepted_cna %>% dplyr::left_join(clock_assignment)

      draws = fit$draws(format = "matrix")
      log_lik_matrix <- fit$draws(variables = "log_lik")
      log_lik_matrix_list[[as.character(K)]] <- log_lik_matrix

      result_single <- list(draws = draws, summarized_results = summarized_results)
      draws_and_summary[[as.character(K)]] = result_single

      # Get Elbo
      output_files <- fit$latent_dynamics_files()
      #print(paste0("output_files ", output_files,"\n"))
      elbo_data <- utils::read.csv(output_files, header = FALSE, comment.char = "#")
      colnames(elbo_data) <- c("iter", "time_in_seconds", "ELBO")
      iterations <- elbo_data$iter  # iteration column
      elbo_values <- elbo_data$ELBO  # ELBO column
      elbo_df <- data.frame(iteration = iterations, elbo = elbo_values)
      elbo_iterations[[as.character(K)]] = elbo_df
      fit
    }, error = function(e) {message(paste("Inference with number of components = ", K, " could not converge. Inference is available up to ", K - 1, "."))

    },finally = {
      #pass
    })

  }

  # Model selection
  get_bic = function(K, n_obs) {
    log_lik_mat = log_lik_matrix_list[[as.character(K)]]
    if (!any(class(log_lik_mat) == "draws")) return(Inf)

    # Total log-lik per draw = sum over segments
    total_ll_per_draw <- rowSums(log_lik_mat)
    mean_ll <- mean(total_ll_per_draw)

    n_params <- 2 * K - 1  # (K-1) for pi + K for tau
    bic <- -2 * mean_ll + n_params * log(n_obs)
    bic
  }
  bics = lapply(range_k, function(K) {get_bic(K = K, n_obs = input_data$N)}) %>% unlist()
  best_K = which.min(bics)

  results_and_data = list(data = accepted_data,
                          draws_and_summary = draws_and_summary,
                          log_lik_matrix_list = log_lik_matrix_list,
                          elbo_iterations = elbo_iterations,
                          bic_values = bics,
                          best_K = best_K)

  x$results_timing = results_and_data
  return(x)
}
