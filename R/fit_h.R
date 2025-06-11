#' Get the fit to obtain the clocks posteriors
#'
#' @description Obtain the approximate posterior of the clocks for each model fit with the number of components up to k_max.
#'
#' @param x list: A CNAqc object.
#' @param max_attempts num: max number of repeated inference for ADVI
#' @param INIT logical: boolean variable to set the initialization phase to TRUE or FALSE
#' @param tolerance num: tolerance in the ELBO optimization procedure
#'
#' @param possible_k chr: "2:1" "2:2" "2:0"
#' @param alpha num: (type double) confidence interval level to choose the data that fall in the expected binomial intervals
#' @param min_mutations_number num: (type double) minimum number of accepted mutations for a segment to be included in the inference
#' @param n_components  number of components specified from user 
#'
#' @param initial_iter description
#' @param grad_samples description
#' @param elbo_samples description
#' @param tmp_file_path path of the directory where to save the temporary files with the info on the elbo evaluations during the VI inference
#' @param cmd_version_old version of cmdstanr for the draws parameter invariational method 
#'
#' @return   results_and_data = list(data = input_data_list, results = results, output_files_list = output_files_list)
#' @export
fit_h = function(x, max_attempts=2, INIT=TRUE, tolerance = 0.0001, possible_k = c("2:1", "2:2", "2:0"), alpha = .05, min_mutations_number = 4, n_components = 0, initial_iter=200, grad_samples=10, elbo_samples=200, tmp_file_path = NULL, cmd_version_old=FALSE)
{
  # stopifnot(inherits(x, 'cnaqc'))
  
  # mutations <- CNAqc::Mutations(x)
  mutations <- Mutations(x)
  
  if(!nrow(mutations)*ncol(mutations)){
    stop("No mutations have been called on this CNAqc object.")
  }
  
  # cna <- CNAqc::CNA(x)
  segments <- CNA(x)
  
  if(!nrow(segments)*ncol(segments)){
    stop("No CNA events have been called on this CNAqc object.")
  }
  
  # extract from the CNAqc object: data = A tibble [N × 9] (S3: tbl_df/tbl/data.frame): (key + attribute type) $mutation: chr, $allele: chr, $type: chr "private", $karyotype:chr "2:2" "2:2" "2:2" "2:2",$segment_id: int, $DP: int, $NV: int, $tau: num (0 for real case analysis or num for validation analysis), $segment_name:chr, $segment_name_real:chr = segment_idx to keep track of the accepted/non accepted segments. The table is ordered with respect to the segment_id attribute (is it necessary?)
  mutations = mutations
  
  # temporarly set the purity here or give it in input before implementing a getter for the purity
  purity = x$metadata$purity
  
  accepted_data <- prepare_input_data(mutations, segments, purity, possible_k = possible_k, alpha = alpha, min_mutations_number = min_mutations_number)
  
  input_data = accepted_data$input_data
  accepted_cna = accepted_data$accepted_cna
  
  
  message(" ", input_data$S, " segments will be included in the inference")
  
  
  if (n_components != 0){
    range_k = n_components
  } else if (input_data$S <= 2){ range_k = (1:input_data$S)
  } else if (input_data$S <= 7){ range_k = (1:(input_data$S-1))
  } else if (input_data$S <= 10) {range_k = (1:((ceiling(input_data$S/2))))
  } else {range_k = (1:5)
  # } else if (input_data$S <= 15) {range_k = (1:((floor(input_data$S/2))+1))
  # } else if (input_data$S <= 25) {range_k = (1:((floor(input_data$S/2))))
  # } else{
  #   optimal_k <- get_k_max_k_means(input_data, purity)
  #     # Define the number of additional K values to try
  #     n_additional <- 8
  #   range_k <- round(seq(optimal_k - floor(input_data$S / 8), optimal_k + floor(input_data$S / 8), length.out = n_additional))
  #   range_k <- c(range_k,1,2)
  #   range_k <- sort(unique(range_k[range_k >= 1 & range_k <= input_data$S]))
  }
  

  message("Performing inference with the following number of components ", range_k, ". Insert a specificset of values in the <range> parameter if a different set of components is desired! ")
  
  draws_and_summary = c()
  elbo_iterations = list()
  log_lik_matrix_list = list()
  
  if (input_data$S==1 & tolerance>=0.0001){
    message("Performing inference with ", range_k, " component. Decreasing tolerance to 0.01")
    tolerance = 0.01
  }
  
  # before inference add K to the list obtained as input_data
  for (K in range_k){
    input_data$K = K
    
    
      # inits_chain <- get_initialization(input_data, purity = purity)
      inits_chain <- NULL
      res <-  tryCatch({res <-fit_variational_h(input_data,
                                                purity = purity,
                                               initialization = inits_chain,
                                               max_attempts = max_attempts,
                                               INIT = INIT,
                                               initial_iter = initial_iter,
                                               grad_samples = grad_samples,
                                               elbo_samples = elbo_samples,
                                               tolerance = tolerance,
                                               tmp_file_path = tmp_file_path,
                                               cmd_version_old = cmd_version_old)
    
    
    S = input_data$S
    
    # extract only what's needed from the fit
    names_tau <- paste("tau[", 1:K, "]", sep = "")
    names_weights <- outer(1:S, 1:K,
                           FUN = function(i, j) paste0("w[", i, ",", j, "]"))
    names_weights <- as.vector(names_weights)
    names_alpha_beta <- apply(
      expand.grid(i = 1:S, j = 1:K, z = 1:2),
      1,
      function(x) paste0("theta[", x[1], ",", x[2], ",", x[3], "]")
    )
    names_alpha_beta <- as.vector(names_alpha_beta)
    vars <- append (names_tau, names_weights)
    vars <- append (vars, names_alpha_beta)
    
    
    ###################################################################
    tau_and_w_draws <- res$draws(variables = vars)
    tau_and_w_summary <- res$summary(variables = vars)
    
    clock_assignment <- compute_clock_assignment(accepted_cna, tau_and_w_draws, tau_and_w_summary, K, data_simulation=NULL)
    
    summarized_results <- accepted_cna %>%
      dplyr::mutate(clock_mean = clock_assignment$predicted_clocks$tau_inferred_median) %>%
      dplyr::mutate(clock_low = clock_assignment$clock_assignment$predicted_clocks$tau_inferred_low) %>%
      dplyr::mutate(clock_high = clock_assignment$clock_assignment$predicted_clocks$tau_inferred_high) %>%
      dplyr::mutate(alpha = clock_assignment$alpha$alpha_median) %>%
      dplyr::mutate(beta = clock_assignment$beta$beta_median)


    
    log_lik_matrix <- res$draws(variables = "log_lik")
    log_lik_matrix_list[[as.character(K)]] <- log_lik_matrix
    result_single <- list(draws = tau_and_w_draws, summary = tau_and_w_summary, summarized_results = summarized_results)
    draws_and_summary[[as.character(K)]] = result_single
    
    output_files <- res$latent_dynamics_files()
    print(paste0("output_files ", output_files,"\n"))
    elbo_data <- utils::read.csv(output_files, header = FALSE, comment.char = "#")
    colnames(elbo_data) <- c("iter", "time_in_seconds", "ELBO")
    iterations <- elbo_data$iter  # iteration column
    elbo_values <- elbo_data$ELBO  # ELBO column
    elbo_df <- data.frame(iteration = iterations, elbo = elbo_values)
    
    elbo_iterations[[as.character(K)]] = elbo_df
    res
    
    # implement a function to get the fit like "get_posterior"
    ####################################################################################
    
      }, error = function(e) {message(paste("Inference with number of components = ", K, " could not converge. Inference is available up to ", K - 1, "."))
        
  },finally = {
    #pass
  })
    
  }
  
  results_and_data = list(data = accepted_data, draws_and_summary = draws_and_summary, log_lik_matrix_list = log_lik_matrix_list, elbo_iterations = elbo_iterations)
  x$results_timing = results_and_data
  return(x)
}
