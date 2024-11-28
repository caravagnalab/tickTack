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
#' @return   results_and_data = list(data = input_data_list, results = results, output_files_list = output_files_list)
#' @export
fit_h = function(x, max_attempts=2, INIT=TRUE, tolerance = 0.01, possible_k = c("2:1", "2:2", "2:0"), alpha = .05, min_mutations_number = 2, n_components = 0)
{
  # stopifnot(inherits(x, 'cnaqc'))
  
  # mutations <- CNAqc::Mutations(x)
  mutations <- Mutations(x)
  
  if(!nrow(mutations)*ncol(mutations)){
    stop("No mutations have been called on this CNAqc object.")
  }
  
  # cna <- CNAqc::CNA(x)
  segments <- CNA(x)
  # just while developing the package to accelerate the inference
  # cna <- cna[1:13,]
  
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
  
  if (n_components != 0){
    k_max = n_components
  } else if (input_data$S <= 2){
    k_max = (input_data$S)
  } else if (input_data$S <= 7){
    k_max = (input_data$S-1)
  } else if (input_data$S <= 15) {
    k_max = ((floor(input_data$S/2))-1)
  } else{
    k_max = ceiling(sqrt(input_data$S))
  }
  
  draws_and_summary = c()
  elbo_iterations = list()
  log_lik_matrix_list = list()
  
  
  # set k_max
  # before inference add K to the list obtained as input_data
  for (K in 1:k_max){
    input_data$K = K
    
    if (INIT==TRUE){
      inits_chain <- get_initialization(input_data, purity = purity)
      res <- fit_variational_h(input_data,
                               initialization = inits_chain,
                               max_attempts = max_attempts,
                               INIT = INIT,
                               initial_iter = 1000,
                               grad_samples = 10,
                               elbo_samples = 100,
                               tolerance = tolerance)
    }
    
    S = input_data$S
    # extract only what's needed from the fit
    names_tau <- paste("tau[", 1:K, "]", sep = "")
    names_weights <- outer(1:S, 1:K,
                           FUN = function(i, j) paste0("w[", i, ",", j, "]"))
    names_weights <- as.vector(names_weights)
    vars <- append (names_tau, names_weights)
    
    tau_and_w_draws <- res$draws(variables = vars)
    tau_and_w_summary <- res$summary(variables = vars)
    
    clock_assignment <- compute_clock_assignment(accepted_cna, tau_and_w_draws, tau_and_w_summary, K, data_simulation=NULL)
    
    summarized_results <- accepted_cna %>%
      dplyr::mutate(clock_mean = clock_assignment$tau_inferred_median) %>%
      dplyr::mutate(clock_low = clock_assignment$tau_inferred_low) %>%
      dplyr::mutate(clock_high = clock_assignment$tau_inferred_high)


    
    # Check log likelihood values  POSSO ESTRARRE LA LOG LIK DIRETTAMENTE DAL MODELLO E NON DALLE GENERATED QUANTITIES?
    
    # log_lik_contributions <- res$draws(variables = "log_lik_matrix")
    log_lik_matrix <- res$draws(variables = "log_lik")
    log_lik_matrix_list[[K]] <- log_lik_matrix
    result_single <- list(draws = tau_and_w_draws, summary = tau_and_w_summary, summarized_results = summarized_results)
    draws_and_summary[[K]] = result_single
    
    output_files <- res$latent_dynamics_files()
    print(paste0("output_files ", output_files,"\n"))
    elbo_data <- utils::read.csv(output_files, header = FALSE, comment.char = "#")
    colnames(elbo_data) <- c("iter", "time_in_seconds", "ELBO")
    iterations <- elbo_data$iter  # iteration column
    elbo_values <- elbo_data$ELBO  # ELBO column
    elbo_df <- data.frame(iteration = iterations, elbo = elbo_values)
    
    elbo_iterations[[K]] = elbo_df
    
    # get summarized results from the draws (estimate of the paramenters and the relative confidence associated to them)
    # implement a function to get the fit like "get_posterior"
    
    
    
    
  }
  
  results_and_data = list(data = accepted_data, draws_and_summary = draws_and_summary, log_lik_matrix_list = log_lik_matrix_list, elbo_iterations = elbo_iterations)
  x$results_timing = results_and_data
  return(x)
}
