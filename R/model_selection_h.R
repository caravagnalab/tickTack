#' model_selection Function
#'
#' Perform model selection among the models fit with varying number of mixture components (number of clocks).
#' @param data list of k_max of lists: $input_data:list: List of 7: $S: int, $N: int, $karyotype: num (0 or 1), $seg_assignment: num, $peaks:List of N of num (1:2), $NV: num, $DP: num
#' @param draws_and_summary  list of lenght k_max of draws form the variational method with the summary statistics of the draws from the approximate posterior
#' @param log_lik_matrix_list list of lenght k_max
#' @param elbo_iterations list of lenght k_max
#' 
#' @keywords fit
#'
#' @return result_model_selection: list $best_fit, $best_K, $model_selection_tibble, $entropy_list)
#'
#' @export
model_selection_h = function(data, draws_and_summary, log_lik_matrix_list, elbo_iterations) {
  karyo <- data$karyotype
  
  if (length(karyo) <= 2){
    k_max = (length(karyo))
  } else if (length(karyo) <= 7){
    k_max = (length(karyo)-1)
  } else if (length(karyo) <= 15) {
    k_max = ((floor(length(karyo)/2))-1)
  } else{
    k_max = ceiling(sqrt(length(karyo)))
  }
  
  model_selection_tibble <- dplyr::tibble()
  S <- length(length(karyo))
  
  entropy_per_segment_matrix = matrix(0, k_max, S) # k_max rows and S columns
  entropy_per_segment_matrix_norm = matrix(0, k_max, S)
  
  for (K in 1:k_max) {
    
    # select it from results
    # elbo_df <- elbo_iterations[[K]]
    
    ############################################
    # SELECT THE BEST MODEL
    draws = 1000
    
    total_number_params <- K+((K-1)*S)+2                    # tau = K, w = K*S, phi, kappa (dirichlet reparametrization)
    N <- length(data$NV)
    
    log_lik_matrix <- log_lik_matrix_list[[K]]
    log_lik_total_per_sample <- rowSums(log_lik_matrix)
    L <- stats::median(log_lik_total_per_sample)
    
    BIC <- ((total_number_params * log(N)) - 2 * L) # %>% unname()
    AIC <- 2 * total_number_params - 2 * L
    
    loo_result <- loo::loo(log_lik_matrix)
    loo_value <- loo_result$estimates[3, "Estimate"]  # LOO-CV estimate                                                    #getter?
    
    # Calculate ICL --> move outside the main function
    # Generate names for w[sim_params_num, K] format
    
    
    # extract w from results instead that from stanfit
    names_weights <- outer(1:S, 1:K,
                           FUN = function(i, j) paste0("w[", i, ",", j, "]"))
    names_weights <- as.vector(names_weights)
    tau_and_w_draws <- draws_and_summary[[K]]$draws
    w_ICL <- tau_and_w_draws[, colnames(tau_and_w_draws) %in% names_weights]
    
    dim(w_ICL) = c(draws,S*K)
    w_ICL <- apply(w_ICL, 2, .data$median) # median check over draws
    dim(w_ICL) = c(S,K) # check by row
    w_ICL = t(w_ICL)
    
    num_mutations_all <- c()
    for (i in seq_along(unique(data$seg_assignment))) {
      segment <- unique(data$seg_assignment)[i]
      num_mutations_single <- length( data$seg_assignment[(data$seg_assignment == segment)])
      num_mutations_all <- c(num_mutations_all, num_mutations_single)
    }
    mut_per_seg = num_mutations_all
    
    res_entropy = 0
    post = w_ICL
    for (k in 1:K ){
      post_k = post[k]
      log_post_k = log(post_k + 0.000001)
      post_k_entr = post_k * log_post_k * mut_per_seg
      post_k_entr = sum(post_k_entr)
      post_k_entr = -1 * (post_k_entr)
      res_entropy = res_entropy + post_k_entr
    }
    entropy = res_entropy
    ICL = BIC + entropy
    
    
    model_selection_tibble <- dplyr::bind_rows(model_selection_tibble, dplyr::tibble(K = K, BIC = BIC, AIC = AIC, LOO = loo_value, Log_lik = L, ICL = ICL))
    
    # ICL PER SEGMENT to see its behaviour with increasing K
    post_Segments = t(w_ICL)
    entropy_per_segment = c()
    entropy_per_segment_norm = c()
    for (s in 1:S ){
      post_s = post_Segments[s]
      log_post_s = log(post_s + 0.000001)
      post_s_entr = post_s * log_post_s * mut_per_seg
      post_s_entr = sum(post_s_entr)
      post_s_entr = -1 * (post_s_entr)
      entropy_per_segment = c(entropy_per_segment, post_s_entr)
      
      post_s_entr_norm = post_s * log_post_s
      post_s_entr_norm = sum(post_s_entr_norm)
      post_s_entr_norm = -1 * (post_s_entr_norm)
      entropy_per_segment_norm = c(entropy_per_segment_norm, post_s_entr_norm)
      
    }
    
    print("entropy per segment: ")
    print(entropy_per_segment)
    
    print("entropy per segment normalized: ")
    print(entropy_per_segment_norm)
    
    entropy_per_segment_matrix_norm[K,] = entropy_per_segment_norm
    entropy_per_segment_matrix[K,] = entropy_per_segment
    
  }
  
  entropy_list <- list(entropy_per_segment_matrix = entropy_per_segment_matrix, entropy_per_segment_matrix_norm = entropy_per_segment_matrix_norm)
  
  model_selection_tibble_temp <- model_selection_tibble[1:2, bycol= TRUE]
  best_K_temp <- model_selection_tibble_temp %>% dplyr::filter(BIC == min(BIC)) %>% dplyr::pull(K)
  
  if (best_K_temp!=1){
    if (k_max==2){
      best_K <- 2
    }else{
      while(mean(entropy_per_segment_matrix_norm[best_K_temp+1,]) - mean(entropy_per_segment_matrix_norm[best_K_temp,]) < 0 & best_K_temp < k_max ){
        best_K_temp = best_K_temp + 1
        if ( best_K_temp == k_max ){
          break
        }}}
  } else {
    best_K <- 1
  }
  best_K <- best_K_temp
  
  if(best_K==k_max){
    cli::cli_alert_info("The algorithm should be run with more Components ")
  }
  
  result_model_selection = list(best_fit = draws_and_summary[[best_K]], best_K = best_K, model_selection_tibble = model_selection_tibble, entropy_list = entropy_list)
  return (result_model_selection)
}

