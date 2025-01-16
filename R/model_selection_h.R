#' model_selection Function
#'
#' Perform model selection among the models fit with varying number of mixture components (number of clocks).
#' @param results list of 4: $data, $draws_and_summary, $log_lik_matrix_list and $elbo_iterations    
#' @param n_components number of components specified from user 
#' 
#' @keywords fit
#'
#' @return result_model_selection: list $best_fit, $best_K, $model_selection_tibble, $entropy_list)
#'
#' @export
model_selection_h = function(results, n_components = 0) {
  
  data <- results$data$input_data 
  draws_and_summary <- results$draws_and_summary 
  log_lik_matrix_list <- results$log_lik_matrix_list 
  elbo_iterations <- results$elbo_iterations 
  
  
  karyo <- data$karyotype
  
  range_k <- as.numeric(names(results$draws_and_summary)) 
  
  # if (n_components != 0){
  #   range_k = n_components
  # } else if (length(karyo) <= 2){ range_k = (1:length(karyo))
  # } else if (length(karyo) <= 7){ range_k = (1:length(karyo)-1)
  # } else if (length(karyo) <= 10) {range_k = (1:(ceiling(length(karyo)/2))+1)
  # } else if (length(karyo)<= 15) {range_k = (1:(floor(length(karyo)))+1)
  # } else if (length(karyo) <= 25) {range_k = (1:(floor(length(karyo)/2)))
  # } else{
  #   optimal_k <- get_k_max_k_means(data, purity)
  #   # Define the number of additional K values to try
  #   n_additional <- 8
  #   range_k <- round(seq(optimal_k - floor(data$S / 8), optimal_k + floor(data$S / 8), length.out = n_additional))
  #   range_k <- c(range_k,1,2)
  #   range_k <- sort(unique(range_k[range_k >= 1 & range_k <= data$S]))
  # }
  
  
  model_selection_tibble <- dplyr::tibble()
  S <- length(length(karyo))
  
  entropy_per_segment_matrix = list() # k_max rows and S columns
  entropy_per_segment_matrix_norm = list()
  
  for (K in range_k) {
    
    # select it from results
    # elbo_df <- elbo_iterations[[K]]
    
    ############################################
    # SELECT THE BEST MODEL
    draws = 1000
    
    total_number_params <- K+((K-1)*S)+2                    # tau = K, w = K*S, phi, kappa (dirichlet reparametrization)
    N <- length(data$NV)
    
    log_lik_matrix <- log_lik_matrix_list[[as.character(K)]]
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
    tau_and_w_draws <- draws_and_summary[[as.character(K)]]$draws
    w_ICL <- tau_and_w_draws[, colnames(tau_and_w_draws) %in% names_weights]
    
    dim(w_ICL) = c(draws,S*K)
    w_ICL <- apply(w_ICL, 2, stats::median) # median check over draws
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
    
    
    model_selection_tibble <- dplyr::bind_rows(model_selection_tibble, dplyr::tibble(K = K, BIC = BIC, Log_lik = L, ICL = ICL)) #  AIC = AIC, LOO = loo_value,
    
    # ICL PER SEGMENT to see its behaviour with increasing K
    post_Segments = t(w_ICL)
    entropy_per_segment = list()
    entropy_per_segment_norm = list()
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
    
    entropy_per_segment_matrix_norm[[as.character(K)]] = entropy_per_segment_norm
    entropy_per_segment_matrix[[as.character(K)]] = entropy_per_segment
    
  }
  
  entropy_list <- list(entropy_per_segment_matrix = entropy_per_segment_matrix, entropy_per_segment_matrix_norm = entropy_per_segment_matrix_norm)

  model_selection_tibble_temp <- model_selection_tibble[1:2, bycol= TRUE]
  best_K_temp <- min(model_selection_tibble_temp %>% dplyr::filter(BIC == min(BIC)) %>% dplyr::pull(K))

  if (best_K_temp!=1){
    if (length(range_k)==2){
      best_K <- 2
      cli::cli_alert_info("The algorithm should be run with more Components ")
    }else{
      
      while ( best_K_temp < range_k[length(range_k)] & as.numeric((entropy_per_segment_matrix_norm[[as.character(range_k[best_K_temp+1])]])) - as.numeric(entropy_per_segment_matrix_norm[[as.character(range_k[best_K_temp])]]) < 0 ){
        best_K_temp = range_k[best_K_temp + 1]
        if ( best_K_temp == range_k[length(range_k)] ){
          break
        }}
  }} else {
    best_K <- 1
  }
  best_K <- best_K_temp
      # model_selection_tibble_temp <- model_selection_tibble[2:k_max, bycol= TRUE]
      # best_K <- model_selection_tibble_temp %>% dplyr::filter(ICL == min(ICL)) %>% dplyr::pull(K)
    # }
  # }else {
    # best_K <- 1}
       
  if(best_K==range_k[length(range_k)]){
    cli::cli_alert_info("The algorithm should be run with more Components ")
  }
  
  result_model_selection = list(best_fit = draws_and_summary[[as.character(best_K)]], best_K = best_K, model_selection_tibble = model_selection_tibble, entropy_list = entropy_list)
  return (result_model_selection)
}

