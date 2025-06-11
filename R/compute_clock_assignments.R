compute_clock_assignment <- function(accepted_cna, tau_and_w_draws, tau_and_w_summary, K, data_simulation=NULL){

  names_tau <- paste("tau[", 1:K, "]", sep = "")
  tau_inferred <- tau_and_w_draws[, colnames(tau_and_w_draws) %in% names_tau]
  tau_inferred_medianap <- lapply(1:ncol(tau_inferred), function(i) {stats::median(tau_inferred[,i])} ) %>% unlist() 
  
  tau_original_list = c()
  tau_inferred_assigned_list = tibble::tibble()
  alpha_median = tibble::tibble()
  beta_median = tibble::tibble()
  identity_matrix_RI = c()
  
  for (i in 1:nrow(accepted_cna)){
    segment <- accepted_cna[i,]
    segment_number <- segment$segment_id
    # 
    # tau_original <- data_simulation$taus[segment$segment_original_indx]
    # tau_original_list <- c(tau_original_list, tau_original)
    # 
    names_weights <- paste("w[",segment_number,",", 1:K, "]", sep = "")  #regex_pars = c("w")
    weights_inferred <- tau_and_w_draws[, colnames(tau_and_w_draws) %in% names_weights]
    
    weights_inferred_median <- lapply(1:ncol(weights_inferred), function(i) {stats::median(weights_inferred[,i])} ) %>% unlist() 
    
    tau_index_assigned <- which.max(weights_inferred_median)
    
    tau_inferred_assigned =  tau_inferred_medianap[tau_index_assigned]
    tau_inferred_low <- tau_and_w_summary$q5[tau_index_assigned]
    tau_inferred_high <- tau_and_w_summary$q95[tau_index_assigned]
    tau_inferred_single <- tibble::tibble (tau_inferred_low = tau_inferred_low,tau_inferred_median = tau_inferred_assigned, tau_inferred_high = tau_inferred_high)
    # tau_inferred_assigned_list <- c(tau_inferred_assigned_list, tau_inferred_single)
    tau_inferred_assigned_list <- dplyr::bind_rows(tau_inferred_assigned_list, tau_inferred_single)
    
    identity_matrix_RI = c(identity_matrix_RI, tau_index_assigned) #extract the vector of taus (unordered) to which the ordered segments are assigned
    
    names_alpha <- paste0("theta[", segment_number, ",", tau_index_assigned, ",1]")
    names_alpha <- as.vector(names_alpha)
    names_beta <- paste0("theta[", segment_number, ",", tau_index_assigned, ",2]")
    names_beta <- as.vector(names_beta)
    
    
    alpha_inferred <- tau_and_w_draws[, colnames(tau_and_w_draws) %in% names_alpha]
    beta_inferred <- tau_and_w_draws[, colnames(tau_and_w_draws) %in% names_beta]
    
    alpha_inferred_median <- lapply(1:ncol(alpha_inferred), function(i) {stats::median(alpha_inferred[,i])} ) %>% unlist() 
    beta_inferred_median <- lapply(1:ncol(beta_inferred), function(i) {stats::median(beta_inferred[,i])} ) %>% unlist() 
    alpha_median <- dplyr::bind_rows(alpha_median, tibble::tibble (alpha_median = alpha_inferred_median))
    beta_median <- dplyr::bind_rows(beta_median, tibble::tibble (beta_median = beta_inferred_median))
    
  }
  
  # result = list(original_clocks = tau_original_list, predicted_clocks = tau_inferred_assigned_list )
  predictions=list(predicted_clocks=tau_inferred_assigned_list, alpha = alpha_median, beta = beta_median)
  return(predictions)
  
}

