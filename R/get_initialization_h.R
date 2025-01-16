#' get_initialization
#'
#' Perform c-means on the proportions of mutations that fall in to the first or second binomial interval (in accepted mutations) and retrieve the centroids and u matrix to be used as initialization parameters.
#' @param input_data list: List of 7: $S: int, $N: int, $karyotype: num (0 or 1), $seg_assignment: num, $peaks:List of N of num (1:2), $NV: num, $DP: num
#' @param purity num: (type = double) sample purity between 0 and 1
#' @param phi parameters of dirichlet reparametrization uniformly initialized
#' @param kappa parameters of dirichlet reparametrization
#' @param alpha num: (type double) confidence interval level to choose the data that fall in the expected binomial intervals
#' @keywords cluster
#'
#' @return inits list: List of 4: $w: num (1:S, 1:3), $tau: num (1:K), $phi: num (1:K), $kappa: num
#'
#' @export

#take as input the data ready to be used for inference
get_initialization = function(input_data, purity, phi=c(), kappa=5, alpha = 0.05){
  K = input_data$K
  
  karyotype <- input_data$karyotype
  
  peaks = input_data$peaks
  
  # Check if mutation is inside CI
  probs <- c(alpha/2, 1 - alpha/2)
  
  DP <- input_data$DP
  NV <- input_data$NV
  
  
  #write it as in prepare_input_data as the process is the same
  
  # print(unique(input_data$seg_assignment))
  alpha_beta_all <- lapply(1:length(DP), function(i) {
    for (j in unique(input_data$seg_assignment)) {
      if (input_data$seg_assignment[i]==j){
        quantiles_1 <- stats::qbinom(probs, DP[i], peaks[[j]][1])
        quantiles_2 <- stats::qbinom(probs, DP[i], peaks[[j]][2])
        
        if ((NV[i] >= quantiles_1[1]) && (NV[i] <= quantiles_1[2])) {
          return("omega1")
        }else if ((NV[i] >= quantiles_2[1]) && (NV[i] <= quantiles_2[2])){
          return("omega2")
          
        }
      }
    }
  }) %>% unlist()
  
  
  # print(paste0("alpha_beta_all ",length(alpha_beta_all)))
  df <- data.frame(alpha_beta_all = alpha_beta_all, segment_id = input_data$seg_assignment)
  
  # add karyotype for each row by respecting the segment id contraint
  for (i in 1:nrow(df)){
    df$karyotype[i] = karyotype[input_data$seg_assignment[i]]
  }
  
  
  # Group by segment_id and calculate the proportions of alpha and beta
  proportions <- df %>%
    dplyr::group_by(.data$segment_id, .data$karyotype) %>%
    dplyr::summarise(
      proportion_alpha = mean(alpha_beta_all == "omega1"),
      proportion_beta = mean(alpha_beta_all == "omega2"),
      .groups = 'drop'                                                    #check dplyr
    )
  
  
  
  myReps <- function(x, y, n) rep(x, (x %in% y) * (n-1) +1)
  phi = myReps(1/K, 1/K, K)
  
  
  # 1) calcolo prima i tau e poi clusterizzo (perché dipendono dal karyotypo)
  # Initialize the tau_posterior column
  proportions$tau_posterior <- NA  # Create an empty column in the dataframe
  
  # Loop through each row to calculate tau_posterior
  for (i in 1:nrow(proportions)) {
    if (proportions$karyotype[i] == '2:1') {
      proportions$tau_posterior[i] <- 3 * (proportions$proportion_beta[i]+0.0001) /
        (2 * (proportions$proportion_beta[i]+0.0001) + (proportions$proportion_alpha[i]+0.0001))
    } else {
      proportions$tau_posterior[i] <- 2 * (proportions$proportion_beta[i]+0.0001) /
        (2 * (proportions$proportion_beta[i]+0.0001) + (proportions$proportion_alpha[i]))
    }
  }
  
  # Check if tau_posterior contains any NA values
  if (any(is.na(proportions$tau_posterior))) {
    stop("tau_posterior contains NA values, check the loop calculations.")
  }
  
  # Apply fuzzy c-means clustering
  # res.fcm <- ppclust::fcm(as.matrix(proportions$tau_posterior), centers = K)
  
  res.fcm <- run_fcm_with_fallback(proportions$tau_posterior, K)
  
  
  
  # print(paste0("dim(proportions$tau_posterior): ", length(proportions$tau_posterior) ))
  if (length(proportions$tau_posterior)==K){
    init_taus <- (proportions$tau_posterior)
  }else{
    init_taus <- c(res.fcm$v)                                                                 #getter?
  }
  cat(paste0("init_taus from clustering  ",init_taus))
  
  # init_taus <- c(res.fcm$v)
  init_taus[init_taus >= 0.88] <- 0.88  # Check for elements greater than 1 and replace them with 1 otherwise fit fails
  init_taus[init_taus == 0] <- 0.07
  init_w <- as.matrix(res.fcm$u)
  
  
  epsilon <- 1e-4
  if (all(init_w > 0.5) ){
    perturbed_probabilities <- init_w - epsilon
  } else {
    perturbed_probabilities <- init_w + epsilon
  }
  
  
  # Rinormalizza i pesi lungo l'asse 1 (per riga)
  normalized_probabilities <- (apply(perturbed_probabilities, 1, function(x) x / sum(x)))
  # Verifica che i pesi siano stati rinormalizzati correttamente
  #print(rowSums(normalized_probabilities))
  # Stampa i pesi normalizzati
  init_w <- (normalized_probabilities)
  if (K==1){
    init_w = t(init_w)
  }
  
  inits <- list(w = t(init_w), tau = init_taus, phi=phi, kappa=kappa)
  
  # cat(paste0("w = ",inits$w,"\n "))
  # cat(paste0("tau = ", inits$tau, "\n "))
  # cat(paste0("phi = ", inits$phi, "\n "))
  # cat(paste0("kappa = ", inits$kappa))
  #
  return(inits)
  
}





# Function to adjust cluster results
adjust_fcm_results <- function(result, original_K) {
  current_K <- nrow(result$v)
  
  # If current K is smaller than original K, expand results
  if (current_K < original_K) {
    message(paste("Adjusting results to match original K =", original_K, "..."))
    
    # Adjust cluster centers (res.fcm$v)
    sampled_v <- result$v[sample(1:current_K, original_K, replace = TRUE), , drop = FALSE]
    
    # Adjust membership matrix (res.fcm$u)
    sampled_u <- result$u[, sample(1:current_K, original_K, replace = TRUE), drop = FALSE]
    
    # Normalize the membership matrix rows to ensure they sum to 1
    sampled_u <- sweep(sampled_u, 1, rowSums(sampled_u), FUN = "/")

    # Update result
    result$v <- sampled_v
    result$u <- sampled_u
  }
  
  return(result)
}


# data <- proportions$tau_posterior

# Main function to run FCM with fallback and adjustment
run_fcm_with_fallback <- function(data, K, min_K = 2) {
  original_K <- K
  current_K <- K
  result <- NULL
  
  if (current_K == 1){
    result <- ppclust::fcm(as.matrix(data), centers = current_K)
  }
  
  while (is.null(result) && current_K > 1) {
    print(current_K)
    
   current_K <- tryCatch({
      # Attempt to run FCM
      result <- ppclust::fcm(as.matrix(data), centers = current_K)
      current_K <- result
    }, error = function(e) {
      if (grepl("too few positive probabilities", e$message)) {
        message(paste("Error encountered with K =", current_K, ". Trying K =", current_K - 1, "..."))
        
        # Update current_K
        current_K <- current_K - 1
        
      } else {
        stop(e)  # Re-throw unexpected errors
      }
    },
    finally = {
      #pass
    })
  }
  
    result <- adjust_fcm_results(result, original_K)

  
  return(result)
}





