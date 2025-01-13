#' model_selection Function
#'
#' Perform model selection among the models fit with varying number of mixture components (number of clocks).
#' @param input_data input data    
#' @param purity sample purity 
#' @param alpha description
#' 
#' @keywords fit
#'
#' @return optimal_k: suggested number from c-means
#'
#' @export
get_k_max_k_means = function(input_data, purity, alpha = 0.05){
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
  
  
  inertia = list()
  list_k = c(1:input_data$S)
  wss_df <- data.frame()
  
  for(k in list_k){
    results <- ppclust::fcm(x=proportions$tau_posterior, centers=k, numseed=1)
    wss_df <- dplyr::bind_rows( wss_df, data.frame(K=k,inertia=results$sumsqrs$tot.within.ss))
  }
  
  total_inertia <- wss_df$inertia[1]  
  wss_df$cumulative_explained <- (total_inertia - wss_df$inertia) / total_inertia  
  explained_threshold <- 0.99 
  optimal_k <- min(wss_df$K[wss_df$cumulative_explained >= explained_threshold], na.rm = TRUE)
  cat("Optimal number of clusters (K):", optimal_k, "\n")
  
  # scree_plot <- ggplot(wss_df, aes(x = list_k, y = inertia, group = 1)) +
  #   geom_point(size = 4)+
  #   geom_line() +
  #   scale_x_continuous(breaks = list_k) +
  #   xlab('Number of clusters')
  # scree_plot
  
  return(optimal_k)
  
}