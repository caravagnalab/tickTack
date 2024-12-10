#' Plot posterior weights distributions obtained from the hierarchical model fit
#'
#' @param results list of 4: $data, $draws_and_summary, $log_lik_matrix_list and $elbo_iterations 
#' @param K index of inference whose results want to be plotted  
#' 
#' @return areas_tau 
#' @export

plot_posterior_weights_h <- function(results, K){
  
  draws <- results$draws_and_summary[[K]]$draws
  S = nrow(results$data$accepted_cna)
  
  intervals_weigths_per_tau <- list()
  for (k in 1:K){
    names_weights <- paste("w[",1:S,",", k, "]", sep = "")
    p <- bayesplot::mcmc_intervals(draws, pars = names_weights, point_est = "median", prob = 0.8, prob_outer = 0.95)+
      ggplot2::labs(
        title =  stringr::str_wrap( paste0("Posterior distributions of the weigths for tau ",k), width = 30 + K + sqrt(S)),
        subtitle = "With median and 80% and 95% intervals"
      )
    intervals_weigths_per_tau[[k]] <- ggplot2::ggplotGrob(p)
  }
  p <- gridExtra::grid.arrange(grobs = intervals_weigths_per_tau, ncol = K) #add global title
  
  return(p)
}


