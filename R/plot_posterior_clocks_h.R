#' Plot posterior clocks distributions obtained from the hierarchical model fit
#'
#' @param results list of 4: $data, $draws_and_summary, $log_lik_matrix_list and $elbo_iterations 
#' @param K index of inference whose results want to be plotted   
#' 
#' @return areas_tau 
#' @export
#' 
plot_posterior_clocks_h <- function(results, K){
  draws <- results$draws_and_summary[[K]]$draws
  
  names <- paste("tau[", 1:K, "]", sep = "")
  
  areas_tau <- bayesplot::mcmc_areas(
    draws,
    pars = names,
    prob = 0.8, # 80% intervals
    prob_outer = 0.95, # 99%
    point_est = "median"
  )+
    ggplot2::labs(
      title = "Approximate Posterior distributions",
      subtitle = "With median and 80% and 95% intervals"
    )+
    ggplot2::xlim(0, 1) 
 
   return(areas_tau)
}


