#' Plot elbo Function
#'
#' @description Plot the ELBO behaviour for each fit fixing the number of cluster.
#' @param elbo_iteration data.frame':	#iterations until convergence of elbo obs. of  2 variables: iteration: int, elbo     : num  
#' one of the element list elbo_iteration obtained from fit_h "results"
#' 
#' @return p 
#' @export
plot_elbo_h <- function(elbo_iteration){
  
  data <- elbo_iteration
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$iteration, y = .data$elbo)) +
    ggplot2::geom_line(color = "blue", size = 0.2) +  
    # geom_point(color = "black", size = 1.5) +  
    ggplot2::labs(
      title = paste0("ELBO Convergence for K = ", nrow(data) - 1),  
      x = "Iteration",
      y = "ELBO"
    ) +
    # theme_classic(base_size = 12) +  
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 14), 
      axis.text = ggplot2::element_text(size = 12), 
      panel.grid = ggplot2::element_line(color = "gray100"),
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(min(data$iteration), max(data$iteration), by = 1)  # Ensure only integers appear
    ) +
    ggplot2::scale_y_continuous(labels = scales::comma)  
  
  
  return(p)
}
