#' Plot the behavior of the model selection scores vs the number of clusters
#'
#' @description Plot the BIC, Log Likelihood, ICL scores obtatined after the fit for each number of cluster.
#'
#' @param model_selection_tibble a tibble with 3 scores and k_max values, one for each inference 
#' @param best_K the integer corresponding to the best number of components
#' 
#' @return model_selection_plot 
#' @export

plot_model_selection_h <- function(model_selection_tibble, best_K) {
  create_plot <- function(data, x, y, best_K, y_label, score_name) {
    ggplot2::ggplot(data, ggplot2::aes(x = !!ggplot2::sym(x), y = !!ggplot2::sym(y))) +
      ggplot2::geom_line(color = "steelblue", linewidth = 1) +
      ggplot2::geom_point(color = "steelblue", size = 3) +
      ggplot2::geom_point(data = data[data[[x]] == best_K, ],
                          ggplot2::aes(x = !!ggplot2::sym(x), y = !!ggplot2::sym(y)), color = "firebrick", size = 4) +
      ggplot2::labs(y = y_label, x = "Number of Clusters (K)", title = score_name) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
        axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
        panel.grid.minor = ggplot2::element_blank()
      )
  }
  
  bic_plot <- create_plot(model_selection_tibble, "K", "BIC", best_K, "BIC", "BIC vs K")
  log_lik_plot <- create_plot(model_selection_tibble, "K", "Log_lik", best_K, "Log-Likelihood", "Log-Likelihood vs K")
  icl_plot <- create_plot(model_selection_tibble, "K", "ICL", best_K, "ICL", "ICL vs K")
  
  model_selection_plot <- (bic_plot) / icl_plot / (log_lik_plot)  +
    patchwork::plot_annotation(
      title = "Model Selection Graphs: Scores vs Number of Clusters",
      caption = "Source: Your Data",
      theme = ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
          plot.caption = ggplot2::element_text(size = 10, hjust = 0.5, margin = ggplot2::margin(t = 10))
        )
    )
  
  return(model_selection_plot)
}
