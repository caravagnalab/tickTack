
#' Plot timing of clonal peaks in cancer genome sequencing data
#'
#' This function generates a plot showing the timing of clonal peaks in cancer genome sequencing data.
#'
#' @param fit_results A list containing the results of the `fit_timing` function, specifically `summarized_results`.
#' @param segments A data frame containing segment information with columns `chr`, `from`, `to`, `Major`, and `minor`.
#' @param colour_by A character string specifying the variable to color the plot by. Default is "karyotype".
#' @param ref Reference genome desired. Either 'GRCh38' or 'hg19'
#' @return A ggplot object showing the timing of clonal peaks.
#' @export
plot_timing = function(fit_results, segments, colour_by = "karyotype", ref = 'GRCh38') {
  reference_genome <- get_reference(ref)

  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  absoulte_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
           to = .data$to + vfrom[.data$chr])

  summarized_results <- fit_results$summarized_results %>%
    dplyr::mutate(from = absoulte_segments[.data$segment,]$from) %>%
    dplyr::mutate(to = absoulte_segments[.data$segment,]$to) %>%
    dplyr::mutate(tau_mean = ifelse(.data$tau_mean < 1, .data$tau_mean, 1)) %>%
    dplyr::mutate(tau_high = ifelse(.data$tau_high < 1, .data$tau_high, 1))

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = summarized_results, ggplot2::aes(xmin=.data$from, xmax=.data$to, ymin=.data$tau_low, ymax=.data$tau_high, fill = as.factor(.data[[colour_by]])), alpha = .5) +
    ggplot2::geom_segment(data = summarized_results, ggplot2::aes(y = .data$tau_mean, yend = .data$tau_mean, x = .data$from, xend = .data$to)) +
    ggplot2::scale_x_continuous(breaks = reference_genome$to, labels = gsub("chr", "", reference_genome$chr)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::lims(y = c(0,1)) +
    ggplot2::labs(x = "chromosome", y = bquote(tau))
  # scale_fill_manual(values = c("forestgreen", "indianred3", "steelblue"), name = "")

  return(p)
}
