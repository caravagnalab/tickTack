#' plot_timing_h Function
#'
#' This function obtains the list of input data to be used in the stan model.
#' @param results   list(data = input_data_list, results = results, output_files_list = output_files_list)
#' @param colour_by chr: default =  "karyotype"
#' @param K mun: number of clocks
#' @param split_contiguous_segments option to plot segments' setalarion lines
#' 
#' @return p : plot of inference results with credibility intervals in the chromosome absolute positions
#' 
#' @keywords plot
#' @export

plot_timing_h = function(results, K, colour_by = "karyotype", split_contiguous_segments = TRUE) {
  
  #segments <- x_segments[ x_segments$chr %in% results$data$accepted_cna$chr, ]
  
  segments <- results$data$accepted_cna
  segments$from = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[2]}) %>% unlist() %>% as.numeric()
  segments$to = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[3]}) %>% unlist() %>% as.numeric()
  
  reference_genome <- tickTack::chr_coordinates_GRCh38
  
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  
  absoulte_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])
  
  summarized_results <- results$draws_and_summary[[K]]$summarized_results %>%
    dplyr::mutate(from = absoulte_segments[.data$segment_id,]$from) %>%
    dplyr::mutate(to = absoulte_segments[.data$segment_id,]$to) %>%
    dplyr::mutate(tau_mean = ifelse(.data$clock_mean  < 1, .data$clock_mean , 1)) %>%
    dplyr::mutate(tau_high = ifelse(.data$clock_high < 1, .data$clock_high, 1)) %>%
    dplyr::mutate(tau_low = .data$clock_low) 
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = summarized_results, ggplot2::aes(xmin=.data$from, xmax=.data$to, ymin=.data$tau_low, ymax=.data$tau_high, fill = as.factor(.data[[colour_by]])), alpha = .5) +
    ggplot2::geom_segment(data = summarized_results, ggplot2::aes(y = .data$tau_mean, yend = .data$tau_mean, x = .data$from, xend = .data$to)) +
    ggplot2::scale_x_continuous(breaks = reference_genome$to, labels = gsub("chr", "", reference_genome$chr)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::lims(y = c(0,1)) +
    ggplot2::labs(x = "Chromosome", y = bquote("Pseudotime"~tau), fill=colour_by)
  # scale_fill_manual(values = c("forestgreen", "indianred3", "steelblue"), name = "")
  
  if (split_contiguous_segments) {
    p <- p +
      ggplot2::geom_segment(data = summarized_results, ggplot2::aes(y = .data$tau_low, yend = .data$tau_high, x = .data$from, xend = .data$from), linetype = "dashed", color = "darkslategray") +
      ggplot2::geom_segment(data = summarized_results, ggplot2::aes(y = .data$tau_low, yend = .data$tau_high, x = .data$to, xend = .data$to), linetype = "dashed", color = "darkslategray")
  }
  
  return(p)
}









