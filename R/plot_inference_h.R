plot_inference_h = function(results, x_segments, input_data, colour_by = "karyotype", K) {
  
  segments <- x_segments[ x_segments$chr %in% results$data$accepted_cna$chr, ]
  
  # reference_genome <- CNAqc::chr_coordinates_GRCh38
  load("../tickTack/data/chr_coordinates_GRCh38.rda")
  reference_genome <- chr_coordinates_GRCh38
  
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
    ggplot2::labs(x = "chromosome", y = bquote(tau))
  # scale_fill_manual(values = c("forestgreen", "indianred3", "steelblue"), name = "")
  
  return(p)
}
