#' prepare_input_data Function
#'
#' This function obtains the list of input data to be used in the stan model.
#' @param results   list(data = input_data_list, results = results, output_files_list = output_files_list)
#' @param x_segments  tibble((S3: tbl_df/tbl/data.frame) chr, from, to, Major, minor, total_cn data$cna) 
#' @param input_data List of 7: $S: int, $N: int, $karyotype: num (0 or 1), $seg_assignment: num, $peaks:List of N of num (1:2), $NV: num, $DP: num
#' @param colour_by chr: default =  "karyotype"
#' @param K mun: number of clocks
#' 
#' @return p : plot of inference results with credibility intervals in the chromosome absolute positions
#' 
#' @keywords plot
#' @export

plot_inference_h = function(results, x_segments, input_data, colour_by = "karyotype", K) {
  
  segments <- x_segments[ x_segments$chr %in% results$data$accepted_cna$chr, ]
  
  # reference_genome <- CNAqc::chr_coordinates_GRCh38
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
    ggplot2::labs(x = "chromosome", y = "clock estimate") #bquote(tau)
  # scale_fill_manual(values = c("forestgreen", "indianred3", "steelblue"), name = "")
  
  return(p)
}
