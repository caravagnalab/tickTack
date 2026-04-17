
plot_vaf_mobsterlike = function(x, K) {
  input_data = x$results_timing$data$input_data
  cna_data = x$results_timing$draws_and_summary[[K]]$summarized_results

  cna_data_with_vaf = lapply(cna_data$segment_id, function(i) {
    idxs = which(input_data$seg_assignment == i)
    DP = input_data$DP[idxs]
    NV = input_data$NV[idxs]
    VAF = NV / DP

    cna_data %>% dplyr::filter(segment_id == i) %>%
      dplyr::bind_cols(dplyr::tibble(VAF = VAF))
  }) %>% do.call(dplyr::bind_rows, .)

  cna_data_with_vaf %>%
    dplyr::mutate(clock_rank = dplyr::dense_rank(clock_mean)) %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = VAF, fill = factor(clock_rank))) +
    ggplot2::geom_histogram(binwidth = 0.01) +
    ggplot2::facet_grid(clock_rank~karyotype) +
    ggplot2::theme_light() +
    ggplot2::scale_x_continuous(limits = c(-0.01,1.01)) +
    labs(fill = "Cluster", col = "Cluster")
}

plot_vaf_ecdf = function(x, K) {
  input_data = x$results_timing$data$input_data
  cna_data = x$results_timing$draws_and_summary[[K]]$summarized_results

  cna_data_with_vaf = lapply(cna_data$segment_id, function(i) {
    idxs = which(input_data$seg_assignment == i)
    DP = input_data$DP[idxs]
    NV = input_data$NV[idxs]
    VAF = NV / DP

    cna_data %>% dplyr::filter(segment_id == i) %>%
      dplyr::bind_cols(dplyr::tibble(VAF = VAF))
  }) %>% do.call(dplyr::bind_rows, .)

  cna_data_with_vaf %>%
    dplyr::mutate(clock_rank = dplyr::dense_rank(clock_mean)) %>%
    ggplot2::ggplot(aes(x = VAF, color = factor(clock_rank), group = segment_name)) +
    ggplot2::stat_ecdf() +
    ggplot2::facet_grid(clock_rank~karyotype) +
    ggplot2::theme_light() +
    ggplot2::scale_x_continuous(limits = c(-0.01,1.01)) +
    labs(fill = "Cluster", col = "Cluster")
}
