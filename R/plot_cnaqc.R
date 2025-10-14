# library(CNAqc)
# library(tickTack)
# library(dplyr)
# library(stringr)
# library(fuzzyjoin)
# library(ggplot2)
# library(patchwork)
# library(data.table)
# library(tidyverse)
# library(mobster)

ALPHA = .8
# color = c(
#   'AmplificationTimeR' = alpha('forestgreen', alpha = ALPHA),
#   'MutationTimeR' = alpha('steelblue', alpha = ALPHA),
#   'tickTack' = alpha('orange', alpha = ALPHA),
#   'tickTackH' = alpha('firebrick3', alpha = ALPHA)
# )

color = c(
  'AmplificationTimeR' = ggplot2::alpha('#868686', alpha = ALPHA),
  'MutationTimeR' = ggplot2::alpha('#EFC000', alpha = ALPHA),
  'tickTack baseline' = ggplot2::alpha('#7AA6DC', alpha = ALPHA),
  'tickTack full' = ggplot2::alpha('#CD534C', alpha = ALPHA)
)

k_colors = list(
  '1:1' = '#228B22CC',
  '1:0' = 'steelblue',
  '2:0' = 'turquoise4',
  '2:1' = ggplot2::alpha('orange', .8),
  '2:2' = 'firebrick3'
)


#' Convert model names to standardized names
#'
#' @param n A character string representing the model name.
#' @return A standardized model name as a character string.
#' @export

convert_name = function(n) {
  if (grepl("AmpTime", n)) return("AmplificationTimeR")
  if (grepl("MutTimeR", n)) return("MutationTimeR")
  if (grepl("tickTack_h", n)) return("tickTack full")
  if (grepl("tickTack", n)) return("tickTack baseline")
  stop("error name not recognized")
}




#' Get reference genome coordinates
#'
#' @param ref A character string specifying the reference genome ('hg19' or 'hg38').
#' @return A data frame with chromosome coordinates for the specified reference genome.
#' @export
get_reference <- function(ref) {
  if (ref %in% c("hg19", "GRCh37")) {
    return(CNAqc::chr_coordinates_hg19)
  }
  
  if (ref %in% c("hg38", "GRCh38")) {
    return(CNAqc::chr_coordinates_GRCh38)
  }
  
  stop("Available references: 'hg19' (or 'GRCh37') and 'hg38' (or 'GRCh38').")
}

#' Extract driver mutations from mutation data
#'
#' @param x A CNAqc object containing mutation data.
#' @param chromosomes A character vector of chromosome names to filter for (default: 'chr1' to 'chr22', 'X', 'Y').
#' @param which A character string specifying whether to extract 'VAF' or 'CCF' values.
#' @return A data frame of driver mutations filtered by chromosome and data type.
#' @export

get_drivers = function(x,
                       chromosomes = paste0('chr', c(1:22, 'X', 'Y')),
                       which = 'VAF') {
  if(!has_driver_data(x)) return(NULL)
  
  # CCF?
  if (all(is.null(x$CCF_estimates)) & which == "CCF")
  {
    warning("Input does not have CCF estimates, see ?compute_CCF to determine CCF values.")
    return(NULL)
  }
  
  drivers_list = NULL
  
  if(which == "VAF")
    drivers_list = x$mutations %>%
    dplyr::filter(.data$is_driver, .data$chr %in% chromosomes)
  
  if(which == "CCF")
    drivers_list = CNAqc::CCF(x) %>%
    dplyr::filter(.data$is_driver, .data$chr %in% chromosomes)
  
  return(drivers_list)
}



#' Check if mutation data has driver annotations
#'
#' @param x A CNAqc object containing mutation data.
#' @return Logical (`TRUE` if driver annotation columns exist, otherwise `FALSE`).
#' @export

has_driver_data = function(x)
{
  
  cn = colnames(x$mutations)
  
  if (all(c("is_driver", "driver_label") %in% cn)) return(TRUE)
  
  if ("is_driver" %in% cn & !("driver_label" %in% cn)){
    cli::cli_warn("Column 'is_driver' is annotated but 'driver_label' no -- did you try to add drivers data?")
  }
  
  if ("driver_label" %in% cn & !("is_driver" %in% cn)){
    cli::cli_warn("Column 'driver_label' is annotated but 'is_driver' no -- did you try to add drivers data?")
  }
  
  return(FALSE)
}



#' Add driver mutation annotations to a segment plot
#'
#' @param x A CNAqc object containing mutation data.
#' @param drivers_list A data frame of driver mutations.
#' @param base_plot A ggplot object representing the base plot.
#' @return A modified ggplot object with driver mutation annotations.
#' @export

add_drivers_to_segment_plot = function(x, drivers_list, base_plot)
{
  if(is.null(drivers_list)) return(base_plot)
  if (nrow(drivers_list) == 0) return(base_plot)
  
  L = ggplot2::ggplot_build(base_plot)$layout$panel_params[[1]]
  
  drivers_list = CNAqc:::relative_to_absolute_coordinates(
    x,
    drivers_list %>% dplyr::filter(.data$is_driver)
  )
  
  drivers_list$y = L$y.range[2] * 0.9
  
  base_plot +
    ggplot2::geom_vline(
      data = drivers_list,
      show.legend = FALSE,
      ggplot2::aes(xintercept = .data$from),
      linetype = 'dashed',
      color = 'black',
      size = .3
    ) +
    ggrepel::geom_label_repel(
      data = drivers_list,
      ggplot2::aes(
        x = .data$from,
        y = .data$y,
        label = .data$driver_label,
      ),
      ylim = c(L$y.range[2] * .9, NA),
      size = 2,
      nudge_y = 0,
      nudge_x = 0,
      show.legend = FALSE
    ) +
    ggplot2::coord_cartesian(clip = 'off')
}



#' Plot 
#'
#' @param x A CNAqc object containing mutation data.
#' @param chromosomes A character vector of chromosome names to filter for (default: 'chr1' to 'chr22', 'X', 'Y').
#' @param add_mobster TRUE or FALSE
#' @param max_Y_height heiht plot
#' @param cn type of cn 
#' @param highlight which karyotype to highlight
#' @param highlight_QC FALSE
#' @return Plot.
#' @export
plot_cnaqc <- function(x, chromosomes = paste0('chr', c(1:22)), add_mobster=FALSE, max_Y_height = 6, cn = 'absolute', highlight = x$most_prevalent_karyotype, highlight_QC = FALSE) {
  
  results_model_selection <- model_selection_h(x$results, n_components = 0)
  x$K = results_model_selection$best_K
  
  K = x$K
  
  cnaqc_x = CNAqc::init(mutations = x$mutations, cna = x$cna, purity = x$metadata$purity, ref = x$reference_genome)
  
  plot_CNA = plot_segments_tick_tack_CN(cnaqc_x, K = K) +
    ggplot2::theme(legend.position='right', panel.spacing = ggplot2::unit(0, "lines")) +
    ggplot2::labs(caption = NULL) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  
  x$reference_genome = "hg38"
  data_plot <- plot_segments_tick_tack_data(x, K = K) +
    ggplot2::theme(legend.position='right',panel.spacing = ggplot2::unit(0, "lines")) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  
  timing_plot <- plot_segments_tick_tack(x, colour_by = "clock_mean", K = K) +
    ggplot2::theme(legend.position='right',panel.spacing = ggplot2::unit(0, "lines"))
  
  if(add_mobster){
      vaf_plot <- plot_vaf(x, K) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = "white", size = 20))
  
  # segment_plot <- plot_segments_h(x, chromosomes, max_Y_height, cn, highlight, highlight_QC) +
  #   ggplot2::theme(axis.title.x = element_blank())  # Keep chromosome labels only on this plot
  
  pA = timing_plot + CNAqc:::my_ggplot_theme() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      legend.position = "left"
    )
  pB = plot_CNA + CNAqc:::my_ggplot_theme() +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      legend.position = "left"
    )
  pC = data_plot + CNAqc:::my_ggplot_theme() +
    ggplot2::theme(legend.position = "left") +
    ggplot2::labs(y = "VAF")
  pD = vaf_plot + CNAqc:::my_ggplot_theme()
  
  des_left = "
  AAAA#
  AAAAE
  AAAAE
  AAAAE
  AAAAE
  BBBBE
  BBBBE
  BBBBE
  CCCCE
  CCCCE
  DDDDD"
  
  pp = pA + pB + pC + patchwork::guide_area() + pD +
    patchwork::plot_layout(design = des_left, guides = "collect") &
    ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")
  pp
  } else {

    # segment_plot <- plot_segments_h(x, chromosomes, max_Y_height, cn, highlight, highlight_QC) +
    #   ggplot2::theme(axis.title.x = element_blank())  # Keep chromosome labels only on this plot
    
    pA = timing_plot + CNAqc:::my_ggplot_theme() +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.position = "left"
      )
    pB = plot_CNA + CNAqc:::my_ggplot_theme() +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.position = "left"
      )
    pC = data_plot + CNAqc:::my_ggplot_theme() +
      ggplot2::theme(legend.position = "left") +
      ggplot2::labs(y = "VAF")
    pD = ggplot2::ggplot() + ggplot2::theme_void() + CNAqc:::my_ggplot_theme()
    
    des_left = "
  AAAA#
  AAAAE
  AAAAE
  AAAAE
  AAAAE
  BBBBE
  BBBBE
  BBBBE
  CCCCE
  CCCCE
  DDDDD"
    
    pp = pA + pB + pC + patchwork::guide_area() + pD +
      patchwork::plot_layout(design = des_left, guides = "collect") &
      ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")
    pp
  }

}






#' Plot segmentation results from tickTack analysis
#'
#' @param x An object containing segmentation results.
#' @param colour_by A character string specifying the variable to color segments by (default: 'clock_mean').
#' @param K Numeric, number of clusters (default: 1).
#' @return A ggplot object visualizing tickTack segmentation.
#' @export

plot_segments_tick_tack <- function(x, colour_by = "clock_mean", K = 1) {
  
  results <- x$results_timing
  reference_genome <- get_reference(x$reference_genome)
  segments <- results$data$accepted_cna
  segments$from = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[2]}) %>% unlist() %>% as.numeric()
  segments$to = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[3]}) %>% unlist() %>% as.numeric()
  
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
  ##############
  
  k_colors = list(
    '1:1' = '#228B22CC',
    '1:0' = 'steelblue',
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )
  
  # six_color_palette <- six_color_palette <- RColorBrewer::brewer.pal(6, "RdYlBu")
  summarized_results <- summarized_results %>%
    dplyr::mutate(tau_cluster = factor(
      .data[[colour_by]],
      levels = sort(unique(.data[[colour_by]])),
      labels = paste(seq_along(sort(unique(.data[[colour_by]]))))
    ))
  
  capt_label = paste0("Tumour type ", x$ttype,
                      " Ploidy ", x$ploidy, "; Purity  ", x$purity,
                      '; n = ', x$results_timing$data$input_data$N, ' accepted mutations in ',
                      x$results_timing$data$input_data$S,
                      ' accepted segments'
  )
  
  my_palette <- c(  "#66a61e",  "#7570b3", "#e7298a", "#1b9e77", "#d95f02")
  
  p <- CNAqc:::blank_genome('hg19',chromosomes = paste0("chr", c(1:22)))+
    ggplot2::geom_segment(
      data = summarized_results,
      ggplot2::aes(
        y = .data$clock_mean,
        yend = .data$clock_mean,
        x = .data$from,
        xend = .data$to
      )
    ) +
    ggplot2::geom_rect(
      data = summarized_results,
      ggplot2::aes(
        xmin = .data$from,
        xmax = .data$to,
        ymin = .data$tau_low,
        ymax = .data$tau_high,
        fill = .data$tau_cluster
      ),
      alpha = 0.5
    ) +
    ggplot2::scale_fill_manual(values = my_palette)+
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 0)
    ) +
    ggplot2::ggtitle(label = "", subtitle = capt_label) +
    ggplot2::labs(
      x = "Chromosome",
      y = bquote("Pseudotime"~tau),
      fill = "Cluster"
    )
  
  drivers_list = get_drivers(x, chromosomes = paste0('chr', c(1:22)))
  base_plot = add_drivers_to_segment_plot(x, drivers_list = drivers_list, p)
  base_plot
}



#' Plot copy number segments from CNAqc
#'
#' @param cnaqc_x A CNAqc object containing CNA data.
#' @param K Numeric, number of clusters.
#' @param max_alleles Numeric, maximum number of alleles to plot (default: 6).
#' @param chromosomes A character vector specifying which chromosomes to plot.
#' @return A ggplot object showing CNA segments.
#' @export

plot_segments_tick_tack_CN <- function(cnaqc_x, K = K, max_alleles = 6, chromosomes = paste0("chr", c(1:22))) {
  
  reference_genome <- get_reference(cnaqc_x$reference_genome)
  
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  
  segments <- cnaqc_x$cna %>%
    dplyr::filter(.data$chr %in% chromosomes)
  
  absolute_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])
  
  absolute_segments$karyotype = paste(absolute_segments$Major, absolute_segments$minor, sep = ":")
  absolute_segments = absolute_segments %>%
    dplyr::filter(.data$Major <= max_alleles & .data$minor <= max_alleles)
  
  k_colors = list(
    '1:1' = '#228B22CC',
    '1:0' = 'steelblue',
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )
  p = CNAqc:::blank_genome(cnaqc_x$reference_genome, chromosomes = paste0("chr", c(1:22)))
  
  # Add shadows
  p = p + ggplot2::geom_rect(
    data = absolute_segments %>% dplyr::filter(.data$karyotype %in% c('2:0', '1:0', '1:1', '2:1', '2:2')),
    ggplot2::aes(
      xmin = .data$from,
      xmax = .data$to,
      ymin = -Inf,
      ymax = Inf,
      fill = factor(.data$karyotype, levels = c('2:0', '1:0', '1:1', '2:1', '2:2'))
    ),
    alpha = .3
  ) +
    ggplot2::scale_fill_manual(values = k_colors) +
    ggplot2::guides(fill = ggplot2::guide_legend('', override.aes = list(alpha = 1)))
  
  # Add segments
  p = p +
    ggplot2::geom_segment(
      data = absolute_segments %>%
        dplyr::mutate(Major = as.numeric(.data$Major) + .1, minor = as.numeric(.data$minor) - .1) %>%
        dplyr::select(.data$karyotype, .data$chr, .data$from, .data$to, .data$Major, .data$minor) %>%
        tidyr::pivot_longer(!c(.data$karyotype, .data$chr, .data$from, .data$to)),
      ggplot2::aes(
        x = .data$from,
        xend = .data$to,
        y = .data$value,
        color = as.factor(.data$name)
      ),
      size=1.5
    ) +
    ggplot2::scale_color_manual(values = c("Major" = "red", "minor" ="steelblue")) +
    ggplot2::guides(color = ggplot2::guide_legend('')) +
    ggplot2::labs(y = "Allel count") +
    ggplot2::ylim(c(0, max_alleles))
  
  p
}



#' Plot mutation data from tickTack analysis
#'
#' @param x An object containing mutation data.
#' @param colour_by A character string specifying the variable to color points by (default: 'clock_mean').
#' @param K Numeric, number of clusters.
#' @return A ggplot object showing mutation data.
#' @export


plot_segments_tick_tack_data <- function(x, colour_by = "clock_mean", K = K) {
  
  mutations <- x$mutations
  results <- x$results_timing
  reference_genome <- get_reference(x$reference_genome)
  
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  
  absolute_mutations <- mutations  %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])
  
  segments <- results$data$accepted_cna
  segments$from = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[2]}) %>% unlist() %>% as.numeric()
  segments$to = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[3]}) %>% unlist() %>% as.numeric()
  
  absolute_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])
  
  summarized_results <- results$draws_and_summary[[K]]$summarized_results %>%
    dplyr::mutate(from = absolute_segments[.data$segment_id,]$from) %>%
    dplyr::mutate(to = absolute_segments[.data$segment_id,]$to) %>%
    dplyr::mutate(tau_mean = ifelse(.data$clock_mean  < 1, .data$clock_mean , 1)) %>%
    dplyr::mutate(tau_high = ifelse(.data$clock_high < 1, .data$clock_high, 1)) %>%
    dplyr::mutate(tau_low = .data$clock_low)
  
  accepted_mutations = data.frame()
  for (segment_idx in 1:nrow(summarized_results)) {
    segment <- summarized_results[segment_idx, ]
    # print(segment$chr)
    segment_mutations <- absolute_mutations %>%
      dplyr::filter(.data$chr == segment$chr, .data$from > segment$from, .data$to < segment$to) %>%
      tidyr::drop_na(.data$DP)
    segment_mutations <- segment_mutations %>% dplyr::mutate(karyotype = segment$karyotype)
    # print(nrow(segment_mutations))
    # if (nrow(segment_mutations)> 40){
    accepted_mutations <- dplyr::bind_rows(accepted_mutations, segment_mutations)
    # }
  }
  
  matched_mutations <- accepted_mutations
  
  k_colors = list(
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )
  
  CNAqc:::blank_genome(ref=x$reference_genome, chromosomes = paste0("chr", c(1:22)))+
    ggplot2::geom_point(
      data = matched_mutations,
      ggplot2::aes(
        x = .data$from,
        y = .data$NV / .data$DP,
        color = as.factor(.data$karyotype)
      ),
      alpha = 0.5,
      size=0.1
    ) +
    ggplot2::scale_color_manual(values = k_colors, name = "CN") +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::labs(
      x = "Chromosome",
      y = bquote("Variant Allele Frequency (VAF)")
    )
}



#' Plot variant allele frequency (VAF) distributions
#'
#' @param x An object containing mutation and segmentation data.
#' @param K Numeric, number of clusters.
#' @return A ggplot object visualizing VAF distributions.
#' @export

plot_vaf = function(x, K){
  
  mutations <- x$mutations
  results <- x$results_timing
  reference_genome <- get_reference(x$reference_genome)
  purity = x$purity
  
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  
  absolute_mutations <- mutations  %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])
  
  segments <- results$data$accepted_cna
  segments$from = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[2]}) %>% unlist() %>% as.numeric()
  segments$to = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[3]}) %>% unlist() %>% as.numeric()
  
  absolute_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])
  
  summarized_results <- results$draws_and_summary[[K]]$summarized_results %>%
    dplyr::mutate(from = absolute_segments[.data$segment_id,]$from) %>%
    dplyr::mutate(to = absolute_segments[.data$segment_id,]$to) %>%
    dplyr::mutate(tau_mean = ifelse(.data$clock_mean  < 1, .data$clock_mean , 1)) %>%
    dplyr::mutate(tau_high = ifelse(.data$clock_high < 1, .data$clock_high, 1)) %>%
    dplyr::mutate(tau_low = .data$clock_low)
  
  
  
  accepted_mutations = data.frame()
  for (segment_idx in 1:nrow(summarized_results)) {
    segment <- summarized_results[segment_idx, ]
    segment_mutations <- absolute_mutations %>%
      dplyr::filter(.data$chr == segment$chr, .data$from > segment$from, .data$to < segment$to) %>%
      tidyr::drop_na(.data$DP)
    # print(nrow(segment_mutations))
    cluster = summarized_results$tau_mean[segment_idx] # tau mean o clock mean?
    peaks = x$results_timing$data$input_data$peaks[[segment_idx]]
    # if (nrow(segment_mutations)> 40){
    segment_mutations <- segment_mutations %>% dplyr::mutate(cluster = cluster)
    
    alpha = .05
    probs <- c(alpha/2, 1 - alpha/2)
    peaks = get_clonal_peaks(unique(segment_mutations$karyotype), purity)
    accepted_idx <- lapply(1:nrow(segment_mutations), function(i) {
      
      for (p in peaks) {
        quantiles <- stats::qbinom(probs, segment_mutations$DP[i], p)
        if ((segment_mutations$NV[i] >= quantiles[1]) && (segment_mutations$NV[i] <= quantiles[2])) {
          return(i)
        }
      }
    }) %>% unlist()
    segment_mutations <- segment_mutations[accepted_idx,]
    
    accepted_mutations <- dplyr::bind_rows(accepted_mutations, segment_mutations)
    # }
  }
  as.factor(accepted_mutations$cluster)
  accepted_mutations <- accepted_mutations %>% dplyr::mutate(VAF = .data$NV/.data$DP)
  
  accepted_mutations$cluster <- as.numeric(factor(accepted_mutations$cluster))
  
  k_colors = list(
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )
  my_palette <- c("#66a61e","#7570b3","#e7298a", "#1b9e77", "#d95f02")
  
  comb = tidyr::expand_grid(unique(accepted_mutations$cluster), unique(accepted_mutations$karyotype))
  
  DENS = dplyr::tibble()
  PLOT_DATA = dplyr::tibble()
  for (i in 1:nrow(comb)) {
    cl = as.numeric(comb[i,1])
    kar = as.character(comb[i,2])
    
    df = accepted_mutations %>%
      dplyr::filter(cluster == cl, .data$karyotype == kar) %>%
      dplyr::select(.data$DP, .data$NV, .data$VAF) %>%
      dplyr::rename(successes = .data$NV, trials = .data$DP)
    
    if (nrow(df) > 0) {
      mobfit = mobster::mobster_fit(df, K = c(1,2), tail = FALSE, auto_setup = "FAST")
      
      binwidth = 0.01
      domain = seq(0, 1, binwidth)
      plot_data = mobster::Clusters(mobfit$best)
      clusters = sort(unique(plot_data$cluster), na.last = TRUE)
      Beta_peaks = mobfit$best$Clusters %>%
        dplyr::filter(.data$type == 'Mean', cluster != 'Tail')
      
      densities = suppressWarnings(mobster:::template_density(
        mobfit$best,
        x.axis = domain,
        binwidth = binwidth,
        reduce = TRUE
      ))
      
      DENS = dplyr::bind_rows(DENS, densities %>% dplyr::mutate(tickTack_cl = cl, karyotype = kar))
      PLOT_DATA = dplyr::bind_rows(PLOT_DATA, plot_data %>% dplyr::mutate(tickTack_cl = cl, karyotype = kar))
    }
  }
  
  COL = c("C1" = "#E41A1C", "C2" = "#377EB8")
  PLOT_DATA = PLOT_DATA %>%
    dplyr::group_by(.data$karyotype, .data$tickTack_cl) %>%
    dplyr::mutate(dens = .data$VAF / sum(.data$VAF))
  
  PLOT_DATA = PLOT_DATA %>% dplyr::mutate(cluster = ifelse(cluster == "C2", "P1", "P2"))
  DENS = DENS %>% dplyr::mutate(cluster = ifelse(cluster == "C2", "P1", "P2"))
  
  ggplot2::ggplot(PLOT_DATA, ggplot2::aes(.data$VAF, fill = factor(.data$tickTack_cl),              y = ggplot2::after_stat(.data$..count.. / tapply(.data$..count.., .data$PANEL, sum)[.data$PANEL]))) +
    ggplot2::geom_histogram(alpha = 0.5, color = NA, position = "identity", binwidth = binwidth) +
    ggplot2::geom_line(data = DENS, ggplot2::aes(y = .data$y, x = x, linetype=factor(cluster)), size = .5) +
    ggplot2::scale_color_manual(values = COL) +
    ggplot2::scale_fill_manual(values = my_palette) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n=3)) +
    ggplot2::labs(x = "VAF", y = "Density", fill = "", col = "Cluster", linetype="Peak") +
    ggh4x::facet_nested(.data$tickTack_cl~"CN"+karyotype) +
    # facet_grid(tickTack_cl ~ karyotype,
    #            labeller = labeller(var1 = label_both, var2 = label_both)) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.2),
      strip.background = ggplot2::element_rect(fill = "gray80", color = "gray30")
    )
}
