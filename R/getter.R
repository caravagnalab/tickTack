
get_model <- function(model_name) {
  all_paths <- list(
    "timing_betabinomial" = "mixture_CNA_timing_betabinomial.stan",
    "timing_binomial" = "mixture_CNA_timing_binomial.stan",
    "timing" = "mixture_CNA_timing.stan",
    "clustering" = "clustering.stan"
  )

  if (!(model_name) %in% names(all_paths)) stop("model_name not recognized")

  model_path <- system.file("cmdstan", all_paths[[model_name]], package = "tickTack", mustWork = T)
  tmp <- utils::capture.output(suppressMessages(model <- cmdstanr::cmdstan_model(model_path)))
  model
}

get_clonal_peaks = function(k, purity) {
  multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
  major <- multiplicities[1]
  n_tot <- sum(multiplicities)

  # get only Major and 1
  multiplicities <- c(1, major)

  peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
  return(sort(peaks))
}

get_tau_posteriors = function(fit, k) {
  omega1 <- fit$draws("omega[1]", format = 'matrix') %>% dplyr::as_tibble() %>% `colnames<-`('value') %>% dplyr::mutate(value = as.numeric(.data$value)) %>% dplyr::pull(.data$value)
  omega2 <- fit$draws("omega[2]", format = 'matrix') %>% dplyr::as_tibble() %>% `colnames<-`('value') %>% dplyr::mutate(value = as.numeric(.data$value)) %>% dplyr::pull(.data$value)
  #omega1 <- rstan::extract(fit, pars="omega[1]") %>% unlist()
  #omega2 <- rstan::extract(fit, pars="omega[2]") %>% unlist()

  if (k == '2:1') {
    tau_posterior <- 3 * omega2 / (2*omega2 + omega1)
  } else {
    tau_posterior <- 2 * omega2 / (2*omega2 + omega1)
  }

  tau_posteriors <- dplyr::tibble(tau = tau_posterior)
  tau_posteriors
}

get_reference = function(ref) {
  if(ref %in% c("hg19", "GRCh37"))
    return(tickTack::chr_coordinates_hg19)

  if(ref %in% c("hg38", "GRCh38"))
    return(tickTack::chr_coordinates_GRCh38)

  stop("Available references: hg19 (or GRCh37), and hg38 (or GRCh38)")
}
