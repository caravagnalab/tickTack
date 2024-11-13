#' get_clonal_peaks Function
#'
#' Obtain the theoretical peaks that would be observed in a VAF spectrum.
#' @param k karyotype
#' @param purity peaks
#' @keywords peaks
#' @export

get_clonal_peaks = function(k, purity) {
  multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
  major <- multiplicities[1]
  n_tot <- sum(multiplicities)
  # get only Major and 1
  multiplicities <- c(1, major)
  peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
  return(sort(peaks))
}


#' karyo_to_int Function
#'
#' This function allows you to change the string representation of the CN events into numeric values 0 - 1.
#' @param k karyotype string
#' @keywords karyotype
#' @export

karyo_to_int <- function(k) {
  if (k == "2:1") return(1)
  return(0)
}



#' Extract mutations.
#'
#' @description Getter to obtain mutation calls from an object.
#'
#' @param x A CNAqc object.
#' @param cna \code{"clonal"} for clonal CNAs, \code{"subclonal"} for subclonal CNAs.
#' @param type \code{"SNV"} for single-nucleotide variants, \code{"indel"} for insertion-deletions.
#'
#' @return A tibble with the data.
#' @export
#'

Mutations = function(x, cna = c("clonal", "subclonal"), type = c("SNV", "indel"))
{
  # stopifnot(inherits(x, 'cnaqc'))
  
  clonal = NULL
  # if("clonal" %in% cna)
  clonal = x$mutations %>% dplyr::mutate(cna = 'clonal')
  
  subclonal = NULL
  # if(("subclonal" %in% cna) & x$has_subclonal_CNA)
  # subclonal = x$cna_subclonal$mutations %>% Reduce(f = dplyr::bind_rows) %>% dplyr::mutate(cna = 'subclonal')
  
  mutations = dplyr::bind_rows(clonal, subclonal) %>%
    dplyr::select(.data$chr, .data$from, .data$to, .data$ref, .data$alt, .data$NV, .data$DP, .data$VAF, dplyr::everything())
  # %>%
  #   dplyr::filter(type %in% !!type)
  
  if((mutations %>% nrow())== 0) cli::cli_alert_danger("No mutations with these parameters: CNA {.field {cna}}, type {.field {type}}.")
  
  return(mutations)
}

#' Extract CNAs.
#'
#' @description Getter to obtain copy number calls from an object.
#'
#' @param x A CNAqc object.
#' @param type \code{"clonal"} for clonal CNAs, \code{"subclonal"} for subclonal CNAs.
#'
#' @return A tibble with the data.
#' @export
#'

CNA = function(x, type = c("clonal", "subclonal"))
{
  # stopifnot(inherits(x, 'cnaqc'))
  
  clonal = NULL
  # if("clonal" %in% type)
  clonal = x$cna
  
  subclonal = NULL
  # if(("subclonal" %in% type) & x$has_subclonal_CNA)
  #   subclonal = x$cna_subclonal
  
  cna = dplyr::bind_rows(clonal, subclonal) %>%
    dplyr::select(.data$chr, .data$from, .data$to, dplyr::starts_with('Major'), dplyr::starts_with('minor'), .data$CCF, dplyr::everything())
  
  if((cna %>% nrow())== 0) cli::cli_alert_danger("No CNAs with these parameters: {.field {cna}}.")
  
  return(cna)
}



