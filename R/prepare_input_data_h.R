#' prepare_input_data Function
#'
#' This function obtains the list of input data to be used in the stan model.
#' @param mutations  list: output from CNAqc Mutations(x), where x is a CNAqc object
#' @param cna list: output from CNAqc CNA(x), where x is a CNAqc object
#' @param purity num: (type = double) sample purity between 0 and 1
#' 
#' @param possible_k chr: "2:1" "2:2" "2:0"
#' @param alpha num: (type double) confidence interval level to choose the data that fall in the expected binomial intervals
#' @param min_mutations_number num: (type double) minimum number of accepted mutations for a segment to be included in the inference
#' 
#' @return accepted_data: list $input_data:list: List of 7: $S: int, $N: int, $karyotype: num (0 or 1), $seg_assignment: num, $peaks:List of N of num [1:2], $NV: num, $DP: num
#'                             $accepted_cna: tibble [S × 5] (S3: tbl_df/tbl/data.frame): $segment_original_indx: int, $segment_name: chr, $segment_id: num, $karyotype: chr, $chr: chr 
#'                              
#' @keywords input
#' @export
#' @examples
#' prepare_input_data()

prepare_input_data = function(mutations, segments, purity, possible_k = c("2:1", "2:2", "2:0"), alpha = .05, min_mutations_number = 2 ){
  
  #  already done in the preprocessing : data <- data[sort(data$segment_id), ]
  
  segments <- segments %>%
    tidyr::drop_na(Major, minor)
  
  n_segments <- nrow(segments)
  accepted_segments_info <- dplyr::tibble()
  # summarized_results <- dplyr::tibble()
  
  accapted_mutations_all <- dplyr::tibble()
  

  
  accepted_segment_idx <- 0
  
  for (segment_idx in 1:n_segments) {
    
    print(segment_idx)
    
    # Segments
    segment <- segments[segment_idx, ]
    chr <- segment$chr                                               #getter ?$ 
    
    segment_id <- paste(chr, segment$from, segment$to, sep = "_")
    
    # Get karyotype
    Major <- segment$Major
    minor <- segment$minor
    
    k <- paste(Major, minor, sep=':')
    
    peaks <- get_clonal_peaks(k, purity)
    
    if (k %in% possible_k) {
      # Get info for mutations
      segment_mutations <- mutations %>%
        filter(chr == segment$chr,from > segment$from, to < segment$to) %>%
        drop_na(DP)
      
      accepted_mutations <- data.frame()
      if (nrow(segment_mutations) > 0) {
        # Check if mutation is inside CI
        probs <- c(alpha/2, 1 - alpha/2)
        
        DP <- segment_mutations$DP
        NV <- segment_mutations$NV
        
        accepted_idx <- lapply(1:length(DP), function(i) {
          #fisso i picchi per il segmento e vedo se tutte le mutazioni ricadono in almeno uno dei due intervalli intorno ai picchiche sono diversi a seconda del valore di DP per la specifica mutazione
          for (p in peaks) {
            quantiles <- qbinom(probs, DP[i], p)
            if ((NV[i] >= quantiles[1]) && (NV[i] <= quantiles[2])) {
              return(i)
            }
          }
        }) %>% unlist()
        
        # Get only good mutations
        accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx])
      
      }
      
      if (nrow(accepted_mutations) >= min_mutations_number) {
        
        accepted_segment_idx <- accepted_segment_idx + 1
        
        cli::cli_alert_info("Adding segment with index {.val {segment_idx}} to segments included in the inference.")
        
        # return the fit result directly together with the information about the segments associated to the clocks
        accepted_segments_info <- dplyr::bind_rows(accepted_segments_info, dplyr::tibble(segment_original_indx = segment_idx , segment_name = segment_id, segment_id = accepted_segment_idx, karyotype = k, chr = chr))
        
        
        accepted_mutations <- dplyr::tibble(DP = DP[accepted_idx], NV = NV[accepted_idx], segment_original_indx = segment_idx , segment_name = segment_id, segment_id = accepted_segment_idx, karyotype = k, chr = chr)
        accapted_mutations_all <- dplyr::bind_rows(accapted_mutations_all, accepted_mutations)
        # accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx], segment_id=data$segment_name[accepted_idx], karyotype=data$karyotype[accepted_idx], tau=data$tau[accepted_idx] , segment_index=data$segment_id[accepted_idx])
        
      }
    }
  }
  
        # check that there is at least one segment 
        if(!nrow(accepted_segments_info)*ncol(accepted_segments_info)){
          stop("No segments respected the constraint to perform the clock inference in this CNAqc object.")
        }
        
        # then obtain the data "all together" and get the input for the variational model inference 
        
        accepted_cna <- accepted_segments_info[sort(accepted_segments_info$segment_id), ]
        
        
        
        input_data <- list(
          S = nrow(accepted_cna),
          N = nrow(accapted_mutations_all),
          karyotype = lapply(accepted_cna$karyotype, karyo_to_int) %>% unlist(),
          seg_assignment = accapted_mutations_all$segment_id, 
          peaks = lapply(accepted_cna$karyotype, get_clonal_peaks, purity=purity ), # get_clonal_peaks(accepted_cna$karyotype, purity),
          NV = accapted_mutations_all$NV,
          DP = accapted_mutations_all$DP
        )

    accepted_data = list(input_data = input_data, accepted_cna = accepted_cna)
  return(accepted_data)
}
