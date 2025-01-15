#' fit_variational Function
#'
#' @description This function performs the inference using the ADVI algorithm. Repeat inference if it fails and repeat to avoid local minima, taking the best run.
#'
#' @param input_data list: List of 7: $S: int, $N: int, $karyotype: num (0 or 1), $seg_assignment: num, $peaks:List of N of num (1:2), $NV: num, $DP: num
#' @param max_attempts num: max number of repeated inference for ADVI
#' @param initialization list: List of 4: $w: num (1:S, 1:3), $tau: num (1:K), $phi: num (1:K), $kappa: num
#' @param INIT logical: boolean variable to set the initialization phase to TRUE or FALSE
#' @param initial_iter description
#' @param grad_samples description
#' @param elbo_samples description
#' @param tolerance num: tolerance in the ELBO optimization procedure
#'
#' @keywords variational
#'
#' @return best_fit
#'
#' @export

fit_variational_h <- function(input_data, max_attempts = 2, initialization = NULL, INIT = TRUE, initial_iter = 100, grad_samples = 200, elbo_samples = 200, tolerance = 0.01) {
  # Load the Stan model
  
  get_model <- function(model_name) {
    all_paths <- list(
      "timing_betabinomial" = "mixture_CNA_timing_betabinomial.stan",
      "timing_binomial" = "mixture_CNA_timing_binomial.stan",
      "timing" = "mixture_CNA_timing.stan",
      "clustering" = "clustering.stan",
      "hierarchical" = "timing_mixed_simple.stan"
    )
    
    if (!(model_name) %in% names(all_paths)) stop("model_name not recognized")
    
    model_path <- system.file("cmdstan", all_paths[[model_name]], package = "tickTack", mustWork = T)
    tmp <- utils::capture.output(suppressMessages(model <- cmdstanr::cmdstan_model(model_path)))
    model
  }
  
  model <- get_model("hierarchical") 
  # cmdstanr::cmdstan_model("../models/timing_mixed_simple.stan")
  best_elbo <- -Inf  # Start with the worst possible ELBO
  best_fit <- NULL  # To store the best model fit
  total_attempts <- 0
  
  for (attempt in 1:max_attempts) {
    message("Attempt ", attempt, " of ", max_attempts)
    
    # retries <- 0  # Initialize retries for each attempt
    fit_successful <- FALSE  # Reset success flag for each attempt
    iter <- initial_iter  # Reset iterations for each new attempt
    
    retries <- 0
    
    tol_rel_obj <- tolerance
    
    while (!fit_successful && retries < 3) {  # Set a max number of retries for each attempt
      message("Retries ", retries, " of ", 3)
      
      # Increment total_attempts at the start of each retry
      total_attempts <- total_attempts + 1
      
     
       result <- tryCatch({
        # Attempt variational inference
        print(list(initialization))
        if (INIT == TRUE) {
          res <- model$variational(
            data = input_data,
            init = list(initialization),  # Use the provided initialization
            iter = iter,
            grad_samples = grad_samples,
            elbo_samples = elbo_samples,
            save_latent_dynamics = TRUE,
            draws = 1000,
            # output_dir = "./",
            eval_elbo = 1,
            tol_rel_obj = tolerance
          )
          # print(res$init())
          
          
          # Shuffle and moderately perturb initialization values
          ###############################################################################################################################################################
          shuffled_taus <- sample(initialization$tau)
          # print(paste0("init_tau: ", initialization$tau ))
          perturbed_taus <- shuffled_taus + stats::rnorm(length(shuffled_taus), 0, 0.1) # Adjust the perturbation scale as needed
          perturbed_taus[perturbed_taus >= 0.88] <- 0.88  # Check for elements greater than 1 and replace them with 1 otherwise fit fails
          perturbed_taus[perturbed_taus <= 0] <- 0.07
          
          initialization$tau = perturbed_taus
          # phi does not change
          
          # Apply the perturbation to w as wells
          # print (paste0("w_init after inference: ",initialization$w," ",ncol(initialization$w) == 1))
          #######################################################################################################################################################################
          
          
        } else {
          res <- model$variational(                                                #getter?
            data = input_data,
            iter = iter,
            grad_samples = grad_samples,
            elbo_samples = elbo_samples,
            save_latent_dynamics = TRUE,
            draws = 1000,
            # output_dir = "./",
            eval_elbo = 1,
            tol_rel_obj = tolerance
          )
        }
        
        
        
        output_files <- res$latent_dynamics_files()                                 #getter?
        elbo_data <- utils::read.csv(output_files, header = FALSE, comment.char = "#")
        colnames(elbo_data) <- c("iter", "time_in_seconds", "ELBO")
        elbo_values <- elbo_data$ELBO  # The ELBO column
        elbo <- elbo_values[length(elbo_data)]
        # Update the best model if this ELBO is better
        if (!is.na(elbo) && elbo > best_elbo) {
          best_elbo <- elbo
          best_fit <- res
        }
        
        message("ELBO for this run: ", elbo)
        
        # Check for invalid log evaluations
        output_files <- res$output_files()                                           #getter? + implement condition
        fit_successful <- TRUE  # Mark this fit as successful
        
      }, error = function(e) {
        message("An error occurred during inference: ", e$message)
        retries <- retries + 1  # Increment retry count
        fit_successful <- FALSE  # Mark fit as unsuccessful
        iter <- iter + 400
        tol_rel_obj <- tol_rel_obj * 5
        grad_samples <- grad_samples * 3
        elbo_samples <- elbo_samples * 3
        
        # if it failed then       retries <- retries + 1 and fit_succ <- false else continue 
        
        NULL  # Ensure NULL is returned so loop can continue
      })
      
      # No return from `tryCatch` - just continue the loop if not successful
    }
    
    if (retries == 3) {
      message("Max retries reached for attempt ", attempt)
    }
  }
  
  if (!is.null(best_fit)) {
    # message("Best ELBO after ", total_attempts, " attempts: ", best_elbo)
    return(best_fit)  # Return the result with the best ELBO
  } else {
    message("Inference could not be completed successfully.")
    return(NULL)  # Return NULL if no fit was successful
  }
}
