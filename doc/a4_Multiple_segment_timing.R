## ----echo=TRUE----------------------------------------------------------------
library(tickTack)

# View example dataset components
mutations <- tickTack::pcawg_example_2$mutations
cna <- tickTack::pcawg_example_2$cna
metadata <- tickTack::pcawg_example_2$metadata

head(mutations)
head(cna)
metadata

## ----echo=FALSE---------------------------------------------------------------
# Extract input data
data <- tickTack::pcawg_example_2
tolerance = 0.1

# Run the fit function
data <- fit_h(
  x = data,
  max_attempts = 2,
  INIT = TRUE,
  tolerance = tolerance
)


## ----echo=TRUE----------------------------------------------------------------
# View summary for a specific K, here K = 2
results <- data$results

## ----echo=TRUE----------------------------------------------------------------
# View summary for a specific K, here K = 2
results$draws_and_summary[[2]]$summary

# View detailed summarized results for a specific K, here K = 2
results$draws_and_summary[[2]]$summarized_results

## -----------------------------------------------------------------------------
results_model_selection <- tickTack::model_selection_h(results, n_components = 0)

best_K <- results_model_selection$best_K
model_selection_tibble <- results_model_selection$model_selection_tibble
entropy <- results_model_selection$entropy_list


## -----------------------------------------------------------------------------
tickTack::plot_timing_h(results, 2)

## -----------------------------------------------------------------------------
posterior_clocks <- tickTack::plot_posterior_clocks_h(results, 2)
posterior_weights <- tickTack::plot_posterior_weights_h(results, 2)


## -----------------------------------------------------------------------------
K = nrow(results_model_selection$model_selection_tibble)

p_elbo <- list()
for (i in 1:K){
  p_elbo[[i]] <- tickTack::plot_elbo_h(results$elbo_iterations[[i]]) + ggplot2::ggtitle(paste0("K = ", i))
}
p_elbo <- gridExtra::grid.arrange(grobs = p_elbo, nrow=K)  #add global title
p_elbo


## -----------------------------------------------------------------------------

plot_model_selection_inference <- list()
for (i in 1:K){
  plot_model_selection_inference[[i]] <- tickTack::plot_timing_h(results, i) + ggplot2::ggtitle(paste0("K = ", i))
}
plot_model_selection_inference <- gridExtra::grid.arrange(grobs = plot_model_selection_inference, nrow = K) #add global title
plot_model_selection_inference


