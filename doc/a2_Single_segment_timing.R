## ----echo=TRUE----------------------------------------------------------------
library(tickTack)

# View example dataset components
mutations <- tickTack::pcawg_example$mutations
cna <- tickTack::pcawg_example$cna
metadata <- tickTack::pcawg_example$metadata

head(mutations)
head(cna)
metadata

## ----echo=FALSE---------------------------------------------------------------
# Extract input data
segments <- tickTack::pcawg_example$cna
mutations <- tickTack::pcawg_example$mutations
purity <- tickTack::pcawg_example$metadata$purity

# Run the fit function
results <- fit(
  segments = segments,
  mutations = mutations,
  purity = purity,
  possible_k = c("2:1", "2:2", "2:0"),
  beta_binomial = TRUE
)

## ----echo=TRUE----------------------------------------------------------------
# View summarized results
results$summarized_results

# View detailed inference results
results$inference_results

## -----------------------------------------------------------------------------
tickTack::plot_timing(results, segments, colour_by = "karyotype")

