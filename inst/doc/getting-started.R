## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(socialSim)

## ----simulate, eval=FALSE-----------------------------------------------------
# sim <- simulate_data(
#   ind = 200,
#   partners = 4,
#   repeats = 1,
#   iterations = 5,
#   B_0 = 1,
#   psi = 0.3,
#   Valpha = 0.2,
#   Vepsilon = 0.1
# )

## ----model, eval=FALSE--------------------------------------------------------
# res <- run_model(sim, model = "Trait.stan", iter = 2000, cores = 4)

## ----summarise, eval=FALSE----------------------------------------------------
# summary <- summarise_results(res)
# print(summary)

## ----example, eval=FALSE------------------------------------------------------
# sim <- simulate_data(ind = 50, partners = 2, iterations = 4, Valpha = 0.2, Vepsilon = 0.1)
# res <- run_model(sim, model = "Trait.stan", iter = 500, cores = 4)
# summary <- summarise_results(res)
# print(summary)

## -----------------------------------------------------------------------------
?simulate_data
?run_model
?summarise_results

