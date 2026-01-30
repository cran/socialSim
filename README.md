
<!-- README.md is generated from README.Rmd. Please edit that file -->

# socialSim <a href="https://github.com/RoriWijnhorst/socialSim"><img src="https://img.shields.io/badge/GitHub-RoriWijnhorst/socialSim-blue?logo=github" alt="GitHub"></a>

**Simulate and Analyse Social Interaction Models**

<!-- badges: start -->

[![R-CMD-check](https://github.com/RoriWijnhorst/socialSim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RoriWijnhorst/socialSim/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The **socialSim R package** provides tools to simulate and analyse
datasets of social interactions between individuals using hierarchical
Bayesian models implemented in Stan. This packages accompanies
[Wijnhorst et al. (2025)](https://doi.org/10.32942/X2F65M) *EcoEvoRxiv*,
which details the underlying statistical models.

It enables users to generate realistic social interaction data, where
individual phenotypes influence and respond to those of their partners.
You can simulate a sampling design by adjusting the **number of
individuals, partners, and repeated dyads**. The simulation framework
allows control over variation in **mean trait values, social
responsiveness, and social impact** and correlation, making it suitable
for research on **direct and indirect genetic effects (DGEs and IGEs)**
and **interacting phenotypes**. See `?simulate_data` for a full list of
adjustable parameters.

The package also provides analysis functions to evaluate model
performance in terms of bias and dispersion, using both established and
novel approaches to modelling social effects, including
impact–responsiveness, variance–partitioning, and trait-based models.

------------------------------------------------------------------------

## 🧭 Installation

You can install the development version from GitHub using:

``` r
# install.packages("remotes")
remotes::install_github("RoriWijnhorst/socialSim")
# Then load the package:
library(socialSim)
```

## ⚙️ Example workflow:

``` r
library(socialSim)

# 1. Simulate data. See ?simulate_data for all adjustable parameters
sim <- simulate_data(
  ind = 400,            # number of unique focal individuals
  partners = 4,         # number of social partners per individual
  repeats = 1,          # number of repeats of dyads   
  iterations = 10,      # number of datasets created    
  B_0 = 1,              # population intercept
  psi = 0.3,            # population-level response
  Valpha = 0.2,         # variance in direct effects
  Vepsilon = 0.1        # variance in residual partner effects
)

# 2. Fit a Stan model. For the analyses, rstan needs to be installed.
res <- run_model(sim, model = "Trait.stan", iter=2000, cores = 6)

# 3. Summarise results
summary <- summarise_results(res)
print(summary)
```

## 🧪 Available IGE models

| Model name      | Description                                                  |
|-----------------|--------------------------------------------------------------|
| I&R.stan        | Full impact–responsiveness model                             |
| VP.stan         | Variance-partitioning model                                  |
| Trait.stan      | Trait-based model with residual partner effects              |
| Trait_only.stan | Simple trait-based model without residual partner effects    |
| Trait_RS.stan   | Random-slope trait model with residual partner effects       |
| Trait_EIV.stan  | Errors-in-variable trait model with residual partner effects |

See [Wijnhorst et al. (2025)](https://doi.org/10.32942/X2F65M)
*EcoEvoRxiv* for detailed explanation of the models.
