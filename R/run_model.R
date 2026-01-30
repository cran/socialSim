#' Fit one of the available Stan models to simulated datasets
#'
#' @param sim Output from \code{simulate_data()}.
#' @param model Name of the Stan model to use (choose from available options).
#' @param iter Number of iterations per chain (default = 1000).
#' @param seed Random seed for reproducibility.
#' @param cores Number of CPU cores used for parallel execution.
#' @examples
#' \donttest{
#' if (requireNamespace("rstan", quietly = TRUE)) {
#'   sim <- simulate_data(ind = 50, Valpha = 0.2, Vepsilon = 0.1, iterations = 2)
#'   res <- run_model(sim, model = "Trait.stan", iter = 100, cores = 2)
#'   summarise_results(res)
#' } else {
#'   message("rstan not available; example skipped.")
#' }
#' }
#' @return A list of fitted model summaries, one per dataset.
#' @export
run_model <- function(sim,
                      model = NULL,
                      iter = 2000,
                      seed = 1234,
                      cores = 1) {

  stopifnot(inherits(sim, "socialSim_data"))

  if (missing(cores) || !is.numeric(cores) || cores < 1)
    stop("Please specify a valid number of cores (e.g., cores = 4).")

  # -------------------------------------------------------------------------
  # 1. Define available models
  # -------------------------------------------------------------------------
  available_models <- c(
    "I&R.stan",
    "VP.stan",
    "Trait.stan",
    "Trait_only.stan",
    "Trait_RS.stan",
    "Trait_EIV.stan"
  )

  if (is.null(model)) {
    message("Please choose one of the following available models:\n",
            paste0("  - ", available_models, collapse = "\n"))
    stop("No model selected. Specify a model, e.g. model = 'Trait.stan'.")
  }

  if (!model %in% available_models)
    stop("Invalid model name. Available options are:\n",
         paste0("  - ", available_models, collapse = "\n"))

  model_path <- system.file("stan", model, package = "socialSim")
  if (model_path == "")
    stop("Model not found in inst/stan/. Please check model name.")

  # -------------------------------------------------------------------------
  # 2. Sanity checks for data/model compatibility
  # -------------------------------------------------------------------------
  params <- sim$params

  if (model %in% c("I&R.stan", "Trait_EIV.stan") && params$Vxe == 0) {
    message(
      "\nThe selected model (", model, ") estimates measurement error (Vxe),",
      "\nBut Vxe = 0 in your simulated data."
    )
    ans <- utils::menu(c("Yes", "No"), title = "Continue anyway?")
    if (ans == 2) stop("Model fitting cancelled by user.")
  }

  if (model %in% c("I&R.stan", "Trait_RS.stan") && params$Vpsi == 0) {
    message(
      "\nThe selected model (", model, ") estimates individual variation in responsiveness (Vpsi),",
      "\nBut Vpsi = 0 in your simulated data."
    )
    ans <- utils::menu(c("Yes", "No"), title = "Continue anyway?")
    if (ans == 2) stop("Model fitting cancelled by user.")
  }

  # -------------------------------------------------------------------------
  # 3. Detect backend: rstan only
  # -------------------------------------------------------------------------
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop(
      "The 'rstan' package is required but not installed.\n",
      "Please install it with:\n",
      "  install.packages('rstan', repos = 'https://cloud.r-project.org')"
    )
  }

  # -------------------------------------------------------------------------
  # 4. Compile model once (outside workers)
  # -------------------------------------------------------------------------
  message("Compiling Stan model with rstan...")
  rstan::rstan_options(auto_write = TRUE)
  mod <- rstan::stan_model(model_path)

  fit_one <- function(dat, id) {
    standata <- list(
      n_obs = dat$n_obs,
      n_ind = dat$n_ind,
      individual = dat$individual,
      opponent = dat$opponent,
      xj = dat$xj,
      z = dat$z
    )
    fit <- rstan::sampling(
      object = mod,
      data = standata,
      iter = iter,
      chains = 1,
      seed = seed + id,
      refresh = 0
    )
    list(
      summary = rstan::summary(fit)$summary,
      params  = params,
      data_id = id
    )
  }

  # -------------------------------------------------------------------------
  # 5. Parallel execution across datasets
  # -------------------------------------------------------------------------
  if (!requireNamespace("future.apply", quietly = TRUE))
    stop("Please install 'future.apply' for parallel execution.")

  message("Running ", length(sim$data), " datasets on ", cores, " cores")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  cores <- min(cores, future::availableCores())
  future::plan(future::multisession, workers = cores)

  timing <- system.time({
    results <- future.apply::future_lapply(
      seq_along(sim$data),
      function(i) fit_one(sim$data[[i]], i),
      future.seed = TRUE,
      future.globals = list(mod = mod)  # share compiled model
    )
  })

  runtime_sec <- round(unname(timing["elapsed"]), 1)

  message("All datasets finished. Runtime: ", runtime_sec, " s")

  structure(results, class = "socialSim_results")
  }


