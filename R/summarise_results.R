#' Summarise bias and dispersion (MADm) across simulated fits
#'
#' @param results Output from \code{run_model()}.
#' @examples
#' \donttest{
#' if (requireNamespace("rstan", quietly = TRUE)){
#'   sim <- simulate_data(ind = 50, Valpha = 0.2, Vepsilon = 0.1, iterations = 2)
#'   res <- run_model(sim, model = "Trait.stan", iter = 100, cores = 2)
#'   summarise_results(res)
#' } else {
#'   message("rstan not available; example skipped.")
#' }
#' }
#' @return Data frame with only parameters that were estimated in the Stan model.
#' @export
summarise_results <- function(results) {
  stopifnot(inherits(results, "socialSim_results"))
  if (length(results) == 0L) stop("No results provided.")

  # Map simulation parameter names -> Stan variable names
  map <- c(
    B_0              = "B_0",
    psi              = "psi",
    Valpha           = "Sigma2_intercept",
    Vepsilon         = "Sigma2_epsilon",
    Vpsi             = "Sigma2_psi",
    Vx               = "Sigma2_x",
    Vphi             = "Sigma2_phi",
    Ve               = "Sigma2_e",
    Vxe              = "Sigma2_ex",
    cov_alpha_psi    = "cov_int_psi",
    cov_alpha_phi    = "cov_int_phi",
    cov_psi_phi      = "cov_psi_phi",
    r_alpha_psi      = "cor_1",
    r_alpha_epsilon  = "cor_2",
    r_epsilon_psi    = "cor_3",
    r_alpha_x        = "cor_4",
    r_psi_x          = "cor_5",
    r_epsilon_x      = "cor_6"
  )

  true_all <- results[[1]]$params
  keep <- intersect(names(map), names(true_all))
  map <- map[keep]
  par_names <- names(map)

  # Extract posterior means for each dataset
  est_mat <- matrix(NA_real_, nrow = length(results), ncol = length(par_names),
                    dimnames = list(NULL, par_names))
  for (i in seq_along(results)) {
    summ <- results[[i]]$summary
    for (p in seq_along(par_names)) {
      stan_var <- unname(map[par_names[p]])
      row <- summ[rownames(summ) == stan_var, , drop = FALSE]
      est_mat[i, p] <- if (nrow(row)) row[1, "mean"] else NA_real_
    }
  }

  mean_est <- colMeans(est_mat, na.rm = TRUE)
  true_vec <- unlist(true_all[par_names], use.names = FALSE)

  abs_bias <- mean_est - true_vec
  rel_bias <- ifelse(abs(true_vec) > .Machine$double.eps, abs_bias / true_vec, NA_real_)
  madm <- apply(est_mat, 2, function(x) mean(abs(x - mean(x, na.rm = TRUE)), na.rm = TRUE))

  out <- data.frame(
    Parameter = par_names,
    True      = as.numeric(true_vec),
    Mean_est  = as.numeric(mean_est),
    RelativeBias_percentage   = as.numeric(rel_bias)*100,
    RelativeDispersion_percentage  = as.numeric(madm)*100,
    check.names = FALSE
  )

  # Filter out parameters not estimated (NaN mean)
  out <- out[!is.nan(out$Mean_est) & !is.na(out$Mean_est), , drop = FALSE]

  # Reset row names
  rownames(out) <- NULL
  out
}

