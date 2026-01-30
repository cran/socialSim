#' Simulate social interaction datasets
#'
#' @description
#' This function generates datasets where individual phenotypes are influenced
#' by both direct and indirect (social) effects, under a specified sampling design.
#'
#' @param ind Number of individuals.
#' @param partners Partners per individual.
#' @param repeats Repeats per unique dyad.
#' @param iterations Number of datasets to simulate.
#' @param B_0 Population intercept.
#' @param psi Population-level responsiveness (social slope).
#' @param Valpha Direct effect (focal variance).
#' @param Vepsilon Indirect effect (partner variance).
#' @param Vpsi Social responsiveness (among individual variance in slopes).
#' @param Vx Partner trait variance.
#' @param Ve Residual variance.
#' @param Vxe Measurement error/within-individual variation in partner trait.
#' @param r_alpha_epsilon Corr(alpha, epsilon).
#' @param r_alpha_psi Corr(alpha, psi).
#' @param r_epsilon_psi Corr(epsilon, psi).
#' @param r_alpha_x Corr(alpha, x).
#' @param r_psi_x Corr(psi, x).
#' @param r_epsilon_x Corr(epsilon, x).
#' @param fix_total_var Logical; if TRUE (default), residual variance is
#'   adjusted so total phenotypic variance is approx. 1.
#'
#' @return A list with:
#' \itemize{
#'   \item data: list of datasets
#'   \item params: named list of effect sizes
#'   \item design: sample design (n_ind, partners, repeats, iterations)
#' }
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm aggregate
#' @examples
#' sim <- simulate_data(ind =50, partners = 4, iterations = 5,
#'                      B_0 = 1, Valpha=0.2, Vepsilon = 0.1)
#'
#' @export
simulate_data <- function(
    ind = 200,
    partners = 4,
    repeats = 1,
    iterations = 100,
    B_0 = 0,                 # new population intercept
    psi = NULL,              # required if no Vepsilon
    Valpha,                  # required
    Vepsilon = NULL,         # optional
    Vpsi = 0,
    Vx = 1,
    Ve = 0.6,
    Vxe = 0,
    r_alpha_epsilon = 0,
    r_alpha_psi = 0,
    r_epsilon_psi = 0,
    r_alpha_x = 0,
    r_psi_x = 0,
    r_epsilon_x = 0,
    fix_total_var = TRUE
) {

  # --- 1. Required checks ---
  if (missing(Valpha) || !is.numeric(Valpha) || Valpha <= 0)
    stop("Specify a focal direct effect 'Valpha'. Must be a positive numeric value")

  # Must specify at least one partner-effect parameter
  if (is.null(Vepsilon) && is.null(psi)) {
    stop(
      "You must specify at least one partner-effect parameter:\n",
      "  - 'Vepsilon' (variance in partner effects), or\n",
      "  - 'psi' (population-level responsiveness).\n",
      "Both cannot be missing."
    )
  }

  # If missing, set to 0 (safe default)
  if (is.null(Vepsilon)) Vepsilon <- 0
  if (is.null(psi)) psi <- 0

  # --- 2. Fix total phenotypic variance (optional) ---
  if (fix_total_var) {
    var_components <- Valpha + Vpsi * Vx + (psi^2) * Vx + Vepsilon + 2 * psi * r_epsilon_x
    Ve_new <- 1 - var_components
    if (Ve_new < 0) {
      warning(sprintf(paste0(
        "Computed residual variance would be negative (%.3f).\n",
        "Replaced with 1e-8. Consider lowering variance components or\n",
        "set fix_total_var = FALSE and set residual variance by 'Ve = '."
      ), Ve_new))
      Ve_new <- 1e-8
    }
    Ve <- Ve_new
  }

  # --- 3. Sampling design ---
  sampling_design_balanced <- function(ind, partners, repeats) {
    ind_seq <- seq_len(ind)
    IDi <- rep(ind_seq, partners)
    IDj <- vector("list", partners)
    for (i in seq_len(partners)) {
      seq_j <- c(ind_seq[-seq_len(i)], ind_seq[seq_len(i)])
      IDj[[i]] <- seq_j
    }
    IDj <- do.call(c, IDj)
    df <- data.frame(IDi = rep(IDi, each = 1), IDj = IDj)
    df[rep(seq_len(nrow(df)), repeats), , drop = FALSE]
  }

  # --- 4. Means and covariance structure ---
  M <- c(0, 0, psi, 0) # Malpha, Mepsilon, psi, Mx (all means fixed to 0 except psi)
  G <- matrix(0, 4, 4)
  diag(G) <- c(Valpha, Vepsilon, Vpsi, Vx)
  G[1,2] <- G[2,1] <- r_alpha_epsilon * sqrt(Valpha * Vepsilon)
  G[1,3] <- G[3,1] <- r_alpha_psi      * sqrt(Valpha * Vpsi)
  G[2,3] <- G[3,2] <- r_epsilon_psi    * sqrt(Vepsilon * Vpsi)
  G[4,1] <- G[1,4] <- r_alpha_x        * sqrt(Valpha * Vx)
  G[4,2] <- G[2,4] <- r_epsilon_x      * sqrt(Vepsilon * Vx)
  G[4,3] <- G[3,4] <- r_psi_x          * sqrt(Vpsi * Vx)

  # --- 5. Simulation for one dataset ---
  sim_one <- function(df) {
    a <- as.data.frame(MASS::mvrnorm(ind, mu = M, Sigma = G))
    names(a) <- c("alpha", "epsilon", "psi", "x")
    a$ID <- seq_len(ind)

    df$alpha_i   <- a$alpha[match(df$IDi, a$ID)]
    df$epsilon_j <- a$epsilon[match(df$IDj, a$ID)]
    df$psi_i     <- a$psi[match(df$IDi, a$ID)]
    df$xj        <- a$x[match(df$IDj, a$ID)]

    df$x_ijk <- df$xj + rnorm(nrow(df), 0, sqrt(Vxe))
    df$e_ijk <- rnorm(nrow(df), 0, sqrt(Ve))
    df$z_i   <- B_0 + df$alpha_i + df$psi_i * df$xj + df$epsilon_j + df$e_ijk

    xi <- aggregate(df$x_ijk, list(IDi = df$IDj), mean)
    list(
      n_obs = nrow(df),
      n_ind = ind,
      individual = df$IDi,
      opponent = df$IDj,
      xj = df$x_ijk,
      xi = xi$x,
      z  = df$z_i
    )
  }

  df_design <- sampling_design_balanced(ind, partners, repeats)
  data_list <- lapply(seq_len(iterations), function(i) sim_one(df_design))

  # --- 6. Derived parameters ---
  Vphi <- (psi^2) * Vx + Vepsilon + 2 * psi * r_epsilon_x * sqrt(Vx * Vepsilon)
  cov_alpha_psi <- r_alpha_psi * sqrt(Valpha * Vpsi)
  cov_alpha_phi <- r_alpha_epsilon * sqrt(Valpha * Vepsilon) +
    psi * r_alpha_x * sqrt(Valpha * Vx)
  cov_psi_phi   <- r_epsilon_psi * sqrt(Vepsilon * Vpsi) +
    psi * r_psi_x * sqrt(Vpsi * Vx)

  params <- list(
    B_0 = B_0, psi = psi,
    Valpha = Valpha, Vepsilon = Vepsilon, Vpsi = Vpsi, Vx = Vx,
    Ve = Ve, Vxe = Vxe,
    r_alpha_epsilon = r_alpha_epsilon, r_alpha_psi = r_alpha_psi,
    r_epsilon_psi = r_epsilon_psi, r_alpha_x = r_alpha_x,
    r_psi_x = r_psi_x, r_epsilon_x = r_epsilon_x,
    Vphi = Vphi,
    cov_alpha_psi = cov_alpha_psi,
    cov_alpha_phi = cov_alpha_phi,
    cov_psi_phi = cov_psi_phi
  )

  design <- list(
    ind = ind,
    partners = partners,
    repeats = repeats,
    iterations = iterations
  )

  structure(list(data = data_list, params = params, design = design),
            class = "socialSim_data")
}
