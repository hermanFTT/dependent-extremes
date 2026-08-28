# =============================================================================
# bias_core.R
#
# Monte-Carlo engine for the bias / RMSE / coverage sensitivity study
# described in bias_sensitivity.md.
#
# The unit of work is ONE GRID CELL:
#
#     (dep_type, dep_par, k, n, r, m)
#
# `m` (the return period) only affects the z_m target -- k and sigma estimates
# are identical across cells that differ only in m, since m enters the fit
# solely through the return-level computation, not the GPD parameter fit.
#
# Within a cell, S dependent series are simulated and ALL requested approaches
# are fitted to the SAME sample, so the comparison is paired.  For each of the
# three targets
#
#     k     the shape
#     sigma the scale ABOVE the threshold,  sigma* = sigma + k*u
#     zm    the m-observation return level,  z_m = F^-1(1 - 1/m)
#
# we report mean estimate, bias, RMSE and (optionally) coverage of a credible
# interval at a nominal level.
#
# Returns a tidy long data.frame: one row per approach x target parameter.
#
# This file defines functions only -- it runs nothing on source().
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

APPROACHES <- c("Conventional", "Weighted GPD", "Cluster maxima (POT)")

# dep_type -> (generator, name of the dependence parameter)
DEP_TYPES <- c("ar1", "log", "nlog")


# -----------------------------------------------------------------------------
# source_gpd_functions()
#
# weighted_gpdfit.R ends with a live scratch block (`sim <- sim_compare(...)`,
# `x11()`, ...), so a plain source() would run a simulation and open a graphics
# device in every pipeline job.  Parse it and keep only function definitions.
# -----------------------------------------------------------------------------
source_gpd_functions <- function(path = "weighted_gpdfit.R",
                                 envir = parent.frame()) {

  if (!file.exists(path))
    stop("cannot find '", path, "' (run from the project root, or pass --gpdfit=)")

  for (e in parse(path)) {
    if (is.call(e) && length(e) >= 3L &&
        as.character(e[[1L]])[1L] %in% c("<-", "=", "<<-") &&
        is.call(e[[3L]]) &&
        identical(as.character(e[[3L]][[1L]])[1L], "function"))
      eval(e, envir = envir)
  }

  needed <- c("sim_series", "auto_run_length", "weights_fun", "gpdfit_dens",
              "geom_weighted_gpdfit_dens", "qgeneralized_pareto", "cred_interval")
  missing <- needed[!vapply(needed, exists, logical(1), envir = envir)]
  if (length(missing))
    stop("weighted_gpdfit.R did not provide: ", paste(missing, collapse = ", "))

  invisible(TRUE)
}


# -----------------------------------------------------------------------------
# sim_dep_series()
#
# Single entry point mapping the task's three dependence structures onto the
# two generators in weighted_gpdfit.R.
#
#   dep_type   parameter   grid (per bias_sensitivity.md)   low end -> high end
#   --------   ---------   -----------------------------   -------------------
#   "ar1"      phi         (0.1, 0.98)     approx. independence -> total dependence
#   "log"      alpha       (0.1, 0.98)     total dependence     -> approx. independence
#   "nlog"     alpha       (0.1, 30)       total independence   -> some dependence
# -----------------------------------------------------------------------------
sim_dep_series <- function(n, dep_type, dep_par, sigma = 1, k = 0) {
  switch(dep_type,
    "ar1"  = sim_series(n, dep = "ar1", phi = dep_par,
                        mu = 0, sigma = sigma, k = k),
    "log"  = sim_series(n, dep = "markov", alpha = dep_par, model = "log",
                        mu = 0, sigma = sigma, k = k),
    "nlog" = sim_series(n, dep = "markov", alpha = dep_par, model = "nlog",
                        mu = 0, sigma = sigma, k = k),
    stop("unknown dep_type: ", dep_type, " (expected one of ",
         paste(DEP_TYPES, collapse = ", "), ")")
  )
}


# -----------------------------------------------------------------------------
# fit_approach()
#
# Fit one approach and return the three point estimates plus, if `density`,
# the posterior grids needed for credible intervals.
#
# All three approaches share the same Zhang-Stephens backend and differ only
# in the weights supplied:
#   Conventional          r = 0  -> every exceedance is its own cluster
#   Cluster maxima (POT)  q = 0  -> all weight on the cluster maximum
#   Weighted GPD                 -> CRPS-optimal geometric weights
# -----------------------------------------------------------------------------
fit_approach <- function(approach, x, u, run_s, rparams,
                         n_grid = 50L, min_grid_pts = 100L, wip = FALSE,
                         density = FALSE) {
  switch(approach,
    "Weighted GPD" =
      geom_weighted_gpdfit_dens(x, u, run_s, n_grid = n_grid, wip = wip,
                                min_grid_pts = min_grid_pts, sort_x = FALSE,
                                rparams = rparams, density = density),
    "Cluster maxima (POT)" = {
      w <- weights_fun(x, u, r = run_s, q = 0)
      gpdfit_dens(w$excesses, wip = wip, min_grid_pts = min_grid_pts,
                  sort_x = TRUE, weights = w$weights,
                  rparams = rparams, density = density)
    },
    "Conventional" = {
      w <- weights_fun(x, u, r = 0L)
      gpdfit_dens(w$excesses, wip = wip, min_grid_pts = min_grid_pts,
                  sort_x = TRUE, weights = w$weights,
                  rparams = rparams, density = density)
    },
    stop("unknown approach: ", approach)
  )
}


# -----------------------------------------------------------------------------
# bias_one_cell()
#
# Runs S replicates for a single (dep_type, dep_par, k, n, r, m) cell.
#
# Reproducibility: one independent L'Ecuyer-CMRG substream per replicate, so
# results do not depend on `ncores` or on scheduling order.
#
# Returns a tidy long data.frame with one row per approach x target parameter.
# -----------------------------------------------------------------------------
bias_one_cell <- function(dep_type,
                          dep_par,
                          k,
                          n,
                          r            = NULL,      # NULL -> auto_run_length()
                          S            = 1000L,
                          sigma        = 1,
                          p_thresh     = 0.95,
                          threshold    = c("theoretical", "empirical"),
                          m            = 1000,      # return period
                          approaches   = APPROACHES,
                          coverage     = TRUE,
                          cred_mass    = 0.95,
                          interval     = c("hpd", "eti"),
                          n_draws      = 10000L,
                          n_grid       = 50L,
                          min_grid_pts = 100L,
                          wip          = FALSE,
                          seed         = 1L,
                          ncores       = 1L,
                          verbose      = TRUE) {

  threshold <- match.arg(threshold)
  interval  <- match.arg(interval)
  stopifnot(dep_type %in% DEP_TYPES)
  if (!all(approaches %in% APPROACHES))
    stop("unknown approach(es): ",
         paste(setdiff(approaches, APPROACHES), collapse = ", "))
  if (dep_type != "ar1")
    stopifnot(requireNamespace("POT", quietly = TRUE))
  if (coverage && interval == "hpd")
    stopifnot(requireNamespace("HDInterval", quietly = TRUE))

  A <- length(approaches)

  # --- cores ---------------------------------------------------------------
  if (is.null(ncores) || is.na(ncores))
    ncores <- max(1L, parallel::detectCores() - 1L)
  ncores <- max(1L, min(as.integer(ncores), S))
  if (ncores > 1L && .Platform$OS.type == "windows") {
    warning("forking unavailable on Windows; running serially", call. = FALSE)
    ncores <- 1L
  }

  # --- fixed truths --------------------------------------------------------
  u_fixed <- qgeneralized_pareto(p_thresh,   mu = 0, sigma = sigma, k = k)
  zm_true <- qgeneralized_pareto(1 - 1 / m,  mu = 0, sigma = sigma, k = k)

  # per replicate we store, for each approach: k, sigma, xp,
  # hit_k, hit_zm, wid_k, wid_zm  (7 slots), plus 3 shared values.
  n_slot <- 7L
  SHARED <- 3L                      # ok, sigma_true, n_exc

  one_rep <- function(stream) {

    assign(".Random.seed", stream, envir = globalenv())
    out <- numeric(SHARED + A * n_slot)

    x <- sim_dep_series(n, dep_type, dep_par, sigma = sigma, k = k)

    u <- if (threshold == "theoretical") u_fixed
         else as.numeric(stats::quantile(x, probs = p_thresh))

    exc <- x > u
    n_exc <- sum(exc)
    if (n_exc < 10L) return(out)             # unusable replicate

    rparams <- list(u = u, m = m, p_u = mean(exc))
    run_s   <- if (is.null(r)) auto_run_length(x, u) else r

    vals <- numeric(A * n_slot)
    for (j in seq_len(A)) {
      f <- tryCatch(
        fit_approach(approaches[j], x, u, run_s, rparams,
                     n_grid = n_grid, min_grid_pts = min_grid_pts,
                     wip = wip, density = coverage),
        error = function(e) NULL)

      if (is.null(f) || !is.finite(f$k) || !is.finite(f$sigma))
        return(out)                          # paired design: drop whole rep

      o <- (j - 1L) * n_slot
      vals[o + 1L] <- f$k
      vals[o + 2L] <- f$sigma
      vals[o + 3L] <- if (is.finite(f$xp)) f$xp else NA_real_

      if (coverage) {
        ci_k <- tryCatch(cred_interval(f$k_grid, f$k_dens, f$w_theta,
                                       cred_mass, interval, n_draws),
                         error = function(e) c(NA_real_, NA_real_))
        ci_z <- tryCatch(cred_interval(f$xp_grid, f$xp_dens, f$w_theta,
                                       cred_mass, interval, n_draws),
                         error = function(e) c(NA_real_, NA_real_))
        vals[o + 4L] <- as.numeric(k       >= ci_k[1L] & k       <= ci_k[2L])
        vals[o + 5L] <- as.numeric(zm_true >= ci_z[1L] & zm_true <= ci_z[2L])
        vals[o + 6L] <- ci_k[2L] - ci_k[1L]
        vals[o + 7L] <- ci_z[2L] - ci_z[1L]
      } else {
        vals[o + (4L:7L)] <- NA_real_
      }
    }

    c(1, sigma + k * u, n_exc, vals)
  }

  # --- independent RNG substreams, one per replicate -----------------------
  old_kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old_kind[1L]), add = TRUE)
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  streams <- vector("list", S)
  streams[[1L]] <- .Random.seed
  for (j in seq_len(S - 1L))
    streams[[j + 1L]] <- parallel::nextRNGStream(streams[[j]])

  t0  <- Sys.time()
  res <- if (ncores > 1L)
    parallel::mclapply(streams, one_rep, mc.cores = ncores, mc.preschedule = TRUE)
  else
    lapply(streams, one_rep)

  want <- SHARED + A * n_slot
  bad  <- !vapply(res, function(z) is.numeric(z) && length(z) == want, logical(1))
  if (any(bad)) res[bad] <- list(numeric(want))

  R  <- matrix(unlist(res, use.names = FALSE), nrow = want)   # want x S
  ok <- R[1L, ] == 1
  n_ok <- sum(ok)

  sig_true <- R[2L, ok]
  n_exc    <- R[3L, ok]

  # --- summarise one approach x one target ---------------------------------
  summarise <- function(app, par, est, truth, hit, wid) {
    fin <- is.finite(est) & is.finite(truth)
    e   <- est[fin]; tv <- truth[fin]
    if (!length(e))
      return(data.frame(approach = app, parameter = par, truth = NA_real_,
                        mean_est = NA_real_, bias = NA_real_, rmse = NA_real_,
                        mae = NA_real_, sd_est = NA_real_,
                        coverage = NA_real_, ci_width = NA_real_,
                        n_ok = 0L))
    d <- e - tv
    data.frame(
      approach = app,
      parameter = par,
      truth     = mean(tv),
      mean_est  = mean(e),
      bias      = mean(d),
      rmse      = sqrt(mean(d^2)),
      mae       = mean(abs(d)),
      sd_est    = stats::sd(e),
      coverage  = if (is.null(hit)) NA_real_ else mean(hit[fin], na.rm = TRUE),
      ci_width  = if (is.null(wid)) NA_real_ else mean(wid[fin], na.rm = TRUE),
      n_ok      = length(e)
    )
  }

  rows <- vector("list", A * 3L)
  i <- 0L
  for (j in seq_len(A)) {
    o   <- SHARED + (j - 1L) * n_slot
    app <- approaches[j]
    est_k <- R[o + 1L, ok]; est_s <- R[o + 2L, ok]; est_z <- R[o + 3L, ok]
    hit_k <- R[o + 4L, ok]; hit_z <- R[o + 5L, ok]
    wid_k <- R[o + 6L, ok]; wid_z <- R[o + 7L, ok]

    rows[[i <- i + 1L]] <- summarise(app, "k",     est_k, rep(k, n_ok),
                                     if (coverage) hit_k, if (coverage) wid_k)
    rows[[i <- i + 1L]] <- summarise(app, "sigma", est_s, sig_true,
                                     NULL, NULL)
    rows[[i <- i + 1L]] <- summarise(app, "zm",    est_z, rep(zm_true, n_ok),
                                     if (coverage) hit_z, if (coverage) wid_z)
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # --- cell descriptors, repeated on every row (tidy long format) ----------
  out$dep_type  <- dep_type
  out$dep_par   <- dep_par
  out$k_true    <- k
  out$n         <- n
  out$r         <- if (is.null(r)) NA_integer_ else as.integer(r)
  out$r_label   <- if (is.null(r)) "auto" else as.character(r)
  out$n_exc     <- if (n_ok > 0L) mean(n_exc) else NA_real_   # realised
  out$n_exc_nom <- n * (1 - p_thresh)                         # nominal
  out$S         <- S
  out$m         <- m
  out$sigma     <- sigma
  out$p_thresh  <- p_thresh
  out$threshold <- threshold
  out$cred_mass <- if (coverage) cred_mass else NA_real_
  out$interval  <- if (coverage) interval else NA_character_
  out$seed      <- seed
  out$secs      <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (verbose)
    cat(sprintf("%s=%g k=%+.2f n=%d r=%s m=%g | %d/%d ok | n_exc~%.0f | %.1fs\n",
                dep_type, dep_par, k, n, out$r_label[1L], m, n_ok, S,
                out$n_exc[1L], out$secs[1L]))

  out
}
