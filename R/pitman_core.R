# =============================================================================
# pitman_core.R
#
# Monte-Carlo engine for the Pitman-closeness / relative-efficiency experiment
# described in experiments.md.
#
# For a single true shape value k on the grid it compares two GPD fitting
# procedures (T1 = challenger, T2 = baseline) on the SAME simulated samples,
# for BOTH parameters:
#
#   * k       -- the shape
#   * sigma*  -- the true scale ABOVE the threshold,  sigma* = sigma + k * u
#
# and returns, per parameter,
#
#   p_hat    = (1/S) sum_s 1{ |T1_s - theta| < |T2_s - theta| }   (ties -> 0.5)
#   se       = sqrt( p_hat (1 - p_hat) / n_ok )
#   rel_eff  = E[(T2 - theta)^2] / E[(T1 - theta)^2]      (> 1  ->  T1 better)
#   rmse_t1, rmse_t2, bias_t1, bias_t2, ...
#
# This file defines functions only -- it runs nothing on source().
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a


# -----------------------------------------------------------------------------
# source_gpd_functions()
#
# Loads the estimator machinery from weighted_gpdfit.R.
#
# NOTE: weighted_gpdfit.R ends with a live scratch block (`sim <- sim_compare(...)`,
# `x11()`, `plot_sim_compare()`), so a plain source() would launch a full
# simulation and try to open an X11 device on every pipeline job.  We therefore
# parse the file and evaluate ONLY the top-level function definitions.
# -----------------------------------------------------------------------------
source_gpd_functions <- function(path = "weighted_gpdfit.R",
                                 envir = parent.frame()) {

  if (!file.exists(path))
    stop("cannot find '", path, "' (run from the project root, or pass --gpdfit=)")

  exprs <- parse(path)
  n_fun <- 0L
  for (e in exprs) {
    if (is.call(e) && length(e) >= 3L &&
        as.character(e[[1L]])[1L] %in% c("<-", "=", "<<-") &&
        is.call(e[[3L]]) &&
        identical(as.character(e[[3L]][[1L]])[1L], "function")) {
      eval(e, envir = envir)
      n_fun <- n_fun + 1L
    }
  }

  needed <- c("sim_series", "auto_run_length", "weights_fun",
              "gpdfit", "geom_weighted_gpdfit", "qgeneralized_pareto")
  missing <- needed[!vapply(needed, exists, logical(1), envir = envir)]
  if (length(missing))
    stop("weighted_gpdfit.R did not provide: ", paste(missing, collapse = ", "))

  invisible(n_fun)
}


# -----------------------------------------------------------------------------
# make_k_grid() -- the single source of truth for the grid, so that the R
# runner and the Snakemake job indices always agree.
# -----------------------------------------------------------------------------
make_k_grid <- function(k_min = -0.5, k_max = 0.5, n_k = 20L)
  seq(k_min, k_max, length.out = n_k)


PROCEDURES <- c("Weighted GPD", "Cluster maxima (POT)", "Conventional")


# -----------------------------------------------------------------------------
# fit_shape_scale()
#
# Fit one procedure to series `x` above threshold `u`; returns c(k_hat, sigma_hat)
# where sigma_hat is the scale of the excesses, i.e. an estimate of sigma*.
# Returns c(NA, NA) if the fit fails.
# -----------------------------------------------------------------------------
fit_shape_scale <- function(proc, x, u, run_s,
                            n_grid = 50L, min_grid_pts = 100L, wip = FALSE) {
  f <- tryCatch(switch(
    proc,
    "Weighted GPD" =
      geom_weighted_gpdfit(x, u, run_s, n_grid = n_grid, wip = wip,
                           min_grid_pts = min_grid_pts, sort_x = FALSE),
    "Cluster maxima (POT)" = {
      w <- weights_fun(x, u, r = run_s, q = 0)
      gpdfit(w$excesses, wip = wip, min_grid_pts = min_grid_pts,
             sort_x = TRUE, weights = w$weights)
    },
    "Conventional" = {
      w <- weights_fun(x, u, r = 0L)
      gpdfit(w$excesses, wip = wip, min_grid_pts = min_grid_pts,
             sort_x = TRUE, weights = w$weights)
    },
    stop("unknown procedure: ", proc)
  ), error = function(e) NULL)

  if (is.null(f)) return(c(NA_real_, NA_real_))
  c(as.numeric(f$k), as.numeric(f$sigma))
}


# -----------------------------------------------------------------------------
# pitman_one_k()
#
# Runs S Monte-Carlo replicates at a single true shape value.
#
# Reproducibility: one independent L'Ecuyer-CMRG substream is created per
# replicate, so the answer is bit-identical whatever `ncores` is (and whatever
# order the scheduler happens to finish jobs in).
#
# Returns a tidy 2-row data.frame (one row per parameter).
# -----------------------------------------------------------------------------
pitman_one_k <- function(k_true,
                         S            = 100L,
                         n            = 2000L,
                         dep          = c("ar1", "markov"),
                         phi          = 0.95,
                         alpha        = 0.5,
                         dep_model    = "log",
                         sigma        = 1,
                         p_thresh     = 0.95,
                         threshold    = c("theoretical", "empirical"),
                         run          = NULL,
                         n_grid       = 50L,
                         min_grid_pts = 100L,
                         wip          = FALSE,
                         t1           = "Weighted GPD",
                         t2           = "Conventional",
                         seed         = 1L,
                         ncores       = 1L,
                         verbose      = TRUE) {

  dep       <- match.arg(dep)
  threshold <- match.arg(threshold)
  if (!t1 %in% PROCEDURES) stop("bad t1: ", t1)
  if (!t2 %in% PROCEDURES) stop("bad t2: ", t2)
  if (identical(t1, t2)) stop("`t1` and `t2` must be two different procedures.")
  if (dep == "markov") stopifnot(requireNamespace("POT", quietly = TRUE))

  # --- cores --------------------------------------------------------------
  if (is.null(ncores) || is.na(ncores))
    ncores <- max(1L, parallel::detectCores() - 1L)
  ncores <- max(1L, min(as.integer(ncores), S))
  if (ncores > 1L && .Platform$OS.type == "windows") {
    warning("forking unavailable on Windows; running serially", call. = FALSE)
    ncores <- 1L
  }

  u_fixed <- qgeneralized_pareto(p_thresh, mu = 0, sigma = sigma, k = k_true)

  # --- one replicate ------------------------------------------------------
  # returns c(ok, k_t1, k_t2, s_t1, s_t2, sigma_star_true)
  one_rep <- function(stream) {

    assign(".Random.seed", stream, envir = globalenv())

    x <- sim_series(n, dep = dep, phi = phi, alpha = alpha, model = dep_model,
                    mu = 0, sigma = sigma, k = k_true)

    u <- if (threshold == "theoretical") u_fixed
         else as.numeric(stats::quantile(x, probs = p_thresh))

    run_s <- if (is.null(run)) auto_run_length(x, u) else run

    f1 <- fit_shape_scale(t1, x, u, run_s, n_grid, min_grid_pts, wip)
    f2 <- fit_shape_scale(t2, x, u, run_s, n_grid, min_grid_pts, wip)

    if (!all(is.finite(c(f1, f2)))) return(c(0, 0, 0, 0, 0, 0))

    c(1, f1[1L], f2[1L], f1[2L], f2[2L], sigma + k_true * u)
  }

  # --- independent RNG substreams, one per replicate ----------------------
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

  bad <- !vapply(res, function(z) is.numeric(z) && length(z) == 6L, logical(1))
  if (any(bad)) res[bad] <- list(c(0, 0, 0, 0, 0, 0))

  R  <- matrix(unlist(res, use.names = FALSE), nrow = 6L)   # 6 x S
  ok <- R[1L, ] == 1
  n_ok <- sum(ok)

  # --- summarise one parameter -------------------------------------------
  # e1, e2 = signed errors of T1 / T2 ; truth may vary by replicate (sigma*)
  summarise_par <- function(par, est1, est2, truth) {
    if (n_ok == 0L)
      return(data.frame(parameter = par, k_true = k_true, truth = NA_real_,
                        p_hat = NA_real_, se = NA_real_, rel_eff = NA_real_,
                        rmse_t1 = NA_real_, rmse_t2 = NA_real_,
                        bias_t1 = NA_real_, bias_t2 = NA_real_,
                        mean_t1 = NA_real_, mean_t2 = NA_real_,
                        n_wins = NA_integer_, n_ties = NA_integer_,
                        n_ok = 0L, S = S))

    e1 <- est1 - truth
    e2 <- est2 - truth
    a1 <- abs(e1); a2 <- abs(e2)

    score <- ifelse(a1 < a2, 1, ifelse(a1 > a2, 0, 0.5))
    p     <- mean(score)

    mse1 <- mean(e1^2); mse2 <- mean(e2^2)

    data.frame(
      parameter = par,
      k_true    = k_true,
      truth     = mean(truth),
      p_hat     = p,
      se        = sqrt(p * (1 - p) / n_ok),
      rel_eff   = mse2 / mse1,                 # > 1  ->  T1 preferable
      rmse_t1   = sqrt(mse1),
      rmse_t2   = sqrt(mse2),
      bias_t1   = mean(e1),
      bias_t2   = mean(e2),
      mean_t1   = mean(est1),
      mean_t2   = mean(est2),
      n_wins    = sum(score == 1),
      n_ties    = sum(score == 0.5),
      n_ok      = n_ok,
      S         = S
    )
  }

  out <- rbind(
    summarise_par("k",     R[2L, ok], R[3L, ok], k_true),
    summarise_par("sigma", R[4L, ok], R[5L, ok], R[6L, ok])
  )
  rownames(out) <- NULL

  # --- metadata carried alongside every row (handy after combining) -------
  out$t1        <- t1
  out$t2        <- t2
  out$dep       <- dep
  out$dep_par   <- if (dep == "ar1") phi else alpha
  out$dep_model <- if (dep == "ar1") NA_character_ else dep_model
  out$n         <- n
  out$p_thresh  <- p_thresh
  out$threshold <- threshold
  out$seed      <- seed
  out$secs      <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (verbose)
    cat(sprintf(paste0("k = %+.4f | n_ok = %d/%d | ",
                       "k: P=%.3f e=%.3f | sigma*: P=%.3f e=%.3f | %.1fs\n"),
                k_true, n_ok, S,
                out$p_hat[1L], out$rel_eff[1L],
                out$p_hat[2L], out$rel_eff[2L], out$secs[1L]))

  out
}
