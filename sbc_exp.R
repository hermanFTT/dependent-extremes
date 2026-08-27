


library(patchwork)
#library(ggplot2)
#library(bayesplot)
#source("dependence_aware_unif_tests.R")
#source("gpd.R")
#if (!exists("find_opt_weights"))
source("weighted_gpdfit.R")
#source("R/dependence_aware_unif_tests.R")

devtools::load_all("../posterior")

rank_stat <- function(element, vec) { sum(vec < element) + 1 }


# SBC 
calculate_ranks_draws_matrix <- function(variables, dm, params = NULL) {

  if(!is.null(params)) {
    warning("The `params` argument is deprecated use `variables` instead.")
    if(is.null(variables)) {
      variables <- params
    }
  }

  stopifnot(posterior::is_draws_matrix(dm))
  stopifnot(posterior::nvariables(dm) == length(variables))

  max_rank <- posterior::ndraws(dm)
  nvars <- posterior::nvariables(dm)

  less_matrix <- sweep(dm, MARGIN = 2, STATS = variables, FUN = "<")
  rank_min <- colSums(less_matrix, na.rm = TRUE)

  # When there are ties (e.g. for discrete variables), the rank is currently drawn stochastically
  # among the ties
  # NA is assumed to be potentially equal to any value (Issue #78)
  equal_or_NA <- function(a,b) {
    is.na(a) | is.na(b) | a == b
  }
  equal_matrix <- sweep(dm, MARGIN = 2, STATS = variables, FUN = equal_or_NA)
  rank_range <- colSums(equal_matrix)

  ranks <- rank_min + rdunif(posterior::nvariables(dm), a = 0, b = rank_range)

  attr(ranks, "max_rank") <- max_rank
  ranks
}





#########

#' Pareto-smoothed probability integral transform
#'
#' Compute PIT values using the empirical CDF, then refine values in
#' the tails by fitting a generalized Pareto distribution (GPD) to
#' the tail draws. This gives smoother, more accurate PIT values in
#' the tails where the ECDF is coarse, and avoids PIT values of 0 and 1.
#' Due to use of generalized Pareto distribution CDF in tails, the
#' PIT values are not anymore rank based and continuous uniformity
#' test is appropriate.
#'
#' @name pareto_pit
#'
#' @param x (draws) A [`draws_matrix`] object or one coercible to a
#'   `draws_matrix` object, or an [`rvar`] object.
#'
#' @param y (observations) A 1D vector, or an array of dim(x), if x is `rvar`.
#'   Each element of `y` corresponds to a variable in `x`.
#'
#' @param weights A matrix of weights for each draw and variable. `weights`
#'   should have one column per variable in `x`, and `ndraws(x)` rows.
#'
#' @param log (logical) Are the weights passed already on the log scale? The
#'   default is `FALSE`, that is, expecting `weights` to be on the standard
#'   (non-log) scale.
#'
#' @param ndraws_tail (integer) Number of tail draws to use for GPD
#'   fitting. If `NULL` (the default), computed using [ps_tail_length()].
#'
#' @template args-methods-dots
#'
#' @details The function first computes raw PIT values identically to
#'   [pit()] (including support for weighted draws). It then fits a
#'   GPD to both tails of the draws (using the same approach as
#'   [pareto_smooth()]) and replaces PIT values for observations falling in
#'   the tail regions:
#'
#'   For a right-tail observation \eqn{y_i > c_R} (where \eqn{c_R} is
#'   the right-tail cutoff):
#'
#'   \deqn{PIT(y_i) = 1 - p_{tail}(1 - F_{GPD}(y_i; c_R, \sigma_R, k_R))}
#'
#'   For a left-tail observation \eqn{y_i < c_L}:
#'
#'   \deqn{PIT(y_i) = p_{tail}(1 - F_{GPD}(-y_i; -c_L, \sigma_L, k_L))}
#'
#'   where \eqn{p_{tail}} is the proportion of (weighted) mass in the
#'   tail.
#'
#'   When (log-)weights in `weights` are provided, they are used for
#'   the raw PIT computation (as in [pit()]) and for GPD fit.
#'
#' @return A numeric vector of length `length(y)` containing the PIT values, or
#'   an array of shape `dim(y)`, if `x` is an `rvar`.
#'
#' @seealso [pit()] for the unsmoothed version, [pareto_smooth()] for
#'   Pareto smoothing of draws.
#'
#' @examples
#' x <- example_draws()
#' y <- rnorm(nvariables(x), 5, 5)
#' posterior::pareto_pit(x, y)


# ============================================================
# Matrix version of pareto_pit
# x: plain matrix (ndraws x nvariables)
# y: numeric vector (length nvariables)
# All draws-object helpers replaced by base R equivalents.
# Returns an unnamed numeric vector of pointwise PIT values.
# ============================================================
pareto_pit_mat <- function(x, y, ndraws_tail = NULL, weights = NULL,
                           log = FALSE) {
  stopifnot(is.matrix(x), is.numeric(y), ncol(x) == length(y))

  ndraws <- nrow(x)
  nvar   <- ncol(x)

  # convert and normalize weights column-wise to log scale
  if (!is.null(weights)) {
    if (!log) weights <- base::log(weights)
    weights <- apply(weights, 2, function(col) col - log_sum_exp(col))
  }

  if (is.null(ndraws_tail)) {
    ndraws_tail <- ps_tail_length(ndraws, 1)
  } else {
    ndraws_tail <- as.integer(ndraws_tail)
  }

  gpd_ok <- !is.na(ndraws_tail) && ndraws_tail >= 5
  if (gpd_ok && ndraws_tail > ndraws / 2) ndraws_tail <- floor(ndraws / 2)
  if (gpd_ok && ndraws_tail >= ndraws)     gpd_ok <- FALSE

  if (gpd_ok) tail_ids <- seq(ndraws - ndraws_tail + 1, ndraws)

  pit_values <- vapply(seq_len(nvar), function(j) {
    draws <- x[, j]

    # --- raw PIT ---
    sel_min <- draws < y[j]
    if (!any(sel_min)) {
      raw_pit <- 0
    } else {
      if (is.null(weights)) {
        raw_pit <- mean(sel_min)
      } else {
        raw_pit <- exp(log_sum_exp(weights[sel_min, j]))
      }
    }

    sel_sup <- draws == y[j]
    if (any(sel_sup)) {
      if (is.null(weights)) {
        pit_sup <- raw_pit + mean(sel_sup)
      } else {
        pit_sup <- raw_pit + exp(log_sum_exp(weights[sel_sup, j]))
      }
      raw_pit <- runif(1, raw_pit, pit_sup)
    }

    if (!gpd_ok) return(raw_pit)

    # sort draws and carry weights along
    ord           <- sort.int(draws, index.return = TRUE)
    sorted        <- ord$x
    log_wt_sorted <- if (!is.null(weights)) weights[ord$ix, j] else NULL

    # tail proportion
    if (!is.null(log_wt_sorted)) {
      tail_proportion <- exp(log_sum_exp(log_wt_sorted[tail_ids]))
    } else {
      tail_proportion <- ndraws_tail / ndraws
    }

    # --- right tail ---
    right_replaced <- FALSE
    right_tail     <- sorted[tail_ids]
    if (!is_constant(right_tail)) {
      right_cutoff <- sorted[min(tail_ids) - 1]
      if (right_cutoff == right_tail[1])
        right_cutoff <- right_cutoff - .Machine$double.eps

      right_tail_wt <- if (!is.null(log_wt_sorted)) exp(log_wt_sorted[tail_ids]) else NULL
      right_fit <- gpdfit(right_tail - right_cutoff, sort_x = FALSE,
                          weights = right_tail_wt)
      if (is.finite(right_fit$k) && !is.na(right_fit$sigma)) {
        if (y[j] > right_cutoff) {
          gpd_cdf <- pgeneralized_pareto(
            y[j], mu = right_cutoff, sigma = right_fit$sigma, k = right_fit$k
          )
          raw_pit        <- 1 - tail_proportion * (1 - gpd_cdf)
          right_replaced <- TRUE
        }
      }
    }

    # --- left tail (negate trick) ---
    if (!right_replaced) {
      left_draws         <- -draws
      left_ord           <- sort.int(left_draws, index.return = TRUE)
      left_sorted        <- left_ord$x
      log_wt_left_sorted <- if (!is.null(weights)) weights[left_ord$ix, j] else NULL

      left_tail <- left_sorted[tail_ids]
      if (!is_constant(left_tail)) {
        left_cutoff <- left_sorted[min(tail_ids) - 1]
        if (left_cutoff == left_tail[1])
          left_cutoff <- left_cutoff - .Machine$double.eps

        left_tail_wt <- if (!is.null(log_wt_left_sorted)) exp(log_wt_left_sorted[tail_ids]) else NULL
        left_fit <- gpdfit(left_tail - left_cutoff, sort_x = FALSE,
                           weights = left_tail_wt)
        if (is.finite(left_fit$k) && !is.na(left_fit$sigma)) {
          if (-y[j] > left_cutoff) {
            gpd_cdf <- pgeneralized_pareto(
              -y[j], mu = left_cutoff, sigma = left_fit$sigma, k = left_fit$k
            )
            left_tail_proportion <- if (!is.null(log_wt_left_sorted)) {
              exp(log_sum_exp(log_wt_left_sorted[tail_ids]))
            } else {
              tail_proportion
            }
            raw_pit <- left_tail_proportion * (1 - gpd_cdf)
          }
        }
      }
    }

    raw_pit
  }, FUN.VALUE = 1.0)

  min_tail_prob <- 1 / ndraws / 1e4
  pmin(pmax(pit_values, min_tail_prob), 1 - min_tail_prob)
}



# ---------------------------------------------------------------------------
# pareto_pit_opt_weights_mat()
# Same as pareto_pit but uses find_opt_weights() (weigthed_gpdfit.R) to determine tail weights.

# Key difference in weight handling:
#   find_opt_weights() returns $exc / $wts already *sorted by excess value*
#   (ascending), so we pass them directly to gpdfit(sort_x = FALSE) without
#   any index-matching.  User-supplied weights (if any) are reordered to the
#   same excess-sorted permutation before combining.
# x are considered replicate (y_rep) and y are the observations. 
# ---------------------------------------------------------------------------
pareto_pit_opt_weights_mat <- function(x, y, run=NULL, ndraws_tail = NULL,
                                       weights = NULL, log = FALSE,
                                       n_grid = 50) {
  stopifnot(is.matrix(x), is.numeric(y), ncol(x) == length(y))

  ndraws <- nrow(x)
  nvar   <- ncol(x)

  # Normalise user weights column-wise to log scale (once, outside the loop)
  if (!is.null(weights)) {
    if (!log) weights <- base::log(weights)
    weights <- apply(weights, 2, function(col) col - log_sum_exp(col))
  }

  if (is.null(ndraws_tail)) {
    ndraws_tail <- ps_tail_length(ndraws, 1)
  } else {
    ndraws_tail <- as.integer(ndraws_tail)
  }

  gpd_ok <- !is.na(ndraws_tail) && ndraws_tail >= 5
  if (gpd_ok && ndraws_tail > ndraws / 2) ndraws_tail <- floor(ndraws / 2)
  if (gpd_ok && ndraws_tail >= ndraws)     gpd_ok <- FALSE

  if (gpd_ok) tail_ids <- seq(ndraws - ndraws_tail + 1L, ndraws)

  # combine_wts_sorted():
  #   log_w_user_sorted  -- user log-weights already reordered to match the
  #                         excess-sorted order that find_opt_weights returns
  #   w_opt              -- numeric weights from find_opt_weights (same order)
  # Returns combined (product) weights on the original scale.
  combine_wts_sorted <- function(log_w_user_sorted, w_opt) {
    if (!is.null(log_w_user_sorted)) exp(log_w_user_sorted) * w_opt else w_opt
  }

  pit_values <- vapply(seq_len(nvar), function(j) {
    draws <- x[, j]

    # ---- raw PIT (no tail smoothing) ----------------------------------------
    sel_min <- draws < y[j]
    raw_pit <- if (!any(sel_min)) {
      0
    } else if (is.null(weights)) {
      mean(sel_min)
    } else {
      exp(log_sum_exp(weights[sel_min, j]))
    }

    sel_sup <- draws == y[j]
    if (any(sel_sup)) {
      pit_sup <- raw_pit + if (is.null(weights)) mean(sel_sup) else exp(log_sum_exp(weights[sel_sup, j]))
      raw_pit <- runif(1, raw_pit, pit_sup)
    }

    if (!gpd_ok) return(raw_pit)

    # Sort draws once; reused by both tail branches.
    ord          <- sort.int(draws, index.return = TRUE)
    sorted       <- ord$x
    # User log-weights in the sorted-draw order (NULL if no user weights)
    log_w_user_s <- if (!is.null(weights)) weights[ord$ix, j] else NULL

    # ---- right tail ----------------------------------------------------------
    right_replaced <- FALSE
    right_tail     <- sorted[tail_ids]

    if (!is_constant(right_tail)) {
      right_cutoff <- sorted[min(tail_ids) - 1L]
      if (right_cutoff == right_tail[1L])
        right_cutoff <- right_cutoff - .Machine$double.eps

      if (y[j] > right_cutoff) {
        # find_opt_weights() returns exc & wts sorted by excess (ascending)
        ep_r <- find_opt_weights(draws, right_cutoff, run, n_grid = n_grid)

        # Tail proportion: weighted mass above the cutoff
        tail_proportion <- if (!is.null(log_w_user_s)) {
          exp(log_sum_exp(log_w_user_s[tail_ids]))
        } else {
          ndraws_tail / ndraws
        }

        # Reorder user weights to match the excess-sorted order of ep_r$exc.
        # ep_r$exc are the sorted excesses; the corresponding original indices
        # are the tail draw indices sorted by their excess value.
        if (!is.null(log_w_user_s)) {
          # indices of tail draws sorted by excess (ascending)
          exc_order_in_tail  <- order(right_tail - right_cutoff)
          log_w_user_exc_r   <- log_w_user_s[tail_ids][exc_order_in_tail]
        } else {
          log_w_user_exc_r <- NULL
        }

        right_fit <- gpdfit(ep_r$exc, sort_x = FALSE,
                            weights = combine_wts_sorted(log_w_user_exc_r,
                                                         ep_r$wts))
        if (is.finite(right_fit$k) && !is.na(right_fit$sigma)) {
          gpd_cdf <- pgeneralized_pareto(y[j], mu = right_cutoff,
                                         sigma = right_fit$sigma,
                                         k     = right_fit$k)
          raw_pit        <- 1 - tail_proportion * (1 - gpd_cdf)
          right_replaced <- TRUE
        }
      }
    }

    # ---- left tail (negation trick) — only when right tail did not fire ------
    if (!right_replaced) {
      left_draws  <- -draws
      left_ord    <- sort.int(left_draws, index.return = TRUE)
      left_sorted <- left_ord$x
      log_w_user_ls <- if (!is.null(weights)) weights[left_ord$ix, j] else NULL

      left_tail <- left_sorted[tail_ids]
      if (!is_constant(left_tail)) {
        left_cutoff <- left_sorted[min(tail_ids) - 1L]
        if (left_cutoff == left_tail[1L])
          left_cutoff <- left_cutoff - .Machine$double.eps

        if (-y[j] > left_cutoff) {
          ep_l <- find_opt_weights(left_draws, left_cutoff, run,
                                   n_grid = n_grid)

          left_tail_proportion <- if (!is.null(log_w_user_ls)) {
            exp(log_sum_exp(log_w_user_ls[tail_ids]))
          } else {
            ndraws_tail / ndraws
          }

          if (!is.null(log_w_user_ls)) {
            exc_order_in_ltail <- order(left_tail - left_cutoff)
            log_w_user_exc_l   <- log_w_user_ls[tail_ids][exc_order_in_ltail]
          } else {
            log_w_user_exc_l <- NULL
          }

          left_fit <- gpdfit(ep_l$exc, sort_x = FALSE,
                             weights = combine_wts_sorted(log_w_user_exc_l,
                                                          ep_l$wts))
          if (is.finite(left_fit$k) && !is.na(left_fit$sigma)) {
            gpd_cdf <- pgeneralized_pareto(-y[j], mu = left_cutoff,
                                            sigma = left_fit$sigma,
                                            k     = left_fit$k)
            raw_pit <- left_tail_proportion * (1 - gpd_cdf)
          }
        }
      }
    }

    raw_pit
  }, FUN.VALUE = 1.0)

  min_tail_prob <- 1 / ndraws / 1e4
  pmin(pmax(pit_values, min_tail_prob), 1 - min_tail_prob)
}



#########
#library(gnorm)

N <-2000 #4096-1
M <-10000
phi=0.9
#y <- rnorm(M, 0, 1) #rgnorm(M, mu = 0, alpha = sqrt(2), beta =2.5) #  # #rt(M, df = 30)  #rt(M, df = 4) # rgnorm(M, mu = 0, alpha = sqrt(2), beta = 4)  # exactly normal#rt(M, df = 100) # rnorm(M, 0, 1) # i.i.d e.g., prior draws  #Simulated autocorrelated e.g., MCMC draws

#X_ar<- replicate(M, arima.sim(n = N, list(ar =c(phi)), sd = sqrt(1 - phi^2)))
#library(POT)

y <- rt(M, df =3)

X_ar <- replicate( M, qt(pnorm(arima.sim(n =N,list(ar = phi), sd = sqrt(1 - phi^2))),  df=3) )


graphics.off()

# original
emp_ranks<- vapply(seq_len(length(y)), function(j) rank_stat(y[j], X_ar[,j]),FUN.VALUE = 1.0)

pit_emp <-emp_ranks/(N+1)

x11()
hist(pit_emp, breaks=N+1)

# weighted GPD-based smoothed tail 
pits1 <- pareto_pit_opt_weights_mat(X_ar, y, run=NULL, ndraws_tail=NULL)

x11()
hist(pits1, breaks=N+1) # we see spikes at the edge ( characterization of auto correlation)

# naive tail smoothing (ignore dependence)
pits2 <-pareto_pit_mat(X_ar,y,ndraws_tail = 100)

x11()
hist(pits2, breaks=N+1)


graphics.off()

uniformity_test(pits1, "POT", truncate = T)$pvalue

uniformity_test(pit_emp, "PRIT",T)$pvalue

#x11()
#graphical_rep(pit_emp,test="POT-C", infl_points_only=FALSE,diff=TRUE) #+
## graphical_rep(pits1,test="POT-C", infl_points_only=FALSE,diff=TRUE)+
## graphical_rep(pits2,test="POT-C", infl_points_only=FALSE,diff=TRUE)

#ppc_pit_ecdf(pit=pit_emp, prob=0.95, plot_diff = TRUE) 



##################################  Check Calibration


S <- 50

phi <- 0.9
N <- 2000
M <- 4000
  
x11()
ranks_mat1<- matrix(numeric(S*M), nrow=S, ncol =M)
ranks_mat2<- matrix(numeric(S*M), nrow=S, ncol =M)
ranks_mat3<- matrix(numeric(S*M), nrow=S, ncol =M)

system.time({
  
for(s in 1:S){
  #print(s)
  
  y <-rnorm(M, 0, 1) #rt(M, df = 4) # rgnorm(M, mu = 0, alpha = sqrt(2), beta = 4)  # exactly normal#rt(M, df = 100) # rnorm(M, 0, 1) # i.i.d e.g., prior draws
  #Simulated autocorrelated e.g., MCMC draws
 
  # compute normalized ranks/pits for each replicate 
  X_ar<- replicate(M, arima.sim(n = N, list(ar =c(phi)), sd = sqrt(1 - phi^2)))
 # ranks_mat1[s,]<-pareto_pit_mat(X_ar,y,ndraws_tail = 100)
  ranks_mat2[s,]<-pareto_pit_opt_weights_mat(X_ar, y, run=NULL, ndraws_tail = NULL)
  ranks_mat3[s, ]<- vapply(seq_len(length(y)), function(i) rank_stat(y[i], X_ar[,i]), FUN.VALUE=1.0)/(N+1)
 }

})


rep_pits <- ranks_mat2 # ranks_mat1,  ranks_mat3

pvals <- apply( rep_pits, 1, function(x) uniformity_test(x, "PRIT")$pvalue )

mean(pvals<0.01)


####################### Application to  SBC ( address bias due to autocorrelation)

#' ## Bias in extreme ranks when using one Markov chain only once
#'
#' The second example address the behavior in SBC, where we
#' compare one $y_j \sim g(y)$ to dependent draws $x \sim p(y)$, but
#' for each $y_j$ we generate new dependent draws $x$. This means that
#' the rank statistics of $y_j$ are independent from each other, but
#' the rank statistics are still biased.
#' 


# ===========================================================================
# SBC tail fixing
#
# Differences from pareto_pit_opt_weights_mat():
#   * draws are NEVER stored in an M x N matrix.  Each y_j gets its own fresh
#     chain, which is consumed and discarded, so memory is O(N) not O(M*N).
#   * no user-supplied weights.
#   * the PIT is the NORMALIZED RANK  r/(N+1),  r = #{x < y} + 1, and is
#     replaced by a (weighted-)GPD PIT only when y falls in a tail.
#
# The weighted fit uses find_opt_weights() (runs declustering + CRPS-optimal
# geometric weights), so the tail fit is not fooled by the autocorrelation.
# weighted = FALSE gives the dependence-naive comparator (plain gpdfit on the
# raw tail excesses), so one code path serves both arms of the comparison.
# ===========================================================================


# ---------------------------------------------------------------------------
# default_tail_length()
# Same rule as posterior:::ps_tail_length(ndraws, r_eff = 1), inlined so this
# section does not depend on posterior internals.
# ---------------------------------------------------------------------------
default_tail_length <- function(ndraws, r_eff = 1) {
  floor(ifelse(ndraws * r_eff > 225, 3 * sqrt(ndraws / r_eff), ndraws / 5))
}


# ---------------------------------------------------------------------------
# sbc_pit_one()
#
# PIT of ONE observation y against ONE chain of dependent draws x.
#
#   bulk  : pit = r/N,  r = #{x < y}                 (ECDF-scale rank)
#   tails : pit from a GPD fitted to the tail draws
#
# Continuity.  With n_t tail draws and cutoffs
#     c_R = x_(N - n_t)     (n_t draws lie strictly above it)
#     c_L = x_(n_t + 1)     (n_t draws lie strictly below it)
# the tail masses are on the SAME scale as the bulk,
#     p_R = p_L = n_t/N
# so that
#     y -> c_R from above :  #{x<y} = N - n_t  ->  pit = 1 - n_t/N = 1 - p_R
#     y -> c_L from below :  #{x<y} = n_t      ->  pit =     n_t/N =     p_L
# i.e. the GPD limit equals the rank value at each cutoff -- no jump at the
# splice.  (The r/(N+1) convention would instead need p_R = n_t/(N+1) and
# p_L = (n_t+1)/(N+1); mixing the two leaves a gap of order 1/N on the left.)
#
#   x         draws IN CHAIN ORDER -- declustering needs the time ordering
#   run       run length; NULL -> auto_run_length() inside find_opt_weights()
#   weighted  TRUE = dependence-aware, FALSE = naive comparator
#   details   TRUE -> list(pit, tail, q_opt, run, rank) instead of a scalar
#
# Any failure of the tail fit falls back to the plain normalized rank.
# ---------------------------------------------------------------------------
sbc_pit_one <- function(y, x, run = NULL, ndraws_tail = NULL, n_grid = 50,
                        weighted = TRUE, wip = TRUE, details = FALSE) {

  stopifnot(length(y) == 1L, is.numeric(x))
  N <- length(x)

  # --- bulk value: raw PIT, exactly as in pareto_pit() but with no user
  #     weights -- the ECDF at y, i.e. #{x < y}/N, with any ties resolved by a
  #     CONTINUOUS uniform draw across the tied probability mass.
  sel_min <- x < y
  pit     <- if (!any(sel_min)) 0 else mean(sel_min)

  sel_sup <- x == y
  if (any(sel_sup)) {
    pit_sup <- pit + mean(sel_sup)
    pit     <- runif(1, pit, pit_sup)
  }
  raw_pit <- pit

  out <- function(p, tail = "bulk", q = NA_real_, rl = NA_integer_) {
    #eps <- 1 / (N + 1) / 1e4
    #p   <- min(max(p, eps), 1 - eps)
    if (details) list(pit = p, tail = tail, q_opt = q, run = rl,
                      raw_pit = raw_pit)
    else p
  }

  n_t <- as.integer(if (is.null(ndraws_tail)) default_tail_length(N)
                    else ndraws_tail)
  if (!is.finite(n_t) || n_t < 5L || n_t >= N / 2) return(out(pit))

  xs  <- sort.int(x)
  c_R <- xs[N - n_t]
  c_L <- xs[n_t + 1L]

  # --- fit one tail; `series` stays in chain order ------------------------
  fit_tail <- function(series, cutoff) {
    if (weighted) {
      ep <- tryCatch(find_opt_weights(series, cutoff, run, n_grid = n_grid),
                     error = function(e) NULL)
      if (is.null(ep) || length(ep$exc) < 5L) return(NULL)
      f <- tryCatch(gpdfit(ep$exc, wip = wip, sort_x = FALSE, weights = ep$wts),
                    error = function(e) NULL)
      if (is.null(f)) return(NULL)
      list(fit = f, q = ep$q_opt, run = ep$run)
    } else {
      exc <- series[series > cutoff] - cutoff
      if (length(exc) < 5L) return(NULL)
      f <- tryCatch(gpdfit(exc, wip = wip, sort_x = TRUE),
                    error = function(e) NULL)
      if (is.null(f)) return(NULL)
      list(fit = f, q = NA_real_, run = NA_integer_)
    }
  }

  ok <- function(f) !is.null(f) && is.finite(f$fit$k) &&
                    !is.na(f$fit$sigma) && f$fit$sigma > 0

  # --- right tail ---------------------------------------------------------
  if (y > c_R) {
    top <- xs[(N - n_t + 1L):N]
    if (diff(range(top)) == 0) return(out(pit))            # constant tail
    if (c_R == top[1L]) c_R <- c_R - .Machine$double.eps   # ensure exc > 0

    ft <- fit_tail(x, c_R)
    if (!ok(ft)) return(out(pit))

    p_R <- n_t / N
    cdf <- pgeneralized_pareto(y, mu = c_R, sigma = ft$fit$sigma, k = ft$fit$k)
    return(out(1 - p_R * (1 - cdf), "right", ft$q, ft$run))
  }

  # --- left tail: negate, so it becomes a right tail ----------------------
  if (y < c_L) {
    bot <- xs[1L:n_t]
    if (diff(range(bot)) == 0) return(out(pit))
    cut_L <- -c_L
    if (cut_L == -bot[n_t]) cut_L <- cut_L - .Machine$double.eps

    ft <- fit_tail(-x, cut_L)
    if (!ok(ft)) return(out(pit))

    p_L <- n_t / N
    cdf <- pgeneralized_pareto(-y, mu = cut_L, sigma = ft$fit$sigma,
                               k = ft$fit$k)
    return(out(p_L * (1 - cdf), "left", ft$q, ft$run))
  }

  out(pit)   # y sits in the bulk: plain normalized rank
}


# ---------------------------------------------------------------------------
# sbc_pit()
#
# Streaming driver over M SBC replicates -- no M x N matrix is ever built.
# `generator(j)` must return list(y = <scalar>, x = <numeric vector>): a fresh
# prior draw and a fresh dependent chain for it.
#
# Reproducibility: one independent L'Ecuyer-CMRG substream per replicate, so
# the result does not depend on ncores.  The generator draws from that stream
# too, so the simulated chains are reproducible as well.
#
# Returns a numeric vector of M PIT values, with attributes "tail", "q_opt"
# and "run" (one entry per replicate) for diagnostics.
# ---------------------------------------------------------------------------
sbc_pit <- function(M, generator, run = NULL, ndraws_tail = NULL, n_grid = 50,
                    weighted = TRUE, wip = FALSE, seed = NULL, ncores = 1L,
                    verbose = TRUE) {

  stopifnot(is.function(generator), M >= 1L)

  if (is.null(ncores) || is.na(ncores))
    ncores <- max(1L, parallel::detectCores() - 1L)
  ncores <- max(1L, min(as.integer(ncores), M))
  if (ncores > 1L && .Platform$OS.type == "windows") {
    warning("forking unavailable on Windows; running serially", call. = FALSE)
    ncores <- 1L
  }

  old_kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(RNGkind(old_kind[1L]), add = TRUE)
  set.seed(if (is.null(seed)) 1L else seed, kind = "L'Ecuyer-CMRG")
  streams <- vector("list", M)
  streams[[1L]] <- .Random.seed
  for (j in seq_len(M - 1L))
    streams[[j + 1L]] <- parallel::nextRNGStream(streams[[j]])

  one <- function(j) {
    assign(".Random.seed", streams[[j]], envir = globalenv())
    rep_j <- generator(j)                    # fresh y and fresh chain
    if (!is.list(rep_j) || is.null(rep_j$y) || is.null(rep_j$x))
      stop("generator(j) must return list(y = scalar, x = numeric vector)")
    d <- sbc_pit_one(rep_j$y, rep_j$x, run = run, ndraws_tail = ndraws_tail,
                     n_grid = n_grid, weighted = weighted, wip = wip,
                     details = TRUE)
    c(d$pit, d$q_opt, d$run, match(d$tail, c("bulk", "right", "left")))
    # rep_j goes out of scope here -> memory stays O(N)
  }

  t0  <- Sys.time()
  res <- if (ncores > 1L)
    parallel::mclapply(seq_len(M), one, mc.cores = ncores, mc.preschedule = TRUE)
  else
    lapply(seq_len(M), one)

  bad <- !vapply(res, function(z) is.numeric(z) && length(z) == 4L, logical(1))
  if (any(bad)) {
    warning(sum(bad), " replicate(s) failed; NA returned for them", call. = FALSE)
    res[bad] <- list(rep(NA_real_, 4L))
  }

  R   <- matrix(unlist(res, use.names = FALSE), nrow = 4L)
  pit <- R[1L, ]
  attr(pit, "tail")  <- c("bulk", "right", "left")[R[4L, ]]
  attr(pit, "q_opt") <- R[2L, ]
  attr(pit, "run")   <- R[3L, ]

  if (verbose) {
    tb <- table(factor(attr(pit, "tail"), levels = c("bulk", "right", "left")))
    cat(sprintf("[sbc_pit] M = %d | %s | bulk %d, right %d, left %d | %.1fs\n",
                M, if (weighted) "weighted (dependence-aware)" else "naive",
                tb[["bulk"]], tb[["right"]], tb[["left"]],
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  pit
}


# ---------------------------------------------------------------------------
# make_ar1_generator()
#
# Generator for the null case: y and the chain share the same marginal, so the
# PIT would be exactly uniform if the draws were independent.  Any deviation
# is precisely the autocorrelation bias we are trying to remove.
# phi = 0 gives i.i.d. draws (useful as a sanity check).
# ---------------------------------------------------------------------------
make_ar1_generator <- function(N, phi, r_fun = rnorm, q_fun = qnorm) {
  force(N); force(phi); force(r_fun); force(q_fun)
  function(j) {
    z <- if (phi == 0) rnorm(N)
         else as.numeric(arima.sim(n = N, list(ar = phi), sd = sqrt(1 - phi^2)))
    list(y = r_fun(1L), x = q_fun(pnorm(z)))
  }
}


## ---- usage -----------------------------------------------------------------
gen <- make_ar1_generator(N = 2000, phi = 0.9)          # null case

pit_w <- sbc_pit(M = 5000, gen, ndraws_tail = NULL, weighted = TRUE,
                 seed = 221, ncores = 15)# dependence-aware

pit_n <- sbc_pit(M = 20000, gen, ndraws_tail = 1, weighted = FALSE,
                 seed = 221, ncores = 15)                 # naive

graphics.off()

x11(); hist(pit_w, breaks = 4000)
x11(); hist(pit_n, breaks = 4000)   # spikes at the edges

# same seed => same chains => the two arms are paired
uniformity_test(pit_w,"POT", TRUE)$pvalue

uniformity_test(pit_n,"PRIT", TRUE)$pvalue

table(attr(pit_w, "tail"))        # how often each branch fired
summary(attr(pit_w, "run"))       # declustering run lengths chosen

## system.time({
## dfs <- sbc_run(20000, gen, ndraws_tail = NULL, thin =16,
##                seed =11, ncores = 16)
## })

## x11()
## hist(dfs[,1],breaks=4001)

## uniformity_test(dfs[,3],"PRIT")$pvalue
# ===========================================================================
# COMPARING THREE SBC PROCEDURES UNDER AUTOCORRELATION
#
#   1. "usual"     normalized ranks on the full chain            (status quo)
#   2. "weighted"  normalized ranks + weighted-GPD tail fixing   (this work)
#   3. "thinned"   normalized ranks on a thinned chain           (common fix)
#
# All three are computed FROM THE SAME CHAIN in each replicate, so the
# comparison is exactly paired: any difference is attributable to the
# procedure, not to simulation noise.
#
# Two things are measured over repeated SBC studies:
#   calibration = P(reject | sampler is correct)   -> should equal the level
#   power       = P(reject | sampler is wrong)     -> should be large
# A procedure is only useful if it is calibrated FIRST; power comparisons
# between procedures with different type-I error rates are meaningless.
# ===========================================================================


# ---------------------------------------------------------------------------
# simmc_uniform() / make_markov_generator()
#
# First-order Markov chain whose consecutive pairs follow a bivariate extreme
# value distribution (Fawcett & Walshaw 2007), with UNIFORM margins, which are
# then pushed through q_fun to get the target marginal.  Unlike the Gaussian
# AR(1) copula this is asymptotically DEPENDENT for the logistic model, so
# clustering persists at arbitrarily high thresholds -- exactly the regime the
# declustering-based tail fix is designed for.
#
#   model = "log"   alpha in (0, 1]   alpha -> 0    total dependence
#                                     alpha -> 1    independence
#   model = "nlog"  alpha > 0         alpha -> 0    independence
#                                     alpha large   strong dependence
#
# (POT::simmc occasionally fails with a root-finding error; retry as in
# sim_markov_gpd() in weighted_gpdfit.R.)
# ---------------------------------------------------------------------------
simmc_uniform <- function(n, alpha, model = c("log", "nlog"),
                          max_try = 20L, ...) {
  model <- match.arg(model)
  stopifnot(requireNamespace("POT", quietly = TRUE))
  for (i in seq_len(max_try)) {
    uu <- tryCatch(
      POT::simmc(n, alpha = alpha, model = model, margins = "uniform", ...),
      error = function(e) {
        if (!grepl("opposite sign", conditionMessage(e))) stop(e)
        NULL
      })
    if (!is.null(uu) && all(is.finite(uu)) && all(uu > 0 & uu < 1))
      return(as.numeric(uu))
  }
  stop(sprintf("POT::simmc() failed %d consecutive times (model='%s', alpha=%g, n=%d)",
               max_try, model, alpha, n))
}

make_markov_generator <- function(N, alpha, model = c("log", "nlog"),
                                  r_fun = rnorm, q_fun = qnorm, ...) {
  model <- match.arg(model)
  force(N); force(alpha); force(r_fun); force(q_fun)
  dots <- list(...)
  function(j) {
    u <- do.call(simmc_uniform,
                 c(list(n = N, alpha = alpha, model = model), dots))
    list(y = r_fun(1L), x = q_fun(u))
  }
}


# ---------------------------------------------------------------------------
# miscalibrate()
#
# Wraps ANY generator to emulate a WRONG sampler, for power calculations:
# the chain is shifted and/or rescaled while y still comes from the correct
# prior.  These are the two classic SBC failure modes --
#   shift != 0  : biased posterior          -> PIT histogram tilts
#   scale  < 1  : over-confident posterior  -> PIT histogram is U-shaped
#   scale  > 1  : under-confident posterior -> PIT histogram is dome-shaped
# The dependence structure is untouched, so calibration and power are probed
# under identical autocorrelation.
# ---------------------------------------------------------------------------
miscalibrate <- function(generator, shift = 0, scale = 1) {
  force(generator); force(shift); force(scale)
  function(j) {
    rep_j   <- generator(j)
    rep_j$x <- shift + scale * rep_j$x
    rep_j
  }
}


# ---------------------------------------------------------------------------
# norm_rank()  --  the plain normalized rank r/(n+1), ties broken at random.
# ---------------------------------------------------------------------------
norm_rank <- function(y, x) {
  n      <- length(x)
  n_less <- sum(x < y)
  n_eq   <- sum(x == y)
  r      <- if (n_eq > 0L) n_less + sample.int(n_eq + 1L, 1L) else n_less + 1L
  r / (n + 1)
}


# ---------------------------------------------------------------------------
# act_time() / auto_thin()
#
# Integrated autocorrelation time  tau = 1 + 2 * sum_k rho_k, truncated at the
# first negative autocorrelation (Geyer's initial positive sequence).  Thinning
# by ceiling(tau) leaves draws that are approximately independent -- this is
# the standard recipe the "thinned" procedure implements.
# For AR(1), tau = (1+phi)/(1-phi), e.g. 39 at phi = 0.95.
# ---------------------------------------------------------------------------
act_time <- function(x, max_lag = NULL) {
  n <- length(x)
  if (is.null(max_lag)) max_lag <- min(n - 1L, max(50L, floor(n / 10)))
  a   <- stats::acf(x, lag.max = max_lag, plot = FALSE)$acf[-1L]
  neg <- which(a < 0)[1L]
  if (!is.na(neg)) a <- if (neg > 1L) a[seq_len(neg - 1L)] else numeric(0)
  max(1, 1 + 2 * sum(a))
}

auto_thin <- function(x) max(1L, as.integer(ceiling(act_time(x))))


# ---------------------------------------------------------------------------
# sbc_pit_all()
#
# All three PITs from ONE (y, chain) pair.  Returns a named numeric vector
#   usual, weighted, thinned, thin (lag used), n_thin (draws kept)
# ---------------------------------------------------------------------------
sbc_pit_all <- function(y, x, ndraws_tail = NULL, thin =16,
                        run = NULL, n_grid = 50, wip = TRUE) {
  n   <- length(x)
  k   <- if (identical(thin, "auto")) auto_thin(x) else max(1L, as.integer(thin))
  idx <- seq.int(1L, n, by = k)
  c(usual    = norm_rank(y, x),
    weighted = sbc_pit_one(y, x, run = run, ndraws_tail = ndraws_tail,
                           n_grid = n_grid, weighted = TRUE, wip = wip),
    thinned  = norm_rank(y, x[idx]),
    thin     = k,
    n_thin   = length(idx))
}


# ---------------------------------------------------------------------------
# sbc_run()
#
# ONE SBC study: M replicates, each contributing all three PITs.
# Returns an M x 5 matrix (usual, weighted, thinned, thin, n_thin).
# Memory stays O(N): chains are generated, used, and discarded.
#
# seed = NULL and ncores = 1 -> uses the ambient RNG and does NOT reseed,
# which is what lets sbc_study() drive many runs off its own streams.
# ---------------------------------------------------------------------------
sbc_run <- function(M, generator, ndraws_tail = NULL, thin =16,
                    run = NULL, n_grid = 50, wip = TRUE,
                    seed = NULL, ncores = 1L, verbose = FALSE) {

  stopifnot(is.function(generator), M >= 1L)

  if (is.null(ncores) || is.na(ncores))
    ncores <- max(1L, parallel::detectCores() - 1L)
  ncores <- max(1L, min(as.integer(ncores), M))
  if (ncores > 1L && .Platform$OS.type == "windows") ncores <- 1L

  use_streams <- !is.null(seed) || ncores > 1L
  if (use_streams) {
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
    old_kind <- RNGkind("L'Ecuyer-CMRG")
    on.exit(RNGkind(old_kind[1L]), add = TRUE)
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    streams <- vector("list", M)
    streams[[1L]] <- .Random.seed
    for (j in seq_len(M - 1L))
      streams[[j + 1L]] <- parallel::nextRNGStream(streams[[j]])
  }

  one <- function(j) {
    if (use_streams) assign(".Random.seed", streams[[j]], envir = globalenv())
    rep_j <- generator(j)
    sbc_pit_all(rep_j$y, rep_j$x, ndraws_tail = ndraws_tail, thin = thin,
                run = run, n_grid = n_grid, wip = wip)
  }

  t0  <- Sys.time()
  res <- if (ncores > 1L)
    parallel::mclapply(seq_len(M), one, mc.cores = ncores, mc.preschedule = TRUE)
  else
    lapply(seq_len(M), one)

  bad <- !vapply(res, function(z) is.numeric(z) && length(z) == 5L, logical(1))
  if (any(bad)) res[bad] <- list(setNames(rep(NA_real_, 5L),
                                          c("usual","weighted","thinned","thin","n_thin")))

  P <- do.call(rbind, res)
  if (verbose)
    cat(sprintf("[sbc_run] M = %d | thin ~ %.0f (keeps ~%.0f draws) | %.1fs\n",
                M, mean(P[, "thin"]), mean(P[, "n_thin"]),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  P
}


# ---------------------------------------------------------------------------
# sbc_study()
#
# Repeat the whole SBC study R times and report, per procedure, the rejection
# rate of uniformity_test(pit, test, truncate) at `level`.
# Requires dependence_aware_unif_tests.R to be sourced.
#
#   under a CORRECT sampler  -> that rate is the type-I error (calibration);
#                               it should be close to `level`
#   under miscalibrate(...)  -> that rate is the power
#
# Plain for-loop over the R studies.  ALL the parallelism lives inside
# sbc_run() (over the M replicates) -- one level only, so `ncores` means
# exactly what it says and the machine is never oversubscribed.
#
# (The earlier version nested mclapply over the R studies with a second
# mclapply inside sbc_run, forking up to ncores x ncores processes.  The loop
# is simpler, reports progress per study, and scales predictably.)
#
# Returns a data.frame with one row per procedure; the full R x 3 matrix of
# p-values is attached as attr(, "pvalues").
# ---------------------------------------------------------------------------
sbc_study <- function(R, M, generator, level = 0.05,
                      test = "PRIT", truncate = TRUE,
                      ndraws_tail = NULL, thin = 16,  #"auto",
                      run = NULL, n_grid = 50, wip = TRUE,
                      seed = 1L, ncores = 1L, label = "", verbose = TRUE) {

  methods <- c("usual", "weighted", "thinned")

  if (!exists("uniformity_test", mode = "function"))
    stop("uniformity_test() not found -- define/source it before calling sbc_study()")

  pval <- function(p) {
    p <- p[is.finite(p)]
    if (length(p) < 2L) return(NA_real_)
    as.numeric(tryCatch(uniformity_test(p, test, truncate)$pvalue,
                        error = function(e) NA_real_))
  }

  PV <- matrix(NA_real_, nrow = R, ncol = 5L,
               dimnames = list(NULL, c(methods, "thin", "n_thin")))

  t0   <- Sys.time()
  tick <- max(1L, floor(R / 10))

  for (i in seq_len(R)) {

    # each study gets its own reproducible per-replicate streams
    P <- sbc_run(M, generator, ndraws_tail = ndraws_tail, thin = thin,
                 run = run, n_grid = n_grid, wip = wip,
                 seed = seed + i - 1L, ncores = ncores)

    PV[i, methods]  <- vapply(methods, function(m) pval(P[, m]), 1.0)
    PV[i, "thin"]   <- mean(P[, "thin"])
    PV[i, "n_thin"] <- mean(P[, "n_thin"])

    if (verbose && (i %% tick == 0L || i == R)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      so_far <- colMeans(PV[seq_len(i), methods, drop = FALSE] < level,
                         na.rm = TRUE)
      cat(sprintf("  study %3d/%d | rej: %s | %.1f min elapsed, ~%.1f min left\n",
                  i, R,
                  paste(sprintf("%s %.3f", methods, so_far), collapse = "  "),
                  el / 60, el / 60 * (R - i) / i))
      utils::flush.console()
    }
  }

  rate <- colMeans(PV[, methods, drop = FALSE] < level, na.rm = TRUE)
  nOK  <- colSums(is.finite(PV[, methods, drop = FALSE]))

  out <- data.frame(
    label      = label,
    method     = methods,
    rej_rate   = as.numeric(rate),
    mc_se      = sqrt(rate * (1 - rate) / nOK),
    R          = R, M = M, level = level,
    mean_thin  = mean(PV[, "thin"],   na.rm = TRUE),
    mean_ndraw = mean(PV[, "n_thin"], na.rm = TRUE),
    row.names  = NULL
  )
  attr(out, "pvalues") <- PV[, methods, drop = FALSE]

  if (verbose) {
    cat(sprintf("\n--- %s | R = %d x M = %d | test = %s%s | level = %.2f | %.1f min ---\n",
                if (nzchar(label)) label else "sbc_study", R, M, test,
                if (isTRUE(truncate)) " (truncated)" else "", level,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    print(out[, c("method", "rej_rate", "mc_se")], digits = 3)
  }
  out
}


## ---- usage -----------------------------------------------------------------
## ## dependence structures
gen_ar1  <- make_ar1_generator(N = 1000, phi = 0.8)
gen_log  <- make_markov_generator(N = 500, alpha = 0.35, model = "log")
gen_nlog <- make_markov_generator(N = 2000, alpha = 5,   model = "nlog")

# uniformity_test() must be available in the session

## CALIBRATION (correct sampler): rejection rate should be ~ level
cal <- sbc_study(R = 200, M = 2000, gen_ar1, ndraws_tail = NULL,
                 test = "PRIT", truncate = FALSE,level=0.02, 
                 seed = 851, ncores = 15, label = "calibration | ar1")

## POWER (over-confident posterior): rejection rate should be large
pow <- sbc_study(R = 100, M = 500, miscalibrate(gen_log, scale = 0.9),
                 ndraws_tail = 100, test = "POT", truncate = TRUE,
                 seed = 1, ncores = 15, label = "power | scale = 0.9")

rbind(cal, pow)

## other tests, e.g.
## sbc_study(..., test = "PRIT", truncate = TRUE)

