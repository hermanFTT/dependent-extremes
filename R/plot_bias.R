#!/usr/bin/env Rscript
# =============================================================================
# plot_bias.R -- figures for the bias / RMSE / coverage sensitivity study.
#
# Implements the seven plot types of bias_sensitivity.md.  All three approaches
# are overlaid on the same axes (colour = approach); use --approaches= to
# compare a subset.
#
#   #  x-axis            y-axis      held fixed        reference line
#   1  dependence level  bias        observation size  0
#   2  observation size  bias        dependence level  0
#   3  dependence level  RMSE        observation size  --
#   4  observation size  RMSE        dependence level  --
#   5  dependence level  coverage    observation size  nominal level
#   6  observation size  coverage    dependence level  nominal level
#   7  true k            coverage    obs size, dep     nominal level
#
# "Observation size" is the number of exceedances above the threshold,
# n * (1 - p_thresh).
#
# The grid has more dimensions than a plot has axes.  Whatever is neither on
# an axis nor the reference line is either FACETED (the "held fixed" column of
# the table: one panel per value) or FILTERED to a single value.  Filters
# default to the middle value present in the data and are always printed in
# the subtitle, so no figure silently hides a dimension.
#
# Usage:
#   Rscript R/plot_bias.R --results=results/bias_sensitivity/bias_results.rds \
#                         --outdir=results/bias_sensitivity/figures
#   Rscript R/plot_bias.R ... --fix_k=0.2 --fix_r=auto \
#                         --approaches="Weighted GPD,Conventional"
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

`%||%` <- function(a, b) if (is.null(a)) b else a

opt <- list(
  results    = "results/bias_sensitivity/bias_results.rds",
  outdir     = "results/bias_sensitivity/figures",
  approaches = "",      # "" = all present
  fix_k      = "",      # "" = middle value in the data
  fix_n      = "",
  fix_dep    = "",
  fix_r      = "",
  width      = 8.2,
  height     = 5.0,
  dpi        = 130,
  device     = "pdf"
)

for (a in commandArgs(trailingOnly = TRUE)) {
  if (!grepl("^--[^=]+=", a)) next
  key <- sub("^--([^=]+)=.*$", "\\1", a)
  val <- sub("^--[^=]+=", "", a)
  if (!key %in% names(opt)) stop("unknown option --", key)
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(val) else val
}

res <- readRDS(opt$results)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# --- approach subset --------------------------------------------------------
if (nzchar(opt$approaches)) {
  keep <- trimws(strsplit(opt$approaches, ",", fixed = TRUE)[[1L]])
  bad  <- setdiff(keep, levels(res$approach))
  if (length(bad)) stop("unknown approach(es): ", paste(bad, collapse = ", "))
  res <- res[res$approach %in% keep, ]
  res$approach <- droplevels(res$approach)
}

# --- labels -----------------------------------------------------------------
DEP_LAB <- c(ar1  = "AR(1) dependence  phi",
             log  = "Markov logistic  alpha",
             nlog = "Markov negative-logistic  alpha")
PAR_LAB <- c(k     = "shape  k",
             sigma = "scale above threshold  sigma*",
             zm    = "return level  z_m")

dep_type  <- as.character(res$dep_type[1L])
dep_axis  <- DEP_LAB[[dep_type]] %||% "dependence parameter"
nominal   <- suppressWarnings(max(res$cred_mass, na.rm = TRUE))
if (!is.finite(nominal)) nominal <- NA_real_

COL <- c("Conventional"         = "#4DAC26",
         "Weighted GPD"         = "#2166AC",
         "Cluster maxima (POT)" = "#D6604D")

base_theme <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 12, face = "bold"),
        plot.subtitle    = element_text(size = 8.5, colour = "grey30"),
        legend.position  = "top")

save_fig <- function(p, name) {
  f <- file.path(opt$outdir, paste0(name, ".", opt$device))
  dev <- if (identical(opt$device, "pdf") && isTRUE(capabilities("cairo")))
           grDevices::cairo_pdf else NULL
  ggsave(f, p, device = dev, width = opt$width, height = opt$height,
         dpi = opt$dpi)
  cat("[plot] ", f, "\n", sep = "")
}

# --- pick the value a dimension is filtered to ------------------------------
# user-supplied if given, else the middle of the sorted unique values
pick <- function(v, user) {
  u <- sort(unique(v))
  if (!nzchar(user)) return(u[ceiling(length(u) / 2)])   # middle value
  if (is.numeric(v)) return(u[which.min(abs(u - as.numeric(user)))])  # nearest
  if (!user %in% u)
    stop("value '", user, "' not present; have: ", paste(u, collapse = ", "))
  user
}

fmt <- function(x) if (is.numeric(x)) formatC(x, format = "g", digits = 4) else as.character(x)

# --- one figure -------------------------------------------------------------
# x_var   : "dep_par" | "n_exc_nom" | "k_true"
# y_var   : "bias" | "rmse" | "coverage"
# facet_by: NULL or a column faceted over (the "held fixed" column)
# filters : named list of column -> value
make_plot <- function(d, x_var, y_var, facet_by, filters, title, hline = NULL) {

  for (nm in names(filters)) d <- d[d[[nm]] == filters[[nm]], ]
  d <- d[is.finite(d[[y_var]]), ]
  if (!nrow(d)) return(NULL)

  x_lab <- switch(x_var,
                  dep_par   = dep_axis,
                  n_exc_nom = "observation size   n (1 - u)   [expected exceedances]",
                  k_true    = "true shape  k")
  y_lab <- switch(y_var,
                  bias     = "bias",
                  rmse     = "RMSE",
                  coverage = sprintf("coverage probability (nominal %.2f)", nominal))

  sub <- paste0(
    sprintf("dependence: %s", dep_type),
    if (length(filters))
      paste0("   |   fixed: ",
             paste(sprintf("%s = %s", names(filters),
                           vapply(filters, fmt, character(1))),
                   collapse = ", ")) else "",
    sprintf("   |   S = %d", max(d$S)))

  p <- ggplot(d, aes(x = .data[[x_var]], y = .data[[y_var]],
                     colour = approach, shape = approach))
  if (!is.null(hline))
    p <- p + geom_hline(yintercept = hline, linetype = "dashed",
                        colour = "grey40")
  p <- p +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.7) +
    scale_colour_manual(values = COL, drop = TRUE) +
    labs(title = title, subtitle = sub, x = x_lab, y = y_lab,
         colour = NULL, shape = NULL) +
    base_theme

  if (!is.null(facet_by) && length(unique(d[[facet_by]])) > 1L) {
    lab <- switch(facet_by,
                  n_exc_nom = function(v) paste0("obs size = ", v),
                  dep_par   = function(v) paste0(dep_type, " = ", v),
                  k_true    = function(v) paste0("k = ", v),
                  function(v) as.character(v))
    p <- p + facet_wrap(stats::as.formula(paste("~", facet_by)),
                        scales = "free_y", labeller = labeller(.default = lab))
  }
  p
}

# --- what varies in this run? ----------------------------------------------
n_dep <- length(unique(res$dep_par))
n_obs <- length(unique(res$n))
n_k   <- length(unique(res$k_true))
n_r   <- length(unique(res$r_label))

cat(sprintf("[plot] grid: dep_par=%d  n=%d  k=%d  r=%d  approaches=%d\n",
            n_dep, n_obs, n_k, n_r, nlevels(droplevels(res$approach))))

fix_k <- pick(res$k_true,  opt$fix_k)
fix_n <- pick(res$n,       opt$fix_n)
fix_d <- pick(res$dep_par, opt$fix_dep)
fix_r <- pick(res$r_label, opt$fix_r)

cat(sprintf("[plot] filters: k=%s  n=%s  dep=%s  r=%s\n",
            fmt(fix_k), fmt(fix_n), fmt(fix_d), fmt(fix_r)))

n_written <- 0L
emit <- function(p, name) {
  if (is.null(p)) { cat("[plot] skip (no data): ", name, "\n", sep = ""); return(invisible()) }
  save_fig(p, name); n_written <<- n_written + 1L
}

for (par in levels(droplevels(res$parameter))) {

  d   <- res[res$parameter == par, ]
  lab <- PAR_LAB[[par]]

  # --- 1 & 3: x = dependence level, faceted by observation size ------------
  emit(make_plot(d, "dep_par", "bias", "n_exc_nom",
                 list(k_true = fix_k, r_label = fix_r),
                 paste0("Bias vs dependence - ", lab), hline = 0),
       paste0("fig1_bias_vs_dep_", par))

  emit(make_plot(d, "dep_par", "rmse", "n_exc_nom",
                 list(k_true = fix_k, r_label = fix_r),
                 paste0("RMSE vs dependence - ", lab)),
       paste0("fig3_rmse_vs_dep_", par))

  # --- 2 & 4: x = observation size, faceted by dependence level ------------
  emit(make_plot(d, "n_exc_nom", "bias", "dep_par",
                 list(k_true = fix_k, r_label = fix_r),
                 paste0("Bias vs observation size - ", lab), hline = 0),
       paste0("fig2_bias_vs_obs_", par))

  emit(make_plot(d, "n_exc_nom", "rmse", "dep_par",
                 list(k_true = fix_k, r_label = fix_r),
                 paste0("RMSE vs observation size - ", lab)),
       paste0("fig4_rmse_vs_obs_", par))

  # --- 5, 6, 7: coverage (k and z_m only; sigma has no posterior here) -----
  if (par != "sigma" && any(is.finite(d$coverage))) {

    emit(make_plot(d, "dep_par", "coverage", "n_exc_nom",
                   list(k_true = fix_k, r_label = fix_r),
                   paste0("Coverage vs dependence - ", lab), hline = nominal),
         paste0("fig5_cov_vs_dep_", par))

    emit(make_plot(d, "n_exc_nom", "coverage", "dep_par",
                   list(k_true = fix_k, r_label = fix_r),
                   paste0("Coverage vs observation size - ", lab), hline = nominal),
         paste0("fig6_cov_vs_obs_", par))

    emit(make_plot(d, "k_true", "coverage", "n_exc_nom",
                   list(dep_par = fix_d, r_label = fix_r),
                   paste0("Coverage vs true k - ", lab), hline = nominal),
         paste0("fig7_cov_vs_k_", par))
  }
}

# --- run-length comparison, if more than one r was simulated ---------------
if (n_r > 1L) {
  d <- res[res$parameter == "k", ]
  d <- d[d$k_true == fix_k & d$dep_par == fix_d, ]
  d <- d[is.finite(d$rmse), ]
  if (nrow(d))
    emit(ggplot(d, aes(n_exc_nom, rmse, colour = approach, shape = approach)) +
           geom_line(linewidth = 0.8) + geom_point(size = 1.7) +
           scale_colour_manual(values = COL, drop = TRUE) +
           facet_wrap(~ r_label, labeller = labeller(.default = function(v)
             paste0("run length = ", v))) +
           labs(title = "RMSE of k by run length",
                subtitle = sprintf("dependence: %s = %s   |   k = %s   |   S = %d",
                                   dep_type, fmt(fix_d), fmt(fix_k), max(d$S)),
                x = "observation size   n (1 - u)", y = "RMSE",
                colour = NULL, shape = NULL) +
           base_theme,
         "fig8_rmse_by_runlength_k")
}

cat("[plot] wrote ", n_written, " figures -> ", opt$outdir, "\n", sep = "")
