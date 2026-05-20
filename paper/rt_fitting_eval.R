# RT-correction method bake-off evaluation script

# ---- Locate project root ------------------------------------------------
PROJ_ROOT <- "~/code/CEMBIO-EURAC/"
STUDY_DIR <- file.path(PROJ_ROOT, "applications", "tutorial")

# ---- Config ---------------------------------------------------------------
PARALLEL <- TRUE # TRUE = mclapply (one process per method); FALSE = lapply
# Set FALSE for fair per-method timing comparisons.
N_CORES <- min(20, max(1L, parallel::detectCores() - 1L))
N_BOOT <- 500 # bootstrap resamples for non-parametric PI methods
LOG_FILE <- file.path(PROJ_ROOT, "paper", "objects", "rtfit_eval.log")
PLOT_LOO <- TRUE # TRUE = write per-fold PNGs to LOO_FIGURES_DIR
LOO_FIGURES_DIR <- file.path(PROJ_ROOT, "paper", "loo_figures")

# ---- Packages -------------------------------------------------------------
suppressPackageStartupMessages({
  library(xcms)
  library(mgcv)
  library(scam)
  library(parallel)
  library(ggplot2)
})

source(file.path(PROJ_ROOT, "R/lipid_helpers.R"))

# =========================================================================
# 0. Logging helper
# =========================================================================
# Writes timestamped entries to LOG_FILE and echoes to console.
log_msg <- function(level, ...) {
  msg <- sprintf(
    "[%s] [%-5s] %s",
    format(Sys.time(), "%H:%M:%S"),
    level,
    paste(..., sep = " ")
  )
  message(msg)
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

# =========================================================================
# 1. load_artifacts
# =========================================================================
load_artifacts <- function(polarity) {
  dir_name <- switch(polarity, pos = "positive", neg = "negative", polarity)
  path <- file.path(
    STUDY_DIR,
    dir_name,
    "objects",
    "rt_fit_artifacts.rds"
  )
  if (!file.exists(path)) {
    stop(
      "Artifact not found: ",
      path,
      "\nRender Annotation_",
      polarity,
      ".qmd first."
    )
  }
  readRDS(path)
}

# =========================================================================
# 2. zone_of
# =========================================================================
zone_of <- function(rt_sec) {
  cut(
    rt_sec,
    breaks = c(-Inf, 210, 600, 660, Inf),
    labels = c("G1", "I2", "G2", "I3"),
    right = TRUE,
    include.lowest = TRUE
  )
}

# =========================================================================
# 3. build_anchor_table
# =========================================================================
build_anchor_table <- function(artifacts, ppm = 10) {
  eic_is <- artifacts$eic_is
  intern_standard <- artifacts$intern_standard

  group <- eic_is$sample_type
  group[group == "QC"] <- NA
  group <- as.factor(group)
  param_group <- PeakDensityParam(
    sampleGroups = group,
    minFraction = 0.5,
    binSize = 0.01,
    ppm = ppm,
    bw = 2
  )

  # Delegate anchor extraction; discard the fit itself.
  rt_fit <- fit_rt_correction(
    eic_is,
    intern_standard,
    param_group,
    poly_degree = 6
  )

  anchors <- data.frame(
    short_name = intern_standard$short_name,
    ref_rt = rt_fit$ref_rt,
    exp_rt = rt_fit$exp_rt,
    stringsAsFactors = FALSE
  )
  anchors$zone <- zone_of(anchors$ref_rt)
  anchors[!is.na(anchors$exp_rt), ]
}

# =========================================================================
# 4. validate_anchors
# =========================================================================
validate_anchors <- function(anchors) {
  ref_fit <- loess(
    exp_rt ~ ref_rt,
    data = anchors,
    span = 0.75,
    family = "symmetric"
  )
  anchors$loess_resid <- residuals(ref_fit)
  flagged <- anchors[abs(anchors$loess_resid) > 30, ]
  if (nrow(flagged)) {
    message(
      "Flagged (|LOESS resid| > 30 s): ",
      paste(flagged$short_name, collapse = ", ")
    )
  } else {
    message("No anchors flagged by robust LOESS.")
  }
  invisible(flagged)
}

# =========================================================================
# 5. poly_method  (parametric PI via lm predict)
# =========================================================================
poly_method <- function(anchors, degree = 6) {
  fit <- lm(exp_rt ~ poly(ref_rt, degree, raw = TRUE), data = anchors)
  function(ref_rt_vec, alpha = 0.05) {
    p <- predict(
      fit,
      newdata = data.frame(ref_rt = ref_rt_vec),
      interval = "prediction",
      level = 1 - alpha
    )
    list(point = p[, "fit"], lower = p[, "lwr"], upper = p[, "upr"])
  }
}

# =========================================================================
# 6. with_bootstrap_pi  (case-bootstrap PI decorator)
# =========================================================================
# Wraps a point-estimate factory with case-bootstrap prediction intervals.
# Degenerate resamples (factory errors) are caught and contribute NA;
# a warning fires if > 5% fail.
with_bootstrap_pi <- function(factory, n_boot = N_BOOT) {
  function(anchors) {
    pred_fn_full <- factory(anchors)
    n_anch <- nrow(anchors)
    function(ref_rt_vec, alpha = 0.05) {
      pt <- pred_fn_full(ref_rt_vec)$point
      k <- length(ref_rt_vec)
      boot_mat <- matrix(
        vapply(
          seq_len(n_boot),
          function(.) {
            idx <- sample.int(n_anch, replace = TRUE)
            tryCatch(
              suppressWarnings(factory(anchors[idx, ])(ref_rt_vec)$point),
              error = function(e) rep(NA_real_, k)
            )
          },
          numeric(k)
        ),
        nrow = k,
        ncol = n_boot
      )
      n_failed <- sum(is.na(boot_mat[1L, ]))
      if (n_failed > 0.05 * n_boot) {
        warning(sprintf(
          "%.0f%% of bootstrap resamples failed; PI may be unreliable.",
          100 * n_failed / n_boot
        ))
      }
      list(
        point = pt,
        lower = apply(boot_mat, 1L, quantile, probs = alpha / 2, na.rm = TRUE),
        upper = apply(
          boot_mat,
          1L,
          quantile,
          probs = 1 - alpha / 2,
          na.rm = TRUE
        )
      )
    }
  }
}

# =========================================================================
# 7. Method factories
# =========================================================================

# --- GAM  (posterior predictive PI — Wood 2017 §6.10; Wood 2020 TEST) --------
# Two-stage posterior simulation: (1) draw n_sim coefficient vectors from
# N(coef, Vp) to propagate parameter uncertainty; (2) add N(0, sig^2)
# observation noise.  For Gaussian/identity this is analytically equivalent to
# ŷ ± z * sqrt(se.fit^2 + sig2), but simulation is used because it is the
# approach endorsed by Wood (2017, 2020) and avoids conflating confidence with
# prediction intervals (Marra & Wood 2012, Scand. J. Statist. 39:53-74).
# vcov(fit) returns fit$Vp (Bayesian posterior covariance) for mgcv gam
# objects; coef(fit) is consistent with predict(type = "lpmatrix").
gam_method <- function(anchors, n_sim = 2000L) {
  fit <- mgcv::gam(exp_rt ~ s(ref_rt), data = anchors, method = "REML")
  beta <- coef(fit)
  V <- vcov(fit) # == fit$Vp: Bayesian posterior covariance
  Cv <- chol(V)
  sig <- sqrt(fit$sig2)
  nb <- length(beta)
  function(ref_rt_vec, alpha = 0.05) {
    Xp <- predict(
      fit,
      newdata = data.frame(ref_rt = ref_rt_vec),
      type = "lpmatrix"
    )
    beta_sim <- beta + t(Cv) %*% matrix(rnorm(nb * n_sim), nb, n_sim)
    mu_sim <- Xp %*% beta_sim
    y_sim <- mu_sim +
      matrix(rnorm(length(mu_sim), 0, sig), nrow(mu_sim), ncol(mu_sim))
    pt <- as.numeric(Xp %*% beta)
    list(
      point = pt,
      lower = apply(y_sim, 1L, quantile, probs = alpha / 2),
      upper = apply(y_sim, 1L, quantile, probs = 1 - alpha / 2)
    )
  }
}

# --- Monotone GAM via scam  (Pya & Wood 2015; PI as for gam_method above) ---
# bs="mpi": monotone increasing P-splines.  Posterior predictive PI uses the
# same two-stage simulation as gam_method.  scam internally reparameterises
# coefficients to enforce monotonicity: predict(type = "lpmatrix") operates in
# the transformed (constrained) space, so fit$coefficients.t and fit$Vp.t
# must be used — coef(fit)/vcov(fit) live in the unconstrained space and give
# wrong predictions when multiplied against the lpmatrix.
scam_method <- function(anchors, n_sim = 2000L) {
  fit <- scam::scam(exp_rt ~ s(ref_rt, bs = "mpi"), data = anchors)
  beta <- fit$coefficients.t # constrained parameterisation (matches lpmatrix)
  V <- fit$Vp.t # Bayesian posterior covariance in that space
  Cv <- chol(V)
  sig <- sqrt(fit$sig2)
  nb <- length(beta)
  function(ref_rt_vec, alpha = 0.05) {
    Xp <- predict(
      fit,
      newdata = data.frame(ref_rt = ref_rt_vec),
      type = "lpmatrix"
    )
    beta_sim <- beta + t(Cv) %*% matrix(rnorm(nb * n_sim), nb, n_sim)
    mu_sim <- Xp %*% beta_sim
    y_sim <- mu_sim +
      matrix(rnorm(length(mu_sim), 0, sig), nrow(mu_sim), ncol(mu_sim))
    pt <- as.numeric(Xp %*% beta)
    list(
      point = pt,
      lower = apply(y_sim, 1L, quantile, probs = alpha / 2),
      upper = apply(y_sim, 1L, quantile, probs = 1 - alpha / 2)
    )
  }
}

# --- PredRet monotone GAM (Stanstrup et al. 2015) ----------------------------
# Faithful reimplementation: initial tp GAM for sigmoid-robust weights +
# smoothing parameter; constrained cr spline via mgcv::pcls + mono.con.
# Initial parameters set to sm$xp (knot positions) which are always feasible.
# Sigmoid down-weights anchors with |resid| > ~10% of max(exp_rt).
# PI via with_bootstrap_pi decorator.
predret_mono_method_raw <- function(anchors) {
  x <- anchors$ref_rt
  y <- anchors$exp_rt
  k_val <- max(4L, min(length(unique(round(x, 2L))), 10L))
  dat <- data.frame(x = x, y = y)

  sigmoid_w <- function(w, a = -30, b = 0.1) 1 / (1 + exp(-a * (w - b)))

  f_init <- mgcv::gam(y ~ s(x, k = k_val, bs = "tp"), data = dat)
  w <- sigmoid_w(abs(residuals(f_init)) / max(y))

  sm <- mgcv::smoothCon(
    mgcv::s(x, k = k_val, bs = "cr"),
    data = dat,
    knots = NULL
  )[[1]]
  con <- mgcv::mono.con(sm$xp, up = TRUE)

  # sm$xp (knot positions) are monotonically increasing x-values, so they
  # always satisfy the monotonicity constraints — a faithful starting point.
  G <- list(
    X = sm$X,
    C = matrix(0, 0, 0),
    sp = f_init$sp,
    p = sm$xp,
    y = y,
    w = w,
    S = sm$S,
    Ain = con$A,
    bin = con$b,
    off = 0
  )
  p_coef <- mgcv::pcls(G)

  # Fallback: retry with uniform weights if pcls returned NAs (Stanstrup 2015)
  if (any(is.na(p_coef))) {
    G$w <- rep(1, length(G$w))
    p_coef <- mgcv::pcls(G)
  }

  function(ref_rt_vec, alpha = 0.05) {
    fv <- mgcv::Predict.matrix(sm, data.frame(x = ref_rt_vec)) %*% p_coef
    list(point = pmax(as.numeric(fv), 0), lower = NA_real_, upper = NA_real_)
  }
}

# =========================================================================
# 8. CV generators
# =========================================================================
loocv <- function(anchors) {
  n <- nrow(anchors)
  lapply(seq_len(n), function(i) {
    list(train_idx = seq_len(n)[-i], test_idx = i)
  })
}

# lozo kept for reference but not included in the default cv_schemes.
lozo <- function(anchors) {
  folds <- list()
  for (z in levels(anchors$zone)) {
    test_idx <- which(anchors$zone == z)
    if (length(test_idx) < 2) {
      warning(paste("Skipping zone", z, "in LOZO: < 2 anchors."))
      next
    }
    folds <- c(
      folds,
      list(list(
        train_idx = which(anchors$zone != z),
        test_idx = test_idx
      ))
    )
  }
  folds
}

random_k_out <- function(anchors, frac, n_repeats = 100) {
  n <- nrow(anchors)
  k <- max(1L, round(frac * n))
  lapply(seq_len(n_repeats), function(.) {
    idx <- sample.int(n, k)
    list(train_idx = seq_len(n)[-idx], test_idx = idx)
  })
}

stratified_random_k_out <- function(
  anchors,
  frac,
  n_repeats = 100,
  strat_by = c("zone", "rt_bins"),
  n_classes = 5,
  bin_type = c("quantile", "equal_width")
) {
  strat_by <- match.arg(strat_by)
  bin_type <- match.arg(bin_type)
  n <- nrow(anchors)

  strata <- if (strat_by == "zone") {
    anchors$zone
  } else {
    breaks <- if (bin_type == "quantile") {
      quantile(anchors$ref_rt, probs = seq(0, 1, length.out = n_classes + 1))
    } else {
      n_classes
    }
    cut(anchors$ref_rt, breaks = breaks, include.lowest = TRUE)
  }

  all_levels <- levels(strata)
  strata <- as.character(strata)
  for (i in seq_along(all_levels)) {
    lv <- all_levels[i]
    if (sum(strata == lv, na.rm = TRUE) < 2L) {
      target <- if (i > 1L) all_levels[i - 1L] else all_levels[i + 1L]
      strata[strata == lv] <- target
      warning(sprintf(
        "Stratum '%s' < 2 anchors; merged into '%s'.",
        lv,
        target
      ))
    }
  }
  stratum_levels <- unique(strata)

  lapply(seq_len(n_repeats), function(.) {
    test_idx <- integer(0)
    for (lv in stratum_levels) {
      s_idx <- which(strata == lv)
      n_s <- length(s_idx)
      k_s <- min(max(1L, round(frac * n_s)), n_s - 1L)
      test_idx <- c(test_idx, sample(s_idx, k_s))
    }
    list(train_idx = setdiff(seq_len(n), test_idx), test_idx = test_idx)
  })
}

# =========================================================================
# 9. Gap-interpolation diagnostic
# =========================================================================
gap_diagnostic <- function(anchors, method_fn, zone_label = "G2") {
  zone_bounds <- list(
    G1 = c(0, 210),
    I2 = c(210, 600),
    G2 = c(600, 660),
    I3 = c(660, 950)
  )
  b <- zone_bounds[[zone_label]]
  pred_fn <- method_fn(anchors)
  grid <- seq(b[1], b[2], by = 1)
  preds <- pred_fn(grid)
  data.frame(
    ref_rt = grid,
    point = preds$point,
    lower = preds$lower,
    upper = preds$upper,
    pi_width = preds$upper - preds$lower,
    spread = 0 # filled by run_eval after all methods are computed
  )
}

# =========================================================================
# 10. LOO diagnostic plots  (called only when PLOT_LOO = TRUE)
# =========================================================================
plot_loo_fits <- function(
  mname,
  method_fn,
  anchors,
  loocv_splits,
  out_dir,
  polarity
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  rt_grid <- seq(
    min(anchors$ref_rt) * 0.98,
    max(anchors$ref_rt) * 1.02,
    length.out = 300
  )
  zone_bands <- data.frame(
    xmin = c(-Inf, 210, 600, 660),
    xmax = c(210, 600, 660, Inf),
    zone = c("G1", "I2", "G2", "I3"),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(loocv_splits)) {
    sp <- loocv_splits[[i]]
    train <- anchors[sp$train_idx, ]
    held <- anchors[sp$test_idx, ]

    pred_fn <- tryCatch(method_fn(train), error = function(e) NULL)
    if (is.null(pred_fn)) {
      next
    }

    curve <- tryCatch(pred_fn(rt_grid), error = function(e) NULL)
    hpred <- tryCatch(pred_fn(held$ref_rt), error = function(e) NULL)
    if (is.null(curve) || is.null(hpred)) {
      next
    }

    curve_df <- data.frame(ref_rt = rt_grid, pred = curve$point)
    has_pi <- !any(is.na(curve$lower))
    if (has_pi) {
      curve_df$lower <- curve$lower
      curve_df$upper <- curve$upper
    }

    resid_val <- held$exp_rt - hpred$point

    p <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = zone_bands,
        ggplot2::aes(
          xmin = xmin,
          xmax = xmax,
          ymin = -Inf,
          ymax = Inf,
          fill = zone
        ),
        alpha = 0.07,
        show.legend = FALSE
      ) +
      {
        if (has_pi) {
          ggplot2::geom_ribbon(
            data = curve_df,
            ggplot2::aes(x = ref_rt, ymin = lower, ymax = upper),
            alpha = 0.15,
            fill = "steelblue"
          )
        } else {
          list()
        }
      } +
      ggplot2::geom_line(
        data = curve_df,
        ggplot2::aes(x = ref_rt, y = pred),
        colour = "steelblue",
        linewidth = 0.8
      ) +
      ggplot2::geom_point(
        data = train,
        ggplot2::aes(x = ref_rt, y = exp_rt),
        shape = 16,
        size = 2.5,
        colour = "grey30"
      ) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = held$ref_rt,
          xend = held$ref_rt,
          y = hpred$point,
          yend = held$exp_rt
        ),
        colour = "firebrick",
        linewidth = 0.9,
        linetype = "dashed"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = held$ref_rt, y = hpred$point),
        shape = 24,
        size = 3,
        colour = "firebrick",
        fill = "white"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = held$ref_rt, y = held$exp_rt),
        shape = 4,
        size = 4,
        colour = "firebrick",
        stroke = 1.5
      ) +
      ggplot2::annotate(
        "text",
        x = held$ref_rt,
        y = (held$exp_rt + hpred$point) / 2,
        label = sprintf("%.1f s", resid_val),
        hjust = -0.15,
        size = 3.2,
        colour = "firebrick"
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          G1 = "#4DAF4A",
          I2 = "#FF7F00",
          G2 = "#E41A1C",
          I3 = "#984EA3"
        )
      ) +
      ggplot2::labs(
        title = sprintf(
          "%s  |  %s  |  fold %d/%d",
          toupper(polarity),
          mname,
          i,
          length(loocv_splits)
        ),
        subtitle = sprintf(
          "held-out: %s   resid = %.1f s",
          held$short_name,
          resid_val
        ),
        x = "Reference RT (s)",
        y = "Experimental RT (s)"
      ) +
      ggplot2::theme_bw(base_size = 11)

    ggplot2::ggsave(
      file.path(out_dir, paste0(i, ".png")),
      plot = p,
      width = 7,
      height = 5,
      dpi = 120
    )
  }
}

# =========================================================================
# 11. Metrics  (formerly 10)
# =========================================================================
eval_splits <- function(anchors, method_fn, splits) {
  rows <- lapply(splits, function(sp) {
    train <- anchors[sp$train_idx, ]
    test <- anchors[sp$test_idx, ]
    pred_fn <- method_fn(train)
    preds <- pred_fn(test$ref_rt)
    resid <- test$exp_rt - preds$point
    covered <- (test$exp_rt >= preds$lower) & (test$exp_rt <= preds$upper)
    data.frame(
      short_name = test$short_name,
      zone = as.character(test$zone),
      ref_rt = test$ref_rt,
      exp_rt = test$exp_rt,
      pred = preds$point,
      lower = preds$lower,
      upper = preds$upper,
      resid = resid,
      covered = covered,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

compute_metrics <- function(rd) {
  one <- function(df) {
    data.frame(
      n = nrow(df),
      rmse = sqrt(mean(df$resid^2, na.rm = TRUE)),
      mae = mean(abs(df$resid), na.rm = TRUE),
      max_ae = max(abs(df$resid), na.rm = TRUE),
      coverage = mean(df$covered, na.rm = TRUE)
    )
  }
  global <- cbind(zone = "ALL", one(rd))
  per_zone <- do.call(
    rbind,
    lapply(
      names(split(rd, rd$zone)),
      function(z) cbind(zone = z, one(split(rd, rd$zone)[[z]]))
    )
  )
  rbind(global, per_zone)
}

monotonicity_violations <- function(anchors, method_fn) {
  ord <- order(anchors$ref_rt)
  preds <- method_fn(anchors)(anchors$ref_rt[ord])$point
  diffs <- diff(preds)
  viol <- which(diffs < 0)
  list(unweighted = length(viol), magnitude_weighted = sum(-diffs[viol]))
}

# =========================================================================
# 11. eval_one_method  (runs one method across all CV schemes)
# =========================================================================
# Designed to be called from lapply or mclapply.  Returns a list with:
#   $resids  — list of residual data frames, keyed by "method:cv_type"
#   $metrics — list of metrics data frames (includes time_secs column)
#   $log     — character vector of log entries (written by caller)
eval_one_method <- function(mname, method_fn, anchors, cv_schemes, polarity) {
  out_resids <- list()
  out_metrics <- list()
  log_entries <- character(0)

  for (cvname in names(cv_schemes)) {
    key <- paste(mname, cvname, sep = ":")
    t0 <- proc.time()[["elapsed"]]

    rd <- withCallingHandlers(
      tryCatch(
        eval_splits(anchors, method_fn, cv_schemes[[cvname]]),
        error = function(e) {
          log_entries <<- c(
            log_entries,
            sprintf(
              "ERROR | %s | %s | %s | %s",
              polarity,
              mname,
              cvname,
              conditionMessage(e)
            )
          )
          NULL
        }
      ),
      warning = function(w) {
        log_entries <<- c(
          log_entries,
          sprintf(
            "WARN  | %s | %s | %s | %s",
            polarity,
            mname,
            cvname,
            conditionMessage(w)
          )
        )
        invokeRestart("muffleWarning")
      }
    )

    elapsed <- proc.time()[["elapsed"]] - t0

    if (!is.null(rd)) {
      rd$method <- mname
      rd$cv_type <- cvname
      out_resids[[key]] <- rd
      m <- compute_metrics(rd)
      m$method <- mname
      m$cv_type <- cvname
      m$time_secs <- elapsed
      out_metrics[[key]] <- m
    }
  }

  list(resids = out_resids, metrics = out_metrics, log = log_entries)
}

# =========================================================================
# 12. run_eval  (consolidated; replaces former run_eval + run_eval_stage3)
# =========================================================================
run_eval <- function(polarity) {
  # Initialise log (truncate any previous run for this session).
  dir.create(
    file.path(PROJ_ROOT, "paper", "objects"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  if (!file.exists(LOG_FILE)) {
    cat("", file = LOG_FILE)
  }

  log_msg("INFO", sprintf("=== RT bake-off: %s ===", toupper(polarity)))

  artifacts <- load_artifacts(polarity)
  anchors <- suppressMessages(build_anchor_table(artifacts))
  log_msg(
    "INFO",
    polarity,
    "anchors per zone:",
    paste(
      names(table(anchors$zone)),
      table(anchors$zone),
      sep = "=",
      collapse = "  "
    )
  )

  flagged <- validate_anchors(anchors)
  if (nrow(flagged) > 0) {
    log_msg(
      "WARN ",
      polarity,
      "flagged anchors:",
      paste(flagged$short_name, collapse = ", ")
    )
  }

  # ------ Method definitions -----------------------------------------------
  poly_factories <- setNames(
    lapply(3:6, function(d) {
      force(d)
      function(a) poly_method(a, degree = d)
    }),
    paste0("poly_d", 3:6)
  )

  methods <- c(
    poly_factories,
    list(
      gam = gam_method,
      scam = scam_method
      # predret_mono = with_bootstrap_pi(predret_mono_method_raw)
    )
  )

  # ------ CV scheme definitions (LOOCV + RKO + strat_bins) ----------------
  cv_schemes <- c(
    list(loocv = loocv(anchors)),
    setNames(
      lapply(c(0.05, 0.10, 0.20, 0.30), function(f) {
        random_k_out(anchors, frac = f, n_repeats = 100)
      }),
      c("rko_05pct", "rko_10pct", "rko_20pct", "rko_30pct")
    ),
    setNames(
      lapply(c(0.05, 0.10, 0.20, 0.30), function(f) {
        stratified_random_k_out(
          anchors,
          frac = f,
          n_repeats = 100,
          strat_by = "rt_bins"
        )
      }),
      c(
        "strat_bins_05pct",
        "strat_bins_10pct",
        "strat_bins_20pct",
        "strat_bins_30pct"
      )
    )
  )

  log_msg(
    "INFO",
    polarity,
    sprintf(
      "%d methods × %d CV schemes = %d combinations",
      length(methods),
      length(cv_schemes),
      length(methods) * length(cv_schemes)
    )
  )

  # ------ Run evaluations (parallel or serial) -----------------------------
  runner <- if (PARALLEL) {
    log_msg("INFO", polarity, sprintf("parallel mode: %d cores", N_CORES))
    function(nms, FUN) parallel::mclapply(nms, FUN, mc.cores = N_CORES)
  } else {
    function(nms, FUN) lapply(nms, FUN)
  }

  t_total <- proc.time()[["elapsed"]]

  results <- runner(names(methods), function(mname) {
    log_msg("INFO", polarity, mname, "starting ...")
    eval_one_method(mname, methods[[mname]], anchors, cv_schemes, polarity)
  })
  names(results) <- names(methods)

  # Write log entries collected inside child processes / lapply iterations.
  for (res in results) {
    if (length(res$log) > 0) {
      cat(paste(res$log, collapse = "\n"), "\n", file = LOG_FILE, append = TRUE)
    }
  }

  elapsed_total <- proc.time()[["elapsed"]] - t_total
  log_msg(
    "INFO",
    polarity,
    sprintf("all evaluations done in %.1f s", elapsed_total)
  )

  # ------ Optional LOO diagnostic plots ------------------------------------
  if (PLOT_LOO) {
    log_msg("INFO", polarity, "Writing LOO diagnostic plots ...")
    for (mname in names(methods)) {
      out_dir <- file.path(LOO_FIGURES_DIR, polarity, mname)
      tryCatch(
        plot_loo_fits(
          mname,
          methods[[mname]],
          anchors,
          cv_schemes$loocv,
          out_dir,
          polarity
        ),
        error = function(e) {
          log_msg(
            "ERROR",
            polarity,
            mname,
            "loo_plot:",
            conditionMessage(e)
          )
        }
      )
    }
    log_msg("INFO", polarity, "LOO plots saved to:", LOO_FIGURES_DIR)
  }

  all_resids <- do.call(c, lapply(results, `[[`, "resids"))
  all_metrics <- do.call(c, lapply(results, `[[`, "metrics"))

  # ------ Monotonicity violations ------------------------------------------
  mono <- lapply(setNames(names(methods), names(methods)), function(mn) {
    tryCatch(
      monotonicity_violations(anchors, methods[[mn]]),
      error = function(e) {
        log_msg("ERROR", polarity, mn, "mono:", conditionMessage(e))
        list(unweighted = NA_integer_, magnitude_weighted = NA_real_)
      }
    )
  })

  # ------ Gap diagnostic (all methods, all empty zones) --------------------
  empty_zones <- setdiff(
    levels(anchors$zone),
    as.character(unique(anchors$zone))
  )
  gap_diag <- if (length(empty_zones) > 0) {
    lapply(setNames(empty_zones, empty_zones), function(z) {
      method_diags <- lapply(
        setNames(names(methods), names(methods)),
        function(mn) {
          tryCatch(
            gap_diagnostic(anchors, methods[[mn]], zone_label = z),
            error = function(e) {
              log_msg(
                "ERROR",
                polarity,
                mn,
                paste0("gap_diag_", z),
                conditionMessage(e)
              )
              NULL
            }
          )
        }
      )
      method_diags <- Filter(Negate(is.null), method_diags)
      if (length(method_diags) == 0L) {
        return(NULL)
      }

      pts <- do.call(cbind, lapply(method_diags, `[[`, "point"))
      spread_vec <- if (ncol(pts) > 1L) {
        apply(pts, 1L, sd, na.rm = TRUE)
      } else {
        rep(0, nrow(pts))
      }

      lapply(names(method_diags), function(mn) {
        d <- method_diags[[mn]]
        d$spread <- spread_vec
        d$method <- mn
        d
      }) |>
        setNames(names(method_diags))
    })
  } else {
    NULL
  }

  # ------ Assemble & save output -------------------------------------------
  resid_df <- do.call(rbind, all_resids)
  metrics_df <- do.call(rbind, all_metrics)

  output <- list(
    polarity = polarity,
    anchors = anchors,
    flagged = flagged,
    resid_df = resid_df,
    metrics_df = metrics_df,
    mono = mono,
    gap_diag = gap_diag
  )

  out_path <- file.path(
    PROJ_ROOT,
    "paper",
    "objects",
    paste0("stage3_", polarity, ".rds")
  )
  saveRDS(output, out_path)
  log_msg("INFO", polarity, "saved:", out_path)

  # ------ One-screen LOOCV summary -----------------------------------------
  summ <- metrics_df[
    metrics_df$cv_type == "loocv" & metrics_df$zone == "ALL",
    c("method", "n", "rmse", "mae", "max_ae", "coverage", "time_secs")
  ]
  summ$time_secs <- round(summ$time_secs, 2)
  message("\nLOOCV global metrics (sorted by RMSE):")
  print(summ[order(summ$rmse), ], row.names = FALSE, digits = 3)

  message("\nMonotonicity violations:")
  for (mn in names(mono)) {
    message(sprintf(
      "  %-14s  unweighted=%d  magnitude=%.1f s",
      mn,
      mono[[mn]]$unweighted,
      mono[[mn]]$magnitude_weighted
    ))
  }

  invisible(output)
}

# =========================================================================
# Entry point
# =========================================================================
for (pol in c("pos", "neg")) {
  run_eval(pol)
}
