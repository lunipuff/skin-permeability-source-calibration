# 84_diagnose_pooled_calibration.R
#
# Run from the project root, after 81:
# Rscript R/84_diagnose_pooled_calibration.R
#
# Uses saved predictions and anchor panels.
# Does not refit models or alter outputs from 80-83.

input_dir <- "R/results"
output_dir <- input_dir

read_result <- function(filename, required) {
  path <- file.path(input_dir, filename)

  if (!file.exists(path)) stop("Missing file: ", path)

  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(required %in% names(x))) {
    stop("Missing columns in ", filename)
  }

  if (!nrow(x)) stop("Empty input: ", filename)

  if (anyNA(x[required])) {
    stop("Missing required values in ", filename)
  }

  x$source_id <- as.character(x$source_id)
  x
}

write_result <- function(x, filename) {
  write.csv(
    x,
    file.path(output_dir, filename),
    row.names = FALSE,
    na = ""
  )
}

close_enough <- function(x, y) {
  all(
    is.finite(x) &
      is.finite(y) &
      abs(x - y) <= 1e-7 * pmax(1, abs(y))
  )
}

# ------------------------------------------------------------
# 1. Read the primary panels and predictions saved by 81
# ------------------------------------------------------------

draws <- read_result(
  "81_pooled_calibration_draws.csv",
  c(
    "source_id", "anchors", "draw", "selected_compounds",
    "offset", "offset_raw", "shrinkage",
    "RMSE_uncalibrated", "RMSE_calibrated"
  )
)

anchors <- read_result(
  "81_pooled_anchor_residuals.csv",
  c(
    "source_id", "compound_id",
    "mean_residual", "n_measurements"
  )
)

tests <- read_result(
  "81_pooled_uncalibrated_test_predictions.csv",
  c(
    "source_id", "record_id", "compound_id",
    "log_kp", "uncalibrated_prediction"
  )
)

primary <- read_result(
  "81_pooled_calibration_by_source.csv",
  c(
    "source_id", "anchors",
    "RMSE_uncalibrated", "RMSE_calibrated",
    "n_available_anchors",
    "n_test_compounds", "n_test_measurements"
  )
)

# ------------------------------------------------------------
# 2. Check input consistency
# ------------------------------------------------------------

for (check in list(
  draws[c("source_id", "anchors", "draw")],
  anchors[c("source_id", "compound_id")],
  tests[c("source_id", "record_id")],
  primary[c("source_id", "anchors")]
)) {
  if (anyDuplicated(check)) {
    stop("Duplicate input identifiers.")
  }
}

for (x in list(anchors, tests, primary)) {
  if (!setequal(draws$source_id, x$source_id)) {
    stop("Input publications do not match.")
  }
}

for (x in list(draws, anchors, tests, primary)) {
  numeric_fields <- vapply(x, is.numeric, logical(1))

  if (any(!is.finite(as.matrix(x[numeric_fields])))) {
    stop("Non-finite numeric inputs.")
  }
}

if (!close_enough(
  draws$offset,
  draws$offset_raw * draws$shrinkage
)) {
  stop("Saved offsets do not match raw offsets and shrinkage.")
}

if (
  any(draws$anchors < 1 |
        draws$anchors != floor(draws$anchors)) ||
  any(draws$shrinkage < 0 | draws$shrinkage > 1)
) {
  stop("Invalid panel settings.")
}

# ------------------------------------------------------------
# 3. Diagnose the primary sampled anchor panels
# ------------------------------------------------------------

# Uncalibrated residual:
#   e = observed - predicted
#
# Calibrated residual:
#   e - offset
#
# The mean test residual is calculated retrospectively.
# It is never used to select anchors or calibrate predictions.

tests$uncalibrated_residual <-
  tests$log_kp - tests$uncalibrated_prediction

summaries <- list()
panels <- list()
single_anchors <- list()

tolerance <- 1e-10

for (source in sort(unique(draws$source_id))) {

  d <- draws[draws$source_id == source, , drop = FALSE]
  a <- anchors[anchors$source_id == source, , drop = FALSE]
  t <- tests[tests$source_id == source, , drop = FALSE]

  e <- t$uncalibrated_residual

  mu <- mean(e)
  mse0 <- mean(e^2)
  rmse0 <- sqrt(mse0)

  # Population variance across the retained test measurements.
  residual_variance <- mean((e - mu)^2)

  if (!any(d$anchors == 1)) {
    stop("Missing one-anchor draws: ", source)
  }

  for (k in sort(unique(d$anchors))) {

    x <- d[d$anchors == k, , drop = FALSE]

    p <- primary[
      primary$source_id == source &
        primary$anchors == k,
      ,
      drop = FALSE
    ]

    if (nrow(p) != 1L) {
      stop("Primary summary mismatch: ", source)
    }

    if (
      p$n_available_anchors != nrow(a) ||
      p$n_test_compounds != length(unique(t$compound_id)) ||
      p$n_test_measurements != nrow(t)
    ) {
      stop("Input counts differ: ", source)
    }

    # Exact calibrated MSE for a constant offset.
    mse1 <- residual_variance + (mu - x$offset)^2

    if (
      !close_enough(x$RMSE_uncalibrated, rmse0) ||
      !close_enough(x$RMSE_calibrated, sqrt(mse1)) ||
      !close_enough(
        mean(x$RMSE_calibrated),
        p$RMSE_calibrated
      )
    ) {
      stop(
        "Saved predictions and primary RMSE disagree: ",
        source
      )
    }

    mean_offset <- mean(x$offset)

    # Denominator n, rather than n - 1, preserves the identity.
    offset_variance <-
      mean((x$offset - mean_offset)^2)

    mismatch_squared <- (mean_offset - mu)^2
    bias_squared <- mu^2
    delta_mse <- mean(mse1) - mse0

    # Mean MSE change:
    # offset variance + squared mismatch - squared test bias.
    if (!close_enough(
      delta_mse,
      offset_variance + mismatch_squared - bias_squared
    )) {
      stop("MSE decomposition failed: ", source)
    }

    x$test_mean_residual <- mu
    x$MSE_change <- mse1 - mse0
    x$RMSE_change <- x$RMSE_calibrated - rmse0

    x$improved <- x$RMSE_change < -tolerance
    x$worsened <- x$RMSE_change > tolerance

    x$offset_wrong_direction <-
      x$offset * mu < -tolerance

    # Same-direction adjustment that increases error.
    x$offset_overshoot <-
      x$worsened &
      !x$offset_wrong_direction &
      abs(mu) > tolerance

    panels[[length(panels) + 1L]] <- x

    summaries[[length(summaries) + 1L]] <- data.frame(
      source_id = source,
      anchors = k,

      available_anchors = nrow(a),
      test_compounds = length(unique(t$compound_id)),
      test_measurements = nrow(t),

      draws = nrow(x),

      # Count panels as sets, irrespective of selection order.
      distinct_panels = length(unique(vapply(
        strsplit(
          x$selected_compounds,
          " | ",
          fixed = TRUE
        ),
        function(z) paste(sort(z), collapse = " | "),
        character(1)
      ))),

      median_anchor_replication =
        median(a$n_measurements),

      median_test_replication =
        median(as.numeric(table(t$compound_id))),

      RMSE_uncalibrated = rmse0,
      RMSE_calibrated = mean(x$RMSE_calibrated),

      RMSE_change =
        mean(x$RMSE_calibrated) - rmse0,

      RMSE_reduction_percent = if (rmse0 > 0) {
        100 * (1 - mean(x$RMSE_calibrated) / rmse0)
      } else {
        NA_real_
      },

      mean_RMSE_worsened =
        mean(x$RMSE_calibrated) > rmse0 + tolerance,

      mean_RMSE_improved =
        mean(x$RMSE_calibrated) < rmse0 - tolerance,

      draws_improved_percent =
        100 * mean(x$improved),

      draws_wrong_direction_percent =
        100 * mean(x$offset_wrong_direction),

      draws_overshoot_percent =
        100 * mean(x$offset_overshoot),

      test_mean_residual = mu,

      test_bias_fraction_of_MSE = if (mse0 > 0) {
        bias_squared / mse0
      } else {
        NA_real_
      },

      offset_mean = mean_offset,

      offset_sd = if (nrow(x) > 1) {
        sd(x$offset)
      } else {
        0
      },

      offset_min = min(x$offset),
      offset_max = max(x$offset),

      mean_offset_wrong_direction =
        mean_offset * mu < -tolerance,

      mean_shrinkage = mean(x$shrinkage),

      MSE_uncalibrated = mse0,
      mean_MSE_calibrated = mean(mse1),
      mean_MSE_change = delta_mse,

      offset_variance = offset_variance,
      offset_mismatch_squared = mismatch_squared,
      correctable_test_bias_squared = bias_squared,

      dominant_penalty = if (
        mismatch_squared >= offset_variance
      ) {
        "mean_offset_mismatch"
      } else {
        "offset_variability"
      },

      retrospective_optimal_offset_RMSE =
        sqrt(residual_variance),

      stringsAsFactors = FALSE
    )
  }

  # ----------------------------------------------------------
  # Supplementary diagnostic: every single anchor evaluated once
  # ----------------------------------------------------------

  # This is distinct from the sampled primary analysis above.
  # Reuse saved residuals and the fitted one-anchor shrinkage.

  one <- d[d$anchors == 1, , drop = FALSE]

  if (!close_enough(
    one$shrinkage,
    one$shrinkage[1]
  )) {
    stop("Inconsistent shrinkage: ", source)
  }

  j <- match(
    one$selected_compounds,
    a$compound_id
  )

  if (
    anyNA(j) ||
    !close_enough(
      one$offset_raw,
      a$mean_residual[j]
    )
  ) {
    stop(
      "Saved single-anchor choices and residuals disagree: ",
      source
    )
  }

  a$offset <-
    a$mean_residual * one$shrinkage[1]

  a$RMSE_uncalibrated <- rmse0

  a$RMSE_calibrated <-
    sqrt(residual_variance + (mu - a$offset)^2)

  a$test_mean_residual <- mu

  a$improved <-
    a$RMSE_calibrated < rmse0 - tolerance

  a$offset_wrong_direction <-
    a$offset * mu < -tolerance

  single_anchors[[length(single_anchors) + 1L]] <- a
}

summary84 <- do.call(rbind, summaries)
selection84 <- do.call(rbind, panels)
single84 <- do.call(rbind, single_anchors)

if (nrow(summary84) != nrow(primary)) {
  stop("Source/panel summaries are incomplete.")
}

# Within each panel size, place the largest deterioration first.
summary84 <- summary84[
  order(summary84$anchors, -summary84$RMSE_change),
]

rownames(summary84) <- NULL
rownames(selection84) <- NULL
rownames(single84) <- NULL

# ------------------------------------------------------------
# 4. One versus five anchors within the same publications
# ------------------------------------------------------------

fields <- c(
  "source_id",
  "RMSE_calibrated",
  "offset_mean",
  "offset_variance",
  "offset_mismatch_squared",
  "draws_improved_percent",
  "mean_shrinkage"
)

paired84 <- merge(
  summary84[
    summary84$anchors == 1,
    fields,
    drop = FALSE
  ],
  summary84[
    summary84$anchors == 5,
    fields,
    drop = FALSE
  ],
  by = "source_id",
  suffixes = c("_one", "_five")
)

paired84$RMSE_decrease <-
  paired84$RMSE_calibrated_one -
  paired84$RMSE_calibrated_five

paired84$additional_anchors_helped <-
  paired84$RMSE_decrease > tolerance

# ------------------------------------------------------------
# 5. Compact overview
# ------------------------------------------------------------

# Eligible publications can differ across anchor counts.
# Use 81's fixed-source estimates for the overall anchor trend.

overview84 <- do.call(rbind, lapply(
  sort(unique(summary84$anchors)),
  function(k) {

    x <- summary84[
      summary84$anchors == k,
      ,
      drop = FALSE
    ]

    data.frame(
      anchors = k,
      publications = nrow(x),

      publications_improved =
        sum(x$mean_RMSE_improved),

      publications_worsened =
        sum(x$mean_RMSE_worsened),

      failures_with_mean_offset_wrong_direction = sum(
        x$mean_RMSE_worsened &
          x$mean_offset_wrong_direction
      ),

      failures_dominated_by_mean_mismatch = sum(
        x$mean_RMSE_worsened &
          x$dominant_penalty == "mean_offset_mismatch"
      ),

      failures_dominated_by_offset_variability = sum(
        x$mean_RMSE_worsened &
          x$dominant_penalty == "offset_variability"
      )
    )
  }
))

# ------------------------------------------------------------
# 6. Write diagnostics
# ------------------------------------------------------------

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_result(
  summary84,
  "84_failure_source_summary.csv"
)

write_result(
  selection84,
  "84_failure_anchor_selections.csv"
)

write_result(
  single84,
  "84_all_single_anchor_diagnostics.csv"
)

write_result(
  paired84,
  "84_one_vs_five_diagnostics.csv"
)

write_result(
  overview84,
  "84_failure_overview.csv"
)

notes <- c(
  "Pooled calibration diagnostics completed without model refitting.",

  paste(
    "Primary summaries use exactly the sampled panels saved by 81,",
    "including repeated selections."
  ),

  paste(
    "Exhaustive single-anchor diagnostics are supplementary and",
    "can differ from primary sampled summaries."
  ),

  paste(
    "Test residuals and optimal-offset RMSE are retrospective",
    "diagnostics, not deployment rules."
  ),

  paste(
    "MSE change = offset variance + squared mean-offset mismatch",
    "- squared mean test bias."
  ),

  paste(
    "Offset variance uses denominator n.",
    "Mean squared RMSE differs from squared mean RMSE."
  ),

  paste(
    "Dominant penalties are algebraic descriptions,",
    "not identified experimental causes."
  ),

  paste(
    "Increasing panel size also changes shrinkage;",
    "observed changes reflect both."
  ),

  paste(
    "No confidence intervals are calculated here.",
    "Panel ranges are not confidence intervals."
  ),

  "Use 81's fixed-source estimates for the overall anchor-number conclusion.",

  "",
  capture.output(print(overview84, row.names = FALSE))
)

writeLines(
  notes,
  file.path(output_dir, "84_diagnostic_status.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "84_sessionInfo.txt")
)

cat(paste(notes, collapse = "\n"), "\n")