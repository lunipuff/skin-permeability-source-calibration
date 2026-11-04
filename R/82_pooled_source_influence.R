# 82_pooled_source_influence.R
#
# Removes one evaluation publication at a time from the summaries.
# Does not refit models or alter the pooled data.

# Optional settings use --name=value; no arguments retain main-analysis defaults.
parse_options <- function(defaults) {
  for (arg in commandArgs(trailingOnly = TRUE)) {
    if (!grepl("^--[^=]+=", arg)) stop("Use --name=value: ", arg)
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    if (!key %in% names(defaults)) stop("Unknown option: --", key)
    defaults[[key]] <- sub("^--[^=]+=", "", arg)
  }
  defaults
}

options_run <- parse_options(list(
  "input-file" = "R/results/81_pooled_calibration_by_source.csv", "output-dir" = "R/results"
))
input <- options_run[["input-file"]]
output_dir <- options_run[["output-dir"]]
if (!nzchar(input) || !nzchar(output_dir)) stop("Input and output paths must be nonempty.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_boot <- 20000L
set.seed(2026)

if (!file.exists(input)) {
  stop("Run 81_validate_pooled_calibration.R first.")
}

dat <- read.csv(input, stringsAsFactors = FALSE)

required <- c(
  "source_id", "anchors",
  "RMSE_uncalibrated", "RMSE_calibrated"
)

if (!all(required %in% names(dat))) {
  stop("Missing required input columns.")
}

if (anyDuplicated(dat[c("source_id", "anchors")])) {
  stop("Expected one row per source and anchor count.")
}

if (any(!is.finite(dat$RMSE_uncalibrated)) ||
    any(!is.finite(dat$RMSE_calibrated))) {
  stop("Non-finite RMSE values.")
}

# ------------------------------------------------------------
# 1. Prepare the two comparisons
# ------------------------------------------------------------

one <- dat[dat$anchors == 1, ]

one_comparison <- data.frame(
  source_id = one$source_id,
  before = one$RMSE_uncalibrated,
  after = one$RMSE_calibrated
)

five <- dat[dat$anchors == 5, ]

paired <- merge(
  one[c("source_id", "RMSE_calibrated")],
  five[c("source_id", "RMSE_calibrated")],
  by = "source_id",
  suffixes = c("_one", "_five")
)

panel_comparison <- data.frame(
  source_id = paired$source_id,
  before = paired$RMSE_calibrated_one,
  after = paired$RMSE_calibrated_five
)

# ------------------------------------------------------------
# 2. Estimate effect and source-bootstrap interval
# ------------------------------------------------------------

estimate <- function(x, percentage) {
  n <- nrow(x)

  if (!n) stop("No sources available.")

  before <- mean(x$before)
  after <- mean(x$after)

  effect <- if (percentage) {
    100 * (1 - after / before)
  } else {
    before - after
  }

  interval <- c(NA_real_, NA_real_)

  if (n >= 2L) {
    indices <- matrix(
      sample.int(n, n * n_boot, replace = TRUE),
      nrow = n_boot
    )

    boot_before <- rowMeans(matrix(
      x$before[indices], nrow = n_boot
    ))
    boot_after <- rowMeans(matrix(
      x$after[indices], nrow = n_boot
    ))

    boot_effect <- if (percentage) {
      100 * (1 - boot_after / boot_before)
    } else {
      boot_before - boot_after
    }

    interval <- as.numeric(quantile(
      boot_effect,
      c(0.025, 0.975),
      names = FALSE
    ))
  }

  data.frame(
    n_sources = n,
    mean_RMSE_before = before,
    mean_RMSE_after = after,
    effect = effect,
    ci_low = interval[1],
    ci_high = interval[2]
  )
}

# ------------------------------------------------------------
# 3. Omit each evaluation publication in turn
# ------------------------------------------------------------

run_influence <- function(x, analysis, percentage) {
  if (nrow(x) < 3L) {
    stop("At least three sources are needed for ", analysis)
  }

  x <- x[order(x$source_id), ]
  full <- estimate(x, percentage)

  results <- lapply(x$source_id, function(source) {
    remaining <- x[x$source_id != source, ]
    result <- estimate(remaining, percentage)

    data.frame(
      analysis = analysis,
      omitted_source = source,
      full_effect = full$effect,
      result,
      effect_change = result$effect - full$effect,
      direction_reversed = result$effect * full$effect < 0,
      interval_includes_zero =
        result$ci_low <= 0 && result$ci_high >= 0,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results)
  rownames(results) <- NULL

  results <- results[
    order(-abs(results$effect_change)),
  ]

  full$analysis <- analysis

  list(full = full, omissions = results)
}

calibration <- run_influence(
  one_comparison,
  "One-anchor RMSE reduction (%)",
  percentage = TRUE
)

panels <- run_influence(
  panel_comparison,
  "RMSE decrease from one to five anchors (log units)",
  percentage = FALSE
)

full_results <- rbind(calibration$full, panels$full)
omission_results <- rbind(
  calibration$omissions,
  panels$omissions
)

write.csv(
  full_results,
  file.path(output_dir, "82_full_sample_effects.csv"),
  row.names = FALSE
)

write.csv(
  omission_results,
  file.path(output_dir, "82_source_omission_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Compact influence summary
# ------------------------------------------------------------

summary_results <- do.call(rbind, lapply(
  unique(omission_results$analysis),
  function(analysis) {
    x <- omission_results[
      omission_results$analysis == analysis,
    ]

    data.frame(
      analysis = analysis,
      full_effect = x$full_effect[1],
      minimum_omission_effect = min(x$effect),
      maximum_omission_effect = max(x$effect),
      largest_absolute_change = max(abs(x$effect_change)),
      most_influential_source =
        x$omitted_source[which.max(abs(x$effect_change))],
      n_direction_reversals = sum(x$direction_reversed),
      n_intervals_including_zero = sum(x$interval_includes_zero),
      n_omissions = nrow(x),
      stringsAsFactors = FALSE
    )
  }
))

write.csv(
  summary_results,
  file.path(output_dir, "82_source_influence_summary.csv"),
  row.names = FALSE
)

print(summary_results, row.names = FALSE)

cat(
  "\nPositive effects indicate improved prediction.",
  "\nIntervals use ", n_boot, " publication-bootstrap replicates.",
  "\nModels were not refitted.",
  "\nDo not exclude publications merely because they weaken the result.\n",
  sep = ""
)
write.csv(data.frame(
  input_file = normalizePath(input, winslash = "/", mustWork = TRUE),
  input_md5 = unname(tools::md5sum(input)), n_boot = n_boot, seed = 2026L
), file.path(output_dir, "82_analysis_settings.csv"), row.names = FALSE)
