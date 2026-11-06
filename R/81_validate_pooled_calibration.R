# 81_validate_pooled_calibration.R
#
# Run from the project root:
# Rscript R/81_validate_pooled_calibration.R
#
# Input:
#   data/processed/80_pooled.csv
#
# Analysis:
#   Publication-held-out calibration.
#   Anchors: compounds shared with training publications.
#   Tests: compounds absent from all training publications.
#
# Repeated measurements:
#   Retained during hierarchical fitting and test evaluation.
#   Anchor residuals averaged within each compound.
#
# Uncertainty:
#   Publication-level bootstrap of source-level performance.
#   Training models are not refitted during this bootstrap.
#   Performance is averaged across all panels within source before resampling.

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Install lme4 before running this script.")
}

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
  "input-file" = "data/processed/80_pooled.csv", "output-dir" = "R/results",
  "reuse-predictions" = "false"
))
input_file <- options_run[["input-file"]]
output_dir <- options_run[["output-dir"]]
if (!options_run[["reuse-predictions"]] %in% c("true", "false")) {
  stop("--reuse-predictions must be true or false.")
}
reuse_predictions <- options_run[["reuse-predictions"]] == "true"
if (!nzchar(input_file) || !nzchar(output_dir)) stop("Input and output paths must be nonempty.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Missing pooled dataset. Complete script 80 first.")
}

anchor_sizes <- 1:5
panel_method <- "exhaustive"
# Legacy settings field. There are no random panel draws in this analysis.
n_draws <- 0L
n_boot <- 2000L
seed <- 2026L

set.seed(seed)

write_csv <- function(x, filename) {
  write.csv(
    x, file.path(output_dir, filename),
    row.names = FALSE, na = ""
  )
}

bind_rows <- function(x) {
  result <- do.call(rbind, x)
  rownames(result) <- NULL
  result
}

metrics <- function(observed, predicted) {
  c(
    RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted))
  )
}

# Each unordered panel is evaluated once. Panel sizes are evaluated separately.
# The legacy 'draw' column is a panel identifier within source and panel size.
evaluate_panels <- function(source, anchor_residuals, test, variance) {
  anchor_residuals <- anchor_residuals[
    order(anchor_residuals$compound_id), , drop = FALSE
  ]
  shared <- as.character(anchor_residuals$compound_id)
  if (!length(shared) || anyDuplicated(shared) || !nrow(test)) {
    stop("Invalid saved anchors or test measurements for ", source)
  }
  if (any(shared %in% test$compound_id)) {
    stop("Anchor and test compounds overlap for ", source)
  }
  if (any(!is.finite(anchor_residuals$mean_residual)) ||
      any(!is.finite(test$log_kp)) ||
      any(!is.finite(test$uncalibrated_prediction)) ||
      any(!is.finite(variance[c("source", "residual")])) ||
      any(variance[c("source", "residual")] < 0) ||
      sum(variance[c("source", "residual")]) <= 0) {
    stop("Invalid predictions or variance estimates for ", source)
  }
  baseline <- metrics(test$log_kp, test$uncalibrated_prediction)
  results <- list()
  for (k in anchor_sizes[anchor_sizes <= length(shared)]) {
    expected <- choose(length(shared), k)
    if (!is.finite(expected) || expected > .Machine$integer.max) {
      stop("Too many panels to enumerate for ", source, " at ", k, " anchors.")
    }
    combinations <- utils::combn(seq_along(shared), k)
    if (ncol(combinations) != expected) stop("Panel enumeration failed.")
    shrinkage <- unname(k * variance["source"] /
      (variance["residual"] + k * variance["source"]))
    size_results <- lapply(seq_len(ncol(combinations)), function(panel) {
      selected <- combinations[, panel]
      offset_raw <- mean(anchor_residuals$mean_residual[selected])
      offset <- shrinkage * offset_raw
      calibrated <- metrics(
        test$log_kp, test$uncalibrated_prediction + offset
      )
      data.frame(
        source_id = source, anchors = k, draw = panel,
        selected_compounds = paste(shared[selected], collapse = " | "),
        n_available_anchors = length(shared),
        n_test_compounds = length(unique(test$compound_id)),
        n_test_measurements = nrow(test),
        offset_raw = offset_raw, shrinkage = shrinkage, offset = offset,
        RMSE_uncalibrated = unname(baseline["RMSE"]),
        RMSE_calibrated = unname(calibrated["RMSE"]),
        MAE_uncalibrated = unname(baseline["MAE"]),
        MAE_calibrated = unname(calibrated["MAE"]),
        panel_method = panel_method,
        stringsAsFactors = FALSE
      )
    })
    results[[length(results) + 1L]] <- bind_rows(size_results)
  }
  bind_rows(results)
}

dat <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c(
  "record_id", "dataset_origin", "publication_id",
  "compound_id", "log_kp", "mw", "logp"
)

if (!all(required %in% names(dat))) {
  stop("Missing input columns: ",
       paste(setdiff(required, names(dat)), collapse = ", "))
}

if (anyDuplicated(dat$record_id)) {
  stop("Duplicate record_id values in the pooled dataset.")
}

if (anyNA(dat$publication_id) || anyNA(dat$compound_id) ||
    any(!nzchar(dat$publication_id)) ||
    any(!nzchar(dat$compound_id))) {
  stop("Missing publication or compound identity.")
}

for (v in c("log_kp", "mw", "logp")) {
  if (!is.numeric(dat[[v]]) || any(!is.finite(dat[[v]]))) {
    stop("Non-finite or nonnumeric values in ", v)
  }
}

dat$source_id <- as.character(dat$publication_id)
dat$compound_id <- as.character(dat$compound_id)

dat <- dat[order(
  dat$source_id, dat$compound_id, dat$record_id
), ]
rownames(dat) <- NULL

# A shared compound identity should have one descriptor profile.
for (v in c("mw", "logp")) {
  spread <- tapply(
    dat[[v]], dat$compound_id,
    function(z) max(z) - min(z)
  )
  if (any(spread > 1e-8)) {
    stop("Inconsistent ", v, " values within compound_id.")
  }
}

# ------------------------------------------------------------
# 1. Training-only descriptor standardization
# ------------------------------------------------------------

standardize <- function(train, test = NULL) {
  profile <- train[
    !duplicated(train$compound_id),
    c("compound_id", "mw", "logp"),
    drop = FALSE
  ]

  center <- vapply(profile[c("mw", "logp")], mean, numeric(1))
  scale <- vapply(profile[c("mw", "logp")], sd, numeric(1))

  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("Insufficient descriptor variation in training data.")
  }

  for (v in c("mw", "logp")) {
    name <- paste0(v, "_z")
    train[[name]] <- (train[[v]] - center[v]) / scale[v]

    if (!is.null(test)) {
      test[[name]] <- (test[[v]] - center[v]) / scale[v]
    }
  }

  list(
    train = train,
    test = test,
    center = center,
    scale = scale
  )
}

control <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 100000)
)

fit_model <- function(x) {
  lme4::lmer(
    log_kp ~ mw_z + logp_z +
      (1 | compound_id) + (1 | source_id),
    data = x,
    REML = TRUE,
    control = control
  )
}

convergence_text <- function(fit) {
  messages <- fit@optinfo$conv$lme4$messages
  optimizer_code <- fit@optinfo$conv$opt

  result <- character()

  if (length(messages)) result <- c(result, as.character(messages))

  if (length(optimizer_code) && any(optimizer_code != 0)) {
    result <- c(
      result,
      paste("Optimizer code:", paste(optimizer_code, collapse = ", "))
    )
  }

  paste(unique(result), collapse = " | ")
}

variance_values <- function(fit) {
  vc <- as.data.frame(lme4::VarCorr(fit))

  c(
    compound = vc$vcov[vc$grp == "compound_id"],
    source = vc$vcov[vc$grp == "source_id"],
    residual = vc$vcov[vc$grp == "Residual"]
  )
}

# ------------------------------------------------------------
# 2. Determine eligibility before model fitting
# ------------------------------------------------------------

sources <- sort(unique(dat$source_id))

eligibility <- bind_rows(lapply(sources, function(source) {
  held <- dat[dat$source_id == source, ]
  train <- dat[dat$source_id != source, ]

  shared <- intersect(
    unique(held$compound_id),
    unique(train$compound_id)
  )
  novel <- setdiff(
    unique(held$compound_id),
    unique(train$compound_id)
  )

  data.frame(
    source_id = source,
    dataset_origins = paste(
      sort(unique(held$dataset_origin)), collapse = "; "
    ),
    n_measurements = nrow(held),
    n_compounds = length(unique(held$compound_id)),
    n_anchors = length(shared),
    n_test_compounds = length(novel),
    n_test_measurements = sum(held$compound_id %in% novel),
    stringsAsFactors = FALSE
  )
}))

for (k in anchor_sizes) {
  eligibility[[paste0("eligible_", k)]] <-
    eligibility$n_anchors >= k &
    eligibility$n_test_compounds >= 1
}

write_csv(eligibility, "81_pooled_source_eligibility.csv")

eligible_sources <- eligibility$source_id[
  eligibility$n_anchors >= 1 &
  eligibility$n_test_compounds >= 1
]

if (!length(eligible_sources)) {
  stop("No publications have both anchors and novel test compounds.")
}

# ------------------------------------------------------------
# 3. Evaluate one held-out publication
# ------------------------------------------------------------

evaluate_source <- function(source) {
  train <- dat[dat$source_id != source, ]
  held <- dat[dat$source_id == source, ]

  shared <- sort(intersect(
    unique(held$compound_id),
    unique(train$compound_id)
  ))
  novel <- sort(setdiff(
    unique(held$compound_id),
    unique(train$compound_id)
  ))

  if (!length(shared) || !length(novel)) {
    stop("Source is not eligible.")
  }

  sc <- standardize(train, held)
  fit <- fit_model(sc$train)

  variance <- variance_values(fit)

  if (length(variance) != 3L ||
      any(!is.finite(variance)) ||
      variance["source"] + variance["residual"] <= 0) {
    stop("Invalid variance estimates.")
  }

  anchor_rows <- sc$test[
    sc$test$compound_id %in% shared, , drop = FALSE
  ]

  # Include known compound effects, but not the held-out source effect.
  anchor_rows$expected <- predict(
    fit,
    newdata = anchor_rows,
    re.form = ~(1 | compound_id),
    allow.new.levels = TRUE
  )

  anchor_residuals <- aggregate(
    anchor_rows$log_kp - anchor_rows$expected,
    list(compound_id = anchor_rows$compound_id),
    mean
  )
  names(anchor_residuals)[2] <- "mean_residual"

  anchor_counts <- table(anchor_rows$compound_id)
  anchor_residuals$n_measurements <- as.integer(
    anchor_counts[anchor_residuals$compound_id]
  )
  anchor_residuals$source_id <- source

  test <- sc$test[
    sc$test$compound_id %in% novel, , drop = FALSE
  ]

  # Novel test compounds receive fixed-effect predictions only.
  test$uncalibrated_prediction <- predict(
    fit,
    newdata = test,
    re.form = NA,
    allow.new.levels = TRUE
  )

  if (any(!is.finite(anchor_rows$expected)) ||
      any(!is.finite(test$uncalibrated_prediction))) {
    stop("Non-finite predictions.")
  }

  list(
    draws = evaluate_panels(source, anchor_residuals, test, variance),
    anchors = anchor_residuals,
    tests = test[c(
      "record_id", "dataset_origin", "source_id",
      "compound_id", "log_kp", "uncalibrated_prediction"
    )],
    singular = lme4::isSingular(fit),
    convergence = convergence_text(fit),
    variance = variance
  )
}

# ------------------------------------------------------------
# 4. Fit eligible source-held-out models and record diagnostics
# ------------------------------------------------------------

draw_results <- list()
anchor_results <- list()
test_results <- list()
diagnostics <- list()

if (reuse_predictions) {
  read_cached <- function(filename, columns) {
    path <- file.path(output_dir, filename)
    if (!file.exists(path)) stop("Missing cached output ", filename, ". Run 81 without --reuse-predictions.")
    x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!nrow(x) || !all(columns %in% names(x))) stop("Invalid cached output ", filename)
    x
  }
  cached_settings <- read_cached("81_analysis_settings.csv", "input_md5")
  if (nrow(cached_settings) != 1L ||
      !identical(as.character(cached_settings$input_md5), unname(tools::md5sum(input_file)))) {
    stop("The pooled input has changed. Run 81 without --reuse-predictions first.")
  }
  diagnostics <- read_cached("81_pooled_fit_diagnostics.csv",
    c("source_id", "success", "singular", "convergence", "warnings", "error",
      "source_variance", "residual_variance"))
  for (field in c("convergence", "warnings", "error")) {
    diagnostics[[field]][is.na(diagnostics[[field]])] <- ""
  }
  cached_anchors <- read_cached("81_pooled_anchor_residuals.csv",
    c("source_id", "compound_id", "mean_residual", "n_measurements"))
  cached_tests <- read_cached("81_pooled_uncalibrated_test_predictions.csv",
    c("source_id", "record_id", "compound_id", "log_kp", "uncalibrated_prediction"))
  if (anyDuplicated(cached_anchors[c("source_id", "compound_id")]) ||
      anyDuplicated(cached_tests$record_id)) {
    stop("Duplicate cached records. Run 81 without --reuse-predictions.")
  }
  j <- match(cached_tests$record_id, dat$record_id)
  if (anyNA(j) || anyNA(cached_tests[c("source_id", "compound_id", "log_kp")]) ||
      any(cached_tests$source_id != dat$source_id[j]) ||
      any(cached_tests$compound_id != dat$compound_id[j]) ||
      any(abs(cached_tests$log_kp - dat$log_kp[j]) > 1e-12)) {
    stop("Cached test records differ from the input. Run 81 without --reuse-predictions.")
  }
  if (anyDuplicated(diagnostics$source_id) || anyNA(diagnostics$success) ||
      !all(diagnostics$success) ||
      !setequal(diagnostics$source_id, eligible_sources) ||
      !setequal(cached_anchors$source_id, eligible_sources) ||
      !setequal(cached_tests$source_id, eligible_sources)) {
    stop("Cached folds are incomplete. Run 81 without --reuse-predictions.")
  }
  for (source in eligible_sources) {
    a <- cached_anchors[cached_anchors$source_id == source, , drop = FALSE]
    t <- cached_tests[cached_tests$source_id == source, , drop = FALSE]
    v <- diagnostics[diagnostics$source_id == source, ]
    held_ids <- unique(dat$compound_id[dat$source_id == source])
    train_ids <- unique(dat$compound_id[dat$source_id != source])
    if (!setequal(a$compound_id, intersect(held_ids, train_ids)) ||
        !setequal(t$compound_id, setdiff(held_ids, train_ids))) {
      stop("Cached compound partitions differ for ", source)
    }
    draw_results[[length(draw_results) + 1L]] <- evaluate_panels(
      source, a, t,
      c(source = v$source_variance, residual = v$residual_variance)
    )
  }
  anchor_results <- list(cached_anchors)
  test_results <- list(cached_tests)
} else {
for (i in seq_along(eligible_sources)) {
  source <- eligible_sources[i]
  warnings_seen <- character()

  message(
    "[", i, "/", length(eligible_sources), "] ",
    source
  )

  result <- tryCatch(
    withCallingHandlers(
      evaluate_source(source),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  failed <- inherits(result, "error")

  diagnostics[[i]] <- data.frame(
    source_id = source,
    success = !failed,
    singular = if (failed) NA else result$singular,
    convergence = if (failed) "" else result$convergence,
    warnings = paste(unique(warnings_seen), collapse = " | "),
    error = if (failed) conditionMessage(result) else "",
    source_variance = if (
      failed
    ) NA_real_ else unname(result$variance["source"]),
    residual_variance = if (
      failed
    ) NA_real_ else unname(result$variance["residual"]),
    stringsAsFactors = FALSE
  )

  if (failed) {
    message("Fit failed: ", conditionMessage(result))
    next
  }

  draw_results[[length(draw_results) + 1L]] <- result$draws
  anchor_results[[length(anchor_results) + 1L]] <- result$anchors
  test_results[[length(test_results) + 1L]] <- result$tests
}

diagnostics <- bind_rows(diagnostics)
}
write_csv(diagnostics, "81_pooled_fit_diagnostics.csv")

if (!length(draw_results)) {
  stop("All eligible fits failed. Inspect 81_pooled_fit_diagnostics.csv.")
}

draws <- bind_rows(draw_results)

write_csv(draws, "81_pooled_calibration_draws.csv")
write_csv(
  bind_rows(anchor_results),
  "81_pooled_anchor_residuals.csv"
)
write_csv(
  bind_rows(test_results),
  "81_pooled_uncalibrated_test_predictions.csv"
)

# ------------------------------------------------------------
# 5. Average all panels within each source and anchor count
# ------------------------------------------------------------

source_results <- aggregate(
  draws[c(
    "RMSE_uncalibrated", "RMSE_calibrated",
    "MAE_uncalibrated", "MAE_calibrated",
    "offset_raw", "offset"
  )],
  draws[c("source_id", "anchors")],
  mean
)

source_counts <- unique(draws[c(
  "source_id", "anchors", "n_available_anchors",
  "n_test_compounds", "n_test_measurements"
)])

source_results <- merge(
  source_results, source_counts,
  by = c("source_id", "anchors"),
  sort = FALSE
)

panel_counts <- aggregate(
  list(n_panels = draws$draw), draws[c("source_id", "anchors")], length
)
source_results <- merge(source_results, panel_counts,
  by = c("source_id", "anchors"), sort = FALSE)

source_results$RMSE_reduction_percent <- with(
  source_results,
  100 * (1 - RMSE_calibrated / RMSE_uncalibrated)
)

source_results$MAE_reduction_percent <- with(
  source_results,
  100 * (1 - MAE_calibrated / MAE_uncalibrated)
)

source_results <- source_results[order(
  source_results$source_id, source_results$anchors
), ]

write_csv(
  source_results,
  "81_pooled_calibration_by_source.csv"
)

# ------------------------------------------------------------
# 6. Publication bootstrap and equal-source summaries
# ------------------------------------------------------------

summarize_sources <- function(x, analysis) {
  output <- list()

  for (k in sort(unique(x$anchors))) {
    z <- x[x$anchors == k, ]
    n <- nrow(z)

    # One row per source at this anchor count.
    if (anyDuplicated(z$source_id)) {
      stop("Duplicate source rows during bootstrap.")
    }

    baseline <- mean(z$RMSE_uncalibrated)
    calibrated <- mean(z$RMSE_calibrated)
    reduction <- 100 * (1 - calibrated / baseline)

    boot_reduction <- replicate(n_boot, {
      ii <- sample.int(n, n, replace = TRUE)
      100 * (
        1 -
          mean(z$RMSE_calibrated[ii]) /
          mean(z$RMSE_uncalibrated[ii])
      )
    })

    # A single publication cannot support a source-level interval.
    interval <- if (n >= 2L) {
      as.numeric(quantile(
        boot_reduction, c(0.025, 0.975),
        names = FALSE
      ))
    } else {
      c(NA_real_, NA_real_)
    }

    output[[length(output) + 1L]] <- data.frame(
      analysis = analysis,
      anchors = k,
      n_sources = n,
      n_sources_improved = sum(
        z$RMSE_calibrated < z$RMSE_uncalibrated
      ),
      mean_source_RMSE_uncalibrated = baseline,
      mean_source_RMSE_calibrated = calibrated,
      RMSE_reduction_percent = reduction,
      reduction_ci_low = interval[1],
      reduction_ci_high = interval[2],
      mean_source_MAE_uncalibrated = mean(z$MAE_uncalibrated),
      mean_source_MAE_calibrated = mean(z$MAE_calibrated),
      n_boot = n_boot,
      stringsAsFactors = FALSE
    )
  }

  bind_rows(output)
}

set.seed(seed + 1L)

all_summary <- summarize_sources(
  source_results,
  "All eligible sources at each anchor count"
)

write_csv(
  all_summary,
  "81_pooled_calibration_summary_all.csv"
)

# Fixed publications supporting all five anchor counts.
matched_sources <- source_results$source_id[
  source_results$anchors == max(anchor_sizes)
]

matched <- source_results[
  source_results$source_id %in% matched_sources, ]

if (length(matched_sources)) {
  expected_counts <- table(matched$source_id)
  if (any(expected_counts != length(anchor_sizes))) {
    stop("Matched sources do not have every requested anchor count.")
  }

  matched_summary <- summarize_sources(
    matched,
    "Fixed sources eligible for all anchor counts"
  )
} else {
  matched_summary <- all_summary[FALSE, ]
}

write_csv(
  matched_summary,
  "81_pooled_calibration_summary_matched.csv"
)

# ------------------------------------------------------------
# 7. Paired one-versus-five-anchor comparison
# ------------------------------------------------------------

if (length(matched_sources)) {
  one <- matched[
    matched$anchors == min(anchor_sizes),
    c("source_id", "RMSE_calibrated")
  ]
  five <- matched[
    matched$anchors == max(anchor_sizes),
    c("source_id", "RMSE_calibrated")
  ]

  comparison <- merge(
    one, five, by = "source_id",
    suffixes = c("_one", "_five")
  )

  # Positive values mean lower error with five anchors.
  gain <- comparison$RMSE_calibrated_one -
    comparison$RMSE_calibrated_five

  set.seed(seed + 2L)

  gain_boot <- replicate(n_boot, {
    ii <- sample.int(length(gain), length(gain), replace = TRUE)
    mean(gain[ii])
  })

  gain_ci <- if (length(gain) >= 2L) {
    as.numeric(quantile(
      gain_boot, c(0.025, 0.975),
      names = FALSE
    ))
  } else {
    c(NA_real_, NA_real_)
  }

  paired_summary <- data.frame(
    n_sources = length(gain),
    mean_RMSE_one_anchor = mean(comparison$RMSE_calibrated_one),
    mean_RMSE_five_anchors = mean(comparison$RMSE_calibrated_five),
    mean_RMSE_decrease = mean(gain),
    decrease_ci_low = gain_ci[1],
    decrease_ci_high = gain_ci[2],
    n_boot = n_boot
  )

  write_csv(
    paired_summary,
    "81_pooled_one_vs_five_anchors.csv"
  )
} else {
  write_csv(
    data.frame(
      n_sources = 0L,
      mean_RMSE_one_anchor = NA_real_,
      mean_RMSE_five_anchors = NA_real_,
      mean_RMSE_decrease = NA_real_,
      decrease_ci_low = NA_real_,
      decrease_ci_high = NA_real_,
      n_boot = n_boot
    ),
    "81_pooled_one_vs_five_anchors.csv"
  )
}

# ------------------------------------------------------------
# 8. Record run settings and completion
# ------------------------------------------------------------

status <- c(
  "Pooled publication-held-out calibration completed.",
  paste("Input measurements:", nrow(dat)),
  paste("Input compounds:", length(unique(dat$compound_id))),
  paste("Input publications:", length(sources)),
  paste("Eligible publications:", length(eligible_sources)),
  paste("Successful publication folds:", sum(diagnostics$success)),
  paste("Failed publication folds:", sum(!diagnostics$success)),
  paste("Singular successful fits:", sum(diagnostics$singular, na.rm = TRUE)),
  paste(
    "Successful fits with convergence messages:",
    sum(diagnostics$success & nzchar(diagnostics$convergence))
  ),
  paste("Fixed-source five-anchor subset:", length(matched_sources)),
  "Anchor panels were evaluated exhaustively at each eligible panel size.",
  paste("Total panels evaluated:", nrow(draws)),
  paste("Saved predictions reused:", reuse_predictions),
  paste("Publication bootstrap replicates:", n_boot),
  paste("Seed:", seed),
  "",
  "Primary descriptors: RDKit molecular weight and Crippen logP.",
  "Temperature was not included in the primary model.",
  "Training descriptor scaling used unique training compounds.",
  "All records of the target publication were excluded from fitting.",
  "Test compounds were absent from all training publications.",
  "Anchor residuals were averaged within each compound.",
  "Every unordered panel was evaluated once at each anchor count.",
  "Panel identifiers do not link or pair panels across anchor counts.",
  "Test errors were evaluated across measurement rows.",
  "Source-level metrics were averaged with equal source weights.",
  "Intervals resample publications after averaging over all anchor panels.",
  "Singular fits and fits with convergence messages were retained and flagged.",
  "Failed folds were omitted from summaries and explicitly reported.",
  "Results are conditional on pooled curation and identity mappings.",
  "",
  "All-eligible results:",
  capture.output(print(all_summary, row.names = FALSE)),
  "",
  "Matched-source results:",
  capture.output(print(matched_summary, row.names = FALSE))
)

writeLines(
  status,
  file.path(output_dir, "81_pooled_calibration_status.txt")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "81_pooled_sessionInfo.txt")
)

cat(paste(status, collapse = "\n"), "\n")
# Written only after successful completion, for input provenance checks.
write_csv(data.frame(
  input_file = normalizePath(input_file, winslash = "/", mustWork = TRUE),
  input_md5 = unname(tools::md5sum(input_file)),
  n_draws = n_draws, n_boot = n_boot, seed = seed,
  panel_method = panel_method, total_panels = nrow(draws),
  reused_predictions = reuse_predictions
), "81_analysis_settings.csv")
