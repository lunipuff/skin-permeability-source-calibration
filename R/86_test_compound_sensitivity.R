# 86_test_compound_sensitivity.R
# Run from the project root after exhaustive script 81.
# Rscript R/86_test_compound_sensitivity.R
# Uses saved source-level performance. No models or panels are changed.

input_dir <- "R/results"
table_dir <- "R/tables"
minimum_test_compounds <- 2L

read_result <- function(filename, required) {
  path <- file.path(input_dir, filename)
  if (!file.exists(path)) stop("Missing input: ", path)
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(x) || !all(required %in% names(x))) {
    stop("Empty input or missing columns in ", filename)
  }
  x
}
write_result <- function(x, filename) {
  write.csv(x, file.path(input_dir, filename), row.names = FALSE, na = "")
}
close_enough <- function(x, y) {
  length(x) == length(y) && all(is.finite(x) & is.finite(y) &
    abs(x - y) <= 1e-8 * pmax(1, abs(y)))
}

settings <- read_result("81_analysis_settings.csv",
  c("input_md5", "n_boot", "seed", "panel_method"))
if (nrow(settings) != 1L || settings$panel_method != "exhaustive") {
  stop("Run exhaustive script 81 first.")
}
input_file <- "data/processed/80_pooled.csv"
if (!file.exists(input_file) ||
    !identical(as.character(settings$input_md5), unname(tools::md5sum(input_file)))) {
  stop("The pooled input does not match 81. Regenerate calibration first.")
}
n_boot <- as.integer(settings$n_boot)
seed <- as.integer(settings$seed)
if (is.na(n_boot) || n_boot < 2L || is.na(seed)) stop("Invalid bootstrap settings.")

performance_file <- file.path(input_dir, "81_pooled_calibration_by_source.csv")
performance_md5 <- unname(tools::md5sum(performance_file))
d <- read_result(basename(performance_file), c(
  "source_id", "anchors", "RMSE_uncalibrated", "RMSE_calibrated",
  "n_test_compounds", "n_test_measurements", "n_available_anchors", "n_panels"
))
numeric_columns <- c("anchors", "RMSE_uncalibrated", "RMSE_calibrated",
  "n_test_compounds", "n_test_measurements", "n_available_anchors", "n_panels")
if (anyNA(d$source_id) || any(!nzchar(d$source_id)) ||
    anyDuplicated(d[c("source_id", "anchors")]) ||
    any(!vapply(d[numeric_columns], is.numeric, logical(1))) ||
    any(!is.finite(as.matrix(d[numeric_columns])))) {
  stop("Invalid source-level performance records.")
}
for (field in c("anchors", "n_test_compounds", "n_test_measurements", "n_available_anchors", "n_panels")) {
  if (any(d[[field]] < 1 | d[[field]] != floor(d[[field]]))) stop("Invalid counts in ", field)
}
if (any(d$n_test_measurements < d$n_test_compounds) ||
    any(d$n_available_anchors < d$anchors) ||
    any(d$n_panels != choose(d$n_available_anchors, d$anchors))) {
  stop("Counts do not match exhaustive calibration.")
}
for (source in unique(d$source_id)) {
  x <- d[d$source_id == source, ]
  if (length(unique(x$n_test_compounds)) != 1L ||
      length(unique(x$n_test_measurements)) != 1L ||
      !close_enough(x$RMSE_uncalibrated, rep(x$RMSE_uncalibrated[1], nrow(x)))) {
    stop("Test sets or baselines differ across anchor counts for ", source)
  }
}

one <- d[d$anchors == 1L, , drop = FALSE]
one <- one[order(one$source_id), ]
five <- d[d$anchors == 5L, , drop = FALSE]
if (!nrow(one) || !nrow(five)) stop("Both principal comparisons require eligible sources.")
j <- match(five$source_id, one$source_id)
if (anyNA(j)) stop("A five-anchor source lacks a one-anchor result.")
pair <- one[j, , drop = FALSE]
pair$RMSE_five <- five$RMSE_calibrated
pair <- pair[order(pair$source_id), ]

summarize <- function(x, endpoint, minimum, bootstrap_seed) {
  x <- x[x$n_test_compounds >= minimum, , drop = FALSE]
  n <- nrow(x)
  if (!n) stop("No sources remain for ", endpoint, " with minimum ", minimum)
  before <- if (endpoint == "one_anchor") x$RMSE_uncalibrated else x$RMSE_calibrated
  after <- if (endpoint == "one_anchor") x$RMSE_calibrated else x$RMSE_five
  effect <- function(ii) {
    if (endpoint == "one_anchor") {
      if (mean(before[ii]) <= 0) stop("Percentage reduction requires positive baseline RMSE.")
      100 * (1 - mean(after[ii]) / mean(before[ii]))
    } else mean(before[ii] - after[ii])
  }
  set.seed(bootstrap_seed)
  interval <- if (n >= 2L) {
    boot <- replicate(n_boot, effect(sample.int(n, n, replace = TRUE)))
    as.numeric(quantile(boot, c(0.025, 0.975), names = FALSE))
  } else c(NA_real_, NA_real_)
  data.frame(
    endpoint = endpoint, minimum_test_compounds = minimum,
    n_sources = n, n_test_compounds = sum(x$n_test_compounds),
    n_test_measurements = sum(x$n_test_measurements),
    mean_RMSE_before = mean(before), mean_RMSE_after = mean(after),
    effect = effect(seq_len(n)), ci_low = interval[1], ci_high = interval[2],
    effect_unit = if (endpoint == "one_anchor") "percent" else "log10 Kp units",
    n_sources_improved = sum(after < before), n_boot = n_boot,
    bootstrap_seed = bootstrap_seed, stringsAsFactors = FALSE
  )
}

results <- rbind(
  summarize(one, "one_anchor", 1L, seed + 1L),
  summarize(one, "one_anchor", minimum_test_compounds, seed + 3L),
  summarize(pair, "one_to_five", 1L, seed + 2L),
  summarize(pair, "one_to_five", minimum_test_compounds, seed + 4L)
)

# Confirm that the unrestricted summaries reproduce the saved principal effects.
all81 <- read_result("81_pooled_calibration_summary_all.csv",
  c("anchors", "n_sources", "RMSE_reduction_percent"))
paired81 <- read_result("81_pooled_one_vs_five_anchors.csv",
  c("n_sources", "mean_RMSE_decrease"))
z <- all81[all81$anchors == 1L, ]
if (nrow(z) != 1L || nrow(paired81) != 1L ||
    results$n_sources[1] != z$n_sources ||
    results$n_sources[3] != paired81$n_sources ||
    !close_enough(results$effect[c(1, 3)],
      c(z$RMSE_reduction_percent, paired81$mean_RMSE_decrease))) {
  stop("The unrestricted summaries disagree with 81. Check the input files.")
}
write_result(results, "86_test_set_sensitivity.csv")

membership <- rbind(
  data.frame(endpoint = "one_anchor", one[c("source_id", "n_test_compounds", "n_test_measurements")]),
  data.frame(endpoint = "one_to_five", pair[c("source_id", "n_test_compounds", "n_test_measurements")])
)
membership$included_in_restricted_summary <- membership$n_test_compounds >= minimum_test_compounds
write_result(membership, "86_evaluation_sources.csv")

fmt <- function(x, digits) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "")
digits <- ifelse(results$endpoint == "one_anchor", 1L, 3L)
effect_ci <- vapply(seq_len(nrow(results)), function(i) {
  paste0(fmt(results$effect[i], digits[i]), " [", fmt(results$ci_low[i], digits[i]),
    ", ", fmt(results$ci_high[i], digits[i]), "]")
}, character(1))
formatted <- data.frame(
  Comparison = ifelse(results$endpoint == "one_anchor", "One-anchor reduction (%)", "One-to-five decrease (log units)"),
  `Minimum test compounds` = results$minimum_test_compounds,
  Sources = results$n_sources,
  `RMSE before` = fmt(results$mean_RMSE_before, 3L),
  `RMSE after` = fmt(results$mean_RMSE_after, 3L),
  `Effect [95% CI]` = effect_ci, check.names = FALSE
)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(formatted, file.path(table_dir, "86_test_set_sensitivity.csv"), row.names = FALSE, na = "")
write_result(data.frame(
  input_file = normalizePath(performance_file, winslash = "/"),
  input_md5 = performance_md5, pooled_md5 = settings$input_md5,
  minimum_test_compounds = minimum_test_compounds, n_boot = n_boot,
  models_refitted = FALSE, panels_changed = FALSE
), "86_analysis_settings.csv")
writeLines(capture.output(sessionInfo()), file.path(input_dir, "86_sessionInfo.txt"))
print(results, row.names = FALSE)
cat("\nTraining models and anchor panels were unchanged.\n")
