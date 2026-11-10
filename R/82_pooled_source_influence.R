# Omit one evaluation publication at a time from the saved calibration results.
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
	"input-file" = "R/results/81_pooled_calibration_by_source.csv",
	"output-dir" = "R/results"
))
input <- options_run[["input-file"]]
output_dir <- options_run[["output-dir"]]
if (!file.exists(input)) stop("Missing calibration results. Run 81_validate_pooled_calibration.R first.")
if (!dir.exists(output_dir)) stop("Required output directory is missing: ", output_dir)

n_boot <- 20000L
seed <- 2026L
set.seed(seed)
dat <- read.csv(input, stringsAsFactors = FALSE)
required <- c("source_id", "anchors", "RMSE_uncalibrated", "RMSE_calibrated")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing input columns: ", paste(missing, collapse = ", "))
if (anyNA(dat[c("source_id", "anchors")]) || anyDuplicated(dat[c("source_id", "anchors")])) {
	stop("Expected one row per source and anchor count.")
}
if (any(!is.finite(dat$RMSE_uncalibrated)) || any(!is.finite(dat$RMSE_calibrated))) {
	stop("Non-finite RMSE values.")
}

# Compare one anchor with no calibration and five anchors with one on matched sources.
one <- dat[dat$anchors == 1L, ]
five <- dat[dat$anchors == 5L, ]
one_comparison <- data.frame(source_id = one$source_id, before = one$RMSE_uncalibrated, after = one$RMSE_calibrated)
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

rmse_effect <- function(before, after, percentage) {
	if (percentage) 100 * (1 - after / before) else before - after
}

estimate <- function(x, percentage) {
	n <- nrow(x)
	effect <- rmse_effect(mean(x$before), mean(x$after), percentage)
	indices <- matrix(sample.int(n, n * n_boot, replace = TRUE), nrow = n_boot)
	boot_before <- rowMeans(matrix(x$before[indices], nrow = n_boot))
	boot_after <- rowMeans(matrix(x$after[indices], nrow = n_boot))
	boot_effect <- rmse_effect(boot_before, boot_after, percentage)
	interval <- quantile(boot_effect, c(0.025, 0.975), names = FALSE)
	c(effect = effect, ci_low = interval[1], ci_high = interval[2])
}

run_influence <- function(x, analysis, percentage) {
	if (nrow(x) < 3L) stop("At least three sources are needed for ", analysis)
	x <- x[order(x$source_id), ]
	# Bootstrap the full sample first to preserve the seeded resampling sequence.
	full <- estimate(x, percentage)[["effect"]]
	omissions <- vapply(seq_len(nrow(x)), function(i) {
		estimate(x[-i, , drop = FALSE], percentage)
	}, c(effect = 0, ci_low = 0, ci_high = 0))
	effects <- omissions["effect", ]
	change <- abs(effects - full)
	data.frame(
		analysis = analysis,
		full_effect = full,
		minimum_omission_effect = min(effects),
		maximum_omission_effect = max(effects),
		largest_absolute_change = max(change),
		most_influential_source = x$source_id[which.max(change)],
		n_direction_reversals = sum(effects * full < 0),
		n_intervals_including_zero = sum(omissions["ci_low", ] <= 0 & omissions["ci_high", ] >= 0),
		n_omissions = nrow(x),
		stringsAsFactors = FALSE
	)
}

summary_results <- rbind(
	run_influence(one_comparison, "One-anchor RMSE reduction (%)", percentage = TRUE),
	run_influence(panel_comparison, "RMSE decrease from one to five anchors (log units)", percentage = FALSE)
)
write.csv(summary_results, file.path(output_dir, "82_source_influence_summary.csv"), row.names = FALSE)
write.csv(data.frame(
	input_file = normalizePath(input, winslash = "/", mustWork = TRUE),
	input_md5 = unname(tools::md5sum(input)),
	n_boot = n_boot,
	seed = seed
), file.path(output_dir, "82_analysis_settings.csv"), row.names = FALSE)
message("Summarized publication omissions for both calibration comparisons.")
