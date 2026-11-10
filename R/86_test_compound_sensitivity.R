# Recalculate calibration benefits for publications with at least two test compounds.
result_dir <- "R/results"
input_file <- "data/processed/80_pooled.csv"
minimum_test_compounds <- 2L

read_result <- function(filename, required) {
	path <- file.path(result_dir, filename)
	if (!file.exists(path)) stop("Missing input: ", path)
	x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
	if (!nrow(x) || !all(required %in% names(x))) stop("Empty input or missing columns in ", filename)
	x
}

close_enough <- function(x, y) {
	length(x) == length(y) && all(is.finite(x) & is.finite(y) & abs(x - y) <= 1e-8 * pmax(1, abs(y)))
}

# Use saved exhaustive results for the current pooled dataset.
settings <- read_result("81_analysis_settings.csv", c("input_md5", "n_boot", "seed", "panel_method"))
if (nrow(settings) != 1L || !isTRUE(settings$panel_method == "exhaustive")) {
	stop("Run exhaustive script 81 first.")
}
if (!file.exists(input_file) || !identical(as.character(settings$input_md5), unname(tools::md5sum(input_file)))) {
	stop("The pooled input does not match 81. Regenerate calibration first.")
}
n_boot <- as.integer(settings$n_boot)
seed <- as.integer(settings$seed)
if (is.na(n_boot) || n_boot < 2L || is.na(seed)) stop("Invalid bootstrap settings.")

count_columns <- c("anchors", "n_test_compounds", "n_test_measurements", "n_available_anchors", "n_panels")
numeric_columns <- c("RMSE_uncalibrated", "RMSE_calibrated", count_columns)
d <- read_result("81_pooled_calibration_by_source.csv", c("source_id", numeric_columns))
if (anyNA(d$source_id) || any(!nzchar(d$source_id)) ||
	anyDuplicated(d[c("source_id", "anchors")]) ||
	any(!vapply(d[numeric_columns], is.numeric, logical(1))) ||
	any(!is.finite(as.matrix(d[numeric_columns])))) {
	stop("Invalid source-level performance records.")
}
for (field in count_columns) {
	if (any(d[[field]] < 1 | d[[field]] != floor(d[[field]]))) stop("Invalid counts in ", field)
}
if (any(d$n_test_measurements < d$n_test_compounds) || any(d$n_available_anchors < d$anchors) ||
	any(d$n_panels != choose(d$n_available_anchors, d$anchors))) {
	stop("Counts do not match exhaustive calibration.")
}
for (source in unique(d$source_id)) {
	x <- d[d$source_id == source, ]
	if (length(unique(x$n_test_compounds)) != 1L || length(unique(x$n_test_measurements)) != 1L ||
		!close_enough(x$RMSE_uncalibrated, rep(x$RMSE_uncalibrated[1], nrow(x)))) {
		stop("Test sets or baselines differ across anchor counts for ", source)
	}
}

# Match publications supporting both one and five anchors.
one <- d[d$anchors == 1L, , drop = FALSE]
one <- one[order(one$source_id), ]
five <- d[d$anchors == 5L, , drop = FALSE]
j <- match(five$source_id, one$source_id)
if (anyNA(j)) stop("A five-anchor source lacks a one-anchor result.")
pair <- one[j, , drop = FALSE]
pair$RMSE_five <- five$RMSE_calibrated
pair <- pair[order(pair$source_id), ]

summarize <- function(x, endpoint, minimum, bootstrap_seed) {
	x <- x[x$n_test_compounds >= minimum, , drop = FALSE]
	n <- nrow(x)
	if (!n) stop("No sources remain for ", endpoint, " with minimum ", minimum)
	one_anchor <- endpoint == "one_anchor"
	before <- if (one_anchor) x$RMSE_uncalibrated else x$RMSE_calibrated
	after <- if (one_anchor) x$RMSE_calibrated else x$RMSE_five
	effect <- function(ii) {
		if (!one_anchor) return(mean(before[ii] - after[ii]))
		baseline <- mean(before[ii])
		if (baseline <= 0) stop("Percentage reduction requires positive baseline RMSE.")
		100 * (1 - mean(after[ii]) / baseline)
	}
	set.seed(bootstrap_seed)
	interval <- if (n >= 2L) {
		boot <- replicate(n_boot, effect(sample.int(n, n, replace = TRUE)))
		as.numeric(quantile(boot, c(0.025, 0.975), names = FALSE))
	} else {
		c(NA_real_, NA_real_)
	}
	data.frame(
		endpoint = endpoint,
		minimum_test_compounds = minimum,
		n_sources = n,
		n_test_compounds = sum(x$n_test_compounds),
		n_test_measurements = sum(x$n_test_measurements),
		mean_RMSE_before = mean(before),
		mean_RMSE_after = mean(after),
		effect = effect(seq_len(n)),
		ci_low = interval[1],
		ci_high = interval[2],
		effect_unit = if (one_anchor) "percent" else "log10 Kp units",
		n_sources_improved = sum(after < before),
		n_boot = n_boot,
		bootstrap_seed = bootstrap_seed,
		stringsAsFactors = FALSE
	)
}

results <- rbind(
	summarize(one, "one_anchor", 1L, seed + 1L),
	summarize(one, "one_anchor", minimum_test_compounds, seed + 3L),
	summarize(pair, "one_to_five", 1L, seed + 2L),
	summarize(pair, "one_to_five", minimum_test_compounds, seed + 4L)
)

# Check that the unrestricted estimates agree with script 81.
all81 <- read_result("81_pooled_calibration_summary_all.csv", c("anchors", "n_sources", "RMSE_reduction_percent"))
paired81 <- read_result("81_pooled_one_vs_five_anchors.csv", c("n_sources", "mean_RMSE_decrease"))
z <- all81[all81$anchors == 1L, , drop = FALSE]
if (nrow(z) != 1L || nrow(paired81) != 1L ||
	results$n_sources[1] != z$n_sources || results$n_sources[3] != paired81$n_sources ||
	!close_enough(results$effect[c(1, 3)], c(z$RMSE_reduction_percent, paired81$mean_RMSE_decrease))) {
	stop("The unrestricted summaries disagree with 81. Check the input files.")
}

write.csv(results, file.path(result_dir, "86_test_set_sensitivity.csv"), row.names = FALSE, na = "")
message("Completed test-set sensitivity summaries for both calibration comparisons.")
