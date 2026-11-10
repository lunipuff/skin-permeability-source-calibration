# Diagnose calibration errors from the predictions and panels saved by script 81.
result_dir <- "R/results"
tolerance <- 1e-10

read_result <- function(filename, required) {
	path <- file.path(result_dir, filename)
	if (!file.exists(path)) stop("Missing file: ", path)
	x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
	missing <- setdiff(required, names(x))
	if (length(missing)) stop("Missing columns in ", filename, ": ", paste(missing, collapse = ", "))
	if (!nrow(x) || anyNA(x[required])) stop("Empty input or missing required values in ", filename)
	numeric_fields <- vapply(x[required], is.numeric, logical(1))
	if (any(!is.finite(as.matrix(x[required[numeric_fields]])))) stop("Non-finite numeric inputs in ", filename)
	x$source_id <- as.character(x$source_id)
	x
}

write_result <- function(x, filename) {
	write.csv(x, file.path(result_dir, filename), row.names = FALSE, na = "")
}

close_enough <- function(x, y) {
	all(is.finite(x) & is.finite(y) & abs(x - y) <= 1e-7 * pmax(1, abs(y)))
}

draws <- read_result("81_pooled_calibration_draws.csv", c(
	"source_id", "anchors", "draw", "selected_compounds", "offset", "offset_raw",
	"shrinkage", "RMSE_uncalibrated", "RMSE_calibrated", "panel_method"
))
anchors <- read_result("81_pooled_anchor_residuals.csv", c(
	"source_id", "compound_id", "mean_residual", "n_measurements"
))
tests <- read_result("81_pooled_uncalibrated_test_predictions.csv", c(
	"source_id", "record_id", "compound_id", "log_kp", "uncalibrated_prediction"
))
primary <- read_result("81_pooled_calibration_by_source.csv", c(
	"source_id", "anchors", "RMSE_calibrated",
	"n_available_anchors", "n_test_compounds", "n_test_measurements"
))

# Check that the saved files describe the same exhaustive analysis.
if (any(draws$panel_method != "exhaustive")) stop("Run exhaustive script 81 before these diagnostics.")
for (keys in list(
	draws[c("source_id", "anchors", "draw")],
	anchors[c("source_id", "compound_id")],
	tests[c("source_id", "record_id")],
	primary[c("source_id", "anchors")]
)) {
	if (anyDuplicated(keys)) stop("Duplicate input identifiers.")
}
for (x in list(anchors, tests, primary)) {
	if (!setequal(draws$source_id, x$source_id)) stop("Input publications do not match.")
}
if (!close_enough(draws$offset, draws$offset_raw * draws$shrinkage)) {
	stop("Saved offsets do not match raw offsets and shrinkage.")
}
if (any(draws$anchors < 1 | draws$anchors != floor(draws$anchors)) ||
	any(draws$shrinkage < 0 | draws$shrinkage > 1)) stop("Invalid panel settings.")

# The mean test residual is the retrospectively optimal constant offset.
tests$uncalibrated_residual <- tests$log_kp - tests$uncalibrated_prediction
summaries <- list()
panels <- list()
for (source in sort(unique(draws$source_id))) {
	d <- draws[draws$source_id == source, , drop = FALSE]
	a <- anchors[anchors$source_id == source, , drop = FALSE]
	t <- tests[tests$source_id == source, , drop = FALSE]
	n_test_compounds <- length(unique(t$compound_id))
	e <- t$uncalibrated_residual
	mu <- mean(e)
	mse0 <- mean(e^2)
	rmse0 <- sqrt(mse0)
	residual_variance <- mean((e - mu)^2)
	bias_squared <- mu^2

	one <- d[d$anchors == 1L, , drop = FALSE]
	if (!nrow(one)) stop("Missing one-anchor panels: ", source)
	j <- match(one$selected_compounds, a$compound_id)
	if (anyNA(j) || !close_enough(one$offset_raw, a$mean_residual[j]) ||
		!close_enough(one$shrinkage, one$shrinkage[1])) {
		stop("Saved single-anchor choices and residuals disagree: ", source)
	}

	for (k in sort(unique(d$anchors))) {
		x <- d[d$anchors == k, , drop = FALSE]
		panel_sets <- strsplit(x$selected_compounds, " | ", fixed = TRUE)
		valid_sets <- vapply(panel_sets, function(ids) {
			length(ids) == k && !anyDuplicated(ids) && all(ids %in% a$compound_id)
		}, logical(1))
		panel_keys <- vapply(panel_sets, function(ids) paste(sort(ids), collapse = " | "), character(1))
		if (!all(valid_sets) || anyDuplicated(panel_keys) || nrow(x) != choose(nrow(a), k)) {
			stop("Incomplete or duplicated panel enumeration for ", source, " at ", k, " anchors. Rerun 81.")
		}

		p <- primary[primary$source_id == source & primary$anchors == k, , drop = FALSE]
		if (nrow(p) != 1L) stop("Primary summary mismatch: ", source)
		if (p$n_available_anchors != nrow(a) || p$n_test_compounds != n_test_compounds ||
			p$n_test_measurements != nrow(t)) stop("Input counts differ: ", source)

		# Calibrated MSE equals test residual variance plus squared offset mismatch.
		mse1 <- residual_variance + (mu - x$offset)^2
		rmse1_mean <- mean(x$RMSE_calibrated)
		if (!close_enough(x$RMSE_uncalibrated, rmse0) ||
			!close_enough(x$RMSE_calibrated, sqrt(mse1)) ||
			!close_enough(rmse1_mean, p$RMSE_calibrated)) {
			stop("Saved predictions and primary RMSE disagree: ", source)
		}

		mean_offset <- mean(x$offset)
		offset_variance <- mean((x$offset - mean_offset)^2)
		mismatch_squared <- (mean_offset - mu)^2
		mse1_mean <- mean(mse1)
		delta_mse <- mse1_mean - mse0
		if (!close_enough(delta_mse, offset_variance + mismatch_squared - bias_squared)) {
			stop("MSE decomposition failed: ", source)
		}

		x$test_mean_residual <- mu
		x$MSE_change <- mse1 - mse0
		x$RMSE_change <- x$RMSE_calibrated - rmse0
		x$improved <- x$RMSE_change < -tolerance
		x$worsened <- x$RMSE_change > tolerance
		x$offset_wrong_direction <- x$offset * mu < -tolerance
		x$offset_overshoot <- x$worsened & !x$offset_wrong_direction & abs(mu) > tolerance
		panels[[length(panels) + 1L]] <- x
		summaries[[length(summaries) + 1L]] <- data.frame(
			source_id = source,
			anchors = k,
			available_anchors = nrow(a),
			test_compounds = n_test_compounds,
			test_measurements = nrow(t),
			draws = nrow(x),
			distinct_panels = nrow(x),
			median_anchor_replication = median(a$n_measurements),
			median_test_replication = median(as.numeric(table(t$compound_id))),
			RMSE_uncalibrated = rmse0,
			RMSE_calibrated = rmse1_mean,
			RMSE_change = rmse1_mean - rmse0,
			RMSE_reduction_percent = if (rmse0 > 0) 100 * (1 - rmse1_mean / rmse0) else NA_real_,
			mean_RMSE_worsened = rmse1_mean > rmse0 + tolerance,
			mean_RMSE_improved = rmse1_mean < rmse0 - tolerance,
			draws_improved_percent = 100 * mean(x$improved),
			draws_wrong_direction_percent = 100 * mean(x$offset_wrong_direction),
			draws_overshoot_percent = 100 * mean(x$offset_overshoot),
			test_mean_residual = mu,
			test_bias_fraction_of_MSE = if (mse0 > 0) bias_squared / mse0 else NA_real_,
			offset_mean = mean_offset,
			offset_sd = sqrt(offset_variance),
			offset_min = min(x$offset),
			offset_max = max(x$offset),
			mean_offset_wrong_direction = mean_offset * mu < -tolerance,
			mean_shrinkage = mean(x$shrinkage),
			MSE_uncalibrated = mse0,
			mean_MSE_calibrated = mse1_mean,
			mean_MSE_change = delta_mse,
			offset_variance = offset_variance,
			offset_mismatch_squared = mismatch_squared,
			correctable_test_bias_squared = bias_squared,
			dominant_penalty = if (mismatch_squared >= offset_variance) "mean_offset_mismatch" else "offset_variability",
			retrospective_optimal_offset_RMSE = sqrt(residual_variance),
			stringsAsFactors = FALSE
		)
	}
}

summary84 <- do.call(rbind, summaries)
if (nrow(summary84) != nrow(primary)) stop("Source and panel summaries are incomplete.")
summary84 <- summary84[order(summary84$anchors, -summary84$RMSE_change), ]

# Compare one and five anchors within the same publications.
fields <- c(
	"source_id", "RMSE_calibrated", "offset_mean", "offset_variance",
	"offset_mismatch_squared", "draws_improved_percent", "mean_shrinkage"
)
paired84 <- merge(
	summary84[summary84$anchors == 1L, fields, drop = FALSE],
	summary84[summary84$anchors == 5L, fields, drop = FALSE],
	by = "source_id",
	suffixes = c("_one", "_five")
)
paired84$RMSE_decrease <- paired84$RMSE_calibrated_one - paired84$RMSE_calibrated_five
paired84$additional_anchors_helped <- paired84$RMSE_decrease > tolerance

# Summarize failure counts at each anchor size.
overview84 <- do.call(rbind, lapply(sort(unique(summary84$anchors)), function(k) {
	x <- summary84[summary84$anchors == k, , drop = FALSE]
	worse <- x$mean_RMSE_worsened
	data.frame(
		anchors = k,
		publications = nrow(x),
		publications_improved = sum(x$mean_RMSE_improved),
		publications_worsened = sum(worse),
		failures_with_mean_offset_wrong_direction = sum(worse & x$mean_offset_wrong_direction),
		failures_dominated_by_mean_mismatch = sum(worse & x$dominant_penalty == "mean_offset_mismatch"),
		failures_dominated_by_offset_variability = sum(worse & x$dominant_penalty == "offset_variability")
	)
}))

write_result(summary84, "84_failure_source_summary.csv")
write_result(do.call(rbind, panels), "84_failure_anchor_selections.csv")
write_result(paired84, "84_one_vs_five_diagnostics.csv")
write_result(overview84, "84_failure_overview.csv")
message("Summarized calibration diagnostics for ", length(unique(summary84$source_id)), " publications.")
