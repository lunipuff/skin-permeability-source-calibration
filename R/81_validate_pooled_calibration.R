# Hold out each publication, fit on the others, and evaluate all eligible anchor panels.

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
	"input-file" = "data/processed/80_pooled.csv",
	"output-dir" = "R/results",
	"reuse-predictions" = "false"
))
input_file <- options_run[["input-file"]]
output_dir <- options_run[["output-dir"]]
reuse_option <- options_run[["reuse-predictions"]]
if (!reuse_option %in% c("true", "false")) stop("--reuse-predictions must be true or false.")
reuse_predictions <- reuse_option == "true"
if (!file.exists(input_file)) stop("Missing pooled dataset. Complete script 80 first.")
if (!dir.exists(output_dir)) stop("Required output directory is missing: ", output_dir)

anchor_sizes <- 1:5
panel_method <- "exhaustive"
n_boot <- 2000L
seed <- 2026L
input_md5 <- unname(tools::md5sum(input_file))
set.seed(seed)

read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
write_csv <- function(x, filename) {
	write.csv(x, file.path(output_dir, filename), row.names = FALSE, na = "")
}
bind_rows <- function(x) {
	result <- do.call(rbind, x)
	rownames(result) <- NULL
	result
}
metrics <- function(observed, predicted) {
	error <- observed - predicted
	c(RMSE = sqrt(mean(error^2)), MAE = mean(abs(error)))
}
bootstrap_ci <- function(values, n_sources) {
	if (n_sources < 2L) return(c(NA_real_, NA_real_))
	as.numeric(quantile(values, c(0.025, 0.975), names = FALSE))
}

# Validate input and define compound partitions for each held-out publication.
dat <- read_csv(input_file)
required <- c("record_id", "dataset_origin", "publication_id", "compound_id", "log_kp", "mw", "logp")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing input columns: ", paste(missing, collapse = ", "))
if (anyNA(dat$record_id) || anyDuplicated(dat$record_id)) stop("Missing or duplicate record_id values.")
if (anyNA(dat[c("publication_id", "compound_id")]) ||
	any(!nzchar(dat$publication_id)) || any(!nzchar(dat$compound_id))) {
	stop("Missing publication or compound identity.")
}
for (v in c("log_kp", "mw", "logp")) {
	if (!is.numeric(dat[[v]]) || any(!is.finite(dat[[v]]))) stop("Non-finite or nonnumeric values in ", v)
}
dat$source_id <- as.character(dat$publication_id)
dat$compound_id <- as.character(dat$compound_id)
dat <- dat[order(dat$source_id, dat$compound_id, dat$record_id), ]
rownames(dat) <- NULL

for (v in c("mw", "logp")) {
	spread <- tapply(dat[[v]], dat$compound_id, function(z) max(z) - min(z))
	if (any(spread > 1e-8)) stop("Inconsistent ", v, " values within compound_id.")
}

partition_source <- function(source) {
	held <- dat[dat$source_id == source, ]
	train <- dat[dat$source_id != source, ]
	held_ids <- unique(held$compound_id)
	train_ids <- unique(train$compound_id)
	list(held = held, train = train, shared = intersect(held_ids, train_ids), novel = setdiff(held_ids, train_ids))
}

sources <- sort(unique(dat$source_id))
if (!length(sources)) stop("The pooled dataset is empty.")
eligibility <- bind_rows(lapply(sources, function(source) {
	part <- partition_source(source)
	data.frame(
		source_id = source,
		dataset_origins = paste(sort(unique(part$held$dataset_origin)), collapse = "; "),
		n_measurements = nrow(part$held),
		n_compounds = length(unique(part$held$compound_id)),
		n_anchors = length(part$shared),
		n_test_compounds = length(part$novel),
		n_test_measurements = sum(part$held$compound_id %in% part$novel),
		stringsAsFactors = FALSE
	)
}))
for (k in anchor_sizes) {
	eligibility[[paste0("eligible_", k)]] <- eligibility$n_anchors >= k & eligibility$n_test_compounds >= 1L
}
write_csv(eligibility, "81_pooled_source_eligibility.csv")
eligible_sources <- eligibility$source_id[eligibility$eligible_1]
if (!length(eligible_sources)) stop("No publications have both anchors and novel test compounds.")

# Each unordered panel is evaluated once at each eligible size.
evaluate_panels <- function(source, anchor_residuals, test, variance) {
	anchor_residuals <- anchor_residuals[order(anchor_residuals$compound_id), , drop = FALSE]
	shared <- as.character(anchor_residuals$compound_id)
	if (!length(shared) || anyDuplicated(shared) || !nrow(test)) {
		stop("Invalid anchors or test measurements for ", source)
	}
	if (any(shared %in% test$compound_id)) stop("Anchor and test compounds overlap for ", source)
	if (any(!is.finite(anchor_residuals$mean_residual)) ||
		any(!is.finite(test$log_kp)) || any(!is.finite(test$uncalibrated_prediction)) ||
		any(!is.finite(variance)) || any(variance < 0) || sum(variance) <= 0) {
		stop("Invalid predictions or variance estimates for ", source)
	}

	baseline <- metrics(test$log_kp, test$uncalibrated_prediction)
	n_test_compounds <- length(unique(test$compound_id))
	results <- lapply(anchor_sizes[anchor_sizes <= length(shared)], function(k) {
		expected <- choose(length(shared), k)
		if (!is.finite(expected) || expected > .Machine$integer.max) {
			stop("Too many panels to enumerate for ", source, " at ", k, " anchors.")
		}
		combinations <- utils::combn(seq_along(shared), k)
		shrinkage <- unname(k * variance["source"] / (variance["residual"] + k * variance["source"]))
		bind_rows(lapply(seq_len(ncol(combinations)), function(panel) {
			selected <- combinations[, panel]
			offset_raw <- mean(anchor_residuals$mean_residual[selected])
			offset <- shrinkage * offset_raw
			calibrated <- metrics(test$log_kp, test$uncalibrated_prediction + offset)
			data.frame(
				source_id = source,
				anchors = k,
				draw = panel,
				selected_compounds = paste(shared[selected], collapse = " | "),
				n_available_anchors = length(shared),
				n_test_compounds = n_test_compounds,
				n_test_measurements = nrow(test),
				offset_raw = offset_raw,
				shrinkage = shrinkage,
				offset = offset,
				RMSE_uncalibrated = unname(baseline["RMSE"]),
				RMSE_calibrated = unname(calibrated["RMSE"]),
				MAE_uncalibrated = unname(baseline["MAE"]),
				MAE_calibrated = unname(calibrated["MAE"]),
				panel_method = panel_method,
				stringsAsFactors = FALSE
			)
		}))
	})
	bind_rows(results)
}

# Scale descriptors using unique training compounds only.
standardize <- function(train, held) {
	profile <- train[!duplicated(train$compound_id), c("mw", "logp"), drop = FALSE]
	center <- vapply(profile, mean, numeric(1))
	scale <- vapply(profile, sd, numeric(1))
	if (any(!is.finite(scale)) || any(scale <= 0)) stop("Insufficient descriptor variation in training data.")
	for (v in names(profile)) {
		name <- paste0(v, "_z")
		train[[name]] <- (train[[v]] - center[v]) / scale[v]
		held[[name]] <- (held[[v]] - center[v]) / scale[v]
	}
	list(train = train, held = held)
}

convergence_text <- function(fit) {
	result <- as.character(fit@optinfo$conv$lme4$messages)
	optimizer_code <- fit@optinfo$conv$opt
	if (length(optimizer_code) && any(optimizer_code != 0)) {
		result <- c(result, paste("Optimizer code:", paste(optimizer_code, collapse = ", ")))
	}
	paste(unique(result), collapse = " | ")
}

evaluate_source <- function(source) {
	part <- partition_source(source)
	sc <- standardize(part$train, part$held)
	fit <- lme4::lmer(
		log_kp ~ mw_z + logp_z + (1 | compound_id) + (1 | source_id),
		data = sc$train, REML = TRUE, control = control
	)
	vc <- as.data.frame(lme4::VarCorr(fit))
	variance <- c(source = vc$vcov[match("source_id", vc$grp)], residual = vc$vcov[match("Residual", vc$grp)])

	# Anchor expectations include known compound effects. Novel predictions use fixed effects.
	anchor_rows <- sc$held[sc$held$compound_id %in% part$shared, , drop = FALSE]
	anchor_rows$expected <- predict(fit, newdata = anchor_rows, re.form = ~(1 | compound_id), allow.new.levels = TRUE)
	anchor_residuals <- aggregate(
		anchor_rows$log_kp - anchor_rows$expected, list(compound_id = anchor_rows$compound_id), mean
	)
	names(anchor_residuals)[2] <- "mean_residual"
	anchor_counts <- table(anchor_rows$compound_id)
	anchor_residuals$n_measurements <- as.integer(anchor_counts[anchor_residuals$compound_id])
	anchor_residuals$source_id <- source

	test <- sc$held[sc$held$compound_id %in% part$novel, , drop = FALSE]
	test$uncalibrated_prediction <- predict(fit, newdata = test, re.form = NA, allow.new.levels = TRUE)
	list(
		draws = evaluate_panels(source, anchor_residuals, test, variance),
		anchors = anchor_residuals,
		tests = test[c("record_id", "dataset_origin", "source_id", "compound_id", "log_kp", "uncalibrated_prediction")],
		singular = lme4::isSingular(fit),
		convergence = convergence_text(fit),
		variance = variance
	)
}

# Reuse saved predictions or fit each eligible fold and record diagnostics.
draw_results <- list()
if (reuse_predictions) {
	read_cached <- function(filename, columns) {
		path <- file.path(output_dir, filename)
		if (!file.exists(path)) stop("Missing cached output ", filename, ". Run 81 without --reuse-predictions.")
		x <- read_csv(path)
		if (!nrow(x) || !all(columns %in% names(x))) stop("Invalid cached output ", filename)
		x
	}
	cached_settings <- read_cached("81_analysis_settings.csv", "input_md5")
	if (nrow(cached_settings) != 1L || !identical(as.character(cached_settings$input_md5), input_md5)) {
		stop("The pooled input has changed. Run 81 without --reuse-predictions first.")
	}
	diagnostics <- read_cached("81_pooled_fit_diagnostics.csv", c(
		"source_id", "success", "singular", "convergence", "warnings", "error", "source_variance", "residual_variance"
	))
	for (field in c("convergence", "warnings", "error")) diagnostics[[field]][is.na(diagnostics[[field]])] <- ""
	cached_anchors <- read_cached("81_pooled_anchor_residuals.csv", c(
		"source_id", "compound_id", "mean_residual", "n_measurements"
	))
	cached_tests <- read_cached("81_pooled_uncalibrated_test_predictions.csv", c(
		"source_id", "record_id", "compound_id", "log_kp", "uncalibrated_prediction"
	))
	if (anyDuplicated(cached_anchors[c("source_id", "compound_id")]) || anyDuplicated(cached_tests$record_id)) {
		stop("Duplicate cached records. Run 81 without --reuse-predictions.")
	}
	j <- match(cached_tests$record_id, dat$record_id)
	if (anyNA(j) || anyNA(cached_tests[c("source_id", "compound_id", "log_kp")]) ||
		any(cached_tests$source_id != dat$source_id[j]) || any(cached_tests$compound_id != dat$compound_id[j]) ||
		any(abs(cached_tests$log_kp - dat$log_kp[j]) > 1e-12)) {
		stop("Cached test records differ from the input. Run 81 without --reuse-predictions.")
	}
	if (anyDuplicated(diagnostics$source_id) || anyNA(diagnostics$success) || !all(diagnostics$success) ||
		!setequal(diagnostics$source_id, eligible_sources) ||
		!setequal(cached_anchors$source_id, eligible_sources) || !setequal(cached_tests$source_id, eligible_sources)) {
		stop("Cached folds are incomplete. Run 81 without --reuse-predictions.")
	}
	draw_results <- lapply(eligible_sources, function(source) {
		a <- cached_anchors[cached_anchors$source_id == source, , drop = FALSE]
		t <- cached_tests[cached_tests$source_id == source, , drop = FALSE]
		v <- diagnostics[diagnostics$source_id == source, ]
		part <- partition_source(source)
		test_ids <- part$held$record_id[part$held$compound_id %in% part$novel]
		if (!setequal(a$compound_id, part$shared) || !setequal(t$record_id, test_ids)) {
			stop("Cached anchor or test records differ for ", source)
		}
		evaluate_panels(source, a, t, c(source = v$source_variance, residual = v$residual_variance))
	})
	anchor_results <- list(cached_anchors)
	test_results <- list(cached_tests)
} else {
	control <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
	diagnostics <- list()
	anchor_results <- list()
	test_results <- list()
	for (i in seq_along(eligible_sources)) {
		source <- eligible_sources[i]
		warnings_seen <- character()
		message("[", i, "/", length(eligible_sources), "] ", source)
		result <- tryCatch(
			withCallingHandlers(evaluate_source(source), warning = function(w) {
				warnings_seen <<- c(warnings_seen, conditionMessage(w))
				invokeRestart("muffleWarning")
			}),
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
			source_variance = if (failed) NA_real_ else unname(result$variance["source"]),
			residual_variance = if (failed) NA_real_ else unname(result$variance["residual"]),
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
if (!length(draw_results)) stop("All eligible fits failed. Inspect 81_pooled_fit_diagnostics.csv.")

draws <- bind_rows(draw_results)
write_csv(draws, "81_pooled_calibration_draws.csv")
write_csv(bind_rows(anchor_results), "81_pooled_anchor_residuals.csv")
write_csv(bind_rows(test_results), "81_pooled_uncalibrated_test_predictions.csv")

# Average all panels within source and anchor count before pooling sources.
group_fields <- c("source_id", "anchors")
source_results <- aggregate(
	draws[c("RMSE_uncalibrated", "RMSE_calibrated", "MAE_uncalibrated", "MAE_calibrated", "offset_raw", "offset")],
	draws[group_fields], mean
)
source_counts <- unique(draws[c(group_fields, "n_available_anchors", "n_test_compounds", "n_test_measurements")])
panel_counts <- aggregate(list(n_panels = draws$draw), draws[group_fields], length)
source_results <- merge(source_results, source_counts, by = group_fields, sort = FALSE)
source_results <- merge(source_results, panel_counts, by = group_fields, sort = FALSE)
source_results$RMSE_reduction_percent <- with(source_results, 100 * (1 - RMSE_calibrated / RMSE_uncalibrated))
source_results$MAE_reduction_percent <- with(source_results, 100 * (1 - MAE_calibrated / MAE_uncalibrated))
source_results <- source_results[order(source_results$source_id, source_results$anchors), ]
write_csv(source_results, "81_pooled_calibration_by_source.csv")

# Bootstrap publications with equal source weights.
summarize_sources <- function(x, analysis) {
	bind_rows(lapply(sort(unique(x$anchors)), function(k) {
		z <- x[x$anchors == k, ]
		n <- nrow(z)
		baseline <- mean(z$RMSE_uncalibrated)
		calibrated <- mean(z$RMSE_calibrated)
		boot_reduction <- replicate(n_boot, {
			ii <- sample.int(n, n, replace = TRUE)
			100 * (1 - mean(z$RMSE_calibrated[ii]) / mean(z$RMSE_uncalibrated[ii]))
		})
		interval <- bootstrap_ci(boot_reduction, n)
		data.frame(
			analysis = analysis,
			anchors = k,
			n_sources = n,
			n_sources_improved = sum(z$RMSE_calibrated < z$RMSE_uncalibrated),
			mean_source_RMSE_uncalibrated = baseline,
			mean_source_RMSE_calibrated = calibrated,
			RMSE_reduction_percent = 100 * (1 - calibrated / baseline),
			reduction_ci_low = interval[1],
			reduction_ci_high = interval[2],
			mean_source_MAE_uncalibrated = mean(z$MAE_uncalibrated),
			mean_source_MAE_calibrated = mean(z$MAE_calibrated),
			n_boot = n_boot,
			stringsAsFactors = FALSE
		)
	}))
}

set.seed(seed + 1L)
all_summary <- summarize_sources(source_results, "All eligible sources at each anchor count")
write_csv(all_summary, "81_pooled_calibration_summary_all.csv")

# Compare one and five anchors within the same set of eligible publications.
matched_sources <- source_results$source_id[source_results$anchors == max(anchor_sizes)]
matched <- source_results[source_results$source_id %in% matched_sources, ]
matched_summary <- if (length(matched_sources)) {
	summarize_sources(matched, "Fixed sources eligible for all anchor counts")
} else {
	all_summary[FALSE, ]
}
write_csv(matched_summary, "81_pooled_calibration_summary_matched.csv")

if (length(matched_sources)) {
	one <- matched[matched$anchors == min(anchor_sizes), c("source_id", "RMSE_calibrated")]
	five <- matched[matched$anchors == max(anchor_sizes), c("source_id", "RMSE_calibrated")]
	comparison <- merge(one, five, by = "source_id", suffixes = c("_one", "_five"))
	gain <- comparison$RMSE_calibrated_one - comparison$RMSE_calibrated_five
	set.seed(seed + 2L)
	gain_boot <- replicate(n_boot, {
		ii <- sample.int(length(gain), length(gain), replace = TRUE)
		mean(gain[ii])
	})
	gain_ci <- bootstrap_ci(gain_boot, length(gain))
	paired_summary <- data.frame(
		n_sources = length(gain),
		mean_RMSE_one_anchor = mean(comparison$RMSE_calibrated_one),
		mean_RMSE_five_anchors = mean(comparison$RMSE_calibrated_five),
		mean_RMSE_decrease = mean(gain),
		decrease_ci_low = gain_ci[1],
		decrease_ci_high = gain_ci[2],
		n_boot = n_boot
	)
} else {
	paired_summary <- data.frame(
		n_sources = 0L,
		mean_RMSE_one_anchor = NA_real_,
		mean_RMSE_five_anchors = NA_real_,
		mean_RMSE_decrease = NA_real_,
		decrease_ci_low = NA_real_,
		decrease_ci_high = NA_real_,
		n_boot = n_boot
	)
}
write_csv(paired_summary, "81_pooled_one_vs_five_anchors.csv")

# Retain n_draws for readers of the settings file. Panels are exhaustive.
write_csv(data.frame(
	input_file = normalizePath(input_file, winslash = "/", mustWork = TRUE),
	input_md5 = input_md5,
	n_draws = 0L,
	n_boot = n_boot,
	seed = seed,
	panel_method = panel_method,
	total_panels = nrow(draws),
	reused_predictions = reuse_predictions
), "81_analysis_settings.csv")
message("Evaluated ", nrow(draws), " panels across ", sum(diagnostics$success), " of ",
	length(eligible_sources), " eligible publications.")
