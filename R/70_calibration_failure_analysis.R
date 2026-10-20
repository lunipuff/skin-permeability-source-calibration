# ============================================================
# 70_calibration_failure_analysis.R
# Run from the project root after scripts 30, 40 and 41.
# Exploratory diagnostics. Original calibration outputs are unchanged.
# ============================================================

library(lme4)

set.seed(2026)
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

clean_metadata <- function(x) {
	x <- tolower(trimws(as.character(x)))
	x <- gsub("[[:space:]]+", " ", x)
	x[x %in% c("", "na", "n/a", "unknown", "not reported", "not specified", "-", "--")] <- NA_character_
	x
}

fit_training <- function(train) {
	p <- aggregate(train[c("mw", "log_kow")], list(compound_id = train$compound_id), median)
	center <- vapply(p[c("mw", "log_kow")], mean, numeric(1))
	spread <- vapply(p[c("mw", "log_kow")], sd, numeric(1))
	if (any(!is.finite(spread) | spread <= 0)) stop("Invalid training descriptor scale.")
	for (v in names(center)) train[[paste0(v, "_z")]] <- (train[[v]] - center[v]) / spread[v]
	train$compound_id <- factor(train$compound_id)
	train$reference_id <- factor(train$reference_id)
	fit <- lmer(log_kp ~ mw_z + log_kow_z + (1 | compound_id) + (1 | reference_id), data = train, REML = TRUE, control = ctrl)
	vc <- as.data.frame(VarCorr(fit))
	list(fit = fit, center = center, spread = spread, source_variance = vc$vcov[vc$grp == "reference_id"], residual_variance = vc$vcov[vc$grp == "Residual"], singular = isSingular(fit), convergence = paste(unlist(fit@optinfo$conv$lme4$messages), collapse = " | "))
}

prepare_predictions <- function(x, model, anchor = FALSE) {
	for (v in names(model$center)) x[[paste0(v, "_z")]] <- (x[[v]] - model$center[v]) / model$spread[v]
	x$compound_id <- x$model_compound_id
	x$expected <- as.numeric(predict(model$fit, newdata = x, re.form = if (anchor) ~(1 | compound_id) else NA, allow.new.levels = TRUE))
	x$residual <- x$log_kp - x$expected
	x
}

# Examine recorded heterogeneity only. Different values do not establish different protocols.
protocol_fields <- c("skin_source_type", "skin_source_site", "used_layer", "skin_preparation", "cell_type", "donor_type", "acceptor_type", "temperature_c", "acceptor_temperature_c", "donor_ph", "acceptor_ph")

analyze_source <- function(train, held, anchors, tests, dataset, model = NULL) {
	if (is.null(model)) model <- fit_training(train)
	a <- prepare_predictions(anchors, model, TRUE)
	t <- prepare_predictions(tests, model)
	ap <- aggregate(a[c("residual", "mw_z", "log_kow_z")], list(anchor_id = a$model_compound_id), mean)
	tp <- aggregate(t[c("mw_z", "log_kow_z")], list(test_id = t$original_compound_id), mean)
	mu <- mean(t$residual)
	rmse0 <- sqrt(mean(t$residual^2))
	residual_rms <- sqrt(mean((t$residual - mu)^2))
	anchor_n <- table(a$model_compound_id)
	test_n <- table(t$original_compound_id)
	within_sd <- vapply(split(held$log_kp, held$original_compound_id), function(v) if (length(v) > 1) sd(v) else NA_real_, numeric(1))
	within_sd <- within_sd[is.finite(within_sd)]
	protocol <- do.call(rbind, lapply(protocol_fields, function(field) {
		v <- if (field %in% names(held)) clean_metadata(held[[field]]) else rep(NA_character_, nrow(held))
		data.frame(dataset = dataset,
				reference_id = held$reference_id[1],
				field = field,
				reported = sum(!is.na(v)),
				records = length(v),
				distinct_reported_values = length(unique(v[!is.na(v)])),
				values = paste(sort(unique(v[!is.na(v)])), collapse = " | "))
	}))
	distance <- vapply(seq_len(nrow(tp)), function(i) sqrt((ap$mw_z - tp$mw_z[i])^2 + (ap$log_kow_z - tp$log_kow_z[i])^2), numeric(nrow(ap)))
	distance <- matrix(distance, nrow = nrow(ap), ncol = nrow(tp))
	selection_rows <- list()
	for (k in seq_len(min(5, nrow(ap)))) {
		# Enumerate all one-anchor choices. Enumerate larger panels if at most 100 exist.
		exhaustive <- lchoose(nrow(ap), k) <= log(100)
		panels <- if (k == 1 || exhaustive) combn(seq_len(nrow(ap)), k, simplify = FALSE) else replicate(100, sample.int(nrow(ap), k), simplify = FALSE)
		for (j in seq_along(panels)) {
			selected <- panels[[j]]
			raw <- mean(ap$residual[selected])
			shrinkage <- k * model$source_variance / (model$residual_variance + k * model$source_variance)
			u <- raw * shrinkage
			rmse1 <- sqrt(mean((t$residual - u)^2))
			# Nearest selected anchor per test compound, averaged equally across test compounds.
			d <- mean(apply(distance[selected, , drop = FALSE], 2, min))
			category <- if (abs(u) < 1e-10) "No adjustment" else if (u * mu <= 0) "Opposite direction" else if (abs(u) >= 2 * abs(mu)) "Overshoot" else "Reduces squared error"
			selection_rows[[length(selection_rows) + 1L]] <- data.frame(
				dataset = dataset,
				reference_id = held$reference_id[1],
				anchors = k, selection = j,
				exhaustive = k == 1 || exhaustive,
				anchor_ids = paste(ap$anchor_id[selected], collapse = " | "),
				offset_raw = raw,
				shrinkage = shrinkage,
				offset = u,
				test_mean_residual = mu,
				test_residual_rms = residual_rms,
				RMSE_uncalibrated = rmse0,
				RMSE_calibrated = rmse1,
				RMSE_reduction = rmse0 - rmse1,
				improved = rmse1 < rmse0 - 1e-10,
				anchor_distance = d,
				category = category)
		}
	}
	selection <- do.call(rbind, selection_rows)
	summary <- do.call(rbind, lapply(sort(unique(selection$anchors)), function(k) {
		x <- selection[selection$anchors == k, ]
		data.frame(dataset = dataset,
				reference_id = held$reference_id[1],
				anchors = k,
				available_anchors = nrow(ap),
				test_compounds = length(test_n),
				anchor_measurements = nrow(a),
				test_measurements = nrow(t),
				median_anchor_replication = median(as.numeric(anchor_n)),
				median_test_replication = median(as.numeric(test_n)),
				replicated_compounds = length(within_sd),
				median_within_compound_sd = if (length(within_sd)) median(within_sd) else NA_real_,
				protocol_fields_reported = sum(protocol$reported > 0),
				protocol_fields_varying = sum(protocol$distinct_reported_values > 1),
				selections = nrow(x),
				exhaustive = all(x$exhaustive),
				RMSE_uncalibrated = rmse0,
				RMSE_calibrated = mean(x$RMSE_calibrated),
				RMSE_reduction_percent = 100 * (1 - mean(x$RMSE_calibrated) / rmse0),
				test_mean_residual = mu,
				test_residual_rms = residual_rms,
				offset_mean = mean(x$offset),
				offset_sd = if (nrow(x) > 1) sd(x$offset) else 0,
				offset_min = min(x$offset),
				offset_max = max(x$offset),
				offset_bias_squared = (mean(x$offset) - mu)^2,
				offset_variance = mean((x$offset - mean(x$offset))^2),
				improved_percent = 100 * mean(x$improved),
				anchor_distance_mean = mean(x$anchor_distance),
				distance_improvement_rho = if (length(unique(x$anchor_distance)) > 1 && length(unique(x$RMSE_reduction)) > 1) cor(x$anchor_distance, x$RMSE_reduction, method = "spearman") else NA_real_,
				fit_singular = model$singular,
				fit_convergence = model$convergence)
	}))
	# Exact squared-error identity, using the same measurement weights as evaluation.
	stopifnot(max(abs(selection$RMSE_calibrated^2 - (residual_rms^2 + (mu - selection$offset)^2))) < 1e-7)
	ap$dataset <- dataset
	ap$reference_id <- held$reference_id[1]
	t$dataset <- dataset
	list(summary = summary, selections = selection, protocols = protocol, anchors = ap, tests = t)
}

huskin <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)
huskin$original_compound_id <- huskin$compound_id
huskin$model_compound_id <- huskin$compound_id
outputs <- list()
for (ref in unique(huskin$reference_id)) {
	train <- huskin[huskin$reference_id != ref, ]
	held <- huskin[huskin$reference_id == ref, ]
	shared <- held$compound_id %in% train$compound_id
	if (!any(shared) || all(shared)) next
	cat("HuSkinDB: ", ref, "\n", sep = "")
	outputs[[length(outputs) + 1L]] <- analyze_source(train, held, held[shared, ], held[!shared, ], "HuSkinDB")
}

hu_rdkit <- read.csv("data/processed/40_huskindb_rdkit.csv", stringsAsFactors = FALSE)
hu_rdkit$mw <- hu_rdkit$mw_rdkit
hu_rdkit$log_kow <- hu_rdkit$logp_rdkit

sp <- read.csv("data/processed/40_skinpix.csv", stringsAsFactors = FALSE)
sp$mw <- sp$mw_rdkit
sp$log_kow <- sp$logp_rdkit
sp$skin_source_type <- sp$skin_donor_type
sp$original_compound_id <- sp$compound_id
sp$model_compound_id <- ifelse(sp$overlap_huskin, sp$huskin_compound_id, paste0("NOVEL_", sp$compound_id))

external_model <- fit_training(hu_rdkit)

for (ref in unique(sp$reference_id)) {
	held <- sp[sp$reference_id == ref, ]
	shared <- held$overlap_huskin %in% TRUE & !is.na(held$huskin_compound_id)
	if (!any(shared) || !any(!held$overlap_huskin) || any(held$doi_overlaps_huskin)) next
	cat("SkinPiX: ", ref, "\n", sep = "")
	outputs[[length(outputs) + 1L]] <- analyze_source(hu_rdkit, held, held[shared, ], held[!held$overlap_huskin, ], "SkinPiX", external_model)
}

summary70 <- do.call(rbind, lapply(outputs, function(x) x$summary))
selection70 <- do.call(rbind, lapply(outputs, function(x) x$selections))
protocol70 <- do.call(rbind, lapply(outputs, function(x) x$protocols))
anchor70 <- do.call(rbind, lapply(outputs, function(x) x$anchors))
test70 <- do.call(rbind, lapply(outputs, function(x) x$tests[c("dataset", "reference_id", "source_row", "original_compound_id", "compound_name", "log_kp", "expected", "residual", "mw_z", "log_kow_z")]))

write.csv(summary70, "R/results/70_failure_source_summary.csv", row.names = FALSE, na = "")
write.csv(selection70, "R/results/70_failure_anchor_selections.csv", row.names = FALSE, na = "")
write.csv(protocol70, "R/results/70_failure_protocol_variation.csv", row.names = FALSE, na = "")
write.csv(anchor70, "R/results/70_failure_anchor_residuals.csv", row.names = FALSE, na = "")
write.csv(test70, "R/results/70_failure_test_residuals.csv", row.names = FALSE, na = "")


