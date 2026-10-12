# ============================================================
# 60_uncertainty_analysis.R
# ============================================================

library(lme4)

ctrl <- lmerControl(
	optimizer = "bobyqa",
	optCtrl = list(maxfun = 100000)
)

set.seed(2026)

n_boot <- 2000
ci_probs <- c(0.025, 0.975)

rmse <- function(obs, pred) {
	sqrt(mean((obs - pred)^2))
}

ci <- function(x) {
	as.numeric(quantile(x, ci_probs, na.rm = TRUE))
}

# ------------------------------------------------------------
# 1. Compound-bootstrap uncertainty for LOCO performance
# ------------------------------------------------------------

pred <- read.csv(
	"R/results/50_source_blind_loco_predictions.csv",
	stringsAsFactors = FALSE
)

approaches <- c(
	"conventional",
	"hierarchical_source_blind",
	"hierarchical_known_source"
)

groups <- unique(pred[c("dataset", "specification")])
loco_out <- list()
z <- 1

for (g in seq_len(nrow(groups))) {
	d <- groups$dataset[g]
	s <- groups$specification[g]
	x <- pred[pred$dataset == d & pred$specification == s, ]
	compounds <- unique(x$compound_id)

	point_rmse <- sapply(approaches, function(a) rmse(x$observed, x[[a]]))
	point_improvement <- 100 * (1 - point_rmse / point_rmse["conventional"])

	boot_rmse <- matrix(NA_real_, n_boot, length(approaches))
	colnames(boot_rmse) <- approaches
	boot_improvement <- matrix(NA_real_, n_boot, length(approaches))
	colnames(boot_improvement) <- approaches

	rows_by_compound <- split(seq_len(nrow(x)), x$compound_id)

	for (b in seq_len(n_boot)) {
		sampled <- sample(compounds, length(compounds), replace = TRUE)
		idx <- unlist(rows_by_compound[sampled], use.names = FALSE)
		xb <- x[idx, ]

		for (a in approaches) boot_rmse[b, a] <- rmse(xb$observed, xb[[a]])
		boot_improvement[b, ] <- 100 * (1 - boot_rmse[b, ] / boot_rmse[b, "conventional"])
	}

	for (a in approaches) {
		rmse_ci <- ci(boot_rmse[, a])
		imp_ci <- ci(boot_improvement[, a])

		loco_out[[z]] <- data.frame(
			dataset = d,
			specification = s,
			approach = a,
			n_compounds = length(compounds),
			RMSE = unname(point_rmse[a]),
			RMSE_ci_low = rmse_ci[1],
			RMSE_ci_high = rmse_ci[2],
			improvement_percent = unname(point_improvement[a]),
			improvement_ci_low = imp_ci[1],
			improvement_ci_high = imp_ci[2]
		)
		z <- z + 1
	}
}

loco_ci <- do.call(rbind, loco_out)

write.csv(
	loco_ci,
	"R/results/60_loco_bootstrap_ci.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# 2. Parametric-bootstrap uncertainty for variance fractions
# ------------------------------------------------------------

variance_fraction <- function(fit) {
	vc <- as.data.frame(VarCorr(fit))
	v_compound <- vc$vcov[vc$grp == "compound_id"]
	v_source <- vc$vcov[vc$grp == "reference_id"]
	v_residual <- vc$vcov[vc$grp == "Residual"]
	total <- v_compound + v_source + v_residual

	c(
		Compound = 100 * v_compound / total,
		Source = 100 * v_source / total,
		Residual = 100 * v_residual / total
	)
}

model_files <- c(
	"Cheruvu" = "R/results/10_primary_hierarchical_model.rds",
	"HuSkinDB" = "R/results/31_primary_huskindb_model.rds"
)

missing_models <- model_files[!file.exists(model_files)]

if (length(missing_models)) {
	stop(
		"Missing model file(s): ",
		paste(missing_models, collapse = ", "),
		"\nSave the final REML HuSkin model from script 31 before running script 60."
	)
}

variance_out <- list()
z <- 1

for (d in names(model_files)) {
	fit <- readRDS(model_files[d])
	point <- variance_fraction(fit)

	boot <- bootMer(
		fit,
		FUN = variance_fraction,
		nsim = n_boot,
		type = "parametric",
		use.u = FALSE,
		seed = 2026
	)

	success <- apply(boot$t, 1, function(x) all(is.finite(x)))
	bt <- boot$t[success, , drop = FALSE]

	for (j in seq_along(point)) {
		vci <- ci(bt[, j])

		variance_out[[z]] <- data.frame(
			dataset = d,
			component = names(point)[j],
			estimate_percent = unname(point[j]),
			ci_low = vci[1],
			ci_high = vci[2],
			n_success = sum(success),
			n_boot = n_boot
		)
		z <- z + 1
	}
}

variance_ci <- do.call(rbind, variance_out)

write.csv(
	variance_ci,
	"R/results/60_variance_component_ci.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# 3. Source-bootstrap uncertainty for matched HuSkin calibration
# ------------------------------------------------------------

matched <- read.csv(
	"R/results/31_calibration_matched_reference.csv",
	stringsAsFactors = FALSE
)

matched_refs <- unique(matched$reference_id)
anchor_sizes <- sort(unique(matched$anchors))

cal_point <- aggregate(
	cbind(RMSE_uncalibrated, RMSE_calibrated) ~ anchors,
	data = matched,
	mean
)

cal_point$improvement_percent <- 100 * (
	1 - cal_point$RMSE_calibrated / cal_point$RMSE_uncalibrated
)

boot_cal <- matrix(NA_real_, n_boot, length(anchor_sizes))
colnames(boot_cal) <- anchor_sizes

for (b in seq_len(n_boot)) {
	sampled <- sample(matched_refs, length(matched_refs), replace = TRUE)

	for (j in seq_along(anchor_sizes)) {
		k <- anchor_sizes[j]
		x <- matched[matched$anchors == k, ]
		rows <- match(sampled, x$reference_id)
		xb <- x[rows, ]

		boot_cal[b, j] <- 100 * (
			1 - mean(xb$RMSE_calibrated) / mean(xb$RMSE_uncalibrated)
		)
	}
}

calibration_ci <- do.call(
	rbind,
	lapply(seq_along(anchor_sizes), function(j) {
		k <- anchor_sizes[j]
		bci <- ci(boot_cal[, j])
		p <- cal_point[cal_point$anchors == k, ]

		data.frame(
			anchors = k,
			n_references = length(matched_refs),
			RMSE_uncalibrated = p$RMSE_uncalibrated,
			RMSE_calibrated = p$RMSE_calibrated,
			improvement_percent = p$improvement_percent,
			improvement_ci_low = bci[1],
			improvement_ci_high = bci[2]
		)
	})
)

write.csv(
	calibration_ci,
	"R/results/60_huskin_calibration_ci.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# 4. Source-bootstrap uncertainty for wrong-source falsification
# ------------------------------------------------------------

wrong <- read.csv(
	"R/results/50_wrong_source_reference.csv",
	stringsAsFactors = FALSE
)

wrong_sizes <- sort(unique(wrong$anchors))
wrong_out <- list()

for (j in seq_along(wrong_sizes)) {
	k <- wrong_sizes[j]
	x <- wrong[wrong$anchors == k, ]
	refs <- unique(x$reference_id)

	point_uncal <- mean(x$RMSE_uncalibrated)
	point_correct <- mean(x$RMSE_correct)
	point_wrong <- mean(x$RMSE_wrong)
	point_correct_imp <- 100 * (1 - point_correct / point_uncal)
	point_wrong_imp <- 100 * (1 - point_wrong / point_uncal)
	point_difference <- point_wrong - point_correct

	boot_correct <- numeric(n_boot)
	boot_wrong <- numeric(n_boot)
	boot_difference <- numeric(n_boot)

	for (b in seq_len(n_boot)) {
		sampled <- sample(refs, length(refs), replace = TRUE)
		rows <- match(sampled, x$reference_id)
		xb <- x[rows, ]

		u <- mean(xb$RMSE_uncalibrated)
		correct <- mean(xb$RMSE_correct)
		wrong_rmse <- mean(xb$RMSE_wrong)

		boot_correct[b] <- 100 * (1 - correct / u)
		boot_wrong[b] <- 100 * (1 - wrong_rmse / u)
		boot_difference[b] <- wrong_rmse - correct
	}

	correct_ci <- ci(boot_correct)
	wrong_ci <- ci(boot_wrong)
	difference_ci <- ci(boot_difference)

	wrong_out[[j]] <- data.frame(
		anchors = k,
		n_references = length(refs),
		correct_improvement_percent = point_correct_imp,
		correct_ci_low = correct_ci[1],
		correct_ci_high = correct_ci[2],
		wrong_improvement_percent = point_wrong_imp,
		wrong_ci_low = wrong_ci[1],
		wrong_ci_high = wrong_ci[2],
		RMSE_wrong_minus_correct = point_difference,
		difference_ci_low = difference_ci[1],
		difference_ci_high = difference_ci[2]
	)
}

wrong_ci <- do.call(rbind, wrong_out)

write.csv(
	wrong_ci,
	"R/results/60_wrong_source_ci.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat(
	"\nUncertainty analysis complete\n\n",
	"LOCO compound bootstrap: ", n_boot, " replicates\n",
	"Variance-component parametric bootstrap: ", n_boot, " replicates\n",
	"HuSkin matched-source bootstrap: ", n_boot, " replicates\n",
	"Wrong-source bootstrap: ", n_boot, " replicates\n\n",
	"LOCO improvement intervals\n",
	paste(
		paste0(
			loco_ci$dataset, " ", loco_ci$specification, " ",
			loco_ci$approach, ": ",
			round(loco_ci$improvement_percent, 1), "% [",
			round(loco_ci$improvement_ci_low, 1), ", ",
			round(loco_ci$improvement_ci_high, 1), "]"
		),
		collapse = "\n"
	),
	"\n\nVariance-component intervals\n",
	paste(
		paste0(
			variance_ci$dataset, " ", variance_ci$component, ": ",
			round(variance_ci$estimate_percent, 1), "% [",
			round(variance_ci$ci_low, 1), ", ",
			round(variance_ci$ci_high, 1), "]"
		),
		collapse = "\n"
	),
	"\n\nMatched HuSkin calibration intervals\n",
	paste(
		paste0(
			calibration_ci$anchors, " anchors: ",
			round(calibration_ci$improvement_percent, 1), "% [",
			round(calibration_ci$improvement_ci_low, 1), ", ",
			round(calibration_ci$improvement_ci_high, 1), "]"
		),
		collapse = "\n"
	),
	"\n\nWrong-source falsification intervals\n",
	paste(
		paste0(
			wrong_ci$anchors, " anchors: matched ",
			round(wrong_ci$correct_improvement_percent, 1), "% [",
			round(wrong_ci$correct_ci_low, 1), ", ",
			round(wrong_ci$correct_ci_high, 1), "]; different-source ",
			round(wrong_ci$wrong_improvement_percent, 1), "% [",
			round(wrong_ci$wrong_ci_low, 1), ", ",
			round(wrong_ci$wrong_ci_high, 1), "]"
		),
		collapse = "\n"
	),
	"\n",
	sep = ""
)