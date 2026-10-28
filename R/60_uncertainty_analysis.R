# 60_uncertainty_analysis.R

library(lme4)

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
set.seed(2026)
n_boot <- 2000
ci_probs <- c(0.025, 0.975)

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))

ci <- function(x) as.numeric(quantile(x, ci_probs, na.rm = TRUE))

# Compound-bootstrap uncertainty for LOCO performance

pred <- read.csv("R/results/50_source_blind_loco_predictions.csv", stringsAsFactors = FALSE)
approaches <- c("conventional", "hierarchical_source_blind", "hierarchical_known_source")
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
	boot_rmse <- matrix(NA_real_, n_boot, length(approaches), dimnames = list(NULL, approaches))
	boot_improvement <- matrix(NA_real_, n_boot, length(approaches), dimnames = list(NULL, approaches))
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

write.csv(loco_ci, "R/results/60_loco_bootstrap_ci.csv", row.names = FALSE)

# Parametric-bootstrap uncertainty for variance fractions

variance_fraction <- function(fit) {
	vc <- as.data.frame(VarCorr(fit))
	v_compound <- vc$vcov[vc$grp == "compound_id"]
	v_source <- vc$vcov[vc$grp == "reference_id"]
	v_residual <- vc$vcov[vc$grp == "Residual"]
	total <- v_compound + v_source + v_residual
	c(Compound = 100 * v_compound / total, Source = 100 * v_source / total, Residual = 100 * v_residual / total)
}

model_files <- c(
	Cheruvu = "R/results/10_primary_hierarchical_model.rds",
	HuSkinDB = "R/results/31_primary_huskindb_model.rds"
)

missing_models <- model_files[!file.exists(model_files)]

if (length(missing_models)) stop("Missing model file(s): ", paste(missing_models, collapse = ", "))

variance_out <- list()
z <- 1

bootstrap_log <- "R/results/60_variance_bootstrap_diagnostics.txt"
writeLines("Variance-component bootstrap diagnostics", bootstrap_log)

for (d in names(model_files)) {
	fit <- readRDS(model_files[d])
	point <- variance_fraction(fit)
	boot <- bootMer(fit, FUN = variance_fraction, nsim = n_boot, type = "parametric", use.u = FALSE, seed = 2026)
	success <- apply(boot$t, 1, function(x) all(is.finite(x)))
	n_success <- sum(success)
	n_failed <- nrow(boot$t) - n_success

	diagnostics <- c(
		paste0("\nDataset: ", d),
		paste("Requested replicates:", n_boot),
		paste("Usable replicates:", n_success),
		paste("Excluded non-finite replicates:", n_failed),
			"Model-fitting errors:",
		capture.output(print(attr(boot, "bootFail"))),
			"Error messages:",
		capture.output(print(attr(boot, "boot.fail.msgs"))),
			"All messages, warnings, and errors:",
		capture.output(print(attr(boot, "boot.all.msgs")))
	)

	cat(
		paste(diagnostics, collapse = "\n"), "\n",
			file = bootstrap_log,
			append = TRUE
	)

	if (n_success < 2L) {
		stop(
			"Too few usable bootstrap replicates for ", d,
			". See ", bootstrap_log
		)
	}

	bt <- boot$t[success, , drop = FALSE]
	for (j in seq_along(point)) {
		vci <- ci(bt[, j])
		variance_out[[z]] <- data.frame(
			dataset = d,
			component = names(point)[j],
			estimate_percent = unname(point[j]),
			ci_low = vci[1],
			ci_high = vci[2],
			n_boot = n_boot,
			n_boot_success = n_success,
			n_boot_failed = n_failed
		)
		z <- z + 1
	}
}

variance_ci <- do.call(rbind, variance_out)

write.csv(variance_ci, "R/results/60_variance_component_ci.csv", row.names = FALSE)

# Source-bootstrap uncertainty for matched HuSkinDB calibration

matched <- read.csv("R/results/31_calibration_matched_reference.csv", stringsAsFactors = FALSE)

matched_refs <- unique(matched$reference_id)
anchor_sizes <- sort(unique(matched$anchors))

cal_point <- aggregate(cbind(RMSE_uncalibrated, RMSE_calibrated) ~ anchors, data = matched, mean)
cal_point$improvement_percent <- 100 * (1 - cal_point$RMSE_calibrated / cal_point$RMSE_uncalibrated)

boot_cal <- matrix(NA_real_, n_boot, length(anchor_sizes), dimnames = list(NULL, anchor_sizes))

for (b in seq_len(n_boot)) {
	sampled <- sample(matched_refs, length(matched_refs), replace = TRUE)
	for (j in seq_along(anchor_sizes)) {
		x <- matched[matched$anchors == anchor_sizes[j], ]
		xb <- x[match(sampled, x$reference_id), ]
		boot_cal[b, j] <- 100 * (1 - mean(xb$RMSE_calibrated) / mean(xb$RMSE_uncalibrated))
	}
}

calibration_ci <- do.call(rbind, lapply(seq_along(anchor_sizes), function(j) {
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
}))

write.csv(calibration_ci, "R/results/60_huskin_calibration_ci.csv", row.names = FALSE)

# Source-bootstrap uncertainty for wrong-source control

wrong <- read.csv("R/results/50_wrong_source_reference.csv", stringsAsFactors = FALSE)

wrong_sizes <- sort(unique(wrong$anchors))
wrong_out <- vector("list", length(wrong_sizes))

for (j in seq_along(wrong_sizes)) {
	k <- wrong_sizes[j]
	x <- wrong[wrong$anchors == k, ]
	refs <- unique(x$reference_id)
	point_uncal <- mean(x$RMSE_uncalibrated)
	point_correct <- mean(x$RMSE_correct)
	point_wrong <- mean(x$RMSE_wrong)
	boot_correct <- numeric(n_boot)
	boot_wrong <- numeric(n_boot)
	boot_difference <- numeric(n_boot)
	for (b in seq_len(n_boot)) {
		sampled <- sample(refs, length(refs), replace = TRUE)
		xb <- x[match(sampled, x$reference_id), ]
		u <- mean(xb$RMSE_uncalibrated)
		correct <- mean(xb$RMSE_correct)
		wrong_rmse <- mean(xb$RMSE_wrong)
		boot_correct[b] <- 100 * (1 - correct / u)
		boot_wrong[b] <- 100 * (1 - wrong_rmse / u)
		boot_difference[b] <- wrong_rmse - correct
	}
	correct_ci <- ci(boot_correct)
	wrong_imp_ci <- ci(boot_wrong)
	difference_ci <- ci(boot_difference)
	wrong_out[[j]] <- data.frame(
		anchors = k,
		n_references = length(refs),
		correct_improvement_percent = 100 * (1 - point_correct / point_uncal),
		correct_ci_low = correct_ci[1],
		correct_ci_high = correct_ci[2],
		wrong_improvement_percent = 100 * (1 - point_wrong / point_uncal),
		wrong_ci_low = wrong_imp_ci[1],
		wrong_ci_high = wrong_imp_ci[2],
		RMSE_wrong_minus_correct = point_wrong - point_correct,
		difference_ci_low = difference_ci[1],
		difference_ci_high = difference_ci[2]
	)
}

wrong_ci <- do.call(rbind, wrong_out)

write.csv(wrong_ci, "R/results/60_wrong_source_ci.csv", row.names = FALSE)

# Paired source bootstrap for conventional versus hierarchical calibration

comp32 <- read.csv("R/results/32_calibration_method_reference.csv", stringsAsFactors = FALSE)

comp1 <- comp32[comp32$anchors == 1, ]

conv1 <- comp1[comp1$method == "Conventional", c("reference_id", "RMSE_uncalibrated", "RMSE_calibrated")]
hier1 <- comp1[comp1$method == "Hierarchical", c("reference_id", "RMSE_uncalibrated", "RMSE_calibrated")]

names(conv1)[2:3] <- c("RMSE_uncalibrated_conventional", "RMSE_calibrated_conventional")
names(hier1)[2:3] <- c("RMSE_uncalibrated_hierarchical", "RMSE_calibrated_hierarchical")

comp1_pair <- merge(conv1, hier1, by = "reference_id", sort = FALSE)
n_comp1 <- nrow(comp1_pair)

comp1_conv_imp <- 100 * (1 - mean(comp1_pair$RMSE_calibrated_conventional) / mean(comp1_pair$RMSE_uncalibrated_conventional))

comp1_hier_imp <- 100 * (1 - mean(comp1_pair$RMSE_calibrated_hierarchical) / mean(comp1_pair$RMSE_uncalibrated_hierarchical))

comp1_diff <- mean(comp1_pair$RMSE_calibrated_conventional) - mean(comp1_pair$RMSE_calibrated_hierarchical)

comp1_boot <- replicate(n_boot, {
	i <- sample(seq_len(n_comp1), n_comp1, replace = TRUE)
	conv_uncal <- mean(comp1_pair$RMSE_uncalibrated_conventional[i])
	conv_cal <- mean(comp1_pair$RMSE_calibrated_conventional[i])
	hier_uncal <- mean(comp1_pair$RMSE_uncalibrated_hierarchical[i])
	hier_cal <- mean(comp1_pair$RMSE_calibrated_hierarchical[i])
	c(
		conventional_improvement = 100 * (1 - conv_cal / conv_uncal),
		hierarchical_improvement = 100 * (1 - hier_cal / hier_uncal),
		RMSE_conventional_minus_hierarchical = conv_cal - hier_cal
	)
})

comp1_ci <- apply(comp1_boot, 1, quantile, probs = ci_probs, names = FALSE)

comp1_out <- data.frame(
	RMSE_conventional_minus_hierarchical = comp1_diff,
	difference_ci_low = comp1_ci[1, "RMSE_conventional_minus_hierarchical"],
	difference_ci_high = comp1_ci[2, "RMSE_conventional_minus_hierarchical"]
)

write.csv(comp1_out, "R/results/60_calibration_method_comparison_ci.csv", row.names = FALSE)

# One-anchor deployment estimate across all eligible HuSkinDB sources

cal31 <- read.csv("R/results/31_calibration_performance.csv", stringsAsFactors = FALSE)

all1 <- cal31[cal31$task == "New compound" & cal31$anchors == 1, ]
all1_ref <- aggregate(cbind(RMSE_uncalibrated, RMSE_calibrated) ~ reference_id, data = all1, mean)

n_all1 <- nrow(all1_ref)
n_all1_improved <- sum(all1_ref$RMSE_calibrated < all1_ref$RMSE_uncalibrated)

all1_rmse_uncalibrated <- mean(all1_ref$RMSE_uncalibrated)
all1_rmse_calibrated <- mean(all1_ref$RMSE_calibrated)
all1_improvement <- 100 * (1 - all1_rmse_calibrated / all1_rmse_uncalibrated)

all1_boot <- replicate(n_boot, {
	i <- sample(seq_len(n_all1), n_all1, replace = TRUE)
	100 * (1 - mean(all1_ref$RMSE_calibrated[i]) / mean(all1_ref$RMSE_uncalibrated[i]))
})

all1_ci <- quantile(all1_boot, probs = ci_probs, names = FALSE)

all1_ci_out <- data.frame(
	anchors = 1,
	n_references = n_all1,
	n_references_improved = n_all1_improved,
	RMSE_uncalibrated = all1_rmse_uncalibrated,
	RMSE_calibrated = all1_rmse_calibrated,
	improvement_percent = all1_improvement,
	improvement_ci_low = all1_ci[1],
	improvement_ci_high = all1_ci[2]
)

write.csv(all1_ci_out, "R/results/60_huskin_one_anchor_all_eligible_ci.csv", row.names = FALSE)

cat(
	"\nUncertainty analysis complete\n",
	"Bootstrap repetitions: ", n_boot, "\n",
	"\nLOCO improvement intervals\n",
	paste(paste0(loco_ci$dataset, " ", loco_ci$specification, " ", loco_ci$approach, ": ", round(loco_ci$improvement_percent, 1), "% [", round(loco_ci$improvement_ci_low, 1), ", ", round(loco_ci$improvement_ci_high, 1), "]"), collapse = "\n"),
	"\n\nVariance-component intervals\n",
	paste(paste0(variance_ci$dataset, " ", variance_ci$component, ": ", round(variance_ci$estimate_percent, 1), "% [", round(variance_ci$ci_low, 1), ", ", round(variance_ci$ci_high, 1), "]"), collapse = "\n"),
	"\n\nMatched HuSkin calibration intervals\n",
	paste(paste0(calibration_ci$anchors, " anchors: ", round(calibration_ci$improvement_percent, 1), "% [", round(calibration_ci$improvement_ci_low, 1), ", ", round(calibration_ci$improvement_ci_high, 1), "]"), collapse = "\n"),
	"\n\nWrong-source control intervals\n",
	paste(paste0(wrong_ci$anchors, " anchors: target ", round(wrong_ci$correct_improvement_percent, 1), "% [", round(wrong_ci$correct_ci_low, 1), ", ", round(wrong_ci$correct_ci_high, 1), "]; other-source ", round(wrong_ci$wrong_improvement_percent, 1), "% [", round(wrong_ci$wrong_ci_low, 1), ", ", round(wrong_ci$wrong_ci_high, 1), "]"), collapse = "\n"),
	"\n\nCalibration-method comparison\n",
	"Sources: ", n_comp1, "\n",
	"Conventional improvement: ", round(comp1_conv_imp, 1), "% [", round(comp1_ci[1, "conventional_improvement"], 1), ", ", round(comp1_ci[2, "conventional_improvement"], 1), "]\n",
	"Hierarchical improvement: ", round(comp1_hier_imp, 1), "% [", round(comp1_ci[1, "hierarchical_improvement"], 1), ", ", round(comp1_ci[2, "hierarchical_improvement"], 1), "]\n",
	"RMSE conventional - hierarchical: ", round(comp1_diff, 3), " [", round(comp1_ci[1, "RMSE_conventional_minus_hierarchical"], 3), ", ", round(comp1_ci[2, "RMSE_conventional_minus_hierarchical"], 3), "]\n",
	"\nHuSkin one-anchor deployment estimate\n",
	"Sources improved: ", n_all1_improved, "/", n_all1, "\n",
	"RMSE uncalibrated: ", round(all1_rmse_uncalibrated, 3), "\n",
	"RMSE calibrated: ", round(all1_rmse_calibrated, 3), "\n",
	"Improvement: ", round(all1_improvement, 1), "% [", round(all1_ci[1], 1), ", ", round(all1_ci[2], 1), "]\n",
	sep = ""
)