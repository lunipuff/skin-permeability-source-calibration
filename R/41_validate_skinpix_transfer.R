# 41_validate_skinpix_transfer.R

library(lme4)

huskin <- read.csv("data/processed/40_huskindb_rdkit.csv", stringsAsFactors = FALSE)
skinpix <- read.csv("data/processed/40_skinpix.csv", stringsAsFactors = FALSE)

huskin$compound_id <- factor(huskin$compound_id)
huskin$reference_id <- factor(huskin$reference_id)

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

set.seed(2026)

metrics <- function(obs, pred) {
	c(RMSE = sqrt(mean((obs - pred)^2)), MAE = mean(abs(obs - pred)))
}

# HuSkinDB training model

profile <- aggregate(huskin[c("mw_rdkit", "logp_rdkit")], list(compound_id = huskin$compound_id), median)
center <- sapply(profile[c("mw_rdkit", "logp_rdkit")], mean)
scale <- sapply(profile[c("mw_rdkit", "logp_rdkit")], sd)

for (v in c("mw_rdkit", "logp_rdkit")) {
	profile[[paste0(v, "_z")]] <- (profile[[v]] - center[v]) / scale[v]
	huskin[[paste0(v, "_z")]] <- (huskin[[v]] - center[v]) / scale[v]
	skinpix[[paste0(v, "_z")]] <- (skinpix[[v]] - center[v]) / scale[v]
}

fit <- lmer(log_kp ~ mw_rdkit_z + logp_rdkit_z + (1 | compound_id) + (1 | reference_id), data = huskin, REML = TRUE, control = ctrl)

vc <- as.data.frame(VarCorr(fit))

v_ref <- vc$vcov[vc$grp == "reference_id"]
v_res <- vc$vcov[vc$grp == "Residual"]

# HuSkinDB-defined applicability domain

ztrain <- as.matrix(profile[c("mw_rdkit_z", "logp_rdkit_z")])

skinpix$mahalanobis <- mahalanobis(as.matrix(skinpix[c("mw_rdkit_z", "logp_rdkit_z")]), colMeans(ztrain), cov(ztrain))

skinpix$within_range <- skinpix$mw_rdkit >= min(profile$mw_rdkit) &
	skinpix$mw_rdkit <= max(profile$mw_rdkit) &
	skinpix$logp_rdkit >= min(profile$logp_rdkit) &
	skinpix$logp_rdkit <= max(profile$logp_rdkit)

skinpix$in_domain <- skinpix$mahalanobis <= qchisq(0.95, 2) & skinpix$within_range

# Source eligibility and temperature availability

huskin_with_temperature <- unique(as.character(huskin$compound_id[!is.na(huskin$temperature_c)]))

eligible <- do.call(rbind, lapply(unique(skinpix$reference_id), function(ref) {
	x <- skinpix[skinpix$reference_id == ref, ]
	x_temperature <- x[!is.na(x$temperature_c), ]
	anchors <- unique(x$huskin_compound_id[x$overlap_huskin & !is.na(x$huskin_compound_id)])
	novel <- unique(x$compound_id[!x$overlap_huskin])
	temperature_anchors <- unique(x_temperature$huskin_compound_id[
		x_temperature$overlap_huskin & !is.na(x_temperature$huskin_compound_id) &
		x_temperature$huskin_compound_id %in% huskin_with_temperature
	])
	temperature_novel <- unique(x_temperature$compound_id[!x_temperature$overlap_huskin])
	data.frame(
		reference_id = ref,
		n_anchors = length(anchors),
		eligible = length(anchors) >= 1 && length(novel) >= 1 && !any(x$doi_overlaps_huskin),
		n_temperature_anchors = length(temperature_anchors),
		n_temperature_test_compounds = length(temperature_novel)
	)
}))

eligible_refs <- eligible$reference_id[eligible$eligible]

# Cross-dataset source calibration

calibration <- list()
prediction_rows <- list()
z <- 1
q <- 1

for (ref in eligible_refs) {
	x <- skinpix[skinpix$reference_id == ref, ]
	anchor_rows <- x[x$overlap_huskin & !is.na(x$huskin_compound_id), ]
	novel_rows <- x[!x$overlap_huskin, ]
	anchor_rows$expected <- predict(
		fit,
		newdata = data.frame(
			mw_rdkit_z = anchor_rows$mw_rdkit_z,
			logp_rdkit_z = anchor_rows$logp_rdkit_z,
			compound_id = factor(anchor_rows$huskin_compound_id, levels = levels(huskin$compound_id)),
			reference_id = factor(NA, levels = levels(huskin$reference_id))
		),
		re.form = ~(1 | compound_id), allow.new.levels = TRUE
	)
	anchor_residual <- aggregate(anchor_rows$log_kp - anchor_rows$expected, list(huskin_compound_id = anchor_rows$huskin_compound_id), mean)
	names(anchor_residual)[2] <- "residual"
	novel_rows$expected <- predict(
		fit,
		newdata = data.frame(
			mw_rdkit_z = novel_rows$mw_rdkit_z,
			logp_rdkit_z = novel_rows$logp_rdkit_z,
			compound_id = factor(NA, levels = levels(huskin$compound_id)),
			reference_id = factor(NA, levels = levels(huskin$reference_id))
		),
		re.form = NA, allow.new.levels = TRUE
	)
	anchors <- unique(anchor_residual$huskin_compound_id)
	uncal <- metrics(novel_rows$log_kp, novel_rows$expected)
	for (k in seq_len(min(5, length(anchors)))) {
		shrinkage <- k * v_ref / (v_res + k * v_ref)
		for (rep in 1:100) {
			selected <- sample(anchors, k)
			offset <- shrinkage * mean(anchor_residual$residual[anchor_residual$huskin_compound_id %in% selected])
			cal <- metrics(novel_rows$log_kp, novel_rows$expected + offset)
			calibration[[z]] <- data.frame(
				reference_id = ref,
				anchors = k,
				RMSE_uncalibrated = uncal["RMSE"],
				RMSE_calibrated = cal["RMSE"],
				MAE_uncalibrated = uncal["MAE"],
				MAE_calibrated = cal["MAE"]
			)
			z <- z + 1
			if (k == 1) {
				prediction_rows[[q]] <- data.frame(
					reference_id = ref,
					source_row = novel_rows$source_row,
					compound_id = novel_rows$compound_id,
					compound_name = novel_rows$compound_name,
					observed = novel_rows$log_kp,
					pred_uncalibrated = novel_rows$expected,
					pred_calibrated = novel_rows$expected + offset,
					mahalanobis = novel_rows$mahalanobis,
					in_domain = novel_rows$in_domain
				)
				q <- q + 1
			}
		}
	}
}

calibration <- do.call(rbind, calibration)
prediction_rows <- do.call(rbind, prediction_rows)

# Average one-anchor predictions within each measurement before calculating errors

prediction_summary <- aggregate(
	pred_calibrated ~ reference_id + source_row + compound_id + compound_name + observed + pred_uncalibrated + mahalanobis + in_domain,
	data = prediction_rows, mean
)

prediction_summary$abs_error_uncalibrated <- abs(prediction_summary$observed - prediction_summary$pred_uncalibrated)

prediction_summary$abs_error_calibrated <- abs(prediction_summary$observed - prediction_summary$pred_calibrated)

# All-eligible summary for console checks

summary_all <- aggregate(cbind(RMSE_uncalibrated, RMSE_calibrated) ~ anchors, data = calibration, mean)

n_refs <- aggregate(reference_id ~ anchors, data = calibration, function(x) length(unique(x)))
names(n_refs)[2] <- "n_references"

summary_all <- merge(summary_all, n_refs, by = "anchors", sort = FALSE)
summary_all$RMSE_improvement_percent <- 100 * (1 - summary_all$RMSE_calibrated / summary_all$RMSE_uncalibrated)
summary_all <- summary_all[order(summary_all$anchors), ]

# Maximum common anchor-size comparison for console checks

eligible_counts <- eligible[eligible$eligible, ]
matched_max <- min(5, min(eligible_counts$n_anchors))
matched_refs <- eligible_counts$reference_id[eligible_counts$n_anchors >= matched_max]

matched <- calibration[calibration$reference_id %in% matched_refs & calibration$anchors <= matched_max, ]

matched_summary <- aggregate(cbind(RMSE_uncalibrated, RMSE_calibrated) ~ anchors, data = matched, mean)
matched_summary$RMSE_improvement_percent <- 100 * (1 - matched_summary$RMSE_calibrated / matched_summary$RMSE_uncalibrated)
matched_summary <- matched_summary[order(matched_summary$anchors), ]

# Fixed subset supporting one through three anchors

matched3_refs <- eligible_counts$reference_id[eligible_counts$n_anchors >= 3]
matched3 <- calibration[calibration$reference_id %in% matched3_refs & calibration$anchors <= 3, ]

matched3_summary <- aggregate(cbind(RMSE_uncalibrated, RMSE_calibrated) ~ anchors, data = matched3, mean)

matched3_summary$n_references <- length(matched3_refs)
matched3_summary$RMSE_improvement_percent <- 100 * (1 - matched3_summary$RMSE_calibrated / matched3_summary$RMSE_uncalibrated)
matched3_summary <- matched3_summary[order(matched3_summary$anchors), ]

# Source-level results

reference_summary <- aggregate(
	cbind(RMSE_uncalibrated, RMSE_calibrated, MAE_uncalibrated, MAE_calibrated) ~ reference_id + anchors,
	data = calibration, mean
)

reference_summary$RMSE_improvement_percent <- 100 * (1 - reference_summary$RMSE_calibrated / reference_summary$RMSE_uncalibrated)

reference_summary$MAE_improvement_percent <- 100 * (1 - reference_summary$MAE_calibrated / reference_summary$MAE_uncalibrated)

novel_domain <- unique(skinpix[
	skinpix$reference_id %in% eligible_refs & !skinpix$overlap_huskin,
	c("compound_id", "compound_name", "mw_rdkit", "logp_rdkit", "mahalanobis", "within_range", "in_domain")
])

write.csv(eligible[c("reference_id", "eligible", "n_temperature_anchors", "n_temperature_test_compounds")], "R/results/41_skinpix_eligible_references.csv", row.names = FALSE)
write.csv(matched3_summary, "R/results/41_skinpix_calibration_matched_3anchor.csv", row.names = FALSE)
write.csv(reference_summary, "R/results/41_skinpix_calibration_reference.csv", row.names = FALSE)
write.csv(prediction_summary[c("compound_id", "abs_error_uncalibrated", "abs_error_calibrated")], "R/results/41_skinpix_one_anchor_predictions.csv", row.names = FALSE)
write.csv(novel_domain[c("compound_id", "mahalanobis", "in_domain")], "R/results/41_skinpix_transfer_domain.csv", row.names = FALSE)

cat(
	"\nSkinPiX cross-dataset transfer validation complete\n",
	"HuSkin training rows: ", nrow(huskin), "\n",
	"HuSkin training compounds: ", length(unique(huskin$compound_id)), "\n",
	"Eligible SkinPiX references: ", length(eligible_refs), "\n",
	"Novel SkinPiX test compounds: ", length(unique(novel_domain$compound_id)), "\n",
	"Novel compounds in HuSkin chemical domain: ", sum(novel_domain$in_domain), "/", nrow(novel_domain), "\n",
	"\nAll eligible references\n",
	paste(paste0(summary_all$anchors, " anchors: ", round(summary_all$RMSE_improvement_percent, 1), "% RMSE improvement across ", summary_all$n_references, " references"), collapse = "\n"),
	"\n\nMatched-reference comparison\n",
	"References: ", length(matched_refs), "\n",
	"Maximum common anchor size: ", matched_max, "\n",
	paste(paste0(matched_summary$anchors, " anchors: ", round(matched_summary$RMSE_improvement_percent, 1), "% RMSE improvement"), collapse = "\n"),
	"\n\nMatched references supporting >=3 anchors\n",
	"References: ", length(matched3_refs), "\n",
	paste(paste0(matched3_summary$anchors, " anchors: ", round(matched3_summary$RMSE_improvement_percent, 1), "% RMSE improvement"), collapse = "\n"),
	"\n\nOne-anchor predictions saved: ", nrow(prediction_summary), " measurement rows\n",
	sep = ""
)