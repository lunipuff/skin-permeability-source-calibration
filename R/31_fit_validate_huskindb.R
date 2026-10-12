# ============================================================
# 31_fit_validate_huskindb.R
# ============================================================

library(lme4)

dat <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)
dat$compound_id <- factor(dat$compound_id)
dat$reference_id <- factor(dat$reference_id)

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
set.seed(2026)

metrics <- function(obs, pred) {
	rmse <- sqrt(mean((obs - pred)^2))
	mae <- mean(abs(obs - pred))
	r2 <- if (length(obs) > 1 && sd(obs) > 0) 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2) else NA_real_
	r <- if (length(obs) > 1 && sd(obs) > 0 && sd(pred) > 0) cor(obs, pred) else NA_real_

	c(RMSE = rmse, MAE = mae, R2 = r2, r = r)
}

profile <- function(x) {
	p <- aggregate(x[c("mw", "log_kow")], list(compound_id = x$compound_id), median)
	y <- aggregate(x$log_kp, list(compound_id = x$compound_id), mean)
	names(y)[2] <- "mean_log_kp"
	merge(p, y, by = "compound_id", sort = FALSE)
}

scale_data <- function(train, test = NULL) {
	p <- profile(train)
	center <- sapply(p[c("mw", "log_kow")], mean)
	scale <- sapply(p[c("mw", "log_kow")], sd)

	for (v in c("mw", "log_kow")) {
		p[[paste0(v, "_z")]] <- (p[[v]] - center[v]) / scale[v]
		train[[paste0(v, "_z")]] <- (train[[v]] - center[v]) / scale[v]
		if (!is.null(test)) test[[paste0(v, "_z")]] <- (test[[v]] - center[v]) / scale[v]
	}

	list(train = train, test = test, profile = p, center = center, scale = scale)
}

# Hierarchical replication

m0 <- lmer(log_kp ~ 1 + (1 | compound_id), data = dat, REML = FALSE, control = ctrl)
m1 <- lmer(log_kp ~ 1 + (1 | compound_id) + (1 | reference_id), data = dat, REML = FALSE, control = ctrl)

sc <- scale_data(dat)
dat_z <- sc$train

q0 <- lmer(log_kp ~ 1 + (1 | compound_id) + (1 | reference_id), data = dat_z, REML = FALSE, control = ctrl)
q1 <- lmer(log_kp ~ mw_z + log_kow_z + (1 | compound_id) + (1 | reference_id), data = dat_z, REML = FALSE, control = ctrl)

q0r <- update(q0, REML = TRUE)
q1r <- update(q1, REML = TRUE)

v0 <- as.data.frame(VarCorr(q0r))
v1 <- as.data.frame(VarCorr(q1r))
compound_v0 <- v0$vcov[v0$grp == "compound_id"]
compound_v1 <- v1$vcov[v1$grp == "compound_id"]

m1r <- update(m1, REML = TRUE)
saveRDS(m1r, "R/results/31_primary_huskindb_model.rds")

vc <- as.data.frame(VarCorr(m1r))
compound_pct <- 100 * vc$vcov[vc$grp == "compound_id"] / sum(vc$vcov)
reference_pct <- 100 * vc$vcov[vc$grp == "reference_id"] / sum(vc$vcov)
residual_pct <- 100 * vc$vcov[vc$grp == "Residual"] / sum(vc$vcov)

temp <- dat[!is.na(dat$temperature_c), ]
temp$temperature_c <- temp$temperature_c - mean(temp$temperature_c)

mt0 <- lmer(log_kp ~ 1 + (1 | compound_id) + (1 | reference_id), data = temp, REML = FALSE, control = ctrl)
mt1 <- lmer(log_kp ~ temperature_c + (1 | compound_id) + (1 | reference_id), data = temp, REML = FALSE, control = ctrl)

temp_chisq <- 2 * (as.numeric(logLik(mt1)) - as.numeric(logLik(mt0)))
temp_p <- pchisq(temp_chisq, 1, lower.tail = FALSE)

replication <- data.frame(
	n_observations = nrow(dat),
	n_compounds = length(unique(dat$compound_id)),
	n_references = length(unique(dat$reference_id)),
	AIC_compound_only = AIC(m0),
	AIC_compound_reference = AIC(m1),
	compound_variance_percent = compound_pct,
	reference_variance_percent = reference_pct,
	residual_variance_percent = residual_pct,
	qspr_compound_variance_explained = 100 * (1 - compound_v1 / compound_v0),
	qspr_delta_AIC = AIC(q1) - AIC(q0),
	temperature_n = nrow(temp),
	temperature_estimate = fixef(mt1)["temperature_c"],
	temperature_p = temp_p,
	singular = isSingular(m1r)
)

# LOCO validation and applicability domain

compounds <- unique(as.character(dat$compound_id))
loco <- vector("list", length(compounds))

for (i in seq_along(compounds)) {
	id <- compounds[i]
	train <- dat[as.character(dat$compound_id) != id, ]
	test <- dat[as.character(dat$compound_id) == id, ]

	sc <- scale_data(train, test)
	train <- sc$train
	test <- sc$test
	p <- sc$profile

	lm_fit <- lm(mean_log_kp ~ mw_z + log_kow_z, data = p)
	mix_fit <- lmer(log_kp ~ mw_z + log_kow_z + (1 | compound_id) + (1 | reference_id), data = train, REML = TRUE, control = ctrl)

	test_profile <- test[1, c("mw", "log_kow", "mw_z", "log_kow_z"), drop = FALSE]
	pred_lm <- as.numeric(predict(lm_fit, newdata = test_profile))
	pred_mix <- as.numeric(predict(mix_fit, newdata = test, re.form = ~(1 | reference_id), allow.new.levels = TRUE))

	ztrain <- as.matrix(p[c("mw_z", "log_kow_z")])
	ztest <- as.numeric(test_profile[c("mw_z", "log_kow_z")])
	md <- mahalanobis(ztest, colMeans(ztrain), cov(ztrain))
	in_range <- test_profile$mw >= min(p$mw) && test_profile$mw <= max(p$mw) && test_profile$log_kow >= min(p$log_kow) && test_profile$log_kow <= max(p$log_kow)

	loco[[i]] <- data.frame(
		compound_id = id,
		source_row = test$source_row,
		observed = test$log_kp,
		conventional = pred_lm,
		hierarchical = pred_mix,
		mahalanobis = md,
		in_domain = md <= qchisq(0.95, 2) && in_range
	)
}

loco <- do.call(rbind, loco)

a <- metrics(loco$observed, loco$conventional)
b <- metrics(loco$observed, loco$hierarchical)

loco_performance <- data.frame(
	approach = c("Conventional", "Hierarchical"),
	RMSE = c(a["RMSE"], b["RMSE"]),
	MAE = c(a["MAE"], b["MAE"]),
	R2 = c(a["R2"], b["R2"]),
	r = c(a["r"], b["r"])
)

cmp <- split(loco, loco$compound_id)
compound_rmse_conventional <- sapply(cmp, function(x) sqrt(mean((x$observed - x$conventional)^2)))
compound_rmse_hierarchical <- sapply(cmp, function(x) sqrt(mean((x$observed - x$hierarchical)^2)))

loco_performance$mean_compound_RMSE <- c(mean(compound_rmse_conventional), mean(compound_rmse_hierarchical))

domain <- aggregate(
	cbind(error_conventional = abs(loco$observed - loco$conventional), error_hierarchical = abs(loco$observed - loco$hierarchical)),
	list(compound_id = loco$compound_id, mahalanobis = loco$mahalanobis, in_domain = loco$in_domain),
	mean
)

# Source calibration

references <- unique(as.character(dat$reference_id))
anchor_sizes <- 1:5
calibration <- list()
z <- 1

for (ref in references) {
	train <- dat[as.character(dat$reference_id) != ref, ]
	held <- dat[as.character(dat$reference_id) == ref, ]

	train_compounds <- unique(as.character(train$compound_id))
	shared <- intersect(unique(as.character(held$compound_id)), train_compounds)
	novel <- setdiff(unique(as.character(held$compound_id)), train_compounds)

	if (!length(shared)) next

	sc <- scale_data(train, held)
	train <- sc$train
	held <- sc$test

	fit <- lmer(log_kp ~ mw_z + log_kow_z + (1 | compound_id) + (1 | reference_id), data = train, REML = TRUE, control = ctrl)

	vc <- as.data.frame(VarCorr(fit))
	v_ref <- vc$vcov[vc$grp == "reference_id"]
	v_res <- vc$vcov[vc$grp == "Residual"]

	shared_rows <- held[as.character(held$compound_id) %in% shared, ]
	shared_rows$expected <- predict(fit, newdata = shared_rows, re.form = ~(1 | compound_id), allow.new.levels = TRUE)

	shared_residual <- aggregate(shared_rows$log_kp - shared_rows$expected, list(compound_id = shared_rows$compound_id), mean)
	names(shared_residual)[2] <- "residual"

	novel_rows <- held[as.character(held$compound_id) %in% novel, ]
	if (nrow(novel_rows)) novel_rows$expected <- predict(fit, newdata = novel_rows, re.form = NA, allow.new.levels = TRUE)

	for (k in anchor_sizes) {
		if (length(shared) < k) next

		for (rep in 1:100) {
			anchors <- sample(shared, k)
			offset_raw <- mean(shared_residual$residual[as.character(shared_residual$compound_id) %in% anchors])
			shrinkage <- k * v_ref / (v_res + k * v_ref)
			offset <- shrinkage * offset_raw

			test_shared <- shared_rows[!as.character(shared_rows$compound_id) %in% anchors, ]

			if (nrow(test_shared)) {
				uncal <- metrics(test_shared$log_kp, test_shared$expected)
				cal <- metrics(test_shared$log_kp, test_shared$expected + offset)

				calibration[[z]] <- data.frame(
					reference_id = ref,
					task = "Known compound",
					anchors = k,
					replicate = rep,
					n_test = nrow(test_shared),
					offset_raw = offset_raw,
					shrinkage = shrinkage,
					offset = offset,
					RMSE_uncalibrated = uncal["RMSE"],
					RMSE_calibrated = cal["RMSE"],
					MAE_uncalibrated = uncal["MAE"],
					MAE_calibrated = cal["MAE"]
				)
				z <- z + 1
			}

			if (nrow(novel_rows)) {
				uncal <- metrics(novel_rows$log_kp, novel_rows$expected)
				cal <- metrics(novel_rows$log_kp, novel_rows$expected + offset)

				calibration[[z]] <- data.frame(
					reference_id = ref,
					task = "New compound",
					anchors = k,
					replicate = rep,
					n_test = nrow(novel_rows),
					offset_raw = offset_raw,
					shrinkage = shrinkage,
					offset = offset,
					RMSE_uncalibrated = uncal["RMSE"],
					RMSE_calibrated = cal["RMSE"],
					MAE_uncalibrated = uncal["MAE"],
					MAE_calibrated = cal["MAE"]
				)
				z <- z + 1
			}
		}
	}
}

calibration <- do.call(rbind, calibration)

# All eligible references at each anchor size

calibration_summary <- aggregate(
	cbind(RMSE_uncalibrated, RMSE_calibrated, MAE_uncalibrated, MAE_calibrated) ~ task + anchors,
	data = calibration,
	mean
)

calibration_summary$RMSE_improvement_percent <- 100 * (1 - calibration_summary$RMSE_calibrated / calibration_summary$RMSE_uncalibrated)
calibration_summary$MAE_improvement_percent <- 100 * (1 - calibration_summary$MAE_calibrated / calibration_summary$MAE_uncalibrated)

eligible <- aggregate(reference_id ~ task + anchors, data = calibration, function(x) length(unique(x)))
names(eligible)[3] <- "n_references"

calibration_summary <- merge(calibration_summary, eligible, by = c("task", "anchors"), sort = FALSE)
calibration_summary <- calibration_summary[order(calibration_summary$task, calibration_summary$anchors), ]

# Matched-reference anchor-size comparison
# Same references contribute to k = 1, 2, 3, 4, 5 and all contain new test compounds

matched_refs <- unique(calibration$reference_id[calibration$task == "New compound" & calibration$anchors == max(anchor_sizes)])
matched <- calibration[calibration$task == "New compound" & calibration$reference_id %in% matched_refs, ]

matched_summary <- aggregate(
	cbind(RMSE_uncalibrated, RMSE_calibrated, MAE_uncalibrated, MAE_calibrated) ~ anchors,
	data = matched,
	mean
)

matched_summary$n_references <- length(matched_refs)
matched_summary$RMSE_improvement_percent <- 100 * (1 - matched_summary$RMSE_calibrated / matched_summary$RMSE_uncalibrated)
matched_summary$MAE_improvement_percent <- 100 * (1 - matched_summary$MAE_calibrated / matched_summary$MAE_uncalibrated)
matched_summary <- matched_summary[order(matched_summary$anchors), ]

# Reference-level matched calibration summary

matched_reference <- aggregate(
	cbind(RMSE_uncalibrated, RMSE_calibrated, MAE_uncalibrated, MAE_calibrated) ~ reference_id + anchors,
	data = matched,
	mean
)

matched_reference$RMSE_improvement_percent <- 100 * (1 - matched_reference$RMSE_calibrated / matched_reference$RMSE_uncalibrated)
matched_reference$MAE_improvement_percent <- 100 * (1 - matched_reference$MAE_calibrated / matched_reference$MAE_uncalibrated)

# Outputs

write.csv(replication, "R/results/31_replication_summary.csv", row.names = FALSE)
write.csv(loco_performance, "R/results/31_loco_performance.csv", row.names = FALSE)
write.csv(domain, "R/results/31_applicability_domain.csv", row.names = FALSE)
write.csv(calibration_summary, "R/results/31_calibration_summary.csv", row.names = FALSE)
write.csv(matched_summary, "R/results/31_calibration_matched_summary.csv", row.names = FALSE)
write.csv(matched_reference, "R/results/31_calibration_matched_reference.csv", row.names = FALSE)
write.csv(calibration, "R/results/31_calibration_performance.csv", row.names = FALSE)

new_cal <- calibration_summary[calibration_summary$task == "New compound", ]

cat(
	"\nHuSkinDB replication and validation complete\n",
	"Reference variance: ", round(replication$reference_variance_percent, 1), "%\n",
	"MW + logKow compound variance explained: ", round(replication$qspr_compound_variance_explained, 1), "%\n",
	"LOCO conventional RMSE: ", round(loco_performance$RMSE[loco_performance$approach == "Conventional"], 3), "\n",
	"LOCO hierarchical RMSE: ", round(loco_performance$RMSE[loco_performance$approach == "Hierarchical"], 3), "\n",
	"LOCO RMSE improvement: ", round(100 * (1 - loco_performance$RMSE[loco_performance$approach == "Hierarchical"] / loco_performance$RMSE[loco_performance$approach == "Conventional"]), 1), "%\n",
	"In-domain compounds: ", sum(domain$in_domain), "/", nrow(domain), "\n",
	"Temperature subset n: ", replication$temperature_n, "\n",
	"Temperature p-value: ", signif(replication$temperature_p, 3), "\n",
	"\nNew-source + new-compound calibration, all eligible references\n",
	paste(paste0(new_cal$anchors, " anchors: ", round(new_cal$RMSE_improvement_percent, 1), "% RMSE improvement across ", new_cal$n_references, " references"), collapse = "\n"),
	"\n\nMatched-reference anchor comparison\n",
	"References eligible for all anchor sizes: ", length(matched_refs), "\n",
	paste(paste0(matched_summary$anchors, " anchors: ", round(matched_summary$RMSE_improvement_percent, 1), "% RMSE improvement"), collapse = "\n"),
	"\n",
	sep = ""
)