# 31_fit_validate_huskindb.R

library(lme4)

dat <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)

dat$compound_id <- factor(dat$compound_id)
dat$reference_id <- factor(dat$reference_id)

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
set.seed(2026)

metrics <- function(obs, pred) {
	c(RMSE = sqrt(mean((obs - pred)^2)), MAE = mean(abs(obs - pred)))
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
	list(train = train, test = test, profile = p)
}

# Hierarchical replication

m0 <- lmer(log_kp ~ 1 + (1 | compound_id), data = dat, REML = FALSE, control = ctrl)
m1 <- lmer(log_kp ~ 1 + (1 | compound_id) + (1 | reference_id), data = dat, REML = FALSE, control = ctrl)
m1r <- update(m1, REML = TRUE)

saveRDS(m1r, "R/results/31_primary_huskindb_model.rds")

sc <- scale_data(dat)

q1 <- lmer(log_kp ~ mw_z + log_kow_z + (1 | compound_id) + (1 | reference_id), data = sc$train, REML = FALSE, control = ctrl)
q1r <- update(q1, REML = TRUE)

vc <- as.data.frame(VarCorr(m1r))
v1 <- as.data.frame(VarCorr(q1r))

compound_v0 <- vc$vcov[vc$grp == "compound_id"]
compound_v1 <- v1$vcov[v1$grp == "compound_id"]

temp <- dat[!is.na(dat$temperature_c), ]

temp$temperature_c <- temp$temperature_c - mean(temp$temperature_c)

mt0 <- lmer(log_kp ~ 1 + (1 | compound_id) + (1 | reference_id), data = temp, REML = FALSE, control = ctrl)
mt1 <- lmer(log_kp ~ temperature_c + (1 | compound_id) + (1 | reference_id), data = temp, REML = FALSE, control = ctrl)

temp_chisq <- 2 * (as.numeric(logLik(mt1)) - as.numeric(logLik(mt0)))

replication <- data.frame(
	n_observations = nrow(dat),
	n_compounds = length(unique(dat$compound_id)),
	n_references = length(unique(dat$reference_id)),
	AIC_compound_only = AIC(m0),
	AIC_compound_reference = AIC(m1),
	compound_variance_percent = 100 * compound_v0 / sum(vc$vcov),
	reference_variance_percent = 100 * vc$vcov[vc$grp == "reference_id"] / sum(vc$vcov),
	residual_variance_percent = 100 * vc$vcov[vc$grp == "Residual"] / sum(vc$vcov),
	qspr_compound_variance_explained = 100 * (1 - compound_v1 / compound_v0),
	qspr_delta_AIC = AIC(q1) - AIC(m1),
	temperature_n = nrow(temp),
	temperature_estimate = unname(fixef(mt1)["temperature_c"]),
	temperature_p = pchisq(temp_chisq, 1, lower.tail = FALSE),
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
		observed = test$log_kp,
		conventional = pred_lm,
		hierarchical = pred_mix,
		mahalanobis = md,
		in_domain = md <= qchisq(0.95, 2) && in_range
	)
}

loco <- do.call(rbind, loco)

loco_rmse <- c(
	Conventional = sqrt(mean((loco$observed - loco$conventional)^2)),
	Hierarchical = sqrt(mean((loco$observed - loco$hierarchical)^2))
)

domain <- aggregate(
	list(error_hierarchical = abs(loco$observed - loco$hierarchical)),
	list(compound_id = loco$compound_id, mahalanobis = loco$mahalanobis, in_domain = loco$in_domain),
	mean
)

# Temperature availability among calibration sources

references <- unique(as.character(dat$reference_id))

anchor_sizes <- 1:5

temperature_eligibility <- do.call(rbind, lapply(references, function(ref) {
	train <- dat[as.character(dat$reference_id) != ref, ]
	held <- dat[as.character(dat$reference_id) == ref, ]
	train_compounds <- unique(as.character(train$compound_id))
	held_compounds <- unique(as.character(held$compound_id))
	shared <- intersect(held_compounds, train_compounds)
	novel <- setdiff(held_compounds, train_compounds)
	train_with_temperature <- unique(as.character(train$compound_id[!is.na(train$temperature_c)]))
	held_with_temperature <- unique(as.character(held$compound_id[!is.na(held$temperature_c)]))
	data.frame(
		originally_eligible = length(shared) >= 1 && length(novel) >= 1,
		n_temperature_anchors = length(intersect(shared, intersect(train_with_temperature, held_with_temperature))),
		n_temperature_test_compounds = length(intersect(novel, held_with_temperature)),
		n_recorded_temperatures = length(unique(held$temperature_c[!is.na(held$temperature_c)]))
	)
}))

temperature_eligible <- with(temperature_eligibility, originally_eligible & n_temperature_anchors >= 1 & n_temperature_test_compounds >= 1)

replication$temperature_calibration_n_sources <- sum(temperature_eligible)
replication$temperature_calibration_max_source_temperatures <- if (any(temperature_eligible)) max(temperature_eligibility$n_recorded_temperatures[temperature_eligible]) else NA_integer_

# Source calibration

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
		shrinkage <- k * v_ref / (v_res + k * v_ref)
		for (rep in 1:100) {
			anchors <- sample(shared, k)
			offset_raw <- mean(shared_residual$residual[as.character(shared_residual$compound_id) %in% anchors])
			offset <- shrinkage * offset_raw
			test_shared <- shared_rows[!as.character(shared_rows$compound_id) %in% anchors, ]
			if (nrow(test_shared)) {
				uncal_rmse <- sqrt(mean((test_shared$log_kp - test_shared$expected)^2))
				cal_rmse <- sqrt(mean((test_shared$log_kp - (test_shared$expected + offset))^2))
				calibration[[z]] <- data.frame(
					reference_id = ref,
					task = "Known compound",
					anchors = k,
					RMSE_uncalibrated = uncal_rmse,
					RMSE_calibrated = cal_rmse,
					MAE_uncalibrated = NA_real_,
					MAE_calibrated = NA_real_
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

# All-eligible summaries for console checks

new_cal <- calibration[calibration$task == "New compound", ]

calibration_summary <- aggregate(cbind(RMSE_uncalibrated, RMSE_calibrated) ~ anchors, data = new_cal, mean)

eligible <- aggregate(reference_id ~ anchors, data = new_cal, function(x) length(unique(x)))
names(eligible)[2] <- "n_references"

calibration_summary <- merge(calibration_summary, eligible, by = "anchors", sort = FALSE)
calibration_summary <- calibration_summary[order(calibration_summary$anchors), ]
calibration_summary$RMSE_improvement_percent <- 100 * (1 - calibration_summary$RMSE_calibrated / calibration_summary$RMSE_uncalibrated)

# Matched-source anchor-number comparison

matched_refs <- unique(new_cal$reference_id[new_cal$anchors == max(anchor_sizes)])

matched <- new_cal[new_cal$reference_id %in% matched_refs, ]

matched_summary <- aggregate(cbind(RMSE_uncalibrated, RMSE_calibrated) ~ anchors, data = matched, mean)
matched_summary$RMSE_improvement_percent <- 100 * (1 - matched_summary$RMSE_calibrated / matched_summary$RMSE_uncalibrated)
matched_summary <- matched_summary[order(matched_summary$anchors), ]

matched_reference <- aggregate(
	cbind(RMSE_uncalibrated, RMSE_calibrated, MAE_uncalibrated, MAE_calibrated) ~ reference_id + anchors,
	data = matched, mean
)

matched_reference$RMSE_improvement_percent <- 100 * (1 - matched_reference$RMSE_calibrated / matched_reference$RMSE_uncalibrated)
matched_reference$MAE_improvement_percent <- 100 * (1 - matched_reference$MAE_calibrated / matched_reference$MAE_uncalibrated)

write.csv(replication, "R/results/31_replication_summary.csv", row.names = FALSE)
write.csv(domain, "R/results/31_applicability_domain.csv", row.names = FALSE)
write.csv(matched_reference, "R/results/31_calibration_matched_reference.csv", row.names = FALSE)
write.csv(calibration[c("reference_id", "task", "anchors", "RMSE_uncalibrated", "RMSE_calibrated")], "R/results/31_calibration_performance.csv", row.names = FALSE)

cat(
	"\nHuSkinDB replication and validation complete\n",
	"Reference variance: ", round(replication$reference_variance_percent, 1), "%\n",
	"MW + logKow compound variance explained: ", round(replication$qspr_compound_variance_explained, 1), "%\n",
	"LOCO conventional RMSE: ", round(loco_rmse["Conventional"], 3), "\n",
	"LOCO hierarchical RMSE: ", round(loco_rmse["Hierarchical"], 3), "\n",
	"LOCO RMSE improvement: ", round(100 * (1 - loco_rmse["Hierarchical"] / loco_rmse["Conventional"]), 1), "%\n",
	"In-domain compounds: ", sum(domain$in_domain), "/", nrow(domain), "\n",
	"Temperature subset n: ", replication$temperature_n, "\n",
	"Temperature p-value: ", signif(replication$temperature_p, 3), "\n",
	"\nNew-source + new-compound calibration, all eligible references\n",
	paste(paste0(calibration_summary$anchors, " anchors: ", round(calibration_summary$RMSE_improvement_percent, 1), "% RMSE improvement across ", calibration_summary$n_references, " references"), collapse = "\n"),
	"\n\nMatched-reference anchor comparison\n",
	"References eligible for all anchor sizes: ", length(matched_refs), "\n",
	paste(paste0(matched_summary$anchors, " anchors: ", round(matched_summary$RMSE_improvement_percent, 1), "% RMSE improvement"), collapse = "\n"),
	"\n", sep = ""
)