# ============================================================
# 21_validate_cheruvu_qspr.R
# ============================================================

library(lme4)

dat <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)
large <- read.csv("R/results/10_large_residuals.csv", stringsAsFactors = FALSE)

vars <- c("mw", "log_kow", "melting_point_k", "log_saq", "h_donor", "h_acceptor")
s1 <- c("mw", "log_kow")
s2 <- vars
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

metrics <- function(obs, pred) {
	c(
		RMSE = sqrt(mean((obs - pred)^2)),
		MAE = mean(abs(obs - pred)),
		R2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
		r = cor(obs, pred)
	)
}

make_profile <- function(x, descriptors, FUN = median) {
	p <- aggregate(x[descriptors], list(compound_id = x$compound_id), FUN)
	y <- aggregate(x$log_kp, list(compound_id = x$compound_id), mean)
	names(y)[2] <- "mean_log_kp"
	merge(p, y, by = "compound_id", sort = FALSE)
}

run_loco <- function(x, descriptors, label, FUN = median, nonlinear = FALSE, domain = FALSE) {
	compounds <- unique(x$compound_id)
	out <- vector("list", length(compounds))

	for (i in seq_along(compounds)) {
		test_id <- compounds[i]
		train <- x[x$compound_id != test_id, ]
		test <- x[x$compound_id == test_id, ]
		profile <- make_profile(train, descriptors, FUN)

		center <- sapply(profile[descriptors], mean)
		scale <- sapply(profile[descriptors], sd)

		for (v in descriptors) {
			profile[[paste0(v, "_z")]] <- (profile[[v]] - center[v]) / scale[v]
			train[[paste0(v, "_z")]] <- (train[[v]] - center[v]) / scale[v]
			test[[paste0(v, "_z")]] <- (test[[v]] - center[v]) / scale[v]
		}

		if (nonlinear) {
			fixed <- "mw_z + log_kow_z + melting_point_k_z + log_saq_z + I(log_saq_z^2) + h_acceptor_z + log_kow_z:h_acceptor_z"
		} else {
			fixed <- paste(paste0(descriptors, "_z"), collapse = " + ")
		}

		lm_fit <- lm(as.formula(paste("mean_log_kp ~", fixed)), data = profile)

		train$compound_id <- factor(train$compound_id)
		train$reference_id <- factor(train$reference_id)
		train$temperature_c <- train$temperature_c - mean(train$temperature_c)
		test$temperature_c <- test$temperature_c - mean(x[x$compound_id != test_id, "temperature_c"])

		mix_fit <- lmer(
			as.formula(paste("log_kp ~ temperature_c +", fixed, "+ (1 | compound_id) + (1 | reference_id)")),
			data = train,
			REML = TRUE,
			control = ctrl
		)

		test_profile <- test[1, descriptors, drop = FALSE]
		for (v in descriptors) test_profile[[paste0(v, "_z")]] <- (test_profile[[v]] - center[v]) / scale[v]

		pred_lm <- as.numeric(predict(lm_fit, newdata = test_profile))
		pred_mix <- as.numeric(predict(mix_fit, newdata = test, re.form = ~(1 | reference_id), allow.new.levels = TRUE))

		md <- NA
		threshold <- NA
		within_range <- NA
		in_domain <- NA

		if (domain) {
			zvars <- paste0(descriptors, "_z")
			ztrain <- as.matrix(profile[zvars])
			ztest <- as.numeric(test_profile[zvars])
			md <- mahalanobis(ztest, colMeans(ztrain), cov(ztrain))
			threshold <- qchisq(0.95, length(descriptors))
			within_range <- all(sapply(descriptors, function(v) test_profile[[v]] >= min(profile[[v]]) & test_profile[[v]] <= max(profile[[v]])))
			in_domain <- md <= threshold & within_range
		}

		out[[i]] <- data.frame(
			specification = label,
			compound_id = test_id,
			source_row = test$source_row,
			observed = test$log_kp,
			conventional = pred_lm,
			hierarchical = pred_mix,
			mahalanobis = md,
			domain_threshold = threshold,
			within_range = within_range,
			in_domain = in_domain
		)
	}

	do.call(rbind, out)
}

summarize_loco <- function(x) {
	a <- metrics(x$observed, x$conventional)
	b <- metrics(x$observed, x$hierarchical)

	data.frame(
		specification = unique(x$specification),
		approach = c("Conventional", "Hierarchical"),
		RMSE = c(a["RMSE"], b["RMSE"]),
		MAE = c(a["MAE"], b["MAE"]),
		R2 = c(a["R2"], b["R2"]),
		r = c(a["r"], b["r"])
	)
}

compound_balanced <- function(x) {
	cmp <- split(x, x$compound_id)

	a <- sapply(cmp, function(z) sqrt(mean((z$observed - z$conventional)^2)))
	b <- sapply(cmp, function(z) sqrt(mean((z$observed - z$hierarchical)^2)))

	data.frame(
		specification = unique(x$specification),
		approach = c("Conventional", "Hierarchical"),
		mean_compound_RMSE = c(mean(a), mean(b))
	)
}

# Primary LOCO validation

p1 <- run_loco(dat, s1, "S1", domain = TRUE)
p2 <- run_loco(dat, s2, "S2", domain = TRUE)
pred <- rbind(p1, p2)

performance <- rbind(summarize_loco(p1), summarize_loco(p2))
balanced <- rbind(compound_balanced(p1), compound_balanced(p2))

# Applicability domain

domain <- aggregate(
	cbind(abs_error_conventional = abs(pred$observed - pred$conventional), abs_error_hierarchical = abs(pred$observed - pred$hierarchical)),
	list(specification = pred$specification, compound_id = pred$compound_id, mahalanobis = pred$mahalanobis, domain_threshold = pred$domain_threshold, within_range = pred$within_range, in_domain = pred$in_domain),
	mean
)

domain_summary <- do.call(rbind, lapply(split(domain, list(domain$specification, domain$in_domain), drop = TRUE), function(x) {
	data.frame(
		specification = x$specification[1],
		in_domain = x$in_domain[1],
		n_compounds = nrow(x),
		MAE_conventional = mean(x$abs_error_conventional),
		MAE_hierarchical = mean(x$abs_error_hierarchical)
	)
}))
rownames(domain_summary) <- NULL

# Sensitivity analyses

multi <- names(which(tapply(dat$reference_id, dat$compound_id, function(x) length(unique(x))) > 1))
sens_multi <- run_loco(dat[dat$compound_id %in% multi, ], s2, "Multi-reference compounds")
sens_outlier <- run_loco(dat[!dat$source_row %in% large$source_row, ], s2, "Exclude |standardized residual| > 3")
sens_mean <- run_loco(dat, s2, "Mean descriptor aggregation", FUN = mean)
sens_nonlinear <- run_loco(dat, s2, "Expanded nonlinear specification", nonlinear = TRUE)

sens <- rbind(
	data.frame(scenario = "Primary S2", summarize_loco(p2)),
	data.frame(scenario = "Multi-reference compounds", summarize_loco(sens_multi)),
	data.frame(scenario = "Exclude |standardized residual| > 3", summarize_loco(sens_outlier)),
	data.frame(scenario = "Mean descriptor aggregation", summarize_loco(sens_mean)),
	data.frame(scenario = "Expanded nonlinear specification", summarize_loco(sens_nonlinear))
)

sens_bal <- rbind(
	data.frame(scenario = "Primary S2", compound_balanced(p2)),
	data.frame(scenario = "Multi-reference compounds", compound_balanced(sens_multi)),
	data.frame(scenario = "Exclude |standardized residual| > 3", compound_balanced(sens_outlier)),
	data.frame(scenario = "Mean descriptor aggregation", compound_balanced(sens_mean)),
	data.frame(scenario = "Expanded nonlinear specification", compound_balanced(sens_nonlinear))
)

# Outputs

write.csv(performance, "R/results/21_loco_performance.csv", row.names = FALSE)
write.csv(balanced, "R/results/21_loco_compound_balanced.csv", row.names = FALSE)
write.csv(pred, "R/results/21_loco_predictions.csv", row.names = FALSE)
write.csv(domain, "R/results/21_applicability_domain.csv", row.names = FALSE)
write.csv(domain_summary, "R/results/21_applicability_summary.csv", row.names = FALSE)
write.csv(sens, "R/results/21_sensitivity_performance.csv", row.names = FALSE)
write.csv(sens_bal, "R/results/21_sensitivity_compound_balanced.csv", row.names = FALSE)

s1_rmse <- performance$RMSE[performance$specification == "S1"]
s2_rmse <- performance$RMSE[performance$specification == "S2"]

cat(
	"\nCheruvu QSPR validation complete\n",
	"S1 conventional RMSE: ", round(s1_rmse[1], 3), "\n",
	"S1 hierarchical RMSE: ", round(s1_rmse[2], 3), "\n",
	"S1 RMSE improvement: ", round(100 * (1 - s1_rmse[2] / s1_rmse[1]), 1), "%\n",
	"S2 conventional RMSE: ", round(s2_rmse[1], 3), "\n",
	"S2 hierarchical RMSE: ", round(s2_rmse[2], 3), "\n",
	"S2 RMSE improvement: ", round(100 * (1 - s2_rmse[2] / s2_rmse[1]), 1), "%\n",
	"S1 in-domain compounds: ", sum(domain$specification == "S1" & domain$in_domain), "/", sum(domain$specification == "S1"), "\n",
	"S2 in-domain compounds: ", sum(domain$specification == "S2" & domain$in_domain), "/", sum(domain$specification == "S2"), "\n",
	sep = ""
)