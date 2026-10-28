# 21_validate_cheruvu_qspr.R

library(lme4)

dat <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)
large <- read.csv("R/results/10_large_residuals.csv", stringsAsFactors = FALSE)

s1 <- c("mw", "log_kow")
s2 <- c("mw", "log_kow", "melting_point_k", "log_saq", "h_donor", "h_acceptor")

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

run_loco <- function(x, descriptors, label, FUN = median, domain = FALSE) {
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
		fixed <- paste(paste0(descriptors, "_z"), collapse = " + ")
		lm_fit <- lm(as.formula(paste("mean_log_kp ~", fixed)), data = profile)
		train$compound_id <- factor(train$compound_id)
		train$reference_id <- factor(train$reference_id)
		temperature_mean <- mean(train$temperature_c)
		train$temperature_c <- train$temperature_c - temperature_mean
		test$temperature_c <- test$temperature_c - temperature_mean
		mix_fit <- lmer(
			as.formula(paste("log_kp ~ temperature_c +", fixed, "+ (1 | compound_id) + (1 | reference_id)")),
			data = train, REML = TRUE, control = ctrl
		)
		test_profile <- aggregate(
			test[descriptors],
			list(compound_id = test$compound_id),
			FUN
		)
		for (v in descriptors) test_profile[[paste0(v, "_z")]] <- (test_profile[[v]] - center[v]) / scale[v]
		pred_lm <- as.numeric(predict(lm_fit, newdata = test_profile))
		pred_mix <- as.numeric(predict(mix_fit, newdata = test, re.form = ~(1 | reference_id), allow.new.levels = TRUE))
		md <- NA_real_
		in_domain <- NA
		if (domain) {
			zvars <- paste0(descriptors, "_z")
			ztrain <- as.matrix(profile[zvars])
			ztest <- as.numeric(test_profile[zvars])
			md <- mahalanobis(ztest, colMeans(ztrain), cov(ztrain))
			within_range <- all(sapply(descriptors, function(v) test_profile[[v]] >= min(profile[[v]]) & test_profile[[v]] <= max(profile[[v]])))
			in_domain <- md <= qchisq(0.95, length(descriptors)) & within_range
		}
		out[[i]] <- data.frame(
			specification = label,
			compound_id = test_id,
			observed = test$log_kp,
			conventional = pred_lm,
			hierarchical = pred_mix,
			mahalanobis = md,
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

# Primary LOCO validation and applicability domain

p1 <- run_loco(dat, s1, "S1", domain = TRUE)
p2 <- run_loco(dat, s2, "S2", domain = TRUE)

pred <- rbind(p1, p2)

domain <- aggregate(
	list(abs_error_hierarchical = abs(pred$observed - pred$hierarchical)),
	list(specification = pred$specification, compound_id = pred$compound_id, mahalanobis = pred$mahalanobis, in_domain = pred$in_domain),
	mean
)

domain_summary <- do.call(rbind, lapply(split(domain, list(domain$specification, domain$in_domain), drop = TRUE), function(x) {
	data.frame(
		specification = x$specification[1],
		in_domain = x$in_domain[1],
		n_compounds = nrow(x),
		MAE_hierarchical = mean(x$abs_error_hierarchical)
	)
}))
rownames(domain_summary) <- NULL

# Sensitivity analyses

multi <- names(which(tapply(dat$reference_id, dat$compound_id, function(x) length(unique(x))) > 1))
sens_multi <- run_loco(dat[dat$compound_id %in% multi, ], s2, "Multi-reference compounds")
sens_outlier <- run_loco(dat[!dat$source_row %in% large$source_row, ], s2, "Exclude |standardized residual| > 3")
sens_mean <- run_loco(dat, s2, "Mean descriptor aggregation", FUN = mean)

sens <- rbind(
	data.frame(scenario = "Primary S2", summarize_loco(p2)),
	data.frame(scenario = "Multi-reference compounds", summarize_loco(sens_multi)),
	data.frame(scenario = "Exclude |standardized residual| > 3", summarize_loco(sens_outlier)),
	data.frame(scenario = "Mean descriptor aggregation", summarize_loco(sens_mean))
)

write.csv(domain[c("specification", "compound_id", "mahalanobis", "in_domain")], "R/results/21_applicability_domain.csv", row.names = FALSE)
write.csv(domain_summary, "R/results/21_applicability_summary.csv", row.names = FALSE)
write.csv(sens, "R/results/21_sensitivity_performance.csv", row.names = FALSE)

s1_rmse <- c(sqrt(mean((p1$observed - p1$conventional)^2)), sqrt(mean((p1$observed - p1$hierarchical)^2)))
s2_rmse <- sens$RMSE[sens$scenario == "Primary S2"]

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