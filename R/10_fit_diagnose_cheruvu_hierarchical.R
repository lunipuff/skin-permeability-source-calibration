# ============================================================
# 10_fit_diagnose_cheruvu_hierarchical.R
# ============================================================

library(lme4)

dat <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)

dat$compound_id <- factor(dat$compound_id)
dat$reference_id <- factor(dat$reference_id)
dat$temperature_c <- dat$temperature_c - mean(dat$temperature_c)
dat$skin_thickness_mm <- dat$skin_thickness_mm - mean(dat$skin_thickness_mm)

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

m0 <- lmer(log_kp ~ 1 + (1 | compound_id), data = dat, REML = FALSE, control = ctrl)
m1 <- lmer(log_kp ~ 1 + (1 | compound_id) + (1 | reference_id), data = dat, REML = FALSE, control = ctrl)
m2 <- lmer(log_kp ~ temperature_c + (1 | compound_id) + (1 | reference_id), data = dat, REML = FALSE, control = ctrl)
m3 <- lmer(log_kp ~ temperature_c + skin_thickness_mm + (1 | compound_id) + (1 | reference_id), data = dat, REML = FALSE, control = ctrl)

mods <- list(M0 = m0, M1 = m1, M2 = m2, M3 = m3)

cmp <- do.call(rbind, lapply(names(mods), function(x) {
	m <- mods[[x]]
	data.frame(
		model = x,
		AIC = AIC(m),
		BIC = BIC(m),
		logLik = as.numeric(logLik(m)),
		df = attr(logLik(m), "df")
	)
}))

lrt <- rbind(
	data.frame(comparison = "M0 vs M1", chisq = 2 * (logLik(m1) - logLik(m0)), df = attr(logLik(m1), "df") - attr(logLik(m0), "df")),
	data.frame(comparison = "M1 vs M2", chisq = 2 * (logLik(m2) - logLik(m1)), df = attr(logLik(m2), "df") - attr(logLik(m1), "df")),
	data.frame(comparison = "M2 vs M3", chisq = 2 * (logLik(m3) - logLik(m2)), df = attr(logLik(m3), "df") - attr(logLik(m2), "df"))
)
lrt$p <- pchisq(lrt$chisq, lrt$df, lower.tail = FALSE)

m2_reml <- update(m2, REML = TRUE)

fix <- summary(m2_reml)$coefficients
fixed <- data.frame(
	term = rownames(fix),
	estimate = fix[, "Estimate"],
	se = fix[, "Std. Error"],
	row.names = NULL
)

vc <- as.data.frame(VarCorr(m2_reml))
variance <- vc[c("grp", "vcov", "sdcor")]
names(variance) <- c("component", "variance", "sd")
variance$component[variance$component == "compound_id"] <- "Compound"
variance$component[variance$component == "reference_id"] <- "Reference"
variance$component[variance$component == "Residual"] <- "Residual"
variance$percent <- 100 * variance$variance / sum(variance$variance)

dat$fitted <- fitted(m2_reml)
dat$residual <- resid(m2_reml)
dat$std_residual <- dat$residual / sigma(m2_reml)

compound_re <- ranef(m2_reml, condVar = TRUE)$compound_id
compound_se <- sqrt(attr(ranef(m2_reml, condVar = TRUE)$compound_id, "postVar")[1, 1, ])
compound_effects <- data.frame(
	compound_id = rownames(compound_re),
	effect = compound_re[, 1],
	se = compound_se,
	row.names = NULL
)

reference_re <- ranef(m2_reml, condVar = TRUE)$reference_id
reference_se <- sqrt(attr(ranef(m2_reml, condVar = TRUE)$reference_id, "postVar")[1, 1, ])
reference_effects <- data.frame(
	reference_id = rownames(reference_re),
	effect = reference_re[, 1],
	se = reference_se,
	row.names = NULL
)

diagnostics <- data.frame(
	metric = c(
		"singular",
		"residual_mean",
		"residual_sd",
		"median_absolute_residual",
		"max_absolute_residual",
		"n_abs_std_residual_gt_2",
		"n_abs_std_residual_gt_3",
		"cor_residual_fitted",
		"cor_abs_residual_fitted",
		"cor_residual_temperature"
	),
	value = c(
		as.numeric(isSingular(m2_reml)),
		mean(dat$residual),
		sd(dat$residual),
		median(abs(dat$residual)),
		max(abs(dat$residual)),
		sum(abs(dat$std_residual) > 2),
		sum(abs(dat$std_residual) > 3),
		cor(dat$residual, dat$fitted),
		cor(abs(dat$residual), dat$fitted),
		cor(dat$residual, dat$temperature_c)
	)
)

write.csv(cmp, "R/results/10_model_comparison.csv", row.names = FALSE)
write.csv(lrt, "R/results/10_likelihood_ratio_tests.csv", row.names = FALSE)
write.csv(fixed, "R/results/10_fixed_effects.csv", row.names = FALSE)
write.csv(variance, "R/results/10_variance_components.csv", row.names = FALSE)
write.csv(diagnostics, "R/results/10_diagnostics.csv", row.names = FALSE)
write.csv(compound_effects, "R/results/10_compound_effects.csv", row.names = FALSE)
write.csv(reference_effects, "R/results/10_reference_effects.csv", row.names = FALSE)

large_residuals <- dat[abs(dat$std_residual) > 3, c("source_row", "compound_id", "reference_id", "log_kp", "fitted", "std_residual")]
write.csv(large_residuals, "R/results/10_large_residuals.csv", row.names = FALSE)

saveRDS(m2_reml, "R/results/10_primary_hierarchical_model.rds")

cat(
	"\nCheruvu hierarchical model complete\n",
	"Selected model: M2\n",
	"M0 AIC: ", round(AIC(m0), 2), "\n",
	"M1 AIC: ", round(AIC(m1), 2), "\n",
	"M2 AIC: ", round(AIC(m2), 2), "\n",
	"M3 AIC: ", round(AIC(m3), 2), "\n",
	"Temperature estimate: ", round(fixed$estimate[fixed$term == "temperature_c"], 4), "\n",
	"Temperature p-value: ", signif(lrt$p[lrt$comparison == "M1 vs M2"], 3), "\n",
	"Thickness p-value: ", signif(lrt$p[lrt$comparison == "M2 vs M3"], 3), "\n",
	"Compound variance: ", round(variance$percent[variance$component == "Compound"], 1), "%\n",
	"Reference variance: ", round(variance$percent[variance$component == "Reference"], 1), "%\n",
	"Residual variance: ", round(variance$percent[variance$component == "Residual"], 1), "%\n",
	"Singular fit: ", isSingular(m2_reml), "\n",
	"|standardized residual| > 3: ", sum(abs(dat$std_residual) > 3), "\n",
	sep = ""
)