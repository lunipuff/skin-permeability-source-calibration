# 10_fit_diagnose_cheruvu_hierarchical.R

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
	data.frame(comparison = "M1 vs M2",
				chisq = 2 * as.numeric(logLik(m2) - logLik(m1)),
				df = attr(logLik(m2), "df") - attr(logLik(m1), "df")
				),
	data.frame(comparison = "M2 vs M3", 
				chisq = 2 * as.numeric(logLik(m3) - logLik(m2)),
				df = attr(logLik(m3), "df") - attr(logLik(m2), "df")
				)
)
lrt$p <- pchisq(lrt$chisq, lrt$df, lower.tail = FALSE)

m2_reml <- update(m2, REML = TRUE)

fixed <- data.frame(term = "temperature_c", estimate = unname(fixef(m2_reml)["temperature_c"]))

vc <- as.data.frame(VarCorr(m2_reml))

variance <- data.frame(
	component = c(compound_id = "Compound", reference_id = "Reference", Residual = "Residual")[vc$grp],
	percent = 100 * vc$vcov / sum(vc$vcov),
	row.names = NULL
)

dat$std_residual <- resid(m2_reml) / sigma(m2_reml)
large_residuals <- dat[abs(dat$std_residual) > 3, "source_row", drop = FALSE]

write.csv(cmp, "R/results/10_model_comparison.csv", row.names = FALSE)
write.csv(lrt[c("comparison", "p")], "R/results/10_likelihood_ratio_tests.csv", row.names = FALSE)
write.csv(fixed, "R/results/10_fixed_effects.csv", row.names = FALSE)
write.csv(variance, "R/results/10_variance_components.csv", row.names = FALSE)
write.csv(large_residuals, "R/results/10_large_residuals.csv", row.names = FALSE)
saveRDS(m2_reml, "R/results/10_primary_hierarchical_model.rds")

cat(
	"\nCheruvu hierarchical model complete\n",
	"Selected model: M2\n",
	"M0 AIC: ", round(AIC(m0), 2), "\n",
	"M1 AIC: ", round(AIC(m1), 2), "\n",
	"M2 AIC: ", round(AIC(m2), 2), "\n",
	"M3 AIC: ", round(AIC(m3), 2), "\n",
	"Temperature estimate: ", round(fixed$estimate, 4), "\n",
	"Temperature p-value: ", signif(lrt$p[lrt$comparison == "M1 vs M2"], 3), "\n",
	"Thickness p-value: ", signif(lrt$p[lrt$comparison == "M2 vs M3"], 3), "\n",
	"Compound variance: ", round(variance$percent[variance$component == "Compound"], 1), "%\n",
	"Reference variance: ", round(variance$percent[variance$component == "Reference"], 1), "%\n",
	"Residual variance: ", round(variance$percent[variance$component == "Residual"], 1), "%\n",
	"Singular fit: ", isSingular(m2_reml), "\n",
	"|standardized residual| > 3: ", nrow(large_residuals), "\n",
	sep = ""
)