# ============================================================
# 20_fit_cheruvu_qspr_models.R
# ============================================================

library(lme4)

dat <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)
dat$compound_id <- factor(dat$compound_id)
dat$reference_id <- factor(dat$reference_id)
dat$temperature_c <- dat$temperature_c - mean(dat$temperature_c)

vars <- c("mw", "log_kow", "melting_point_k", "log_saq", "h_donor", "h_acceptor")

profile <- aggregate(dat[vars], list(compound_id = dat$compound_id), median)
endpoint <- aggregate(dat$log_kp, list(compound_id = dat$compound_id), mean)
names(endpoint)[2] <- "mean_log_kp"
profile <- merge(profile, endpoint, by = "compound_id", sort = FALSE)

centers <- sapply(profile[vars], mean)
scales <- sapply(profile[vars], sd)

for (x in vars) {
	profile[[paste0(x, "_z")]] <- (profile[[x]] - centers[x]) / scales[x]
	dat[[paste0(x, "_z")]] <- (dat[[x]] - centers[x]) / scales[x]
}

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

f1 <- c("mw_z", "log_kow_z")
f2 <- c("mw_z", "log_kow_z", "melting_point_k_z", "log_saq_z", "h_donor_z", "h_acceptor_z")

fit_spec <- function(terms, label) {
	fixed <- paste(terms, collapse = " + ")
	form_joint <- as.formula(paste("log_kp ~ temperature_c +", fixed, "+ (1 | compound_id) + (1 | reference_id)"))
	form_lm <- as.formula(paste("mean_log_kp ~", fixed))

	baseline_ml <- lmer(log_kp ~ temperature_c + (1 | compound_id) + (1 | reference_id), data = dat, REML = FALSE, control = ctrl)
	joint_ml <- lmer(form_joint, data = dat, REML = FALSE, control = ctrl)
	baseline_reml <- update(baseline_ml, REML = TRUE)
	joint_reml <- update(joint_ml, REML = TRUE)
	conventional <- lm(form_lm, data = profile)

	v0 <- as.data.frame(VarCorr(baseline_reml))
	v1 <- as.data.frame(VarCorr(joint_reml))
	v0 <- v0$vcov[v0$grp == "compound_id"]
	v1 <- v1$vcov[v1$grp == "compound_id"]

	ll0 <- logLik(baseline_ml)
	ll1 <- logLik(joint_ml)
	chisq <- 2 * (as.numeric(ll1) - as.numeric(ll0))
	df <- attr(ll1, "df") - attr(ll0, "df")

	list(
		label = label,
		conventional = conventional,
		joint = joint_reml,
		summary = data.frame(
			specification = label,
			n_compounds = nrow(profile),
			n_observations = nrow(dat),
			compound_variance_baseline = v0,
			compound_variance_joint = v1,
			compound_variance_explained = 100 * (1 - v1 / v0),
			AIC_baseline = AIC(baseline_ml),
			AIC_joint = AIC(joint_ml),
			BIC_baseline = BIC(baseline_ml),
			BIC_joint = BIC(joint_ml),
			LRT_chisq = chisq,
			LRT_df = df,
			LRT_p = pchisq(chisq, df, lower.tail = FALSE)
		)
	)
}

s1 <- fit_spec(f1, "S1")
s2 <- fit_spec(f2, "S2")

summary_out <- rbind(s1$summary, s2$summary)

coef_out <- do.call(rbind, lapply(list(S1 = s1, S2 = s2), function(x) {
	a <- coef(summary(x$conventional))
	b <- coef(summary(x$joint))
	rbind(
		data.frame(specification = x$label, model = "Conventional", term = rownames(a), estimate = a[, "Estimate"], se = a[, "Std. Error"]),
		data.frame(specification = x$label, model = "Hierarchical", term = rownames(b), estimate = b[, "Estimate"], se = b[, "Std. Error"])
	)
}))
rownames(coef_out) <- NULL

write.csv(summary_out, "R/results/20_model_summary.csv", row.names = FALSE)
write.csv(coef_out, "R/results/20_coefficients.csv", row.names = FALSE)

saveRDS(
	list(
		S1 = s1,
		S2 = s2,
		centers = centers,
		scales = scales,
		compound_profile = profile
	),
	"R/results/20_qspr_models.rds"
)

cat(
	"\nCheruvu QSPR models fitted\n",
	"S1 compound variance explained: ", round(summary_out$compound_variance_explained[summary_out$specification == "S1"], 1), "%\n",
	"S1 delta AIC: ", round(summary_out$AIC_joint[summary_out$specification == "S1"] - summary_out$AIC_baseline[summary_out$specification == "S1"], 2), "\n",
	"S2 compound variance explained: ", round(summary_out$compound_variance_explained[summary_out$specification == "S2"], 1), "%\n",
	"S2 delta AIC: ", round(summary_out$AIC_joint[summary_out$specification == "S2"] - summary_out$AIC_baseline[summary_out$specification == "S2"], 2), "\n",
	sep = ""
)