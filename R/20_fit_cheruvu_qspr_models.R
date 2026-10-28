# 20_fit_cheruvu_qspr_models.R

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

baseline_ml <- lmer(log_kp ~ temperature_c + (1 | compound_id) + (1 | reference_id), 
				data = dat, REML = FALSE, control = ctrl)
baseline_reml <- update(baseline_ml, REML = TRUE)

vc <- as.data.frame(VarCorr(baseline_reml))
v0 <- vc$vcov[vc$grp == "compound_id"]
ll0 <- logLik(baseline_ml)
aic0 <- AIC(baseline_ml)
bic0 <- BIC(baseline_ml)

fit_spec <- function(terms, label) {
	fixed <- paste(terms, collapse = " + ")
	form_joint <- as.formula(paste("log_kp ~ temperature_c +", fixed, "+ (1 | compound_id) + (1 | reference_id)"))
	form_lm <- as.formula(paste("mean_log_kp ~", fixed))
	joint_ml <- lmer(form_joint, data = dat, REML = FALSE, control = ctrl)
	joint_reml <- update(joint_ml, REML = TRUE)
	conventional <- lm(form_lm, data = profile)
	vc <- as.data.frame(VarCorr(joint_reml))
	v1 <- vc$vcov[vc$grp == "compound_id"]
	ll1 <- logLik(joint_ml)
	chisq <- 2 * (as.numeric(ll1) - as.numeric(ll0))
	df <- attr(ll1, "df") - attr(ll0, "df")
	a <- coef(summary(conventional))
	b <- coef(summary(joint_reml))
	list(
		summary = data.frame(
			specification = label,
			n_compounds = nrow(profile),
			n_observations = nrow(dat),
			compound_variance_baseline = v0,
			compound_variance_joint = v1,
			compound_variance_explained = 100 * (1 - v1 / v0),
			AIC_baseline = aic0,
			AIC_joint = AIC(joint_ml),
			BIC_baseline = bic0,
			BIC_joint = BIC(joint_ml),
			LRT_chisq = chisq,
			LRT_df = df,
			LRT_p = pchisq(chisq, df, lower.tail = FALSE)
		),
		coefficients = rbind(
			data.frame(specification = label, model = "Conventional", term = rownames(a), estimate = a[, "Estimate"], se = a[, "Std. Error"]),
			data.frame(specification = label, model = "Hierarchical", term = rownames(b), estimate = b[, "Estimate"], se = b[, "Std. Error"])
		)
	)
}

s1 <- fit_spec(c("mw_z", "log_kow_z"), "S1")
s2 <- fit_spec(paste0(vars, "_z"), "S2")

summary_out <- rbind(s1$summary, s2$summary)

coef_out <- rbind(s1$coefficients, s2$coefficients)
rownames(coef_out) <- NULL

write.csv(summary_out, "R/results/20_model_summary.csv", row.names = FALSE)
write.csv(coef_out, "R/results/20_coefficients.csv", row.names = FALSE)

cat(
	"\nCheruvu QSPR models fitted\n",
	"S1 compound variance explained: ", round(s1$summary$compound_variance_explained, 1), "%\n",
	"S1 delta AIC: ", round(s1$summary$AIC_joint - aic0, 2), "\n",
	"S2 compound variance explained: ", round(s2$summary$compound_variance_explained, 1), "%\n",
	"S2 delta AIC: ", round(s2$summary$AIC_joint - aic0, 2), "\n",
	sep = ""
)