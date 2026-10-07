# ============================================================
# 10_fit_hierarchical_models.R
#
# Purpose:
# 	Fit the core hierarchical models for repeated human
# 	skin-permeability measurements.
#
# Models:
#
# 	M0:
# 		logKp ~ 1 + (1 | compound)
#
# 	M1:
# 		logKp ~ 1 + (1 | compound) + (1 | reference)
#
# 	M2:
# 		logKp ~ temperature +
# 			(1 | compound) + (1 | reference)
#
# 	M3:
# 		logKp ~ temperature + skin thickness +
# 			(1 | compound) + (1 | reference)
#
# This script:
# 	- centers experimental covariates
# 	- fits nested hierarchical models
# 	- compares models using maximum likelihood
# 	- reports fixed effects from REML fits
# 	- estimates variance components
# 	- calculates variance proportions
# 	- saves fitted models for downstream scripts
#
# This script DOES NOT:
# 	- perform detailed model diagnostics
# 	- remove influential observations
# 	- estimate standardized compound logKp
# 	- perform sensitivity analyses
#
# Run from the project root:
# 	source("R/10_fit_hierarchical_models.R")
# ============================================================


source("R/00_config.R")


# ============================================================
# 1. Import prepared data
# ============================================================

dat <- readr::read_csv(
	"data/processed/01_permeability_data.csv",
	show_col_types = FALSE
)


# ============================================================
# 2. Check required variables
# ============================================================

required_variables <- c(
	"compound_id",
	"reference_id",
	"logKp",
	"temperature_K",
	"skin_thickness_mm"
)

missing_variables <- setdiff(
	required_variables,
	names(dat)
)

if (length(missing_variables) > 0) {
	stop(
		"Required variables missing: ",
		paste(
			missing_variables,
			collapse = ", "
		)
	)
}


# ============================================================
# 3. Prepare modeling variables
# ============================================================

dat <- dat %>%
	mutate(
		compound_id = factor(compound_id),
		reference_id = factor(reference_id)
	)


# Center temperature and skin thickness around their
# observed means.
#
# This makes the model intercept correspond approximately
# to an observation under average experimental conditions.

temperature_reference <- mean(
	dat$temperature_K,
	na.rm = TRUE
)

thickness_reference <- mean(
	dat$skin_thickness_mm,
	na.rm = TRUE
)

dat <- dat %>%
	mutate(
		temperature_c =
			temperature_K -
			temperature_reference,

		skin_thickness_c =
			skin_thickness_mm -
			thickness_reference
	)


# ============================================================
# 4. Modeling metadata
# ============================================================

model_metadata <- tibble(
	quantity = c(
		"Number of observations",
		"Number of compounds",
		"Number of references",
		"Temperature centering value (K)",
		"Skin thickness centering value (mm)"
	),

	value = c(
		nrow(dat),
		nlevels(dat$compound_id),
		nlevels(dat$reference_id),
		temperature_reference,
		thickness_reference
	)
)

readr::write_csv(
	model_metadata,
	"R/results/10_model_metadata.csv"
)


# ============================================================
# 5. Model-control settings
# ============================================================

model_control <- lme4::lmerControl(
	optimizer = "bobyqa",
	optCtrl = list(
		maxfun = 100000
	)
)


# ============================================================
# 6. Fit REML models
#
# REML fits are used for parameter estimation and variance
# decomposition.
# ============================================================

m0 <- lmerTest::lmer(
	logKp ~
		1 +
		(1 | compound_id),
	data = dat,
	REML = TRUE,
	control = model_control
)


m1 <- lmerTest::lmer(
	logKp ~
		1 +
		(1 | compound_id) +
		(1 | reference_id),
	data = dat,
	REML = TRUE,
	control = model_control
)


m2 <- lmerTest::lmer(
	logKp ~
		temperature_c +
		(1 | compound_id) +
		(1 | reference_id),
	data = dat,
	REML = TRUE,
	control = model_control
)


m3 <- lmerTest::lmer(
	logKp ~
		temperature_c +
		skin_thickness_c +
		(1 | compound_id) +
		(1 | reference_id),
	data = dat,
	REML = TRUE,
	control = model_control
)


# ============================================================
# 7. Fit ML versions for model comparison
#
# Models with different fixed effects should be compared
# using maximum likelihood rather than REML.
# ============================================================

m0_ml <- update(
	m0,
	REML = FALSE
)

m1_ml <- update(
	m1,
	REML = FALSE
)

m2_ml <- update(
	m2,
	REML = FALSE
)

m3_ml <- update(
	m3,
	REML = FALSE
)


# ============================================================
# 8. Model-comparison table
# ============================================================

model_comparison <- tibble(
	model = c(
		"M0",
		"M1",
		"M2",
		"M3"
	),

	description = c(
		"Compound random effect",
		"Compound + reference random effects",
		"Compound + reference + temperature",
		"Compound + reference + temperature + skin thickness"
	),

	n_parameters = c(
		attr(
			logLik(m0_ml),
			"df"
		),

		attr(
			logLik(m1_ml),
			"df"
		),

		attr(
			logLik(m2_ml),
			"df"
		),

		attr(
			logLik(m3_ml),
			"df"
		)
	),

	logLik = c(
		as.numeric(
			logLik(m0_ml)
		),

		as.numeric(
			logLik(m1_ml)
		),

		as.numeric(
			logLik(m2_ml)
		),

		as.numeric(
			logLik(m3_ml)
		)
	),

	AIC = c(
		AIC(m0_ml),
		AIC(m1_ml),
		AIC(m2_ml),
		AIC(m3_ml)
	),

	BIC = c(
		BIC(m0_ml),
		BIC(m1_ml),
		BIC(m2_ml),
		BIC(m3_ml)
	)
) %>%
	mutate(
		delta_AIC =
			AIC -
			min(AIC),

		delta_BIC =
			BIC -
			min(BIC)
	)

readr::write_csv(
	model_comparison,
	"R/results/10_model_comparison.csv"
)


# ============================================================
# 9. Sequential likelihood-ratio comparisons
#
# These comparisons are useful alongside AIC/BIC.
#
# The M0 vs M1 comparison introduces a random-effect variance
# parameter at the boundary of the parameter space, so its
# conventional chi-square p-value should be interpreted
# cautiously.
# ============================================================

lrt_01 <- anova(
	m0_ml,
	m1_ml
)

lrt_12 <- anova(
	m1_ml,
	m2_ml
)

lrt_23 <- anova(
	m2_ml,
	m3_ml
)


extract_lrt <- function(x, comparison) {

	row <- x[nrow(x), , drop = FALSE]

	tibble(
		comparison = comparison,

		chisq = if (
			"Chisq" %in% names(row)
		) {
			as.numeric(row$Chisq)
		} else {
			NA_real_
		},

		df = if (
			"Chi Df" %in% names(row)
		) {
			as.numeric(row[["Chi Df"]])
		} else {
			NA_real_
		},

		p_value = if (
			"Pr(>Chisq)" %in% names(row)
		) {
			as.numeric(row[["Pr(>Chisq)"]])
		} else {
			NA_real_
		}
	)
}


likelihood_ratio_tests <- bind_rows(
	extract_lrt(
		lrt_01,
		"M0 vs M1"
	),

	extract_lrt(
		lrt_12,
		"M1 vs M2"
	),

	extract_lrt(
		lrt_23,
		"M2 vs M3"
	)
)

readr::write_csv(
	likelihood_ratio_tests,
	"R/results/10_likelihood_ratio_tests.csv"
)


# ============================================================
# 10. Fixed effects from REML models
# ============================================================

fixed_effects <- bind_rows(

	broom.mixed::tidy(
		m0,
		effects = "fixed",
		conf.int = TRUE
	) %>%
		mutate(
			model = "M0"
		),

	broom.mixed::tidy(
		m1,
		effects = "fixed",
		conf.int = TRUE
	) %>%
		mutate(
			model = "M1"
		),

	broom.mixed::tidy(
		m2,
		effects = "fixed",
		conf.int = TRUE
	) %>%
		mutate(
			model = "M2"
		),

	broom.mixed::tidy(
		m3,
		effects = "fixed",
		conf.int = TRUE
	) %>%
		mutate(
			model = "M3"
		)
) %>%
	select(
		model,
		term,
		estimate,
		std.error,
		statistic,
		df,
		p.value,
		conf.low,
		conf.high
	)

readr::write_csv(
	fixed_effects,
	"R/results/10_fixed_effects.csv"
)


# ============================================================
# 11. Variance components
# ============================================================

extract_variance <- function(model, model_name) {

	as.data.frame(
		VarCorr(model)
	) %>%
		transmute(
			model = model_name,

			component = case_when(
				grp == "compound_id" ~
					"Compound",

				grp == "reference_id" ~
					"Reference",

				grp == "Residual" ~
					"Residual",

				TRUE ~
					grp
			),

			variance = vcov,

			sd = sdcor
		)
}


variance_components <- bind_rows(
	extract_variance(
		m0,
		"M0"
	),

	extract_variance(
		m1,
		"M1"
	),

	extract_variance(
		m2,
		"M2"
	),

	extract_variance(
		m3,
		"M3"
	)
)

readr::write_csv(
	variance_components,
	"R/results/10_variance_components.csv"
)


# ============================================================
# 12. Variance proportions
#
# These describe the proportion of random/residual variance
# attributed to each component within each model.
#
# Fixed-effect variance is not included in this denominator.
# ============================================================

variance_proportions <- variance_components %>%
	group_by(model) %>%
	mutate(
		total_random_residual_variance =
			sum(
				variance,
				na.rm = TRUE
			),

		proportion =
			variance /
			total_random_residual_variance,

		percent =
			100 *
			proportion
	) %>%
	ungroup()

readr::write_csv(
	variance_proportions,
	"R/results/10_variance_proportions.csv"
)


# ============================================================
# 13. Compact M3 summary
#
# M3 is currently the full prespecified hierarchical model.
# Whether it should remain the primary model depends on the
# results of this comparison and later diagnostics.
# ============================================================

m3_fixed <- fixed_effects %>%
	filter(
		model == "M3"
	)

m3_variance <- variance_proportions %>%
	filter(
		model == "M3"
	)

readr::write_csv(
	m3_fixed,
	"R/results/10_full_model_fixed_effects.csv"
)

readr::write_csv(
	m3_variance,
	"R/results/10_full_model_variance.csv"
)


# ============================================================
# 14. Save models
# ============================================================

saveRDS(
	list(
		M0 = m0,
		M1 = m1,
		M2 = m2,
		M3 = m3
	),
	"R/results/10_hierarchical_models.rds"
)

saveRDS(
	list(
		M0 = m0_ml,
		M1 = m1_ml,
		M2 = m2_ml,
		M3 = m3_ml
	),
	"R/results/10_hierarchical_models_ml.rds"
)


# ============================================================
# 15. Save modeling data
#
# This preserves the exact centered variables used by the
# fitted models for downstream diagnostics.
# ============================================================

readr::write_csv(
	dat,
	"R/results/10_model_data.csv"
)


# ============================================================
# 16. Console summary
# ============================================================

cat("\n")
cat("========================================\n")
cat("10_fit_hierarchical_models.R complete\n")
cat("========================================\n\n")

cat("Model comparison:\n\n")

print(
	model_comparison %>%
		select(
			model,
			AIC,
			BIC,
			delta_AIC,
			delta_BIC
		)
)

cat("\nSequential likelihood-ratio tests:\n\n")

print(
	likelihood_ratio_tests
)

cat("\nFull-model fixed effects:\n\n")

print(
	m3_fixed %>%
		select(
			term,
			estimate,
			std.error,
			p.value,
			conf.low,
			conf.high
		)
)

cat("\nFull-model variance decomposition:\n\n")

print(
	m3_variance %>%
		select(
			component,
			variance,
			sd,
			percent
		)
)

cat("\nOutputs written to:\n")

cat(
	"  R/results/10_model_metadata.csv\n"
)

cat(
	"  R/results/10_model_comparison.csv\n"
)

cat(
	"  R/results/10_likelihood_ratio_tests.csv\n"
)

cat(
	"  R/results/10_fixed_effects.csv\n"
)

cat(
	"  R/results/10_variance_components.csv\n"
)

cat(
	"  R/results/10_variance_proportions.csv\n"
)

cat(
	"  R/results/10_full_model_fixed_effects.csv\n"
)

cat(
	"  R/results/10_full_model_variance.csv\n"
)

cat(
	"  R/results/10_hierarchical_models.rds\n"
)

cat(
	"  R/results/10_hierarchical_models_ml.rds\n"
)

cat(
	"  R/results/10_model_data.csv\n"
)