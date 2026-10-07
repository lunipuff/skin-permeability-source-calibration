# ============================================================
# 11_hierarchical_diagnostics.R
#
# Purpose:
# 	Assess diagnostics of the primary hierarchical model.
#
# Primary model:
#
# 	M2:
# 		logKp ~ temperature +
# 			(1 | compound) +
# 			(1 | reference)
#
# This script evaluates:
# 	- model convergence
# 	- singularity
# 	- residual distribution
# 	- residuals versus fitted values
# 	- residuals versus temperature
# 	- residual behavior by compound and reference
# 	- compound random effects
# 	- reference random effects
# 	- unusually large residuals
#
# This script DOES NOT:
# 	- remove observations
# 	- refit the model after exclusions
# 	- estimate standardized compound logKp
# 	- perform sensitivity analyses
#
# Run from the project root:
# 	source("R/11_hierarchical_diagnostics.R")
# ============================================================


source("R/00_config.R")


# ============================================================
# 1. Import fitted model and modeling data
# ============================================================

models <- readRDS(
	"R/results/10_hierarchical_models.rds"
)

dat <- readr::read_csv(
	"R/results/10_model_data.csv",
	show_col_types = FALSE
)

model <- models$M2


# ============================================================
# 2. Check model and data alignment
# ============================================================

if (nrow(dat) != nobs(model)) {
	stop(
		"Number of rows in model data does not match model observations."
	)
}


# ============================================================
# 3. Convergence and singularity
# ============================================================

optimizer_messages <- model@optinfo$conv$lme4$messages

if (is.null(optimizer_messages)) {
	optimizer_messages <- "None"
} else {
	optimizer_messages <- paste(
		optimizer_messages,
		collapse = "; "
	)
}

singular_model <- lme4::isSingular(
	model,
	tol = 1e-4
)

diagnostic_status <- tibble(
	check = c(
		"Convergence warning",
		"Singular fit"
	),

	result = c(
		optimizer_messages,
		as.character(singular_model)
	)
)

readr::write_csv(
	diagnostic_status,
	"R/results/11_model_status.csv"
)


# ============================================================
# 4. Observation-level diagnostics
# ============================================================

dat_diag <- dat %>%
	mutate(
		fitted = fitted(model),

		residual = resid(model),

		residual_sd = sigma(model),

		standardized_residual =
			residual /
			residual_sd,

		absolute_residual =
			abs(residual),

		absolute_standardized_residual =
			abs(standardized_residual)
	)


# ============================================================
# 5. Residual summary
# ============================================================

residual_summary <- tibble(
	metric = c(
		"Number of observations",
		"Mean residual",
		"Median residual",
		"Residual SD",
		"Median absolute residual",
		"Mean absolute residual",
		"Maximum absolute residual",
		"Observations with |standardized residual| >2",
		"Observations with |standardized residual| >3"
	),

	value = c(
		nrow(dat_diag),

		mean(
			dat_diag$residual,
			na.rm = TRUE
		),

		median(
			dat_diag$residual,
			na.rm = TRUE
		),

		sd(
			dat_diag$residual,
			na.rm = TRUE
		),

		median(
			dat_diag$absolute_residual,
			na.rm = TRUE
		),

		mean(
			dat_diag$absolute_residual,
			na.rm = TRUE
		),

		max(
			dat_diag$absolute_residual,
			na.rm = TRUE
		),

		sum(
			dat_diag$absolute_standardized_residual > 2,
			na.rm = TRUE
		),

		sum(
			dat_diag$absolute_standardized_residual > 3,
			na.rm = TRUE
		)
	)
)

readr::write_csv(
	residual_summary,
	"R/results/11_residual_summary.csv"
)


# ============================================================
# 6. Residual distribution statistics
# ============================================================

residual_mean <- mean(
	dat_diag$residual,
	na.rm = TRUE
)

residual_sd <- sd(
	dat_diag$residual,
	na.rm = TRUE
)

residual_skewness <- mean(
	(
		(
			dat_diag$residual -
				residual_mean
		) /
			residual_sd
	)^3,
	na.rm = TRUE
)

residual_kurtosis <- mean(
	(
		(
			dat_diag$residual -
				residual_mean
		) /
			residual_sd
	)^4,
	na.rm = TRUE
) - 3

residual_distribution <- tibble(
	metric = c(
		"Skewness",
		"Excess kurtosis"
	),

	value = c(
		residual_skewness,
		residual_kurtosis
	)
)

readr::write_csv(
	residual_distribution,
	"R/results/11_residual_distribution.csv"
)


# ============================================================
# 7. Residual trend checks
#
# These correlations are descriptive diagnostics only.
# They are not hypothesis tests.
# ============================================================

residual_trends <- tibble(
	comparison = c(
		"Residual vs fitted",
		"Absolute residual vs fitted",
		"Residual vs temperature",
		"Absolute residual vs temperature"
	),

	correlation = c(
		cor(
			dat_diag$residual,
			dat_diag$fitted,
			use = "complete.obs"
		),

		cor(
			dat_diag$absolute_residual,
			dat_diag$fitted,
			use = "complete.obs"
		),

		cor(
			dat_diag$residual,
			dat_diag$temperature_K,
			use = "complete.obs"
		),

		cor(
			dat_diag$absolute_residual,
			dat_diag$temperature_K,
			use = "complete.obs"
		)
	)
)

readr::write_csv(
	residual_trends,
	"R/results/11_residual_trends.csv"
)


# ============================================================
# 8. Flag observations with large residuals
#
# These observations are retained.
# This table is for inspection only.
# ============================================================

large_residuals <- dat_diag %>%
	filter(
		absolute_standardized_residual > 2
	) %>%
	select(
		source_row,
		original_no,
		compound_id,
		compound_name,
		reference_id,
		reference_raw,
		temperature_K,
		skin_thickness_mm,
		logKp,
		fitted,
		residual,
		standardized_residual,
		absolute_standardized_residual
	) %>%
	arrange(
		desc(
			absolute_standardized_residual
		)
	)

readr::write_csv(
	large_residuals,
	"R/results/11_large_residuals.csv"
)

# ============================================================
# 9. Residual summaries by compound
# ============================================================

compound_residuals <- dat_diag %>%
	group_by(
		compound_id,
		compound_name
	) %>%
	summarise(
		n_observations = n(),

		mean_residual = mean(
			residual,
			na.rm = TRUE
		),

		median_residual = median(
			residual,
			na.rm = TRUE
		),

		residual_sd = if_else(
			n() >= 2,
			sd(
				residual,
				na.rm = TRUE
			),
			NA_real_
		),

		mean_absolute_residual = mean(
			absolute_residual,
			na.rm = TRUE
		),

		max_absolute_residual = max(
			absolute_residual,
			na.rm = TRUE
		),

		.groups = "drop"
	) %>%
	arrange(
		desc(mean_absolute_residual)
	)

readr::write_csv(
	compound_residuals,
	"R/results/11_compound_residuals.csv"
)


# ============================================================
# 10. Residual summaries by reference
# ============================================================

reference_residuals <- dat_diag %>%
	group_by(
		reference_id,
		reference_raw
	) %>%
	summarise(
		n_observations = n(),

		mean_residual = mean(
			residual,
			na.rm = TRUE
		),

		median_residual = median(
			residual,
			na.rm = TRUE
		),

		residual_sd = if_else(
			n() >= 2,
			sd(
				residual,
				na.rm = TRUE
			),
			NA_real_
		),

		mean_absolute_residual = mean(
			absolute_residual,
			na.rm = TRUE
		),

		max_absolute_residual = max(
			absolute_residual,
			na.rm = TRUE
		),

		.groups = "drop"
	) %>%
	arrange(
		desc(mean_absolute_residual)
	)

readr::write_csv(
	reference_residuals,
	"R/results/11_reference_residuals.csv"
)


# ============================================================
# 11. Extract compound random effects
# ============================================================

compound_effects <- ranef(
	model,
	condVar = TRUE
)$compound_id

compound_post_var <- attr(
	compound_effects,
	"postVar"
)

compound_effects <- tibble(
	compound_id = rownames(
		compound_effects
	),

	random_effect = compound_effects[, 1],

	random_effect_se = sqrt(
		compound_post_var[
			1,
			1,
		]
	)
) %>%
	left_join(
		dat %>%
			distinct(
				compound_id,
				compound_name,
				cas
			),
		by = "compound_id"
	) %>%
	mutate(
		conf_low =
			random_effect -
			1.96 *
			random_effect_se,

		conf_high =
			random_effect +
			1.96 *
			random_effect_se
	) %>%
	arrange(
		random_effect
	)

readr::write_csv(
	compound_effects,
	"R/results/11_compound_random_effects.csv"
)


# ============================================================
# 12. Extract reference random effects
# ============================================================

reference_effects <- ranef(
	model,
	condVar = TRUE
)$reference_id

reference_post_var <- attr(
	reference_effects,
	"postVar"
)

reference_effects <- tibble(
	reference_id = rownames(
		reference_effects
	),

	random_effect = reference_effects[, 1],

	random_effect_se = sqrt(
		reference_post_var[
			1,
			1,
		]
	)
) %>%
	left_join(
		dat %>%
			distinct(
				reference_id,
				reference_raw
			),
		by = "reference_id"
	) %>%
	mutate(
		conf_low =
			random_effect -
			1.96 *
			random_effect_se,

		conf_high =
			random_effect +
			1.96 *
			random_effect_se
	) %>%
	arrange(
		random_effect
	)

readr::write_csv(
	reference_effects,
	"R/results/11_reference_random_effects.csv"
)


# ============================================================
# 13. Random-effect distribution summaries
# ============================================================

random_effect_summary <- bind_rows(

	tibble(
		level = "Compound",

		n_levels = nrow(
			compound_effects
		),

		mean = mean(
			compound_effects$random_effect
		),

		sd = sd(
			compound_effects$random_effect
		),

		min = min(
			compound_effects$random_effect
		),

		max = max(
			compound_effects$random_effect
		)
	),

	tibble(
		level = "Reference",

		n_levels = nrow(
			reference_effects
		),

		mean = mean(
			reference_effects$random_effect
		),

		sd = sd(
			reference_effects$random_effect
		),

		min = min(
			reference_effects$random_effect
		),

		max = max(
			reference_effects$random_effect
		)
	)
)

readr::write_csv(
	random_effect_summary,
	"R/results/11_random_effect_summary.csv"
)


# ============================================================
# 14. Q-Q data for later diagnostic figures
#
# Figure generation itself is reserved for 99b.
# ============================================================

residual_qq <- tibble(
	observed = sort(
		dat_diag$standardized_residual
	),

	expected = qnorm(
		ppoints(
			nrow(dat_diag)
		)
	)
)

readr::write_csv(
	residual_qq,
	"R/results/11_residual_qq.csv"
)


compound_qq <- tibble(
	observed = sort(
		compound_effects$random_effect
	),

	expected = qnorm(
		ppoints(
			nrow(
				compound_effects
			)
		)
	)
)

readr::write_csv(
	compound_qq,
	"R/results/11_compound_qq.csv"
)


reference_qq <- tibble(
	observed = sort(
		reference_effects$random_effect
	),

	expected = qnorm(
		ppoints(
			nrow(
				reference_effects
			)
		)
	)
)

readr::write_csv(
	reference_qq,
	"R/results/11_reference_qq.csv"
)


# ============================================================
# 15. Export observation-level diagnostic data
# ============================================================

readr::write_csv(
	dat_diag,
	"R/results/11_observation_diagnostics.csv"
)


# ============================================================
# 16. Console summary
# ============================================================

cat("\n")
cat("========================================\n")
cat("11_hierarchical_diagnostics.R complete\n")
cat("========================================\n\n")

cat("Model status:\n\n")

print(
	diagnostic_status
)

cat("\nResidual summary:\n\n")

print(
	residual_summary
)

cat("\nResidual distribution:\n\n")

print(
	residual_distribution
)

cat("\nResidual trend checks:\n\n")

print(
	residual_trends
)

cat("\nRandom-effect summary:\n\n")

print(
	random_effect_summary
)

cat("\nOutputs written to:\n")

cat(
	"  R/results/11_model_status.csv\n"
)

cat(
	"  R/results/11_residual_summary.csv\n"
)

cat(
	"  R/results/11_residual_distribution.csv\n"
)

cat(
	"  R/results/11_residual_trends.csv\n"
)

cat(
	"  R/results/11_large_residuals.csv\n"
)

cat(
	"  R/results/11_compound_residuals.csv\n"
)

cat(
	"  R/results/11_reference_residuals.csv\n"
)

cat(
	"  R/results/11_compound_random_effects.csv\n"
)

cat(
	"  R/results/11_reference_random_effects.csv\n"
)

cat(
	"  R/results/11_random_effect_summary.csv\n"
)

cat(
	"  R/results/11_residual_qq.csv\n"
)

cat(
	"  R/results/11_compound_qq.csv\n"
)

cat(
	"  R/results/11_reference_qq.csv\n"
)

cat(
	"  R/results/11_observation_diagnostics.csv\n"
)