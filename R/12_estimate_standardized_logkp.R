# ============================================================
# 12_estimate_standardized_logkp.R
#
# Purpose:
# 	Estimate reference- and temperature-standardized
# 	compound-level logKp from the primary hierarchical model.
#
# Primary model:
#
# 	M2:
# 		logKp ~ temperature +
# 			(1 | compound) +
# 			(1 | reference)
#
# Standardized compound logKp is defined as:
#
# 	fixed intercept + compound random effect
#
# This corresponds to:
# 	- temperature at the model centering value
# 	- reference random effect = 0
#
# The resulting quantity represents the estimated
# compound-level logKp after removing systematic reference
# effects and standardizing temperature.
#
# This script:
# 	- extracts compound-level model estimates
# 	- calculates approximate uncertainty
# 	- compares standardized estimates with raw compound means
# 	- quantifies shrinkage and adjustment
#
# This script DOES NOT:
# 	- refit the hierarchical model
# 	- remove observations
# 	- fit QSPR models
# 	- perform sensitivity analyses
#
# Run from the project root:
# 	source("R/12_estimate_standardized_logkp.R")
# ============================================================


source("R/00_config.R")


# ============================================================
# 1. Import model and data
# ============================================================

models <- readRDS(
	"R/results/10_hierarchical_models.rds"
)

model <- models$M2


dat <- readr::read_csv(
	"R/results/10_model_data.csv",
	show_col_types = FALSE
)


model_metadata <- readr::read_csv(
	"R/results/10_model_metadata.csv",
	show_col_types = FALSE
)


# ============================================================
# 2. Recover model centering temperature
# ============================================================

temperature_reference <- model_metadata %>%
	filter(
		quantity ==
			"Temperature centering value (K)"
	) %>%
	pull(value)


if (length(temperature_reference) != 1) {
	stop(
		"Could not recover a unique temperature centering value."
	)
}


# ============================================================
# 3. Extract fixed effects
# ============================================================

fixed_effects <- lme4::fixef(
	model
)

intercept <- unname(
	fixed_effects["(Intercept)"]
)

temperature_effect <- unname(
	fixed_effects["temperature_c"]
)


# ============================================================
# 4. Extract compound random effects and conditional variance
# ============================================================

compound_ranef <- ranef(
	model,
	condVar = TRUE
)$compound_id

compound_post_var <- attr(
	compound_ranef,
	"postVar"
)


compound_effects <- tibble(
	compound_id = rownames(
		compound_ranef
	),

	compound_random_effect =
		compound_ranef[, 1],

	compound_random_effect_se =
		sqrt(
			compound_post_var[
				1,
				1,
			]
		)
)


# ============================================================
# 5. Approximate uncertainty of standardized logKp
#
# standardized logKp =
# 	fixed intercept + compound random effect
#
# The uncertainty below combines:
# 	- fixed-intercept uncertainty
# 	- conditional uncertainty of the compound random effect
#
# Covariance between the two is not included, so these are
# approximate intervals rather than exact prediction intervals.
# ============================================================

fixed_vcov <- as.matrix(
	vcov(model)
)

intercept_se <- sqrt(
	fixed_vcov[
		"(Intercept)",
		"(Intercept)"
	]
)


compound_effects <- compound_effects %>%
	mutate(
		standardized_logKp =
			intercept +
			compound_random_effect,

		standardized_logKp_se =
			sqrt(
				intercept_se^2 +
				compound_random_effect_se^2
			),

		standardized_logKp_low =
			standardized_logKp -
			1.96 *
			standardized_logKp_se,

		standardized_logKp_high =
			standardized_logKp +
			1.96 *
			standardized_logKp_se
	)


# ============================================================
# 6. Compound-level descriptive information
# ============================================================

compound_observed <- dat %>%
	group_by(
		compound_id,
		compound_name,
		cas
	) %>%
	summarise(
		n_observations = n(),

		n_references = n_distinct(
			reference_id
		),

		n_temperatures = n_distinct(
			temperature_K
		),

		raw_mean_logKp = mean(
			logKp,
			na.rm = TRUE
		),

		raw_median_logKp = median(
			logKp,
			na.rm = TRUE
		),

		raw_sd_logKp = if_else(
			n() >= 2,
			sd(
				logKp,
				na.rm = TRUE
			),
			NA_real_
		),

		raw_min_logKp = min(
			logKp,
			na.rm = TRUE
		),

		raw_max_logKp = max(
			logKp,
			na.rm = TRUE
		),

		mean_temperature_K = mean(
			temperature_K,
			na.rm = TRUE
		),

		.groups = "drop"
	)


# ============================================================
# 7. Construct standardized compound endpoint
# ============================================================

standardized_logkp <- compound_effects %>%
	left_join(
		compound_observed,
		by = "compound_id"
	) %>%
	mutate(
		standardization_temperature_K =
			temperature_reference,

		standardization_temperature_C =
			temperature_reference -
			273.15,

		raw_vs_standardized_difference =
			standardized_logKp -
			raw_mean_logKp,

		absolute_raw_vs_standardized_difference =
			abs(
				raw_vs_standardized_difference
			),

		multi_reference =
			n_references >= 2
	) %>%
	select(
		compound_id,
		compound_name,
		cas,
		n_observations,
		n_references,
		n_temperatures,
		multi_reference,
		mean_temperature_K,
		standardization_temperature_K,
		standardization_temperature_C,
		raw_mean_logKp,
		raw_median_logKp,
		raw_sd_logKp,
		raw_min_logKp,
		raw_max_logKp,
		compound_random_effect,
		compound_random_effect_se,
		standardized_logKp,
		standardized_logKp_se,
		standardized_logKp_low,
		standardized_logKp_high,
		raw_vs_standardized_difference,
		absolute_raw_vs_standardized_difference
	) %>%
	arrange(
		standardized_logKp
	)


readr::write_csv(
	standardized_logkp,
	"R/results/12_standardized_logkp.csv"
)


# ============================================================
# 8. Compare raw and standardized compound estimates
# ============================================================

raw_vs_standardized_summary <- tibble(
	metric = c(
		"Number of compounds",
		"Correlation: raw mean vs standardized logKp",
		"Median absolute adjustment",
		"Mean absolute adjustment",
		"Maximum absolute adjustment",
		"Compounds adjusted by >0.25 log units",
		"Compounds adjusted by >0.50 log units",
		"Compounds adjusted by >1.00 log units"
	),

	value = c(
		nrow(
			standardized_logkp
		),

		cor(
			standardized_logkp$raw_mean_logKp,
			standardized_logkp$standardized_logKp,
			use = "complete.obs"
		),

		median(
			standardized_logkp$
				absolute_raw_vs_standardized_difference,
			na.rm = TRUE
		),

		mean(
			standardized_logkp$
				absolute_raw_vs_standardized_difference,
			na.rm = TRUE
		),

		max(
			standardized_logkp$
				absolute_raw_vs_standardized_difference,
			na.rm = TRUE
		),

		sum(
			standardized_logkp$
				absolute_raw_vs_standardized_difference >
				0.25,
			na.rm = TRUE
		),

		sum(
			standardized_logkp$
				absolute_raw_vs_standardized_difference >
				0.50,
			na.rm = TRUE
		),

		sum(
			standardized_logkp$
				absolute_raw_vs_standardized_difference >
				1.00,
			na.rm = TRUE
		)
	)
)


readr::write_csv(
	raw_vs_standardized_summary,
	"R/results/12_raw_vs_standardized_summary.csv"
)


# ============================================================
# 9. Largest adjustments
# ============================================================

largest_adjustments <- standardized_logkp %>%
	arrange(
		desc(
			absolute_raw_vs_standardized_difference
		)
	) %>%
	select(
		compound_id,
		compound_name,
		cas,
		n_observations,
		n_references,
		mean_temperature_K,
		raw_mean_logKp,
		standardized_logKp,
		raw_vs_standardized_difference,
		absolute_raw_vs_standardized_difference
	)


readr::write_csv(
	largest_adjustments,
	"R/results/12_largest_adjustments.csv"
)


# ============================================================
# 10. Quantify shrinkage
#
# Random-effect estimates are partially pooled toward zero.
# The amount of shrinkage is expected to be greater for
# compounds supported by fewer observations.
#
# This table is descriptive only.
# ============================================================

shrinkage_summary <- standardized_logkp %>%
	mutate(
		raw_centered_effect =
			raw_mean_logKp -
			intercept,

		shrinkage_difference =
			compound_random_effect -
			raw_centered_effect,

		absolute_shrinkage_difference =
			abs(
				shrinkage_difference
			)
	) %>%
	select(
		compound_id,
		compound_name,
		n_observations,
		n_references,
		raw_centered_effect,
		compound_random_effect,
		shrinkage_difference,
		absolute_shrinkage_difference,
		compound_random_effect_se
	) %>%
	arrange(
		desc(
			absolute_shrinkage_difference
		)
	)


readr::write_csv(
	shrinkage_summary,
	"R/results/12_shrinkage_summary.csv"
)


# ============================================================
# 11. Standardized endpoint uncertainty summary
# ============================================================

uncertainty_summary <- tibble(
	metric = c(
		"Median standardized logKp SE",
		"Mean standardized logKp SE",
		"Minimum standardized logKp SE",
		"Maximum standardized logKp SE",
		"Median 95% interval width"
	),

	value = c(
		median(
			standardized_logkp$
				standardized_logKp_se,
			na.rm = TRUE
		),

		mean(
			standardized_logkp$
				standardized_logKp_se,
			na.rm = TRUE
		),

		min(
			standardized_logkp$
				standardized_logKp_se,
			na.rm = TRUE
		),

		max(
			standardized_logkp$
				standardized_logKp_se,
			na.rm = TRUE
		),

		median(
			standardized_logkp$
				standardized_logKp_high -
				standardized_logkp$
					standardized_logKp_low,
			na.rm = TRUE
		)
	)
)


readr::write_csv(
	uncertainty_summary,
	"R/results/12_uncertainty_summary.csv"
)


# ============================================================
# 12. Export definition of standardized endpoint
# ============================================================

endpoint_definition <- tibble(
	quantity = c(
		"Primary hierarchical model",
		"Standardization temperature (K)",
		"Standardization temperature (C)",
		"Temperature coefficient",
		"Reference effect used for standardization",
		"Standardized endpoint definition"
	),

	value = c(
		"M2: logKp ~ temperature_c + (1|compound_id) + (1|reference_id)",

		as.character(
			temperature_reference
		),

		as.character(
			temperature_reference -
				273.15
		),

		as.character(
			temperature_effect
		),

		"0",

		"fixed intercept + compound random effect"
	)
)


readr::write_csv(
	endpoint_definition,
	"R/results/12_endpoint_definition.csv"
)


# ============================================================
# 13. Console summary
# ============================================================

cat("\n")
cat("========================================\n")
cat("12_estimate_standardized_logkp.R complete\n")
cat("========================================\n\n")

cat(
	"Standardization temperature: ",
	round(
		temperature_reference,
		2
	),
	" K (",
	round(
		temperature_reference -
			273.15,
		2
	),
	" C)\n\n",
	sep = ""
)

cat("Raw vs standardized endpoint:\n\n")

print(
	raw_vs_standardized_summary
)

cat("\nUncertainty summary:\n\n")

print(
	uncertainty_summary
)

cat("\nOutputs written to:\n")

cat(
	"  R/results/12_standardized_logkp.csv\n"
)

cat(
	"  R/results/12_raw_vs_standardized_summary.csv\n"
)

cat(
	"  R/results/12_largest_adjustments.csv\n"
)

cat(
	"  R/results/12_shrinkage_summary.csv\n"
)

cat(
	"  R/results/12_uncertainty_summary.csv\n"
)

cat(
	"  R/results/12_endpoint_definition.csv\n"
)