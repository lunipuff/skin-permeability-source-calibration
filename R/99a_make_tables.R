# ============================================================
# 99a_make_tables.R
#
# Purpose:
# 	Generate manuscript-ready tables from finalized analysis
# 	outputs.
#
# This script DOES NOT:
# 	- fit models
# 	- recalculate analyses
# 	- modify processed data
#
# Outputs:
# 	manuscript/tables/
#
# Run from project root:
# 	source("R/99a_make_tables.R")
# ============================================================


source("R/00_config.R")


# ============================================================
# 1. Create output directory
# ============================================================

dir.create(
	"manuscript/tables",
	recursive = TRUE,
	showWarnings = FALSE
)


# ============================================================
# 2. Helper functions
# ============================================================

format_number <- function(x, digits = 3) {
	ifelse(
		is.na(x),
		"",
		formatC(
			x,
			format = "f",
			digits = digits
		)
	)
}


format_p <- function(x) {
	ifelse(
		is.na(x),
		"",
		ifelse(
			x < 0.001,
			"<0.001",
			formatC(
				x,
				format = "f",
				digits = 3
			)
		)
	)
}


# ============================================================
# 3. Table 1: Dataset structure
# ============================================================

dataset_summary <- readr::read_csv(
	"R/results/02_dataset_summary.csv",
	show_col_types = FALSE
)


experimental_conditions <- readr::read_csv(
	"R/results/02_experimental_conditions.csv",
	show_col_types = FALSE
)


condition_variation <- readr::read_csv(
	"R/results/02_condition_variation.csv",
	show_col_types = FALSE
)


within_compound_summary <- readr::read_csv(
	"R/results/02_within_compound_summary.csv",
	show_col_types = FALSE
)


table_1a <- dataset_summary %>%
	transmute(
		section = "Dataset structure",
		metric,
		value = format_number(
			value,
			digits = 2
		)
	)


table_1b <- experimental_conditions %>%
	transmute(
		section = "Experimental conditions",
		metric = case_when(
			variable == "temperature_K" ~
				"Temperature, K",

			variable == "skin_thickness_mm" ~
				"Skin thickness, mm",

			TRUE ~
				variable
		),

		value = paste0(
			format_number(
				median,
				2
			),
			" median; ",
			format_number(
				min,
				2
			),
			"–",
			format_number(
				max,
				2
			),
			" range"
		)
	)


table_1c <- condition_variation %>%
	transmute(
		section = "Within-compound condition variation",
		metric = stringr::str_replace_all(
			metric,
			"_",
			" "
		),
		value = format_number(
			value,
			0
		)
	)


table_1d <- within_compound_summary %>%
	transmute(
		section = "Within-compound permeability variation",
		metric,
		value = if_else(
			stringr::str_detect(
				metric,
				"Compounds"
			),
			format_number(
				value,
				0
			),
			format_number(
				value,
				3
			)
		)
	)


table_1 <- bind_rows(
	table_1a,
	table_1b,
	table_1c,
	table_1d
)


readr::write_csv(
	table_1,
	"manuscript/tables/table_1_dataset_structure.csv"
)


# ============================================================
# 4. Table 2: Hierarchical model comparison
# ============================================================

hierarchical_comparison <- readr::read_csv(
	"R/results/10_model_comparison.csv",
	show_col_types = FALSE
)


hierarchical_lrt <- readr::read_csv(
	"R/results/10_likelihood_ratio_tests.csv",
	show_col_types = FALSE
)


table_2 <- hierarchical_comparison %>%
	transmute(
		model,
		description,
		AIC = round(
			AIC,
			1
		),
		BIC = round(
			BIC,
			1
		),
		delta_AIC = round(
			delta_AIC,
			1
		),
		delta_BIC = round(
			delta_BIC,
			1
		)
	)


readr::write_csv(
	table_2,
	"manuscript/tables/table_2_hierarchical_model_comparison.csv"
)


table_s1 <- hierarchical_lrt %>%
	transmute(
		comparison,
		chisq = round(
			chisq,
			2
		),
		df,
		p_value = format_p(
			p_value
		)
	)


readr::write_csv(
	table_s1,
	"manuscript/tables/table_s1_hierarchical_likelihood_tests.csv"
)


# ============================================================
# 5. Table 3: Primary hierarchical model
#
# M2 is the primary hierarchical model:
#
# 	logKp ~ temperature +
# 		(1 | compound) +
# 		(1 | reference)
# ============================================================

fixed_effects <- readr::read_csv(
	"R/results/10_fixed_effects.csv",
	show_col_types = FALSE
)


variance_proportions <- readr::read_csv(
	"R/results/10_variance_proportions.csv",
	show_col_types = FALSE
)


m2_fixed <- fixed_effects %>%
	filter(
		model == "M2"
	) %>%
	transmute(
		section = "Fixed effects",

		term = case_when(
			term == "(Intercept)" ~
				"Intercept",

			term == "temperature_c" ~
				"Temperature, per K",

			TRUE ~
				term
		),

		estimate = format_number(
			estimate,
			4
		),

		SE = format_number(
			std.error,
			4
		),

		CI_95 = paste0(
			format_number(
				conf.low,
				4
			),
			" to ",
			format_number(
				conf.high,
				4
			)
		),

		p_value = format_p(
			p.value
		),

		variance = "",
		SD = "",
		percent_variance = ""
	)


m2_variance <- variance_proportions %>%
	filter(
		model == "M2"
	) %>%
	transmute(
		section = "Variance components",
		term = component,
		estimate = "",
		SE = "",
		CI_95 = "",
		p_value = "",
		variance = format_number(
			variance,
			4
		),
		SD = format_number(
			sd,
			4
		),
		percent_variance = paste0(
			format_number(
				percent,
				1
			),
			"%"
		)
	)


table_3 <- bind_rows(
	m2_fixed,
	m2_variance
)


readr::write_csv(
	table_3,
	"manuscript/tables/table_3_primary_hierarchical_model.csv"
)


# ============================================================
# 6. Table 4: QSPR comparison
#
# Directly compares conventional compound-mean QSPR with
# hierarchical QSPR under the two fixed molecular
# specifications.
# ============================================================

qspr_sample <- readr::read_csv(
	"R/results/20_qspr_sample_summary.csv",
	show_col_types = FALSE
)


qspr_variance <- readr::read_csv(
	"R/results/20_compound_variance_reduction.csv",
	show_col_types = FALSE
)


hierarchy_advantage <- readr::read_csv(
	"R/results/20_hierarchy_advantage.csv",
	show_col_types = FALSE
)


table_4 <- hierarchy_advantage %>%
	left_join(
		qspr_sample %>%
			select(
				specification,
				n_compounds,
				n_observations
			),
		by = "specification"
	) %>%
	left_join(
		qspr_variance %>%
			select(
				specification,
				percent_compound_variance_explained
			),
		by = "specification"
	) %>%
	transmute(
		specification,
		model = specification_label,
		n_compounds,
		n_observations,

		compound_variance_explained =
			paste0(
				format_number(
					percent_compound_variance_explained,
					1
				),
				"%"
			),

		compound_mean_RMSE =
			format_number(
				`RMSE_Compound-mean QSPR`,
				3
			),

		hierarchical_RMSE =
			format_number(
				`RMSE_Joint hierarchical QSPR`,
				3
			),

		RMSE_improvement =
			paste0(
				format_number(
					RMSE_improvement_percent,
					1
				),
				"%"
			),

		compound_mean_MAE =
			format_number(
				`MAE_Compound-mean QSPR`,
				3
			),

		hierarchical_MAE =
			format_number(
				`MAE_Joint hierarchical QSPR`,
				3
			),

		MAE_improvement =
			paste0(
				format_number(
					MAE_improvement_percent,
					1
				),
				"%"
			),

		compound_mean_R2 =
			format_number(
				`R2_pred_Compound-mean QSPR`,
				3
			),

		hierarchical_R2 =
			format_number(
				`R2_pred_Joint hierarchical QSPR`,
				3
			)
	)


readr::write_csv(
	table_4,
	"manuscript/tables/table_4_qspr_performance.csv"
)


# ============================================================
# 7. Table 5: Sensitivity analyses
# ============================================================

sensitivity_summary <- readr::read_csv(
	"R/results/30_sensitivity_summary.csv",
	show_col_types = FALSE
)


sensitivity_balanced <- readr::read_csv(
	"R/results/30_sensitivity_compound_balanced.csv",
	show_col_types = FALSE
)


table_5 <- sensitivity_summary %>%
	left_join(
		sensitivity_balanced %>%
			select(
				scenario,
				compound_mean_mean_RMSE,
				joint_mean_RMSE,
				RMSE_improvement_percent
			) %>%
			rename(
				compound_balanced_raw_RMSE =
					compound_mean_mean_RMSE,

				compound_balanced_joint_RMSE =
					joint_mean_RMSE,

				compound_balanced_RMSE_improvement =
					RMSE_improvement_percent
			),
		by = "scenario"
	) %>%
	transmute(
		scenario,

		compound_variance_explained =
			paste0(
				format_number(
					percent_compound_variance_explained,
					1
				),
				"%"
			),

		raw_RMSE =
			format_number(
				raw_RMSE,
				3
			),

		joint_RMSE =
			format_number(
				joint_RMSE,
				3
			),

		RMSE_improvement =
			paste0(
				format_number(
					RMSE_improvement_percent,
					1
				),
				"%"
			),

		compound_balanced_raw_RMSE =
			format_number(
				compound_balanced_raw_RMSE,
				3
			),

		compound_balanced_joint_RMSE =
			format_number(
				compound_balanced_joint_RMSE,
				3
			),

		compound_balanced_improvement =
			paste0(
				format_number(
					compound_balanced_RMSE_improvement,
					1
				),
				"%"
			)
	)


readr::write_csv(
	table_5,
	"manuscript/tables/table_5_sensitivity_analyses.csv"
)


# ============================================================
# 8. Supplementary table: QSPR coefficients
# ============================================================

qspr_coefficients <- readr::read_csv(
	"R/results/20_qspr_coefficients.csv",
	show_col_types = FALSE
)


table_s2 <- qspr_coefficients %>%
	transmute(
		specification,
		model = specification_label,
		approach,
		term,
		estimate = format_number(
			estimate,
			4
		),
		SE = format_number(
			std.error,
			4
		),
		CI_95 = paste0(
			format_number(
				conf.low,
				4
			),
			" to ",
			format_number(
				conf.high,
				4
			)
		),
		p_value = format_p(
			p.value
		)
	)


readr::write_csv(
	table_s2,
	"manuscript/tables/table_s2_qspr_coefficients.csv"
)


# ============================================================
# 9. Supplementary table: Reference random effects
# ============================================================

reference_effects <- readr::read_csv(
	"R/results/11_reference_random_effects.csv",
	show_col_types = FALSE
)


table_s3 <- reference_effects %>%
	arrange(
		desc(
			abs(
				random_effect
			)
		)
	) %>%
	transmute(
		reference = reference_raw,
		effect = format_number(
			random_effect,
			3
		),
		SE = format_number(
			random_effect_se,
			3
		),
		CI_95 = paste0(
			format_number(
				conf_low,
				3
			),
			" to ",
			format_number(
				conf_high,
				3
			)
		)
	)


readr::write_csv(
	table_s3,
	"manuscript/tables/table_s3_reference_effects.csv"
)


# ============================================================
# 10. Supplementary table: Large residuals
# ============================================================

large_residuals <- readr::read_csv(
	"R/results/11_large_residuals.csv",
	show_col_types = FALSE
)


table_s4 <- large_residuals %>%
	transmute(
		compound = compound_name,
		reference = reference_raw,
		temperature_K =
			format_number(
				temperature_K,
				1
			),
		observed_logKp =
			format_number(
				logKp,
				3
			),
		fitted_logKp =
			format_number(
				fitted,
				3
			),
		residual =
			format_number(
				residual,
				3
			),
		standardized_residual =
			format_number(
				standardized_residual,
				2
			)
	)


readr::write_csv(
	table_s4,
	"manuscript/tables/table_s4_large_residuals.csv"
)


# ============================================================
# 11. Table manifest
# ============================================================

table_manifest <- tibble(
	file = c(
		"table_1_dataset_structure.csv",
		"table_2_hierarchical_model_comparison.csv",
		"table_3_primary_hierarchical_model.csv",
		"table_4_qspr_performance.csv",
		"table_5_sensitivity_analyses.csv",
		"table_s1_hierarchical_likelihood_tests.csv",
		"table_s2_qspr_coefficients.csv",
		"table_s3_reference_effects.csv",
		"table_s4_large_residuals.csv"
	),

	title = c(
		"Dataset structure and repeated-measurement characteristics",
		"Comparison of hierarchical permeability models",
		"Fixed effects and variance decomposition of the primary hierarchical model",
		"Compound-disjoint predictive performance of fixed QSPR specifications",
		"Sensitivity analyses",
		"Sequential likelihood-ratio tests for hierarchical models",
		"QSPR coefficient estimates",
		"Reference-level random effects",
		"Observations with large standardized residuals"
	)
)


readr::write_csv(
	table_manifest,
	"manuscript/tables/table_manifest.csv"
)


# ============================================================
# 12. Console summary
# ============================================================

cat("\n")
cat("========================================\n")
cat("99a_make_tables.R complete\n")
cat("========================================\n\n")

cat(
	"Tables written to manuscript/tables/\n"
)

print(
	table_manifest
)