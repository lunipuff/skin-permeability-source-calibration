############################################################
# 13_export_final_model_outputs.R
# Export final manuscript-ready model outputs
############################################################

source("R/00_config.R")

############################################################
# Create output directories if needed
############################################################

dir.create("results/model_comparison", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

############################################################
# Load cleaned dataset and selected model
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

selected_model_text <- readLines(
	path_loco_cv_selected_model,
	warn = FALSE
)

selected_model_text <- selected_model_text[1]

selected_formula <- as.formula(
	selected_model_text
)

############################################################
# User settings
############################################################

outcome_col <- "logkpl"
id_col <- "compound_id"

selected_predictors <- get_formula_predictors(
	selected_formula,
	outcome_col = outcome_col
)

############################################################
# Check required columns
############################################################

required_cols <- unique(c(
	id_col,
	"CAS.No",
	"Compound",
	outcome_col,
	selected_predictors
))

missing_cols <- required_cols[
	!(required_cols %in% names(df))
]

if (length(missing_cols) > 0) {
	stop(
		paste(
			"Missing required columns:",
			paste(missing_cols, collapse = ", ")
		)
	)
}

############################################################
# Prepare selected-model complete-case dataset
############################################################

model_df <- df[
	complete.cases(df[, required_cols, drop = FALSE]),
	,
	drop = FALSE
]

model_df$row_id <- seq_len(
	nrow(model_df)
)

############################################################
# Fit final selected model on full complete-case dataset
############################################################

final_fit <- lm(
	selected_formula,
	data = model_df
)

final_fit_summary <- summary(
	final_fit
)

############################################################
# Helper functions
############################################################

safe_read_csv <- function(path) {
	if (file.exists(path)) {
		return(
			read.csv(
				path,
				stringsAsFactors = FALSE
			)
		)
	}

	NULL
}

safe_value <- function(df,
					   row_filter,
					   value_col) {

	if (is.null(df)) {
		return(NA)
	}

	if (!(value_col %in% names(df))) {
		return(NA)
	}

	row_index <- which(row_filter)

	if (length(row_index) == 0) {
		return(NA)
	}

	df[[value_col]][row_index[1]]
}

format_numeric <- function(x, digits = 3) {
	if (length(x) == 0 || is.na(x)) {
		return(NA)
	}

	formatC(
		as.numeric(x),
		digits = digits,
		format = "f"
	)
}

safe_rbind <- function(x) {
	x <- x[!sapply(x, is.null)]

	if (length(x) == 0) {
		return(data.frame())
	}

	all_cols <- unique(
		unlist(
			lapply(x, names)
		)
	)

	x2 <- lapply(
		x,
		function(d) {
			missing_cols <- setdiff(
				all_cols,
				names(d)
			)

			if (length(missing_cols) > 0) {
				for (m in missing_cols) {
					d[[m]] <- NA
				}
			}

			d[
				,
				all_cols,
				drop = FALSE
			]
		}
	)

	do.call(
		rbind,
		x2
	)
}

############################################################
# Export coefficient table
############################################################

coefficient_matrix <- coef(
	final_fit_summary
)

coefficient_table <- data.frame(
	term = rownames(coefficient_matrix),
	estimate = coefficient_matrix[, "Estimate"],
	standard_error = coefficient_matrix[, "Std. Error"],
	t_value = coefficient_matrix[, "t value"],
	p_value = coefficient_matrix[, "Pr(>|t|)"],
	stringsAsFactors = FALSE,
	row.names = NULL
)

coefficient_table$abs_estimate <- abs(
	coefficient_table$estimate
)

write.csv(
	coefficient_table,
	path_final_model_coefficients,
	row.names = FALSE
)

############################################################
# Export final equation text
############################################################

equation_lines <- c(
	"Selected interpretable QSPR model",
	"",
	paste0("Formula: ", selected_model_text),
	"",
	"Fitted coefficients:",
	paste0(
		coefficient_table$term,
		" = ",
		signif(coefficient_table$estimate, 6)
	)
)

writeLines(
	equation_lines,
	path_final_model_equation
)

############################################################
# Export final model summary
############################################################

final_model_summary <- data.frame(
	item = c(
		"selected_model_formula",
		"n_observations",
		"n_compounds",
		"n_predictors",
		"predictors",
		"n_terms_excluding_intercept",
		"residual_standard_error",
		"multiple_R2",
		"adjusted_R2",
		"AIC",
		"BIC"
	),
	value = c(
		selected_model_text,
		nrow(model_df),
		length(unique(model_df[[id_col]])),
		length(selected_predictors),
		paste(selected_predictors, collapse = ", "),
		length(attr(terms(selected_formula), "term.labels")),
		final_fit_summary$sigma,
		final_fit_summary$r.squared,
		final_fit_summary$adj.r.squared,
		AIC(final_fit),
		BIC(final_fit)
	),
	stringsAsFactors = FALSE
)

write.csv(
	final_model_summary,
	path_final_model_summary,
	row.names = FALSE
)

############################################################
# Load downstream analysis tables
############################################################

benchmark_main <- safe_read_csv(
	path_table_benchmark_main
)

benchmark_all <- safe_read_csv(
	path_table_benchmark_all
)

ablation_summary <- safe_read_csv(
	path_table_ablation_summary
)

validation_scheme_summary <- safe_read_csv(
	path_validation_scheme_summary
)

ad_error_summary <- safe_read_csv(
	path_table_ad_error_by_domain
)

ad_descriptor_summary <- safe_read_csv(
	path_table_ad_descriptor_summary
)

hetero_summary <- safe_read_csv(
	path_table_hetero_cv_summary
)

rdkit_summary <- safe_read_csv(
	path_table_rdkit_benchmark_summary
)

cleaning_summary <- safe_read_csv(
	path_cleaning_summary
)

solubility_redundancy <- safe_read_csv(
	path_table_solubility_redundancy
)

############################################################
# Extract key results
############################################################

selected_benchmark_row <- NULL

if (!is.null(benchmark_main) && "model" %in% names(benchmark_main)) {
	selected_benchmark_row <- benchmark_main[
		benchmark_main$model == "selected_model",
	]
}

if (is.null(selected_benchmark_row) || nrow(selected_benchmark_row) == 0) {
	if (!is.null(benchmark_all) && "model" %in% names(benchmark_all)) {
		selected_benchmark_row <- benchmark_all[
			benchmark_all$model == "selected_model",
		]
	}
}

selected_rmse <- safe_value(
	selected_benchmark_row,
	rep(TRUE, nrow(selected_benchmark_row)),
	"RMSE"
)

selected_mae <- safe_value(
	selected_benchmark_row,
	rep(TRUE, nrow(selected_benchmark_row)),
	"MAE"
)

selected_r2 <- safe_value(
	selected_benchmark_row,
	rep(TRUE, nrow(selected_benchmark_row)),
	"R2_pred"
)

rf_selected_rmse <- safe_value(
	benchmark_main,
	benchmark_main$model == "random_forest_selected_predictors",
	"RMSE"
)

rdkit_rf_exp_rmse <- safe_value(
	benchmark_main,
	benchmark_main$model == "rdkit_rf_with_experimental_as_reported",
	"RMSE"
)

rdkit_rf_mol_rmse <- safe_value(
	benchmark_main,
	benchmark_main$model == "rdkit_rf_molecular_only_as_reported",
	"RMSE"
)

loco_rmse <- safe_value(
	validation_scheme_summary,
	validation_scheme_summary$validation_scheme == "leave_one_compound_out",
	"RMSE"
)

rowwise_rmse <- safe_value(
	validation_scheme_summary,
	validation_scheme_summary$validation_scheme == "rowwise_cv_repeats",
	"RMSE"
)

rowwise_rmse_sd <- safe_value(
	validation_scheme_summary,
	validation_scheme_summary$validation_scheme == "rowwise_cv_repeats",
	"RMSE_sd"
)

rowwise_leakage <- safe_value(
	validation_scheme_summary,
	validation_scheme_summary$validation_scheme == "rowwise_cv_repeats",
	"proportion_compound_seen_in_training"
)

central_rmse <- safe_value(
	ad_error_summary,
	ad_error_summary$group == "central",
	"RMSE"
)

broad_rmse <- safe_value(
	ad_error_summary,
	ad_error_summary$group == "broad",
	"RMSE"
)

outside_rmse <- safe_value(
	ad_error_summary,
	ad_error_summary$group == "outside",
	"RMSE"
)

central_high_error <- safe_value(
	ad_error_summary,
	ad_error_summary$group == "central",
	"proportion_abs_error_gt_1"
)

outside_high_error <- safe_value(
	ad_error_summary,
	ad_error_summary$group == "outside",
	"proportion_abs_error_gt_1"
)

ablation_logsaqd_delta <- safe_value(
	ablation_summary,
	ablation_summary$model == "remove_LogSaqd_group",
	"delta_RMSE"
)

ablation_logsoce_delta <- safe_value(
	ablation_summary,
	ablation_summary$model == "remove_LogSoce_group",
	"delta_RMSE"
)

ablation_mwa_delta <- safe_value(
	ablation_summary,
	ablation_summary$model == "remove_MWa_group",
	"delta_RMSE"
)

ablation_mptc_delta <- safe_value(
	ablation_summary,
	ablation_summary$model == "remove_Mptc_group",
	"delta_RMSE"
)

ablation_interaction_delta <- safe_value(
	ablation_summary,
	ablation_summary$model == "remove_interaction_terms",
	"delta_RMSE"
)

n_initial <- safe_value(
	cleaning_summary,
	cleaning_summary$item == "initial_observations",
	"value"
)

n_final <- safe_value(
	cleaning_summary,
	cleaning_summary$item == "final_observations",
	"value"
)

n_compounds <- safe_value(
	cleaning_summary,
	cleaning_summary$item == "unique_compounds",
	"value"
)

r_logsaqd_gse <- safe_value(
	solubility_redundancy,
	solubility_redundancy$analysis == "Correlation between LogSaqd and GSE_logS",
	"value"
)

r2_logsaqd_gse <- safe_value(
	solubility_redundancy,
	solubility_redundancy$analysis == "Regression: LogSaqd ~ GSE_logS",
	"value"
)

############################################################
# Export manuscript key-value table
############################################################

manuscript_key_results <- data.frame(
	item = c(
		"selected_model_formula",
		"selected_model_predictors",
		"final_model_complete_case_observations",
		"final_model_complete_case_compounds",
		"initial_observations",
		"final_cleaned_observations",
		"final_cleaned_unique_compounds",
		"selected_model_LOCO_RMSE",
		"selected_model_LOCO_MAE",
		"selected_model_LOCO_R2_pred",
		"random_forest_selected_predictors_RMSE",
		"rdkit_rf_molecular_only_as_reported_RMSE",
		"rdkit_rf_with_experimental_as_reported_RMSE",
		"rowwise_cv_RMSE_mean",
		"rowwise_cv_RMSE_sd",
		"rowwise_cv_proportion_compound_seen_in_training",
		"central_domain_RMSE",
		"broad_domain_RMSE",
		"outside_domain_RMSE",
		"central_domain_high_error_proportion",
		"outside_domain_high_error_proportion",
		"ablation_remove_LogSaqd_delta_RMSE",
		"ablation_remove_LogSoce_delta_RMSE",
		"ablation_remove_MWa_delta_RMSE",
		"ablation_remove_Mptc_delta_RMSE",
		"ablation_remove_interaction_delta_RMSE",
		"correlation_LogSaqd_GSE_logS",
		"R2_LogSaqd_GSE_logS"
	),
	value = c(
		selected_model_text,
		paste(selected_predictors, collapse = ", "),
		nrow(model_df),
		length(unique(model_df[[id_col]])),
		n_initial,
		n_final,
		n_compounds,
		selected_rmse,
		selected_mae,
		selected_r2,
		rf_selected_rmse,
		rdkit_rf_mol_rmse,
		rdkit_rf_exp_rmse,
		rowwise_rmse,
		rowwise_rmse_sd,
		rowwise_leakage,
		central_rmse,
		broad_rmse,
		outside_rmse,
		central_high_error,
		outside_high_error,
		ablation_logsaqd_delta,
		ablation_logsoce_delta,
		ablation_mwa_delta,
		ablation_mptc_delta,
		ablation_interaction_delta,
		r_logsaqd_gse,
		r2_logsaqd_gse
	),
	stringsAsFactors = FALSE
)

write.csv(
	manuscript_key_results,
	path_final_model_manuscript_values,
	row.names = FALSE
)

############################################################
# Export final performance summary table
############################################################

final_performance_summary <- data.frame(
	analysis = c(
		"Selected model LOCO-CV",
		"Random forest using selected predictors",
		"RDKit RF molecular only, as-reported",
		"RDKit RF with experimental variables, as-reported",
		"Row-wise CV sensitivity",
		"Central applicability domain",
		"Broad applicability domain",
		"Outside applicability domain"
	),
	RMSE = c(
		selected_rmse,
		rf_selected_rmse,
		rdkit_rf_mol_rmse,
		rdkit_rf_exp_rmse,
		rowwise_rmse,
		central_rmse,
		broad_rmse,
		outside_rmse
	),
	MAE = c(
		selected_mae,
		safe_value(
			benchmark_main,
			benchmark_main$model == "random_forest_selected_predictors",
			"MAE"
		),
		safe_value(
			benchmark_main,
			benchmark_main$model == "rdkit_rf_molecular_only_as_reported",
			"MAE"
		),
		safe_value(
			benchmark_main,
			benchmark_main$model == "rdkit_rf_with_experimental_as_reported",
			"MAE"
		),
		safe_value(
			validation_scheme_summary,
			validation_scheme_summary$validation_scheme == "rowwise_cv_repeats",
			"MAE"
		),
		safe_value(
			ad_error_summary,
			ad_error_summary$group == "central",
			"MAE"
		),
		safe_value(
			ad_error_summary,
			ad_error_summary$group == "broad",
			"MAE"
		),
		safe_value(
			ad_error_summary,
			ad_error_summary$group == "outside",
			"MAE"
		)
	),
	R2_pred = c(
		selected_r2,
		safe_value(
			benchmark_main,
			benchmark_main$model == "random_forest_selected_predictors",
			"R2_pred"
		),
		safe_value(
			benchmark_main,
			benchmark_main$model == "rdkit_rf_molecular_only_as_reported",
			"R2_pred"
		),
		safe_value(
			benchmark_main,
			benchmark_main$model == "rdkit_rf_with_experimental_as_reported",
			"R2_pred"
		),
		safe_value(
			validation_scheme_summary,
			validation_scheme_summary$validation_scheme == "rowwise_cv_repeats",
			"R2_pred"
		),
		safe_value(
			ad_error_summary,
			ad_error_summary$group == "central",
			"R2_pred"
		),
		safe_value(
			ad_error_summary,
			ad_error_summary$group == "broad",
			"R2_pred"
		),
		safe_value(
			ad_error_summary,
			ad_error_summary$group == "outside",
			"R2_pred"
		)
	),
	stringsAsFactors = FALSE
)

write.csv(
	final_performance_summary,
	path_table_final_model_performance,
	row.names = FALSE
)

############################################################
# Export selected output inventory
############################################################

selected_outputs <- data.frame(
	output_type = c(
		"Selected model equation",
		"Final model coefficients",
		"Manuscript key values",
		"Final model performance summary",
		"Main benchmark table",
		"All benchmark table",
		"Ablation summary",
		"Validation scheme sensitivity",
		"Applicability-domain error summary",
		"Applicability-domain descriptor summary",
		"Partial-effect dataset",
		"Interaction-effect dataset",
		"RDKit benchmark summary"
	),
	path = c(
		path_final_model_equation,
		path_final_model_coefficients,
		path_final_model_manuscript_values,
		path_table_final_model_performance,
		path_table_benchmark_main,
		path_table_benchmark_all,
		path_table_ablation_summary,
		path_validation_scheme_summary,
		path_table_ad_error_by_domain,
		path_table_ad_descriptor_summary,
		path_partial_effect_dataset,
		path_interaction_effect_dataset,
		path_table_rdkit_benchmark_summary
	),
	exists = c(
		file.exists(path_final_model_equation),
		file.exists(path_final_model_coefficients),
		file.exists(path_final_model_manuscript_values),
		file.exists(path_table_final_model_performance),
		file.exists(path_table_benchmark_main),
		file.exists(path_table_benchmark_all),
		file.exists(path_table_ablation_summary),
		file.exists(path_validation_scheme_summary),
		file.exists(path_table_ad_error_by_domain),
		file.exists(path_table_ad_descriptor_summary),
		file.exists(path_partial_effect_dataset),
		file.exists(path_interaction_effect_dataset),
		file.exists(path_table_rdkit_benchmark_summary)
	),
	stringsAsFactors = FALSE
)

write.csv(
	selected_outputs,
	path_table_final_model_selected_outputs,
	row.names = FALSE
)

############################################################
# Figure: final model coefficients
############################################################

make_final_model_coefficient_plot <- function() {

	plot_df <- coefficient_table[
		coefficient_table$term != "(Intercept)",
	]

	plot_df <- plot_df[
		order(
			plot_df$estimate
		),
	]

	par(
		mar = c(5, 12, 3, 2)
	)

	x_range <- range(
		c(
			plot_df$estimate,
			0
		),
		na.rm = TRUE
	)

	plot(
		plot_df$estimate,
		seq_len(nrow(plot_df)),
		pch = 16,
		xlim = x_range,
		yaxt = "n",
		xlab = "Coefficient estimate",
		ylab = "",
		main = "Final selected model coefficients"
	)

	axis(
		side = 2,
		at = seq_len(nrow(plot_df)),
		labels = plot_df$term,
		las = 1
	)

	abline(
		v = 0,
		lty = 2
	)
}

pdf(
	path_fig_final_model_coefficients_pdf,
	width = 7,
	height = 5
)

make_final_model_coefficient_plot()

dev.off()

png(
	path_fig_final_model_coefficients_png,
	width = 1800,
	height = 1200,
	res = 250
)

make_final_model_coefficient_plot()

dev.off()

############################################################
# Console summary
############################################################

cat("\nFinal model export complete.\n\n")

cat("Selected model:\n")
cat(selected_model_text, "\n\n")

cat("Selected predictors:\n")
cat(paste(selected_predictors, collapse = ", "), "\n\n")

cat("Final model complete-case observations:", nrow(model_df), "\n")
cat("Final model complete-case compounds:", length(unique(model_df[[id_col]])), "\n\n")

cat("Selected model performance:\n")
cat("RMSE:", selected_rmse, "\n")
cat("MAE:", selected_mae, "\n")
cat("R2_pred:", selected_r2, "\n\n")

cat("Key outputs written:\n")
print(selected_outputs)

cat("\nManuscript key results written to:\n")
cat(path_final_model_manuscript_values, "\n")