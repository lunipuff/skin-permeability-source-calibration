source("R/00_config.R")

############################################################
# Table 1: Data cleaning and dataset attrition
############################################################

dir.create(
	"tables",
	showWarnings = FALSE,
	recursive = TRUE
)

dir.create(
	"manuscript/tables",
	showWarnings = FALSE,
	recursive = TRUE
)

path_table1_data_cleaning <- "tables/table1_data_cleaning_attrition.csv"
path_manuscript_table1_data_cleaning <- "manuscript/tables/table1_data_cleaning_attrition.csv"

cleaning_flow <- read.csv(
	path_cleaning_flow,
	stringsAsFactors = FALSE
)

required_cleaning_cols <- c(
	"step",
	"observations",
	"unique_compounds"
)

missing_cleaning_cols <- required_cleaning_cols[
	!(required_cleaning_cols %in% names(cleaning_flow))
]

if (length(missing_cleaning_cols) > 0) {
	stop(
		paste(
			"Missing required columns in cleaning_flow:",
			paste(missing_cleaning_cols, collapse = ", ")
		),
		call. = FALSE
	)
}

get_step_row <- function(step_name) {
	out <- cleaning_flow[
		cleaning_flow$step == step_name,
		,
		drop = FALSE
	]

	if (nrow(out) != 1) {
		stop(
			paste0(
				"Expected exactly one row for cleaning step '",
				step_name,
				"', found ",
				nrow(out),
				"."
			),
			call. = FALSE
		)
	}

	out
}

raw_row <- get_step_row("Raw imported dataset")
required_field_row <- get_step_row("After required-field filtering")
abnormal_row <- get_step_row("After abnormal descriptor filtering")
collapsed_row <- get_step_row("After collapsing repeated profiles")
final_row <- get_step_row("Final modeling dataset")

table1_data_cleaning <- data.frame(
	Step = c(
		"Raw imported dataset",
		"After required-field filtering",
		"After abnormal descriptor filtering",
		"After collapsing repeated profiles",
		"Final modeling dataset"
	),
	Observations = c(
		raw_row$observations,
		required_field_row$observations,
		abnormal_row$observations,
		collapsed_row$observations,
		final_row$observations
	),
	"Unique compounds" = c(
		raw_row$unique_compounds,
		required_field_row$unique_compounds,
		abnormal_row$unique_compounds,
		collapsed_row$unique_compounds,
		final_row$unique_compounds
	),
	check.names = FALSE,
	stringsAsFactors = FALSE
)

write.csv(
	table1_data_cleaning,
	path_table1_data_cleaning,
	row.names = FALSE
)

write.csv(
	table1_data_cleaning,
	path_manuscript_table1_data_cleaning,
	row.names = FALSE
)

cat("\nTable 1 written to:\n")
cat("  ", path_table1_data_cleaning, "\n", sep = "")
cat("  ", path_manuscript_table1_data_cleaning, "\n\n", sep = "")

print(table1_data_cleaning)
	
############################################################
# Table 2: Concise LOCO-CV model-selection summary
############################################################

path_loco_cv_candidate_search <- "results/cross_validation/03_loco_cv_candidate_model_search.csv"
path_selected_model <- "results/cross_validation/03_selected_model.txt"

path_table2_model_selection <- "tables/table2_loco_cv_model_selection_summary.csv"
path_manuscript_table2_model_selection <- "manuscript/tables/table2_loco_cv_model_selection_summary.csv"

candidate_search <- read.csv(
	path_loco_cv_candidate_search,
	stringsAsFactors = FALSE
)

selected_model_formula <- readLines(
	path_selected_model,
	warn = FALSE
)

selected_model_formula <- trimws(
	gsub(
		"\\s+",
		" ",
		paste(selected_model_formula, collapse = " ")
	)
)

required_cols <- c(
	"formula",
	"RMSE",
	"MAE",
	"R2_pred",
	"R",
	"n_base_predictors",
	"n_terms"
)

missing_cols <- required_cols[
	!(required_cols %in% names(candidate_search))
]

if (length(missing_cols) > 0) {
	stop(
		"Missing required columns in candidate-model search table: ",
		paste(missing_cols, collapse = ", "),
		call. = FALSE
	)
}

best_rmse_row <- candidate_search[
	order(
		candidate_search$RMSE,
		candidate_search$n_terms
	),
	,
	drop = FALSE
][1, , drop = FALSE]

selected_model_row <- candidate_search[
	candidate_search$formula == selected_model_formula,
	,
	drop = FALSE
]

if (nrow(selected_model_row) != 1) {
	stop(
		"Expected exactly one selected model row in candidate-model search table.",
		call. = FALSE
	)
}

table2_model_selection <- rbind(
	best_rmse_row,
	selected_model_row
)

table2_model_selection$Model <- c(
	"Best-RMSE candidate",
	"Final selected model"
)

table2_model_selection <- table2_model_selection[
	,
	c(
		"Model",
		"RMSE",
		"MAE",
		"R2_pred",
		"R",
		"n_base_predictors",
		"n_terms"
	),
	drop = FALSE
]

names(table2_model_selection) <- c(
	"Model",
	"RMSE",
	"MAE",
	"Predictive R2",
	"Pearson r",
	"Base predictors",
	"Total terms"
)

numeric_cols <- c(
	"RMSE",
	"MAE",
	"Predictive R2",
	"Pearson r"
)

for (col in numeric_cols) {
	table2_model_selection[[col]] <- sprintf(
		"%.3f",
		as.numeric(table2_model_selection[[col]])
	)
}

write.csv(
	table2_model_selection,
	path_table2_model_selection,
	row.names = FALSE
)

write.csv(
	table2_model_selection,
	path_manuscript_table2_model_selection,
	row.names = FALSE
)

print(table2_model_selection)

############################################################
# Table S1: Descriptor redundancy summary
############################################################

dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("manuscript/tables", showWarnings = FALSE, recursive = TRUE)

path_tableS1_descriptor_redundancy <- "tables/tableS1_descriptor_redundancy_summary.csv"
path_manuscript_tableS1_descriptor_redundancy <- "manuscript/tables/tableS1_descriptor_redundancy_summary.csv"

############################################################
# Input paths
############################################################

path_all_descriptor_vif <- "results/descriptor_redundancy/all_descriptor_vif.csv"
path_descriptor_red_flag_sets <- "results/descriptor_redundancy/descriptor_red_flag_sets.csv"
path_descriptor_soft_warning_sets <- "results/descriptor_redundancy/descriptor_soft_warning_sets.csv"
path_solubility_gse_correlation_matrix <- "results/descriptor_redundancy/solubility_gse_correlation_matrix.csv"
path_solubility_gse_regression_summary <- "results/descriptor_redundancy/solubility_gse_regression_summary.csv"

############################################################
# Helper functions
############################################################

check_file_exists <- function(path) {
	if (!file.exists(path)) {
		stop(
			"Missing required file: ",
			path,
			call. = FALSE
		)
	}
}

calculate_vif <- function(data, predictor_cols) {
	complete_data <- data[
		stats::complete.cases(data[, predictor_cols, drop = FALSE]),
		predictor_cols,
		drop = FALSE
	]

	if (nrow(complete_data) < length(predictor_cols) + 2) {
		stop(
			"Too few complete observations to calculate VIF.",
			call. = FALSE
		)
	}

	vif_values <- vapply(
		predictor_cols,
		function(response_col) {
			other_cols <- setdiff(predictor_cols, response_col)

			if (length(other_cols) == 0) {
				return(NA_real_)
			}

			formula_text <- paste(
				response_col,
				"~",
				paste(other_cols, collapse = " + ")
			)

			fit <- stats::lm(
				stats::as.formula(formula_text),
				data = complete_data
			)

			r_squared <- summary(fit)$r.squared

			if (is.na(r_squared)) {
				return(NA_real_)
			}

			if (r_squared >= 1) {
				return(Inf)
			}

			1 / (1 - r_squared)
		},
		numeric(1)
	)

	data.frame(
		descriptor = names(vif_values),
		VIF = as.numeric(vif_values),
		stringsAsFactors = FALSE
	)
}

classify_vif <- function(vif) {
	ifelse(
		vif >= 100,
		"Extreme redundancy",
		ifelse(
			vif >= 10,
			"High redundancy",
			ifelse(
				vif >= 5,
				"Moderate redundancy",
				"Low to acceptable redundancy"
			)
		)
	)
}

############################################################
# Check inputs
############################################################

input_files <- c(
	path_cleaned_dataset,
	path_all_descriptor_vif,
	path_descriptor_red_flag_sets,
	path_descriptor_soft_warning_sets,
	path_solubility_gse_correlation_matrix,
	path_solubility_gse_regression_summary
)

invisible(
	lapply(
		input_files,
		check_file_exists
	)
)

############################################################
# Read inputs
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

all_descriptor_vif <- read.csv(
	path_all_descriptor_vif,
	stringsAsFactors = FALSE
)

descriptor_red_flag_sets <- read.csv(
	path_descriptor_red_flag_sets,
	stringsAsFactors = FALSE
)

descriptor_soft_warning_sets <- read.csv(
	path_descriptor_soft_warning_sets,
	stringsAsFactors = FALSE
)

solubility_gse_correlation <- read.csv(
	path_solubility_gse_correlation_matrix,
	row.names = 1,
	check.names = FALSE
)

solubility_gse_regression_summary <- read.csv(
	path_solubility_gse_regression_summary,
	stringsAsFactors = FALSE
)

############################################################
# Calculate final-model base-predictor VIF
############################################################
# Update this vector if the final selected model changes after
# red-flag model-search filtering is implemented.

final_base_predictors <- c(
	"MWa",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Texpi"
)

missing_final_predictors <- final_base_predictors[
	!(final_base_predictors %in% names(df))
]

if (length(missing_final_predictors) > 0) {
	stop(
		"Missing final model base predictors in cleaned dataset: ",
		paste(missing_final_predictors, collapse = ", "),
		call. = FALSE
	)
}

for (col in final_base_predictors) {
	df[[col]] <- as.numeric(df[[col]])
}

final_base_predictor_vif <- calculate_vif(
	data = df,
	predictor_cols = final_base_predictors
)

############################################################
# Section 1: Full candidate descriptor VIF
############################################################

table_full_vif <- data.frame(
	Section = "Full candidate descriptor VIF",
	Variable_or_test = all_descriptor_vif$descriptor,
	Value = round(all_descriptor_vif$VIF, 3),
	Interpretation = classify_vif(all_descriptor_vif$VIF),
	stringsAsFactors = FALSE
)

############################################################
# Section 2: Final selected model base-predictor VIF
############################################################

table_final_vif <- data.frame(
	Section = "Final selected model base-predictor VIF",
	Variable_or_test = final_base_predictor_vif$descriptor,
	Value = round(final_base_predictor_vif$VIF, 3),
	Interpretation = classify_vif(final_base_predictor_vif$VIF),
	stringsAsFactors = FALSE
)

############################################################
# Section 3: Red-flag and soft-warning sets
############################################################

table_red_flags <- data.frame(
	Section = "Generated red-flag predictor sets",
	Variable_or_test = descriptor_red_flag_sets$predictor_set,
	Value = descriptor_red_flag_sets$severity,
	Interpretation = descriptor_red_flag_sets$reason,
	stringsAsFactors = FALSE
)

table_soft_warnings <- data.frame(
	Section = "Generated soft-warning predictor sets",
	Variable_or_test = descriptor_soft_warning_sets$predictor_set,
	Value = descriptor_soft_warning_sets$severity,
	Interpretation = descriptor_soft_warning_sets$reason,
	stringsAsFactors = FALSE
)

############################################################
# Section 4: Targeted solubility/GSE checks
############################################################

get_regression_r2 <- function(model_name) {
	out <- solubility_gse_regression_summary$r_squared[
		solubility_gse_regression_summary$model == model_name
	]

	if (length(out) != 1) {
		return(NA_real_)
	}

	out
}

table_solubility_gse <- data.frame(
	Section = "Targeted solubility/GSE redundancy check",
	Variable_or_test = c(
		"Correlation: LogSaqd vs GSE_logS",
		"Correlation: LogSoce vs GSE_logS",
		"Regression R2: LogSaqd ~ logKowb + Mptc",
		"Regression R2: LogSoce ~ logKowb + Mptc"
	),
	Value = round(
		c(
			solubility_gse_correlation["LogSaqd", "GSE_logS"],
			solubility_gse_correlation["LogSoce", "GSE_logS"],
			get_regression_r2("LogSaqd ~ logKowb + Mptc"),
			get_regression_r2("LogSoce ~ logKowb + Mptc")
		),
		3
	),
	Interpretation = c(
		"Agreement between reported aqueous solubility and GSE-derived solubility estimate.",
		"Agreement between reported octanol solubility and GSE-derived aqueous solubility estimate.",
		"Fraction of variation in aqueous solubility explained by lipophilicity and melting point.",
		"Fraction of variation in octanol solubility explained by lipophilicity and melting point."
	),
	stringsAsFactors = FALSE
)

############################################################
# Combine and write table
############################################################

tableS1_descriptor_redundancy <- rbind(
	table_full_vif,
	table_final_vif,
	table_red_flags,
	table_soft_warnings,
	table_solubility_gse
)

write.csv(
	tableS1_descriptor_redundancy,
	path_tableS1_descriptor_redundancy,
	row.names = FALSE
)

write.csv(
	tableS1_descriptor_redundancy,
	path_manuscript_tableS1_descriptor_redundancy,
	row.names = FALSE
)

cat("\nTable S1 written to:\n")
cat("  ", path_tableS1_descriptor_redundancy, "\n", sep = "")
cat("  ", path_manuscript_tableS1_descriptor_redundancy, "\n\n", sep = "")

print(tableS1_descriptor_redundancy)

############################################################
# Table 4: Benchmark model performance
############################################################

path_main_benchmark_models <- file.path(
	"tables",
	"table_main_benchmark_models.csv"
)

if (!file.exists(path_main_benchmark_models)) {
	path_main_benchmark_models <- file.path(
		"results",
		"benchmarks",
		"table_main_benchmark_models.csv"
	)
}

if (!file.exists(path_main_benchmark_models)) {
	stop(
		"Could not find table_main_benchmark_models.csv.",
		call. = FALSE
	)
}

benchmark_models <- read.csv(
	path_main_benchmark_models,
	stringsAsFactors = FALSE,
	check.names = FALSE
)

required_cols <- c(
	"model",
	"model_label",
	"n_observations",
	"RMSE",
	"MAE",
	"R2_pred",
	"proportion_abs_error_gt_1"
)

missing_cols <- required_cols[
	!(required_cols %in% names(benchmark_models))
]

if (length(missing_cols) > 0) {
	stop(
		"Missing required columns in benchmark table: ",
		paste(missing_cols, collapse = ", "),
		call. = FALSE
	)
}

model_labels <- c(
	"null_mean" = "Null mean",
	"potts_guy" = "Potts--Guy-style",
	"linear_selected_predictors" = "Linear selected-predictor model",
	"selected_model" = "Selected interpretable QSPR model",
	"random_forest_selected_predictors" = "Random forest, selected predictors",
	"rdkit_rf_molecular_only_as_reported" = "RDKit random forest, molecular only",
	"rdkit_rf_with_experimental_as_reported" = "RDKit random forest, molecular + experimental"
)

input_sets <- c(
	"null_mean" = "Endpoint mean only",
	"potts_guy" = "MWa + logKowb",
	"linear_selected_predictors" = "Selected dataset descriptors",
	"selected_model" = "Selected dataset descriptors",
	"random_forest_selected_predictors" = "Selected dataset descriptors",
	"rdkit_rf_molecular_only_as_reported" = "RDKit molecular descriptors",
	"rdkit_rf_with_experimental_as_reported" = "RDKit molecular descriptors + experimental variables"
)

keep_models <- names(model_labels)

table4_benchmark <- benchmark_models[
	benchmark_models$model %in% keep_models,
	,
	drop = FALSE
]

table4_benchmark$model <- factor(
	table4_benchmark$model,
	levels = keep_models
)

table4_benchmark <- table4_benchmark[
	order(table4_benchmark$model),
	,
	drop = FALSE
]

table4_benchmark <- data.frame(
	"Model" = unname(model_labels[as.character(table4_benchmark$model)]),
	"n" = table4_benchmark$n_observations,
	"RMSE" = sprintf("%.3f", table4_benchmark$RMSE),
	"MAE" = sprintf("%.3f", table4_benchmark$MAE),
	"Predictive R2" = sprintf("%.3f", table4_benchmark$R2_pred),
	"|Error| > 1 log unit (%)" = sprintf(
		"%.1f",
		100 * table4_benchmark$proportion_abs_error_gt_1
	),
	check.names = FALSE,
	stringsAsFactors = FALSE
)

dir.create(
	"tables",
	showWarnings = FALSE,
	recursive = TRUE
)

dir.create(
	file.path("manuscript", "tables"),
	showWarnings = FALSE,
	recursive = TRUE
)

write.csv(
	table4_benchmark,
	file.path(
		"tables",
		"table4_benchmark_model_performance.csv"
	),
	row.names = FALSE
)

write.csv(
	table4_benchmark,
	file.path(
		"manuscript",
		"tables",
		"table4_benchmark_model_performance.csv"
	),
	row.names = FALSE
)

cat("\nTable 4 written to:\n")
cat("  tables/table4_benchmark_model_performance.csv\n")
cat("  manuscript/tables/table4_benchmark_model_performance.csv\n\n")

print(table4_benchmark)

############################################################
# Table S2: Ablation analysis summary
############################################################

path_ablation_summary <- file.path(
	"results",
	"ablation",
	"09_ablation_summary.csv"
)

if (!file.exists(path_ablation_summary)) {
	path_ablation_summary <- file.path(
	"tables",
	"tableS_ablation_summary.csv"
	)
}

if (!file.exists(path_ablation_summary)) {
	stop(
		"Could not find ablation summary table.",
		call. = FALSE
	)
}

ablation_summary <- read.csv(
	path_ablation_summary,
	stringsAsFactors = FALSE,
	check.names = FALSE
)

required_cols <- c(
	"ablation_description",
	"RMSE",
	"MAE",
	"R2_pred",
	"delta_RMSE",
	"delta_MAE",
	"delta_R2_pred"
)

missing_cols <- required_cols[
	!(required_cols %in% names(ablation_summary))
]

if (length(missing_cols) > 0) {
	stop(
		"Missing required columns in ablation summary table: ",
		paste(missing_cols, collapse = ", "),
		call. = FALSE
	)
}

tableS2_ablation <- data.frame(
	"Model variant" = ablation_summary$ablation_description,
	"RMSE" = sprintf(
		"%.3f",
		ablation_summary$RMSE
	),
	"Delta RMSE" = sprintf(
		"%+.3f",
		ablation_summary$delta_RMSE
	),
	"MAE" = sprintf(
		"%.3f",
		ablation_summary$MAE
	),
	"Delta MAE" = sprintf(
		"%+.3f",
		ablation_summary$delta_MAE
	),
	"Predictive R2" = sprintf(
		"%.3f",
		ablation_summary$R2_pred
	),
	"Delta R2" = sprintf(
		"%+.3f",
		ablation_summary$delta_R2_pred
	),
	check.names = FALSE,
	stringsAsFactors = FALSE
)

dir.create(
	"tables",
	showWarnings = FALSE,
	recursive = TRUE
)

dir.create(
	file.path("manuscript", "tables"),
	showWarnings = FALSE,
	recursive = TRUE
)

write.csv(
	tableS2_ablation,
	file.path(
		"tables",
		"tableS2_ablation_analysis_summary.csv"
	),
	row.names = FALSE
)

write.csv(
	tableS2_ablation,
	file.path(
		"manuscript",
		"tables",
		"tableS2_ablation_analysis_summary.csv"
	),
	row.names = FALSE
)

cat("\nTable S2 written to:\n")
cat("  tables/tableS2_ablation_analysis_summary.csv\n")
cat("  manuscript/tables/tableS2_ablation_analysis_summary.csv\n\n")

print(tableS2_ablation)

############################################################
# Table S3: Coefficient-stability summary
############################################################

path_coefficient_stability <- file.path(
	"tables",
	"tableS_selected_model_coefficient_stability.csv"
)

if (!file.exists(path_coefficient_stability)) {
	stop(
		"Could not find coefficient-stability table.",
		call. = FALSE
	)
}

coefficient_stability <- read.csv(
	path_coefficient_stability,
	stringsAsFactors = FALSE,
	check.names = FALSE
)

required_cols <- c(
	"term",
	"mean",
	"sd",
	"rsd_percent",
	"n_folds",
	"n_nonmissing"
)

missing_cols <- required_cols[
	!(required_cols %in% names(coefficient_stability))
]

if (length(missing_cols) > 0) {
	stop(
		"Missing required columns in coefficient-stability table: ",
		paste(missing_cols, collapse = ", "),
		call. = FALSE
	)
}

format_term <- function(x) {
	x <- gsub("I\\(Mptc\\^2\\)", "Mptc^2", x)
	x <- gsub("I\\(LogSaqd\\^2\\)", "LogSaqd^2", x)
	x <- gsub("MWa:LogSaqd", "MWa x LogSaqd", x)
	x
}

tableS3_coefficient_stability <- data.frame(
	"Term" = format_term(coefficient_stability$term),
	"Mean coefficient" = sprintf(
		"%.4g",
		coefficient_stability$mean
	),
	"SD" = sprintf(
		"%.4g",
		coefficient_stability$sd
	),
	"RSD (%)" = sprintf(
		"%.1f",
		coefficient_stability$rsd_percent
	),
	"Folds" = coefficient_stability$n_folds,
	"Nonmissing estimates" = coefficient_stability$n_nonmissing,
	check.names = FALSE,
	stringsAsFactors = FALSE
)

dir.create(
	"tables",
	showWarnings = FALSE,
	recursive = TRUE
)

dir.create(
	file.path("manuscript", "tables"),
	showWarnings = FALSE,
	recursive = TRUE
)

write.csv(
	tableS3_coefficient_stability,
	file.path(
		"tables",
		"tableS3_coefficient_stability_summary.csv"
	),
	row.names = FALSE
)

write.csv(
	tableS3_coefficient_stability,
	file.path(
		"manuscript",
		"tables",
		"tableS3_coefficient_stability_summary.csv"
	),
	row.names = FALSE
)

cat("\nTable S3 written to:\n")
cat("  tables/tableS3_coefficient_stability_summary.csv\n")
cat("  manuscript/tables/tableS3_coefficient_stability_summary.csv\n\n")

print(tableS3_coefficient_stability)

############################################################
# Table S4: Final selected-model coefficients
############################################################

if (!exists("path_cleaned_dataset")) {
	source("R/00_config.R")
}

path_table_final_model_coefficients <- file.path(
	"tables",
	"tableS4_final_model_coefficients.csv"
)

path_manuscript_table_final_model_coefficients <- file.path(
	"manuscript",
	"tables",
	"tableS4_final_model_coefficients.csv"
)

modeling_data <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE,
	check.names = FALSE
)

selected_formula <- logkpl ~ MWa +
	log(Mptc) +
	LogSaqd +
	LogSoce +
	log(Texpi) +
	I(Mptc^2) +
	I(LogSaqd^2) +
	MWa:LogSaqd

model_variables <- all.vars(
	selected_formula
)

final_model_data <- modeling_data[
	complete.cases(
		modeling_data[, model_variables]
	),
	,
	drop = FALSE
]

final_selected_model <- lm(
	selected_formula,
	data = final_model_data
)

coefficient_summary <- summary(
	final_selected_model
)$coefficients

format_term_label <- function(x) {
	x <- gsub(
		"\\(Intercept\\)",
		"Intercept",
		x
	)

	x <- gsub(
		"log\\(Mptc\\)",
		"log(Mptc)",
		x
	)

	x <- gsub(
		"log\\(Texpi\\)",
		"log(Texpi)",
		x
	)

	x <- gsub(
		"I\\(Mptc\\^2\\)",
		"Mptc^2",
		x
	)

	x <- gsub(
		"I\\(LogSaqd\\^2\\)",
		"LogSaqd^2",
		x
	)

	x <- gsub(
		"MWa:LogSaqd",
		"MWa x LogSaqd",
		x
	)

	x
}

tableS4_final_model_coefficients <- data.frame(
	"Term" = format_term_label(
		rownames(
			coefficient_summary
		)
	),
	"Estimate" = sprintf(
		"%.4g",
		coefficient_summary[, "Estimate"]
	),
	"Standard error" = sprintf(
		"%.4g",
		coefficient_summary[, "Std. Error"]
	),
	"t statistic" = sprintf(
		"%.3f",
		coefficient_summary[, "t value"]
	),
	"p value" = ifelse(
		coefficient_summary[, "Pr(>|t|)"] < 0.001,
		"<0.001",
		sprintf(
			"%.3f",
			coefficient_summary[, "Pr(>|t|)"]
		)
	),
	check.names = FALSE,
	stringsAsFactors = FALSE
) %>% select("Term", "Estimate", "Standard error")

write.csv(
	tableS4_final_model_coefficients,
	path_table_final_model_coefficients,
	row.names = FALSE
)

write.csv(
	tableS4_final_model_coefficients,
	path_manuscript_table_final_model_coefficients,
	row.names = FALSE
)

message(
	"Generated final selected-model coefficient table: ",
	path_table_final_model_coefficients
)

############################################################
# Table 6: Applicability-domain error summary
############################################################

if (!exists("p")) {
	source("R/00_config.R")
}

path_applicability_domain_error_summary <- file.path(
	"results",
	"applicability_domain",
	"12_error_by_domain_class.csv"
)

if (!file.exists(path_applicability_domain_error_summary)) {
	path_applicability_domain_error_summary <- file.path(
		"tables",
		"table_applicability_domain_error_summary.csv"
	)
}

if (!file.exists(path_applicability_domain_error_summary)) {
	stop(
		"Could not find applicability-domain error summary.",
		call. = FALSE
	)
}

applicability_domain_error_summary <- read.csv(
	path_applicability_domain_error_summary,
	stringsAsFactors = FALSE,
	check.names = FALSE
)

required_cols <- c(
	"group",
	"RMSE",
	"MAE",
	"median_abs_error",
	"proportion_abs_error_gt_1"
)

missing_cols <- required_cols[
	!(required_cols %in% names(applicability_domain_error_summary))
]

if (length(missing_cols) > 0) {
	stop(
		"Missing required columns in applicability-domain error summary: ",
		paste(missing_cols, collapse = ", "),
		call. = FALSE
	)
}

applicability_domain_error_summary <- applicability_domain_error_summary[
	applicability_domain_error_summary$group %in% c(
		"overall",
		"central",
		"outside"
	),
	,
	drop = FALSE
]

domain_order <- c(
	"overall",
	"central",
	"outside"
)

applicability_domain_error_summary <- applicability_domain_error_summary[
	match(
		domain_order,
		applicability_domain_error_summary$group
	),
	,
	drop = FALSE
]

domain_labels <- c(
	"overall" = "Full dataset",
	"central" = "Central descriptor domain",
	"outside" = "Outside 2.5th-97.5th percentile range"
)

table_applicability_domain_error_summary <- data.frame(
	"Domain category" = domain_labels[
		applicability_domain_error_summary$group
	],
	"RMSE" = sprintf(
		"%.3f",
		applicability_domain_error_summary$RMSE
	),
	"MAE" = sprintf(
		"%.3f",
		applicability_domain_error_summary$MAE
	),
	"Median absolute error" = sprintf(
		"%.3f",
		applicability_domain_error_summary$median_abs_error
	),
	">1 log unit (%)" = sprintf(
		"%.1f",
		100 * applicability_domain_error_summary$proportion_abs_error_gt_1
	),
	check.names = FALSE,
	stringsAsFactors = FALSE
)

write.csv(
	table_applicability_domain_error_summary,
	file.path(
		"tables",
		"table6_applicability_domain_error_summary.csv"
	),
	row.names = FALSE
)

write.csv(
	table_applicability_domain_error_summary,
	file.path(
		"manuscript",
		"tables",
		"table6_applicability_domain_error_summary.csv"
	),
	row.names = FALSE
)

message(
	"Generated applicability-domain error summary table."
)