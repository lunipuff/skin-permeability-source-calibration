source("R/00_config.R")

############################################################
# Table 1: Data cleaning and dataset attrition
############################################################

	path_table_data_cleaning <- "manuscript/tables/table2_data_cleaning_attrition.csv"

	############################################################
	# Read cleaning output
	############################################################

	cleaning_flow <- read.csv(
		path_cleaning_flow,
		stringsAsFactors = FALSE
	)

	############################################################
	# Basic checks
	############################################################

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

	############################################################
	# Helper function
	############################################################

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

	############################################################
	# Pull cleaning-step rows
	############################################################

	raw_row <- get_step_row("Raw imported dataset")
	required_field_row <- get_step_row("After required-field filtering")
	abnormal_row <- get_step_row("After abnormal descriptor filtering")
	collapsed_row <- get_step_row("After collapsing repeated profiles")
	final_row <- get_step_row("Final modeling dataset")

	############################################################
	# Build manuscript-ready table
	############################################################

	table_data_cleaning <- data.frame(
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

	############################################################
	# Save
	############################################################

	write.csv(
		table_data_cleaning,
		path_table_data_cleaning,
		row.names = FALSE
	)

	print(table_data_cleaning)
	
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