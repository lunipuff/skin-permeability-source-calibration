############################################################
# 08c_compile_benchmark_tables.R
# Compile main and supplementary benchmark tables
############################################################

source("R/00_config.R")

############################################################
# Load benchmark summaries
############################################################

benchmark_summary <- read.csv(
	path_benchmark_model_summary,
	stringsAsFactors = FALSE
)

rdkit_summary <- read.csv(
	path_rdkit_benchmark_summary,
	stringsAsFactors = FALSE
)

############################################################
# Helper functions
############################################################

add_missing_columns <- function(df, required_cols) {
	missing_cols <- required_cols[
		!(required_cols %in% names(df))
	]

	if (length(missing_cols) > 0) {
		for (col in missing_cols) {
			df[[col]] <- NA
		}
	}

	df
}

format_formula_text <- function(x) {
	x <- gsub(
		"  +",
		" ",
		x
	)

	x <- trimws(x)

	x
}

############################################################
# Standardize benchmark summary columns
############################################################

standard_cols <- c(
	"model",
	"model_class",
	"representation",
	"descriptor_set",
	"formula",
	"n_observations",
	"n_compounds",
	"RMSE",
	"MAE",
	"R2_pred",
	"R",
	"mean_error",
	"median_abs_error",
	"n_abs_error_gt_1",
	"proportion_abs_error_gt_1",
	"ntree",
	"mtry",
	"nodesize"
)

benchmark_summary <- add_missing_columns(
	benchmark_summary,
	standard_cols
)

rdkit_summary <- add_missing_columns(
	rdkit_summary,
	standard_cols
)

benchmark_summary$model_class <- "dataset_descriptor_model"

rdkit_summary$model_class <- "rdkit_descriptor_model"

benchmark_summary$representation[
	is.na(benchmark_summary$representation)
] <- "not_applicable"

benchmark_summary$descriptor_set[
	is.na(benchmark_summary$descriptor_set)
] <- "dataset_descriptors"

rdkit_summary$representation[
	is.na(rdkit_summary$representation)
] <- "not_reported"

rdkit_summary$descriptor_set[
	is.na(rdkit_summary$descriptor_set)
] <- "rdkit_descriptors"

benchmark_summary$formula <- format_formula_text(
	benchmark_summary$formula
)

rdkit_summary$formula <- format_formula_text(
	rdkit_summary$formula
)

benchmark_summary <- benchmark_summary[
	,
	standard_cols
]

rdkit_summary <- rdkit_summary[
	,
	standard_cols
]

############################################################
# Combine all benchmark models
############################################################

benchmark_all <- rbind(
	benchmark_summary,
	rdkit_summary
)

benchmark_all <- benchmark_all[
	order(
		benchmark_all$RMSE,
		benchmark_all$MAE
	),
]

benchmark_all$model_rank <- seq_len(
	nrow(benchmark_all)
)

############################################################
# Add readable model labels
############################################################

benchmark_all$model_label <- benchmark_all$model

benchmark_all$model_label[
	benchmark_all$model == "null_mean"
] <- "Null mean model"

benchmark_all$model_label[
	benchmark_all$model == "potts_guy"
] <- "Potts-Guy-style model"

benchmark_all$model_label[
	benchmark_all$model == "linear_selected_predictors"
] <- "Linear model using selected predictors"

benchmark_all$model_label[
	benchmark_all$model == "selected_model"
] <- "Selected interpretable QSPR model"

benchmark_all$model_label[
	benchmark_all$model == "extended_descriptor"
] <- "Extended dataset-descriptor model"

benchmark_all$model_label[
	benchmark_all$model == "random_forest_selected_predictors"
] <- "Random forest using selected predictors"

benchmark_all$model_label[
	benchmark_all$model == "rdkit_rf_molecular_only_as_reported"
] <- "RDKit random forest, molecular descriptors only, as-reported structures"

benchmark_all$model_label[
	benchmark_all$model == "rdkit_rf_with_experimental_as_reported"
] <- "RDKit random forest, molecular descriptors plus experimental variables, as-reported structures"

benchmark_all$model_label[
	benchmark_all$model == "rdkit_rf_molecular_only_parent"
] <- "RDKit random forest, molecular descriptors only, parent structures"

benchmark_all$model_label[
	benchmark_all$model == "rdkit_rf_with_experimental_parent"
] <- "RDKit random forest, molecular descriptors plus experimental variables, parent structures"

benchmark_all$model_label[
	benchmark_all$model == "rdkit_lm_molecular_only_as_reported"
] <- "RDKit linear model, molecular descriptors only, as-reported structures"

benchmark_all$model_label[
	benchmark_all$model == "rdkit_lm_with_experimental_as_reported"
] <- "RDKit linear model, molecular descriptors plus experimental variables, as-reported structures"

benchmark_all$model_label[
	benchmark_all$model == "rdkit_lm_molecular_only_parent"
] <- "RDKit linear model, molecular descriptors only, parent structures"

benchmark_all$model_label[
	benchmark_all$model == "rdkit_lm_with_experimental_parent"
] <- "RDKit linear model, molecular descriptors plus experimental variables, parent structures"

############################################################
# Create supplementary benchmark table
############################################################

benchmark_all_table <- benchmark_all[
	,
	c(
		"model_rank",
		"model",
		"model_label",
		"model_class",
		"representation",
		"descriptor_set",
		"n_observations",
		"n_compounds",
		"RMSE",
		"MAE",
		"R2_pred",
		"R",
		"mean_error",
		"median_abs_error",
		"n_abs_error_gt_1",
		"proportion_abs_error_gt_1",
		"ntree",
		"mtry",
		"nodesize",
		"formula"
	)
]

write.csv(
	benchmark_all_table,
	path_table_benchmark_all,
	row.names = FALSE
)

############################################################
# Select main-text benchmark models
############################################################

main_model_names <- c(
	"null_mean",
	"potts_guy",
	"linear_selected_predictors",
	"selected_model",
	"random_forest_selected_predictors",
	"rdkit_rf_molecular_only_as_reported",
	"rdkit_rf_with_experimental_as_reported"
)

benchmark_main <- benchmark_all[
	benchmark_all$model %in% main_model_names,
]

benchmark_main$main_order <- match(
	benchmark_main$model,
	main_model_names
)

benchmark_main <- benchmark_main[
	order(
		benchmark_main$main_order
	),
]

benchmark_main_table <- benchmark_main[
	,
	c(
		"main_order",
		"model",
		"model_label",
		"model_class",
		"representation",
		"descriptor_set",
		"n_observations",
		"n_compounds",
		"RMSE",
		"MAE",
		"R2_pred",
		"median_abs_error",
		"n_abs_error_gt_1",
		"proportion_abs_error_gt_1"
	)
]

names(benchmark_main_table)[
	names(benchmark_main_table) == "main_order"
] <- "table_order"

write.csv(
	benchmark_main_table,
	path_table_benchmark_main,
	row.names = FALSE
)

############################################################
# Console summary
############################################################

cat("Benchmark tables compiled.\n")
cat("All benchmark models:", nrow(benchmark_all_table), "\n")
cat("Main-text benchmark models:", nrow(benchmark_main_table), "\n")
cat("Supplementary table written to:", path_table_benchmark_all, "\n")
cat("Main table written to:", path_table_benchmark_main, "\n\n")

cat("Main-text benchmark table:\n")
print(
	benchmark_main_table[
		,
		c(
			"model_label",
			"n_observations",
			"n_compounds",
			"RMSE",
			"MAE",
			"R2_pred"
		)
	]
)