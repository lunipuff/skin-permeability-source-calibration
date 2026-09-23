############################################################
# 99c_make_supplementary_tables_and_figures.R
# Manuscript-ready supplementary tables and figures
############################################################

source("R/00_config.R")

############################################################
# Directories
############################################################

dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("manuscript/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("manuscript/figures", showWarnings = FALSE, recursive = TRUE)

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

format_vif <- function(x) {
	ifelse(
		is.na(x),
		"\u2014",
		ifelse(
			x >= 1000,
			format(
				round(x, 2),
				big.mark = ",",
				scientific = FALSE,
				trim = TRUE
			),
			sprintf("%.2f", x)
		)
	)
}

############################################################
# Table S1: VIF summary
############################################################

path_tableS1_vif <- "tables/tableS1_vif_summary.csv"
path_manuscript_tableS1_vif <- "manuscript/tables/tableS1_vif_summary.csv"

input_files <- c(
	path_cleaned_dataset,
	path_all_descriptor_vif
)

invisible(lapply(input_files, check_file_exists))

cleaned_data <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

full_vif <- read.csv(
	path_all_descriptor_vif,
	stringsAsFactors = FALSE
)

required_full_vif_cols <- c("descriptor", "VIF")
missing_full_vif_cols <- required_full_vif_cols[
	!(required_full_vif_cols %in% names(full_vif))
]

if (length(missing_full_vif_cols) > 0) {
	stop(
		"Missing required columns in full VIF table: ",
		paste(missing_full_vif_cols, collapse = ", "),
		call. = FALSE
	)
}

############################################################
# Final model base predictors
############################################################
# Update this vector if the selected final model changes.

final_base_predictors <- c(
	"MWa",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Texpi"
)

missing_final_predictors <- final_base_predictors[
	!(final_base_predictors %in% names(cleaned_data))
]

if (length(missing_final_predictors) > 0) {
	stop(
		"Missing final model base predictors in cleaned dataset: ",
		paste(missing_final_predictors, collapse = ", "),
		call. = FALSE
	)
}

for (col in final_base_predictors) {
	cleaned_data[[col]] <- as.numeric(cleaned_data[[col]])
}

final_vif <- calculate_vif(
	data = cleaned_data,
	predictor_cols = final_base_predictors
)

############################################################
# Combine full-pool and final-model VIF values
############################################################

descriptor_order <- c(
	"logKowb",
	"LogSaqd",
	"LogSoce",
	"MWa",
	"MVh",
	"Hag",
	"Hdf",
	"Mptc",
	"Texpi",
	"Skin.thicknessj"
)

tableS1_vif <- merge(
	full_vif,
	final_vif,
	by = "descriptor",
	all.x = TRUE,
	suffixes = c("_full_pool", "_final_model")
)

tableS1_vif <- tableS1_vif[
	match(descriptor_order, tableS1_vif$descriptor),
	,
	drop = FALSE
]

names(tableS1_vif) <- c(
	"Descriptor",
	"Full candidate-pool VIF",
	"Final model base-predictor VIF"
)

tableS1_vif$`Full candidate-pool VIF` <- format_vif(
	tableS1_vif$`Full candidate-pool VIF`
)

tableS1_vif$`Final model base-predictor VIF` <- format_vif(
	tableS1_vif$`Final model base-predictor VIF`
)

write.csv(
	tableS1_vif,
	path_tableS1_vif,
	row.names = FALSE
)

write.csv(
	tableS1_vif,
	path_manuscript_tableS1_vif,
	row.names = FALSE
)

cat("\nTable S1 written to:\n")
cat("  ", path_tableS1_vif, "\n", sep = "")
cat("  ", path_manuscript_tableS1_vif, "\n\n", sep = "")

print(tableS1_vif)