############################################################
# 01b_cleaning_flow_table.R
# Manuscript-ready cleaning and modeling flow table
############################################################

source("R/00_config.R")

############################################################
# Create output directories
############################################################

dir.create(
	"tables",
	recursive = TRUE,
	showWarnings = FALSE
)

dir.create(
	"results/model_comparison",
	recursive = TRUE,
	showWarnings = FALSE
)

############################################################
# Load raw dataset
############################################################

rawxl <- readxl::read_excel(
	path_raw_excel
)

raw <- as.data.frame(
	rawxl
)

############################################################
# Settings
############################################################

outcome_col <- "logkpl"
compound_col <- "compound_id"

############################################################
# Helper functions
############################################################

standardize_column_names <- function(data) {

	names(data) <- trimws(
		names(data)
	)

	names(data) <- gsub(
		" +",
		".",
		names(data)
	)

	names(data) <- gsub(
		"\\.+",
		".",
		names(data)
	)

	names(data) <- gsub(
		"\\.$",
		"",
		names(data)
	)

	data
}

standardize_cas <- function(x) {

	x <- as.character(
		x
	)

	x <- gsub(
		"–",
		"-",
		x
	)

	x <- gsub(
		"—",
		"-",
		x
	)

	x <- gsub(
		"--",
		"-",
		x
	)

	x <- trimws(
		x
	)

	x
}

standardize_compound <- function(x) {

	x <- as.character(
		x
	)

	x <- tolower(
		x
	)

	x <- trimws(
		x
	)

	x
}

count_unique_nonmissing <- function(x) {

	x <- as.character(
		x
	)

	x <- trimws(
		x
	)

	x <- x[
		!is.na(x) &
			x != "" &
			x != "N/A"
	]

	length(
		unique(x)
	)
}

count_compounds <- function(data, id_col = "CAS.No") {

	if (!(id_col %in% names(data))) {
		return(NA_integer_)
	}

	count_unique_nonmissing(
		data[[id_col]]
	)
}

add_flow_row <- function(step,
						 data,
						 id_col = "CAS.No",
						 note = "") {

	data.frame(
		Step = step,
		Observations = nrow(data),
		Unique_compounds = count_compounds(
			data,
			id_col = id_col
		),
		Note = note,
		stringsAsFactors = FALSE
	)
}

extract_log_variables <- function(formula_text) {

	matches <- gregexpr(
		"log\\(([^\\)]+)\\)",
		formula_text
	)

	regmatches_result <- regmatches(
		formula_text,
		matches
	)[[1]]

	if (length(regmatches_result) == 0) {
		return(character(0))
	}

	log_vars <- gsub(
		"log\\(([^\\)]+)\\)",
		"\\1",
		regmatches_result
	)

	unique(
		log_vars
	)
}

safe_rbind <- function(x) {

	x <- x[
		!sapply(
			x,
			is.null
		)
	]

	if (length(x) == 0) {
		return(data.frame())
	}

	do.call(
		rbind,
		x
	)
}

check_required_columns <- function(data,
								   required_cols,
								   step_name) {

	missing_cols <- required_cols[
		!(required_cols %in% names(data))
	]

	if (length(missing_cols) > 0) {
		stop(
			paste(
				"Missing columns in",
				step_name,
				":",
				paste(missing_cols, collapse = ", ")
			)
		)
	}
}

############################################################
# Step 1: raw imported data
############################################################

raw_imported <- raw

############################################################
# Step 2: standardize column names, compound names, and CAS
############################################################

df <- raw

df <- standardize_column_names(
	df
)

check_required_columns(
	data = df,
	required_cols = c(
		"Compound",
		"CAS.No"
	),
	step_name = "raw data after column-name standardization"
)

df$Compound <- standardize_compound(
	df$Compound
)

df$CAS.No <- standardize_cas(
	df$CAS.No
)

df <- df[
	order(df$Compound),
]

initial_formatted <- df

############################################################
# Step 3: convert numeric columns and select relevant columns
############################################################

numeric_cols <- c(
	"MWa",
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Hdf",
	"Hag",
	"MVh",
	"Texpi",
	"Skin.thicknessj",
	"logkpl"
)

check_required_columns(
	data = df,
	required_cols = numeric_cols,
	step_name = "numeric conversion"
)

for (col in numeric_cols) {
	df[[col]] <- as.numeric(
		df[[col]]
	)
}

selected_cols <- c(
	"Compound",
	"CAS.No",
	"MWa",
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Hdf",
	"Hag",
	"MVh",
	"Texpi",
	"Skin.thicknessj",
	"logkpl",
	"Reference"
)

check_required_columns(
	data = df,
	required_cols = selected_cols,
	step_name = "column selection"
)

df <- df[
	,
	selected_cols
]

selected_columns <- df

############################################################
# Step 4: endpoint availability audit
############################################################

with_endpoint <- df[
	!is.na(df$logkpl),
	,
	drop = FALSE
]

############################################################
# Step 5: remove rows missing molecular volume
############################################################

missing_mvh <- is.na(
	df$MVh
)

after_missing_mvh_removal <- df[
	!missing_mvh,
	,
	drop = FALSE
]

df <- after_missing_mvh_removal

############################################################
# Step 6: remove rows missing compound identifier
############################################################

missing_cas <- is.na(df$CAS.No) |
	df$CAS.No == "" |
	df$CAS.No == "N/A"

after_missing_cas_removal <- df[
	!missing_cas,
	,
	drop = FALSE
]

df <- after_missing_cas_removal

############################################################
# Step 7: remove known abnormal caffeine descriptor entry
############################################################

abnormal_caffeine <- df$CAS.No == "58-08-2" &
	df$logKowb == -0.63

after_abnormal_caffeine_removal <- df[
	!abnormal_caffeine,
	,
	drop = FALSE
]

df <- after_abnormal_caffeine_removal

############################################################
# Step 8: descriptor inconsistency audit
############################################################

condition_keys <- c(
	"CAS.No",
	"Compound",
	"Reference",
	"Texpi",
	"Skin.thicknessj"
)

descriptor_cols <- c(
	"MWa",
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Hdf",
	"Hag",
	"MVh"
)

merge_keys <- c(
	condition_keys,
	descriptor_cols
)

check_required_columns(
	data = df,
	required_cols = merge_keys,
	step_name = "descriptor inconsistency audit"
)

descriptor_inconsistency_summary <- df %>%
	group_by(
		CAS.No,
		Compound,
		Reference,
		Texpi,
		Skin.thicknessj
	) %>%
	summarise(
		n_descriptor_profiles = n_distinct(
			paste(
				MWa,
				logKowb,
				Mptc,
				LogSaqd,
				LogSoce,
				Hdf,
				Hag,
				MVh,
				sep = "|"
			)
		),
		.groups = "drop"
	) %>%
	filter(
		n_descriptor_profiles > 1
	)

n_descriptor_inconsistent_groups <- nrow(
	descriptor_inconsistency_summary
)

############################################################
# Step 9: collapse repeated identical profiles
############################################################

collapsed_output <- df %>%
	group_by(
		across(
			all_of(merge_keys)
		)
	) %>%
	summarise(
		logkpl = mean(logkpl, na.rm = TRUE),
		n_merged = n(),
		logkpl_sd = ifelse(
			n() > 1,
			sd(logkpl, na.rm = TRUE),
			NA_real_
		),
		.groups = "drop"
	)

collapsed_output$compound_id <- collapsed_output$CAS.No

collapsed_output <- as.data.frame(
	collapsed_output
)

############################################################
# Step 10: selected-model complete-case dataset
############################################################

selected_model_text <- readLines(
	path_loco_cv_selected_model,
	warn = FALSE
)

selected_model_text <- selected_model_text[1]

selected_formula <- as.formula(
	selected_model_text
)

selected_predictors <- get_formula_predictors(
	selected_formula,
	outcome_col = outcome_col
)

selected_model_required_cols <- unique(c(
	"compound_id",
	"CAS.No",
	"Compound",
	outcome_col,
	selected_predictors
))

check_required_columns(
	data = collapsed_output,
	required_cols = selected_model_required_cols,
	step_name = "selected-model complete-case filtering"
)

selected_model_dataset <- collapsed_output[
	complete.cases(
		collapsed_output[
			,
			selected_model_required_cols,
			drop = FALSE
		]
	),
	,
	drop = FALSE
]

log_variables <- extract_log_variables(
	selected_model_text
)

if (length(log_variables) > 0) {

	for (v in log_variables) {

		if (!(v %in% names(selected_model_dataset))) {
			stop(
				paste(
					"Log-transformed variable not found:",
					v
				)
			)
		}

		selected_model_dataset <- selected_model_dataset[
			!is.na(selected_model_dataset[[v]]) &
				selected_model_dataset[[v]] > 0,
			,
			drop = FALSE
		]
	}
}

############################################################
# Compile manuscript-ready flow table
############################################################

cleaning_flow_table <- safe_rbind(
	list(
		add_flow_row(
			step = "Raw imported data",
			data = raw_imported,
			id_col = "CAS No",
			note = "Rows read from the source Excel file."
		),
		add_flow_row(
			step = "After name and CAS standardization",
			data = initial_formatted,
			id_col = "CAS.No",
			note = "Compound names were lower-cased and trimmed; CAS separators were standardized."
		),
		add_flow_row(
			step = "After selecting relevant columns",
			data = selected_columns,
			id_col = "CAS.No",
			note = "Endpoint, compound identifiers, descriptors, experimental conditions, and reference fields were retained."
		),
		add_flow_row(
			step = "With valid endpoint",
			data = with_endpoint,
			id_col = "CAS.No",
			note = "Rows with non-missing logkpl."
		),
		add_flow_row(
			step = "After removing rows missing molecular volume",
			data = after_missing_mvh_removal,
			id_col = "CAS.No",
			note = "Rows missing MVh were removed."
		),
		add_flow_row(
			step = "After removing rows missing compound identifier",
			data = after_missing_cas_removal,
			id_col = "CAS.No",
			note = "Rows missing CAS number were removed."
		),
		add_flow_row(
			step = "After removing abnormal caffeine descriptor entry",
			data = after_abnormal_caffeine_removal,
			id_col = "CAS.No",
			note = "One known abnormal caffeine logKowb entry was removed."
		),
		add_flow_row(
			step = "After collapsing repeated identical profiles",
			data = collapsed_output,
			id_col = "compound_id",
			note = "Rows with identical compound, reference, temperature, skin thickness, and descriptor profiles were averaged."
		),
		add_flow_row(
			step = "Final selected-model complete-case dataset",
			data = selected_model_dataset,
			id_col = "compound_id",
			note = "Complete cases for the selected model formula, including transformation-compatible variables."
		)
	)
)

names(cleaning_flow_table) <- c(
	"Step",
	"Observations",
	"Unique compounds",
	"Note"
)

############################################################
# Add diagnostic columns
############################################################

cleaning_flow_table$Rows_removed_from_previous_step <- c(
	NA,
	head(cleaning_flow_table$Observations, -1) -
		tail(cleaning_flow_table$Observations, -1)
)

cleaning_flow_table$Compounds_removed_from_previous_step <- c(
	NA,
	head(cleaning_flow_table$`Unique compounds`, -1) -
		tail(cleaning_flow_table$`Unique compounds`, -1)
)

cleaning_flow_table <- cleaning_flow_table[
	,
	c(
		"Step",
		"Observations",
		"Unique compounds",
		"Rows_removed_from_previous_step",
		"Compounds_removed_from_previous_step",
		"Note"
	)
]

############################################################
# Save outputs
############################################################

write.csv(
	cleaning_flow_table,
	path_cleaning_flow_table,
	row.names = FALSE
)

write.csv(
	cleaning_flow_table,
	path_cleaning_flow_results,
	row.names = FALSE
)

############################################################
# Console summary
############################################################

cat("\nCleaning flow table complete.\n\n")

cat("Selected model:\n")
cat(selected_model_text, "\n\n")

cat("Selected predictors:\n")
cat(paste(selected_predictors, collapse = ", "), "\n\n")

cat("Log-transformed variables requiring positive values:\n")
if (length(log_variables) == 0) {
	cat("None\n\n")
} else {
	cat(paste(log_variables, collapse = ", "), "\n\n")
}

cat("Descriptor-inconsistent compound-condition groups retained as separate profiles:\n")
cat(n_descriptor_inconsistent_groups, "\n\n")

print(
	cleaning_flow_table
)

cat("\nOutput written to:\n")
cat(path_cleaning_flow_table, "\n")