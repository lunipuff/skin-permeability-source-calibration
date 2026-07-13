############################################################
# 01_clean_data.R
# Clean human skin permeability dataset
############################################################

source("R/00_config.R")

############################################################
# Helper functions
############################################################

log_removed <- function(data, reason) {
	if (is.null(data) || nrow(data) == 0) {
		return(NULL)
	}

	data$removal_reason <- reason
	data
}

bind_removed <- function(removed_rows, new_rows) {
	dplyr::bind_rows(
		removed_rows,
		new_rows
	)
}

missing_string <- function(x) {
	x_chr <- trimws(as.character(x))

	is.na(x) |
		x_chr == "" |
		x_chr == "N/A" |
		x_chr == "NA" |
		x_chr == "na"
}

missing_numeric <- function(x) {
	is.na(x)
}

count_unique_nonmissing <- function(x) {
	x <- x[
		!is.na(x) &
			trimws(as.character(x)) != "" &
			trimws(as.character(x)) != "N/A"
	]

	length(unique(x))
}

make_count_row <- function(step, data, previous_n = NA_integer_) {
	n_current <- nrow(data)

	change <- ifelse(
		is.na(previous_n),
		NA_integer_,
		n_current - previous_n
	)

	data.frame(
		step = step,
		observations = n_current,
		unique_compounds = count_unique_nonmissing(data$CAS.No),
		change_from_previous_step = change,
		stringsAsFactors = FALSE
	)
}

############################################################
# Read raw dataset
############################################################

rawxl <- read_excel(path_raw_excel)

write.csv(
	rawxl,
	path_interim_raw_csv,
	row.names = FALSE
)

raw <- read.csv(
	path_interim_raw_csv,
	stringsAsFactors = FALSE
)

############################################################
# Initial formatting
############################################################

df <- raw

if (!("Compound" %in% names(df))) {
	stop("Required column missing from raw dataset: Compound", call. = FALSE)
}

if (!("CAS.No" %in% names(df))) {
	stop("Required column missing from raw dataset: CAS.No", call. = FALSE)
}

df$Compound <- tolower(df$Compound)
df$Compound <- trimws(df$Compound)

df$CAS.No <- gsub("–", "-", df$CAS.No)
df$CAS.No <- gsub("--", "-", df$CAS.No)
df$CAS.No <- trimws(df$CAS.No)

df <- df[order(df$Compound), ]

############################################################
# Initialize logs
############################################################

removed_rows <- NULL

cleaning_flow <- make_count_row(
	step = "Raw imported dataset",
	data = df,
	previous_n = NA_integer_
)

############################################################
# Convert numeric columns
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

missing_numeric_cols <- numeric_cols[
	!(numeric_cols %in% names(df))
]

if (length(missing_numeric_cols) > 0) {
	stop(
		paste(
			"Required numeric columns missing from raw dataset:",
			paste(missing_numeric_cols, collapse = ", ")
		),
		call. = FALSE
	)
}

for (col in numeric_cols) {
	df[[col]] <- as.numeric(df[[col]])
}

############################################################
# Select relevant columns
############################################################

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

missing_selected_cols <- selected_cols[
	!(selected_cols %in% names(df))
]

if (length(missing_selected_cols) > 0) {
	stop(
		paste(
			"Required selected columns missing from raw dataset:",
			paste(missing_selected_cols, collapse = ", ")
		),
		call. = FALSE
	)
}

df <- df[, selected_cols]

cleaning_flow <- dplyr::bind_rows(
	cleaning_flow,
	make_count_row(
		step = "After formatting and column selection",
		data = df,
		previous_n = cleaning_flow$observations[nrow(cleaning_flow)]
	)
)

############################################################
# Define required fields
############################################################

required_identifier_cols <- c(
	"Compound",
	"CAS.No",
	"Reference"
)

required_predictor_cols <- c(
	"MWa",
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Hdf",
	"Hag",
	"MVh",
	"Texpi",
	"Skin.thicknessj"
)

required_endpoint_cols <- c(
	"logkpl"
)

required_cols <- c(
	required_identifier_cols,
	required_predictor_cols,
	required_endpoint_cols
)

############################################################
# Inspect missingness in required fields
############################################################

missing_matrix <- data.frame(
	row_index = seq_len(nrow(df)),
	stringsAsFactors = FALSE
)

for (col in required_cols) {
	if (col %in% numeric_cols) {
		missing_matrix[[col]] <- missing_numeric(df[[col]])
	} else {
		missing_matrix[[col]] <- missing_string(df[[col]])
	}
}

missing_summary <- data.frame(
	variable = required_cols,
	n_missing = vapply(
		required_cols,
		function(col) {
			sum(missing_matrix[[col]], na.rm = TRUE)
		},
		numeric(1)
	),
	stringsAsFactors = FALSE
)

write.csv(
	missing_summary,
	path_required_missingness_summary,
	row.names = FALSE
)

############################################################
# Remove rows with missing required fields
############################################################

missing_required <- apply(
	missing_matrix[, required_cols, drop = FALSE],
	1,
	any
)

if (any(missing_required)) {
	missing_reason <- apply(
		missing_matrix[missing_required, required_cols, drop = FALSE],
		1,
		function(x) {
			missing_cols <- required_cols[as.logical(x)]

			paste(
				missing_cols,
				collapse = "; "
			)
		}
	)

	missing_required_rows <- df[missing_required, , drop = FALSE]
	missing_required_rows$missing_required_fields <- missing_reason

	removed_rows <- bind_removed(
		removed_rows,
		log_removed(
			missing_required_rows,
			"Missing required modeling field"
		)
	)
}

df <- df[!missing_required, , drop = FALSE]

cleaning_flow <- dplyr::bind_rows(
	cleaning_flow,
	make_count_row(
		step = "After required-field filtering",
		data = df,
		previous_n = cleaning_flow$observations[nrow(cleaning_flow)]
	)
)

############################################################
# Explicit CAS identifier-quality check
############################################################

invalid_cas <- missing_string(df$CAS.No)

removed_rows <- bind_removed(
	removed_rows,
	log_removed(
		df[invalid_cas, , drop = FALSE],
		"Missing or unusable CAS number"
	)
)

df <- df[!invalid_cas, , drop = FALSE]

cleaning_flow <- dplyr::bind_rows(
	cleaning_flow,
	make_count_row(
		step = "After CAS identifier filtering",
		data = df,
		previous_n = cleaning_flow$observations[nrow(cleaning_flow)]
	)
)

############################################################
# Remove known abnormal entry after manual inspection
############################################################

abnormal_entry <- df$CAS.No == "58-08-2" &
	df$logKowb == -0.63

removed_rows <- bind_removed(
	removed_rows,
	log_removed(
		df[abnormal_entry, , drop = FALSE],
		"Removed manually confirmed abnormal logKowb entry"
	)
)

df <- df[!abnormal_entry, , drop = FALSE]

cleaning_flow <- dplyr::bind_rows(
	cleaning_flow,
	make_count_row(
		step = "After abnormal descriptor filtering",
		data = df,
		previous_n = cleaning_flow$observations[nrow(cleaning_flow)]
	)
)

############################################################
# Collapse repeated entries from same compound, reference,
# temperature, skin thickness, and descriptor values
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

############################################################
# Check descriptor inconsistencies within compound-condition groups
############################################################

descriptor_inconsistency_log <- df %>%
	dplyr::group_by(dplyr::across(dplyr::all_of(condition_keys))) %>%
	dplyr::summarise(
		n_rows = dplyr::n(),
		n_descriptor_profiles = dplyr::n_distinct(
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
			),
			na.rm = TRUE
		),
		dplyr::across(
			dplyr::all_of(descriptor_cols),
			~ dplyr::n_distinct(.x, na.rm = TRUE),
			.names = "n_unique_{.col}"
		),
		.groups = "drop"
	)

descriptor_inconsistency_log <- descriptor_inconsistency_log[
	descriptor_inconsistency_log$n_descriptor_profiles > 1,
]

write.csv(
	descriptor_inconsistency_log,
	path_descriptor_inconsistency_log,
	row.names = FALSE
)

############################################################
# Merge only rows with identical compound-condition-descriptor profiles
############################################################

output <- df %>%
	dplyr::group_by(dplyr::across(dplyr::all_of(merge_keys))) %>%
	dplyr::summarise(
		logkpl = mean(logkpl, na.rm = TRUE),
		n_merged = dplyr::n(),
		logkpl_sd = ifelse(
			dplyr::n() > 1,
			stats::sd(logkpl, na.rm = TRUE),
			NA_real_
		),
		.groups = "drop"
	)

############################################################
# Create compound ID
############################################################

output$compound_id <- output$CAS.No

############################################################
# Reorder columns
############################################################

output <- output[, c(
	"compound_id",
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
	"Reference",
	"n_merged",
	"logkpl_sd"
)]

cleaning_flow <- dplyr::bind_rows(
	cleaning_flow,
	data.frame(
		step = "After collapsing repeated profiles",
		observations = nrow(output),
		unique_compounds = length(unique(output$compound_id)),
		change_from_previous_step = nrow(output) -
			cleaning_flow$observations[nrow(cleaning_flow)],
		stringsAsFactors = FALSE
	)
)

cleaning_flow <- dplyr::bind_rows(
	cleaning_flow,
	data.frame(
		step = "Final modeling dataset",
		observations = nrow(output),
		unique_compounds = length(unique(output$compound_id)),
		change_from_previous_step = NA_integer_,
		stringsAsFactors = FALSE
	)
)

############################################################
# Save cleaned dataset and logs
############################################################

write.csv(
	output,
	path_cleaned_dataset,
	row.names = FALSE
)

if (is.null(removed_rows)) {
	removed_rows <- data.frame()
}

write.csv(
	removed_rows,
	path_removed_rows,
	row.names = FALSE
)

write.csv(
	cleaning_flow,
	path_cleaning_flow,
	row.names = FALSE
)

############################################################
# Save cleaning summary
############################################################

n_collapsed <- sum(output$n_merged > 1, na.rm = TRUE)
n_records_collapsed <- sum(
	output$n_merged[output$n_merged > 1],
	na.rm = TRUE
)

get_cleaning_n <- function(step_name) {
	cleaning_flow$observations[
		cleaning_flow$step == step_name
	]
}

cleaning_summary <- data.frame(
	item = c(
		"initial_observations",
		"rows_after_formatting_and_column_selection",
		"rows_after_required_field_filtering",
		"rows_after_cas_filtering",
		"rows_after_abnormal_entry_removal",
		"final_observations",
		"unique_compounds",
		"removed_rows",
		"collapsed_observations",
		"original_records_in_collapsed_observations",
		"compound_condition_groups_with_multiple_descriptor_profiles"
	),
	value = c(
		get_cleaning_n("Raw imported dataset"),
		get_cleaning_n("After formatting and column selection"),
		get_cleaning_n("After required-field filtering"),
		get_cleaning_n("After CAS identifier filtering"),
		get_cleaning_n("After abnormal descriptor filtering"),
		nrow(output),
		length(unique(output$compound_id)),
		nrow(removed_rows),
		n_collapsed,
		n_records_collapsed,
		nrow(descriptor_inconsistency_log)
	)
)

write.csv(
	cleaning_summary,
	path_cleaning_summary,
	row.names = FALSE
)

############################################################
# Summary messages
############################################################

cat("Cleaning complete.\n")
cat(
	"Initial rows:",
	cleaning_summary$value[cleaning_summary$item == "initial_observations"],
	"\n"
)
cat(
	"Rows after required-field filtering:",
	cleaning_summary$value[cleaning_summary$item == "rows_after_required_field_filtering"],
	"\n"
)
cat(
	"Rows after CAS filtering:",
	cleaning_summary$value[cleaning_summary$item == "rows_after_cas_filtering"],
	"\n"
)
cat(
	"Rows after abnormal descriptor filtering:",
	cleaning_summary$value[cleaning_summary$item == "rows_after_abnormal_entry_removal"],
	"\n"
)
cat("Final rows:", nrow(output), "\n")
cat("Unique compounds:", length(unique(output$compound_id)), "\n")
cat("Removed rows:", nrow(removed_rows), "\n")
cat("Collapsed final observations:", n_collapsed, "\n")
cat("Original records in collapsed observations:", n_records_collapsed, "\n")
cat("Descriptor inconsistency groups:", nrow(descriptor_inconsistency_log), "\n")