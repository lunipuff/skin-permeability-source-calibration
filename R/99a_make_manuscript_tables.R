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