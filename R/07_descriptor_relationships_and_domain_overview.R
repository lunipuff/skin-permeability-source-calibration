############################################################
# 07_descriptor_relationships_and_domain_overview.R
# Explore selected descriptor relationships with observed logKp
############################################################

source("R/00_config.R")

############################################################
# Load cleaned dataset
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

df$Compound <- trimws(df$Compound)

############################################################
# User settings
############################################################

outcome_col <- "logkpl"

selected_model <- readLines(
	path_loco_cv_selected_model,
	warn = FALSE
)

selected_model <- selected_model[1]

selected_predictors <- get_formula_predictors(
	selected_model,
	outcome_col = outcome_col
)

predictor_labels <- c(
	MWa = "Molecular weight",
	logKowb = "Octanol-water partition coefficient",
	Mptc = "Melting point",
	LogSaqd = "Aqueous solubility",
	LogSoce = "Octanol solubility",
	Texpi = "Experimental temperature",
	Skin.thicknessj = "Skin thickness",
	Hdf = "Hydrogen-bond donor count",
	Hag = "Hydrogen-bond acceptor count",
	MVh = "Molecular volume"
)

predictor_units <- c(
	MWa = "Da",
	logKowb = "logKow",
	Mptc = "K",
	LogSaqd = "log mol/mL",
	LogSoce = "log mol/mL",
	Texpi = "K",
	Skin.thicknessj = "mm",
	Hdf = "count",
	Hag = "count",
	MVh = "cm3/mol"
)

############################################################
# Check required columns
############################################################

required_cols <- c(
	"compound_id",
	"CAS.No",
	"Compound",
	outcome_col,
	selected_predictors
)

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
# Prepare analysis dataset
############################################################

plot_df <- df[
	,
	required_cols,
	drop = FALSE
]

plot_df <- plot_df[
	complete.cases(plot_df),
	,
	drop = FALSE
]

write.csv(
	plot_df,
	path_descriptor_relationship_dataset,
	row.names = FALSE
)

############################################################
# Define central descriptor-domain subset
############################################################

# This subset is not used for model training.
# It is an exploratory central-domain subset for diagnostics.
# Thresholds should be reported as user-defined descriptor ranges.

central_domain <- plot_df[
	plot_df$MWa >= 100 &
		plot_df$MWa <= 250 &
		plot_df$logKowb >= 0 &
		plot_df$logKowb <= 4 &
		plot_df$Mptc >= 200 &
		plot_df$Mptc <= 400 &
		plot_df$LogSaqd >= -7 &
		plot_df$LogSaqd <= -4,
]

write.csv(
	central_domain,
	path_descriptor_central_domain_dataset,
	row.names = FALSE
)

central_domain_summary <- data.frame(
	item = c(
		"n_total_observations",
		"n_total_compounds",
		"n_central_domain_observations",
		"n_central_domain_compounds",
		"proportion_observations_in_central_domain",
		"proportion_compounds_in_central_domain",
		"MWa_range",
		"logKowb_range",
		"Mptc_range",
		"LogSaqd_range"
	),
	value = c(
		nrow(plot_df),
		length(unique(plot_df$compound_id)),
		nrow(central_domain),
		length(unique(central_domain$compound_id)),
		nrow(central_domain) / nrow(plot_df),
		length(unique(central_domain$compound_id)) / length(unique(plot_df$compound_id)),
		"100-250",
		"0-4",
		"200-400",
		"-7 to -4"
	)
)

write.csv(
	central_domain_summary,
	path_descriptor_central_domain_summary,
	row.names = FALSE
)

############################################################
# Console summary
############################################################

cat("Descriptor relationship and domain overview complete.\n")
cat("Total complete observations:", nrow(plot_df), "\n")
cat("Total complete compounds:", length(unique(plot_df$compound_id)), "\n")
cat("Central-domain observations:", nrow(central_domain), "\n")
cat("Central-domain compounds:", length(unique(central_domain$compound_id)), "\n")