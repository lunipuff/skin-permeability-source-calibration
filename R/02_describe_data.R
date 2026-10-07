# ============================================================
# 02_describe_data.R
#
# Purpose:
# 	Describe the structure of the prepared skin-permeability
# 	dataset before hierarchical modeling.
#
# This script describes:
# 	- overall dataset composition
# 	- compound repetition
# 	- compound-reference overlap
# 	- reference composition
# 	- experimental-condition variation
# 	- within-compound variation in logKp
#
# This script DOES NOT:
# 	- fit statistical models
# 	- test significance
# 	- estimate reference effects
# 	- estimate standardized permeability
# 	- remove observations
#
# Run from the project root:
# 	source("R/02_describe_data.R")
# ============================================================


source("R/00_config.R")


# ============================================================
# 1. Import prepared data
# ============================================================

dat <- readr::read_csv(
	"data/processed/01_permeability_data.csv",
	show_col_types = FALSE
)

message("Observations: ", nrow(dat))
message("Compounds: ", dplyr::n_distinct(dat$compound_id))
message("References: ", dplyr::n_distinct(dat$reference_id))


# ============================================================
# 2. Basic integrity checks
# ============================================================

required_variables <- c(
	"source_row",
	"compound_id",
	"compound_name",
	"reference_id",
	"logKp",
	"temperature_K",
	"skin_thickness_mm"
)

missing_variables <- setdiff(
	required_variables,
	names(dat)
)

if (length(missing_variables) > 0) {
	stop(
		"Required variables missing: ",
		paste(
			missing_variables,
			collapse = ", "
		)
	)
}


# ============================================================
# 3. Compound-level structure
#
# One row per compound.
# ============================================================

compound_structure <- dat %>%
	group_by(
		compound_id,
		compound_name,
		cas
	) %>%
	summarise(
		n_observations = n(),

		n_references = n_distinct(
			reference_id
		),

		n_temperatures = n_distinct(
			temperature_K
		),

		n_skin_thicknesses = n_distinct(
			skin_thickness_mm
		),

		min_logKp = min(
			logKp,
			na.rm = TRUE
		),

		max_logKp = max(
			logKp,
			na.rm = TRUE
		),

		logKp_range =
			max_logKp -
			min_logKp,

		mean_logKp = mean(
			logKp,
			na.rm = TRUE
		),

		sd_logKp = if_else(
			n() > 1,
			sd(
				logKp,
				na.rm = TRUE
			),
			NA_real_
		),

		.groups = "drop"
	) %>%
	mutate(
		repeated_compound =
			n_observations >= 2,

		multi_reference =
			n_references >= 2,

		multi_temperature =
			n_temperatures >= 2,

		multi_thickness =
			n_skin_thicknesses >= 2
	) %>%
	arrange(
		desc(n_observations),
		compound_name
	)

readr::write_csv(
	compound_structure,
	"R/results/02_compound_structure.csv"
)


# ============================================================
# 4. Reference-level structure
#
# One row per literature reference.
# ============================================================

reference_structure <- dat %>%
	group_by(
		reference_id,
		reference_raw
	) %>%
	summarise(
		n_observations = n(),

		n_compounds = n_distinct(
			compound_id
		),

		n_temperatures = n_distinct(
			temperature_K
		),

		n_skin_thicknesses = n_distinct(
			skin_thickness_mm
		),

		mean_logKp = mean(
			logKp,
			na.rm = TRUE
		),

		sd_logKp = if_else(
			n() > 1,
			sd(
				logKp,
				na.rm = TRUE
			),
			NA_real_
		),

		.groups = "drop"
	) %>%
	arrange(
		desc(n_observations)
	)

readr::write_csv(
	reference_structure,
	"R/results/02_reference_structure.csv"
)


# ============================================================
# 5. Compound repetition summary
#
# Summarizes how much repeated information exists at the
# compound level.
# ============================================================

compound_repetition_summary <- tibble(
	metric = c(
		"Total compounds",
		"Compounds with >=2 observations",
		"Compounds with >=3 observations",
		"Compounds with >=5 observations",
		"Compounds with >=10 observations",
		"Maximum observations for one compound",
		"Median observations per compound",
		"Mean observations per compound"
	),

	value = c(
		nrow(
			compound_structure
		),

		sum(
			compound_structure$n_observations >= 2
		),

		sum(
			compound_structure$n_observations >= 3
		),

		sum(
			compound_structure$n_observations >= 5
		),

		sum(
			compound_structure$n_observations >= 10
		),

		max(
			compound_structure$n_observations
		),

		median(
			compound_structure$n_observations
		),

		mean(
			compound_structure$n_observations
		)
	)
)

readr::write_csv(
	compound_repetition_summary,
	"R/results/02_compound_repetition.csv"
)


# ============================================================
# 6. Compound-reference overlap
#
# This is important for assessing whether compound and
# reference effects can be separated.
# ============================================================

reference_overlap <- compound_structure %>%
	select(
		compound_id,
		compound_name,
		cas,
		n_observations,
		n_references
	) %>%
	filter(
		n_references >= 2
	) %>%
	arrange(
		desc(n_references),
		desc(n_observations),
		compound_name
	)

readr::write_csv(
	reference_overlap,
	"R/results/02_reference_overlap.csv"
)


reference_overlap_summary <- tibble(
	metric = c(
		"Total compounds",
		"Compounds in >=2 references",
		"Compounds in >=3 references",
		"Compounds in >=4 references",
		"Compounds in >=5 references",
		"Maximum references for one compound",
		"Observations from multi-reference compounds"
	),

	value = c(
		nrow(
			compound_structure
		),

		sum(
			compound_structure$n_references >= 2
		),

		sum(
			compound_structure$n_references >= 3
		),

		sum(
			compound_structure$n_references >= 4
		),

		sum(
			compound_structure$n_references >= 5
		),

		max(
			compound_structure$n_references
		),

		dat %>%
			filter(
				compound_id %in%
					reference_overlap$compound_id
			) %>%
			nrow()
	)
)

readr::write_csv(
	reference_overlap_summary,
	"R/results/02_reference_overlap_summary.csv"
)


# ============================================================
# 7. Experimental-condition summary
# ============================================================

experimental_conditions <- tibble(
	variable = c(
		"temperature_K",
		"skin_thickness_mm"
	),

	n = c(
		sum(
			!is.na(
				dat$temperature_K
			)
		),

		sum(
			!is.na(
				dat$skin_thickness_mm
			)
		)
	),

	n_unique = c(
		n_distinct(
			dat$temperature_K,
			na.rm = TRUE
		),

		n_distinct(
			dat$skin_thickness_mm,
			na.rm = TRUE
		)
	),

	min = c(
		min(
			dat$temperature_K,
			na.rm = TRUE
		),

		min(
			dat$skin_thickness_mm,
			na.rm = TRUE
		)
	),

	median = c(
		median(
			dat$temperature_K,
			na.rm = TRUE
		),

		median(
			dat$skin_thickness_mm,
			na.rm = TRUE
		)
	),

	mean = c(
		mean(
			dat$temperature_K,
			na.rm = TRUE
		),

		mean(
			dat$skin_thickness_mm,
			na.rm = TRUE
		)
	),

	sd = c(
		sd(
			dat$temperature_K,
			na.rm = TRUE
		),

		sd(
			dat$skin_thickness_mm,
			na.rm = TRUE
		)
	),

	max = c(
		max(
			dat$temperature_K,
			na.rm = TRUE
		),

		max(
			dat$skin_thickness_mm,
			na.rm = TRUE
		)
	)
)

readr::write_csv(
	experimental_conditions,
	"R/results/02_experimental_conditions.csv"
)


# ============================================================
# 8. Variation in experimental conditions within compounds
# ============================================================

condition_variation <- compound_structure %>%
	summarise(
		total_compounds = n(),

		compounds_multiple_temperatures = sum(
			n_temperatures >= 2
		),

		compounds_multiple_thicknesses = sum(
			n_skin_thicknesses >= 2
		),

		compounds_varying_temperature_or_thickness = sum(
			n_temperatures >= 2 |
				n_skin_thicknesses >= 2
		),

		compounds_varying_temperature_and_thickness = sum(
			n_temperatures >= 2 &
				n_skin_thicknesses >= 2
		)
	) %>%
	pivot_longer(
		everything(),
		names_to = "metric",
		values_to = "value"
	)

readr::write_csv(
	condition_variation,
	"R/results/02_condition_variation.csv"
)


# ============================================================
# 9. Within-compound variability
#
# Restrict summary interpretation to compounds with at least
# two observations.
# ============================================================

within_compound_variability <- compound_structure %>%
	filter(
		n_observations >= 2
	) %>%
	select(
		compound_id,
		compound_name,
		cas,
		n_observations,
		n_references,
		n_temperatures,
		n_skin_thicknesses,
		min_logKp,
		max_logKp,
		logKp_range,
		mean_logKp,
		sd_logKp
	) %>%
	arrange(
		desc(logKp_range)
	)

readr::write_csv(
	within_compound_variability,
	"R/results/02_within_compound_variability.csv"
)


# ============================================================
# 10. Summary of within-compound variability
# ============================================================

within_compound_summary <- tibble(
	metric = c(
		"Compounds with >=2 observations",
		"Median within-compound logKp range",
		"Mean within-compound logKp range",
		"25th percentile logKp range",
		"75th percentile logKp range",
		"Maximum within-compound logKp range",
		"Compounds with logKp range >0.5",
		"Compounds with logKp range >1.0"
	),

	value = c(
		nrow(
			within_compound_variability
		),

		median(
			within_compound_variability$logKp_range,
			na.rm = TRUE
		),

		mean(
			within_compound_variability$logKp_range,
			na.rm = TRUE
		),

		quantile(
			within_compound_variability$logKp_range,
			0.25,
			na.rm = TRUE
		),

		quantile(
			within_compound_variability$logKp_range,
			0.75,
			na.rm = TRUE
		),

		max(
			within_compound_variability$logKp_range,
			na.rm = TRUE
		),

		sum(
			within_compound_variability$logKp_range > 0.5,
			na.rm = TRUE
		),

		sum(
			within_compound_variability$logKp_range > 1,
			na.rm = TRUE
		)
	)
)

readr::write_csv(
	within_compound_summary,
	"R/results/02_within_compound_summary.csv"
)


# ============================================================
# 11. Overall dataset summary
#
# Compact set of values likely to be useful in the manuscript.
# ============================================================

dataset_summary <- tibble(
	metric = c(
		"Observations",
		"Compounds",
		"References",
		"Compound-reference combinations",
		"Compounds with >=2 observations",
		"Compounds in >=2 references",
		"Compounds in >=3 references",
		"Compounds varying in temperature or skin thickness",
		"Median observations per compound",
		"Median within-compound logKp range"
	),

	value = c(
		nrow(
			dat
		),

		n_distinct(
			dat$compound_id
		),

		n_distinct(
			dat$reference_id
		),

		dat %>%
			distinct(
				compound_id,
				reference_id
			) %>%
			nrow(),

		sum(
			compound_structure$n_observations >= 2
		),

		sum(
			compound_structure$n_references >= 2
		),

		sum(
			compound_structure$n_references >= 3
		),

		sum(
			compound_structure$n_temperatures >= 2 |
				compound_structure$n_skin_thicknesses >= 2
		),

		median(
			compound_structure$n_observations
		),

		median(
			within_compound_variability$logKp_range,
			na.rm = TRUE
		)
	)
)

readr::write_csv(
	dataset_summary,
	"R/results/02_dataset_summary.csv"
)


# ============================================================
# 12. Console summary
# ============================================================

cat("\n")
cat("========================================\n")
cat("02_describe_data.R complete\n")
cat("========================================\n\n")

print(dataset_summary)

cat("\nOutputs written to:\n")

cat(
	"  R/results/02_dataset_summary.csv\n"
)

cat(
	"  R/results/02_compound_structure.csv\n"
)

cat(
	"  R/results/02_reference_structure.csv\n"
)

cat(
	"  R/results/02_compound_repetition.csv\n"
)

cat(
	"  R/results/02_reference_overlap.csv\n"
)

cat(
	"  R/results/02_reference_overlap_summary.csv\n"
)

cat(
	"  R/results/02_experimental_conditions.csv\n"
)

cat(
	"  R/results/02_condition_variation.csv\n"
)

cat(
	"  R/results/02_within_compound_variability.csv\n"
)

cat(
	"  R/results/02_within_compound_summary.csv\n"
)