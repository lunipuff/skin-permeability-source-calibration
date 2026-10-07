# ============================================================
# 01_prepare_data.R
#
# Purpose:
# 	Import and prepare the raw Cheruvu skin-permeability data.
#
# This script:
# 	- preserves row-level permeability measurements
# 	- standardizes variable names and data types
# 	- creates canonical compound identities from names
# 	- uses CAS as supporting metadata rather than sole identity
# 	- audits compound-name and CAS inconsistencies
# 	- audits non-numeric values in numeric source columns
# 	- flags repeated experimental profiles
# 	- creates the analysis-ready row-level dataset
#
# This script DOES NOT:
# 	- average repeated permeability measurements
# 	- remove statistical outliers
# 	- fit models
# 	- estimate standardized permeability
# 	- select QSPR predictors
#
# Run from the project root:
# 	source("R/01_prepare_data.R")
# ============================================================


source("R/00_config.R")


# ============================================================
# 1. Import raw data
# ============================================================

raw <- readxl::read_excel(
	"data/raw/Table 1.xlsx",
	sheet = 1
)


# Remove completely blank spreadsheet rows

raw <- raw %>%
	filter(
		if_any(
			everything(),
			~ !is.na(.) &
				stringr::str_squish(as.character(.)) != ""
		)
	)

n_raw <- nrow(raw)

message(
	"Raw observations after blank-row removal: ",
	n_raw
)


# ============================================================
# 2. Audit numeric source columns before coercion
# ============================================================

numeric_source_columns <- c(
	"MWa",
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Hdf",
	"Hag",
	"MVh",
	"Texpi",
	"Skin thicknessj",
	"logkpl",
	"logJmaxm"
)


numeric_coercion_audit <- purrr::map_dfr(
	numeric_source_columns,
	function(var) {

		x <- raw[[var]]

		tibble(
			variable = var,
			raw_value = as.character(x),
			numeric_value =
				suppressWarnings(
					as.numeric(x)
				)
		) %>%
			filter(
				!is.na(raw_value),
				stringr::str_squish(raw_value) != "",
				is.na(numeric_value)
			)
	}
) %>%
	distinct() %>%
	arrange(
		variable,
		raw_value
	)


readr::write_csv(
	numeric_coercion_audit,
	"R/results/01_numeric_coercion_audit.csv"
)


if (nrow(numeric_coercion_audit) > 0) {
	message(
		"Non-numeric values were found in numeric source columns. ",
		"See R/results/01_numeric_coercion_audit.csv"
	)
}


# ============================================================
# 3. Standardize variable names and data types
# ============================================================

dat <- raw %>%
	transmute(

		# Row identity
		source_row = row_number(),

		original_no =
			suppressWarnings(
				as.integer(No)
			),

		# Compound identity
		compound_raw =
			clean_text(
				Compound
			),

		cas =
			clean_cas(
				`CAS No`
			),

		# Dataset grouping variable from source
		set =
			clean_text(
				set
			),

		# Molecular / physicochemical descriptors
		mw =
			suppressWarnings(
				as.numeric(MWa)
			),

		log_kow =
			suppressWarnings(
				as.numeric(logKowb)
			),

		melting_point_K =
			suppressWarnings(
				as.numeric(Mptc)
			),

		log_saq =
			suppressWarnings(
				as.numeric(LogSaqd)
			),

		log_soc =
			suppressWarnings(
				as.numeric(LogSoce)
			),

		h_donor =
			suppressWarnings(
				as.numeric(Hdf)
			),

		h_acceptor =
			suppressWarnings(
				as.numeric(Hag)
			),

		molecular_volume =
			suppressWarnings(
				as.numeric(MVh)
			),

		# Experimental context
		temperature_K =
			suppressWarnings(
				as.numeric(Texpi)
			),

		skin_thickness_mm =
			suppressWarnings(
				as.numeric(
					`Skin thicknessj`
				)
			),

		skin_integrity_test =
			clean_text(
				`Skin Integrity testk`
			),

		skin_selection_criteria =
			clean_text(
				`Skin Selection criteria`
			),

		# Endpoints
		logKp =
			suppressWarnings(
				as.numeric(logkpl)
			),

		logJmax =
			suppressWarnings(
				as.numeric(logJmaxm)
			),

		# Literature source
		reference_raw =
			clean_text(
				Reference
			)
	)


# ============================================================
# 4. Basic source-data integrity checks
# ============================================================

if (anyDuplicated(dat$source_row) > 0) {
	stop(
		"source_row is not unique."
	)
}


if (
	any(!is.na(dat$original_no)) &&
	anyDuplicated(
		dat$original_no[
			!is.na(dat$original_no)
		]
	) > 0
) {
	warning(
		"Duplicate values were found in the original 'No' column."
	)
}


# ============================================================
# 5. Standardize compound identity
#
# Compound identity is based primarily on canonicalized
# compound names.
#
# CAS is retained as supporting metadata and used for audit,
# but is not used as the sole compound identifier because
# the source data contain at least one CAS collision between
# chemically distinct compounds.
# ============================================================

dat <- dat %>%
	mutate(
		compound_name_key =
			normalize_compound_name(
				compound_raw
			)
	)


# ============================================================
# 5a. Explicit canonical-name mapping
#
# Only clear spelling, capitalization, or synonym variants
# are harmonized here.
#
# Everything not listed falls back to its normalized name.
# ============================================================

compound_name_map <- tribble(
	~compound_name_key, ~canonical_name,

	# Estradiol
	"estradiol", "estradiol",
	"beta_estradiol", "estradiol",
	"betaestradiol", "estradiol",

	# Resorcinol
	"resorcinol", "resorcinol",
	"resorcinol_1_3_benzenediol", "resorcinol",

	# Straight-chain alcohols
	"decanol", "decanol",
	"n_decanol", "decanol",

	"hexanol", "hexanol",
	"n_hexanol", "hexanol",

	"octanol", "octanol",
	"n_octanol", "octanol",

	"nonanol", "nonanol",
	"n_nonanol", "nonanol",

	"heptanol", "heptanol",
	"n_heptanol", "heptanol",

	"propanol", "propanol",
	"n_propanol", "propanol",

	"butanol", "butanol",
	"n_butanol", "butanol",

	"pentanol", "pentanol",
	"n_pentanol", "pentanol",

	# Other explicit naming variants
	"5_fluorouracil", "5_fluorouracil",
	"5_fluorouracil", "5_fluorouracil"
)


dat <- dat %>%
	left_join(
		compound_name_map,
		by = "compound_name_key"
	) %>%
	mutate(
		canonical_name =
			coalesce(
				canonical_name,
				compound_name_key
			)
	)


# ============================================================
# 5b. Explicitly separate known chemically distinct entities
#
# These are deliberately kept separate even if the source
# assigns them the same CAS number.
# ============================================================

dat <- dat %>%
	mutate(
		canonical_name = case_when(

			stringr::str_to_lower(
				compound_raw
			) ==
				"pregnenolone" ~
				"pregnenolone",

			stringr::str_to_lower(
				compound_raw
			) ==
				"hydroxypregnenolone" ~
				"hydroxypregnenolone",

			compound_raw ==
				"Tetraethylammonium" ~
				"tetraethylammonium",

			compound_raw ==
				"Tetraethylammonium bromide (+)" ~
				"tetraethylammonium_bromide",

			TRUE ~
				canonical_name
		)
	)


# ============================================================
# 5c. Remove experimental annotations from known same-chemical
# variants
# ============================================================

dat <- dat %>%
	mutate(
		canonical_name = case_when(

			compound_raw ==
				"5-Fluorouracil (+ - + -)" ~
				"5_fluorouracil",

			compound_raw ==
				"5-Fluorouracil" ~
				"5_fluorouracil",

			TRUE ~
				canonical_name
		)
	)


# ============================================================
# 5d. Create stable compound identifier
# ============================================================

dat <- dat %>%
	mutate(
		compound_id = case_when(

			!is.na(canonical_name) &
				canonical_name != "" ~
				paste0(
					"CMP_",
					canonical_name
				),

			TRUE ~
				NA_character_
		)
	)


# ============================================================
# 6. Create canonical display name
#
# compound_id defines identity.
#
# compound_name is a readable display label. The most
# frequently occurring raw name is used.
# ============================================================

compound_name_lookup <- dat %>%
	filter(
		!is.na(compound_id),
		!is.na(compound_raw)
	) %>%
	count(
		compound_id,
		canonical_name,
		compound_raw,
		name = "name_frequency"
	) %>%
	group_by(
		compound_id,
		canonical_name
	) %>%
	arrange(
		desc(name_frequency),
		compound_raw,
		.by_group = TRUE
	) %>%
	slice(1) %>%
	ungroup() %>%
	transmute(
		compound_id,
		canonical_name,
		compound_name =
			compound_raw
	)


dat <- dat %>%
	left_join(
		compound_name_lookup,
		by = c(
			"compound_id",
			"canonical_name"
		)
	)


# ============================================================
# 7. Audit compound identity
#
# Shows all raw names and CAS numbers contributing to each
# canonical compound.
# ============================================================

compound_id_audit <- dat %>%
	filter(
		!is.na(compound_id)
	) %>%
	distinct(
		compound_id,
		canonical_name,
		compound_raw,
		cas
	) %>%
	group_by(
		compound_id,
		canonical_name
	) %>%
	summarise(
		n_raw_names =
			n_distinct(
				compound_raw
			),

		raw_names =
			paste(
				sort(
					unique(
						compound_raw[
							!is.na(compound_raw)
						]
					)
				),
				collapse = " | "
			),

		n_cas =
			n_distinct(
				cas,
				na.rm = TRUE
			),

		cas_values =
			paste(
				sort(
					unique(
						cas[
							!is.na(cas)
						]
					)
				),
				collapse = " | "
			),

		.groups = "drop"
	) %>%
	arrange(
		desc(n_raw_names),
		desc(n_cas),
		canonical_name
	)


readr::write_csv(
	compound_id_audit,
	"R/results/01_compound_id_audit.csv"
)


# ============================================================
# 8. Audit multiple CAS values within one canonical compound
#
# These cases may indicate:
# 	- synonym records
# 	- salts/forms
# 	- source-data errors
# 	- identity problems
# ============================================================

compound_cas_audit <- compound_id_audit %>%
	filter(
		n_cas > 1
	)


readr::write_csv(
	compound_cas_audit,
	"R/results/01_compound_cas_audit.csv"
)


# ============================================================
# 9. Audit CAS collisions
#
# Identifies CAS numbers assigned to more than one canonical
# compound.
# ============================================================

cas_collision_audit <- dat %>%
	filter(
		!is.na(cas),
		!is.na(compound_id)
	) %>%
	distinct(
		cas,
		compound_id,
		canonical_name,
		compound_name
	) %>%
	group_by(
		cas
	) %>%
	summarise(
		n_compounds =
			n_distinct(
				compound_id
			),

		canonical_names =
			paste(
				sort(
					unique(
						canonical_name
					)
				),
				collapse = " | "
			),

		display_names =
			paste(
				sort(
					unique(
						compound_name
					)
				),
				collapse = " | "
			),

		.groups = "drop"
	) %>%
	filter(
		n_compounds > 1
	) %>%
	arrange(
		desc(n_compounds),
		cas
	)


readr::write_csv(
	cas_collision_audit,
	"R/results/01_cas_collision_audit.csv"
)


# ============================================================
# 10. Audit compounds without CAS
#
# Missing CAS does not prevent compound identification because
# canonical compound names define identity.
# ============================================================

missing_cas_audit <- dat %>%
	filter(
		is.na(cas)
	) %>%
	count(
		compound_id,
		canonical_name,
		compound_name,
		name = "n_observations"
	) %>%
	arrange(
		desc(n_observations),
		compound_name
	)


readr::write_csv(
	missing_cas_audit,
	"R/results/01_missing_cas_audit.csv"
)


# ============================================================
# 11. Standardize reference identity
# ============================================================

reference_lookup <- dat %>%
	distinct(
		reference_raw
	) %>%
	filter(
		!is.na(reference_raw)
	) %>%
	arrange(
		reference_raw
	) %>%
	mutate(
		reference_id =
			sprintf(
				"REF_%03d",
				row_number()
			)
	)


dat <- dat %>%
	left_join(
		reference_lookup,
		by = "reference_raw"
	)


readr::write_csv(
	reference_lookup,
	"R/results/01_reference_lookup.csv"
)


# ============================================================
# 12. Flag repeated experimental profiles
#
# These observations are NOT collapsed.
#
# A repeated profile means multiple rows share:
# 	- compound
# 	- reference
# 	- temperature
# 	- skin thickness
# ============================================================

dat <- dat %>%
	group_by(
		compound_id,
		reference_id,
		temperature_K,
		skin_thickness_mm
	) %>%
	mutate(
		n_same_profile =
			n(),

		repeated_profile =
			n_same_profile > 1
	) %>%
	ungroup()


repeated_profile_audit <- dat %>%
	filter(
		repeated_profile,
		!is.na(compound_id),
		!is.na(reference_id)
	) %>%
	group_by(
		compound_id,
		compound_name,
		canonical_name,
		cas,
		reference_id,
		reference_raw,
		temperature_K,
		skin_thickness_mm
	) %>%
	summarise(
		n_observations =
			n(),

		logKp_min = if (
			all(is.na(logKp))
		) {
			NA_real_
		} else {
			min(
				logKp,
				na.rm = TRUE
			)
		},

		logKp_max = if (
			all(is.na(logKp))
		) {
			NA_real_
		} else {
			max(
				logKp,
				na.rm = TRUE
			)
		},

		.groups = "drop"
	) %>%
	mutate(
		logKp_range = case_when(

			is.na(logKp_min) |
				is.na(logKp_max) ~
				NA_real_,

			TRUE ~
				logKp_max -
					logKp_min
		)
	) %>%
	arrange(
		desc(n_observations),
		desc(logKp_range)
	)


readr::write_csv(
	repeated_profile_audit,
	"R/results/01_repeated_profile_audit.csv"
)


# ============================================================
# 13. Audit missingness
# ============================================================

variables_to_audit <- c(
	"compound_id",
	"reference_id",
	"logKp",
	"temperature_K",
	"skin_thickness_mm",
	"mw",
	"log_kow",
	"melting_point_K",
	"log_saq",
	"log_soc",
	"h_donor",
	"h_acceptor",
	"molecular_volume",
	"skin_integrity_test",
	"skin_selection_criteria"
)


missingness_summary <- tibble(
	variable =
		variables_to_audit,

	n_missing =
		purrr::map_int(
			variables_to_audit,
			~ sum(
				is.na(
					dat[[.x]]
				)
			)
		)
) %>%
	mutate(
		n_total =
			nrow(dat),

		percent_missing =
			100 *
				n_missing /
				n_total
	)


readr::write_csv(
	missingness_summary,
	"R/results/01_missingness_summary.csv"
)


# ============================================================
# 14. Define eligibility for primary hierarchical analysis
#
# Required:
# 	- compound identity
# 	- reference identity
# 	- logKp
# 	- temperature
# 	- skin thickness
#
# Descriptor completeness is NOT required here.
# ============================================================

dat <- dat %>%
	mutate(
		missing_compound_id =
			is.na(compound_id),

		missing_reference_id =
			is.na(reference_id),

		missing_logKp =
			is.na(logKp),

		missing_temperature =
			is.na(temperature_K),

		missing_skin_thickness =
			is.na(skin_thickness_mm),

		hierarchical_model_eligible =
			!missing_compound_id &
			!missing_reference_id &
			!missing_logKp &
			!missing_temperature &
			!missing_skin_thickness
	)


# ============================================================
# 15. Export excluded rows with explicit reasons
# ============================================================

excluded_rows <- dat %>%
	filter(
		!hierarchical_model_eligible
	) %>%
	mutate(
		exclusion_reason =
			purrr::pmap_chr(
				list(
					missing_compound_id,
					missing_reference_id,
					missing_logKp,
					missing_temperature,
					missing_skin_thickness
				),
				function(
					missing_compound_id,
					missing_reference_id,
					missing_logKp,
					missing_temperature,
					missing_skin_thickness
				) {

					reasons <- character(0)

					if (missing_compound_id) {
						reasons <- c(
							reasons,
							"missing compound identity"
						)
					}

					if (missing_reference_id) {
						reasons <- c(
							reasons,
							"missing reference"
						)
					}

					if (missing_logKp) {
						reasons <- c(
							reasons,
							"missing logKp"
						)
					}

					if (missing_temperature) {
						reasons <- c(
							reasons,
							"missing temperature"
						)
					}

					if (missing_skin_thickness) {
						reasons <- c(
							reasons,
							"missing skin thickness"
						)
					}

					paste(
						reasons,
						collapse = "; "
					)
				}
			)
	)


readr::write_csv(
	excluded_rows,
	"R/results/01_excluded_rows.csv"
)


# ============================================================
# 16. Add compound-level bookkeeping
# ============================================================

dat <- dat %>%
	group_by(
		compound_id
	) %>%
	mutate(
		n_obs_compound =
			n(),

		n_refs_compound =
			n_distinct(
				reference_id,
				na.rm = TRUE
			),

		multi_reference_compound =
			n_refs_compound >= 2
	) %>%
	ungroup()


# ============================================================
# 17. Create analysis-ready hierarchical dataset
#
# Every retained row remains an individual permeability
# observation.
# ============================================================

analysis_data <- dat %>%
	filter(
		hierarchical_model_eligible
	) %>%
	arrange(
		compound_id,
		reference_id,
		source_row
	)


# ============================================================
# 18. Final integrity checks
# ============================================================

if (any(is.na(
	analysis_data$compound_id
))) {
	stop(
		"Analysis dataset contains missing compound_id."
	)
}


if (any(is.na(
	analysis_data$reference_id
))) {
	stop(
		"Analysis dataset contains missing reference_id."
	)
}


if (any(is.na(
	analysis_data$logKp
))) {
	stop(
		"Analysis dataset contains missing logKp."
	)
}


if (any(is.na(
	analysis_data$temperature_K
))) {
	stop(
		"Analysis dataset contains missing temperature_K."
	)
}


if (any(is.na(
	analysis_data$skin_thickness_mm
))) {
	stop(
		"Analysis dataset contains missing skin_thickness_mm."
	)
}


# ============================================================
# 19. Cleaning summary
# ============================================================

compound_counts <- analysis_data %>%
	distinct(
		compound_id,
		n_obs_compound,
		n_refs_compound
	)


cleaning_summary <- tibble(
	metric = c(
		"Raw observations after blank-row removal",
		"Hierarchical-model observations",
		"Excluded observations",
		"Unique canonical compounds",
		"Unique CAS numbers",
		"Compounds without CAS",
		"CAS numbers assigned to >1 compound",
		"Canonical compounds assigned >1 CAS",
		"Unique references",
		"Compounds with >=2 observations",
		"Compounds appearing in >=2 references",
		"Rows belonging to repeated experimental profiles"
	),

	value = c(
		n_raw,

		nrow(
			analysis_data
		),

		nrow(
			excluded_rows
		),

		n_distinct(
			analysis_data$compound_id
		),

		n_distinct(
			analysis_data$cas[
				!is.na(
					analysis_data$cas
				)
			]
		),

		n_distinct(
			analysis_data$compound_id[
				is.na(
					analysis_data$cas
				)
			]
		),

		nrow(
			cas_collision_audit
		),

		nrow(
			compound_cas_audit
		),

		n_distinct(
			analysis_data$reference_id
		),

		sum(
			compound_counts$n_obs_compound >= 2
		),

		sum(
			compound_counts$n_refs_compound >= 2
		),

		sum(
			analysis_data$repeated_profile,
			na.rm = TRUE
		)
	)
)


readr::write_csv(
	cleaning_summary,
	"R/results/01_cleaning_summary.csv"
)


# ============================================================
# 20. Export prepared datasets
# ============================================================

readr::write_csv(
	dat,
	"data/interim/01_standardized_raw_data.csv"
)


readr::write_csv(
	analysis_data,
	"data/processed/01_permeability_data.csv"
)


compound_lookup <- analysis_data %>%
	distinct(
		compound_id,
		canonical_name,
		compound_name,
		cas
	) %>%
	arrange(
		canonical_name
	)


readr::write_csv(
	compound_lookup,
	"data/processed/01_compound_lookup.csv"
)


# ============================================================
# 21. Console summary
# ============================================================

cat("\n")
cat("========================================\n")
cat("01_prepare_data.R complete\n")
cat("========================================\n\n")

print(
	cleaning_summary
)

cat("\nOutputs written to:\n")

cat(
	"  data/interim/01_standardized_raw_data.csv\n"
)

cat(
	"  data/processed/01_permeability_data.csv\n"
)

cat(
	"  data/processed/01_compound_lookup.csv\n"
)

cat(
	"  R/results/01_cleaning_summary.csv\n"
)

cat(
	"  R/results/01_numeric_coercion_audit.csv\n"
)

cat(
	"  R/results/01_compound_id_audit.csv\n"
)

cat(
	"  R/results/01_compound_cas_audit.csv\n"
)

cat(
	"  R/results/01_cas_collision_audit.csv\n"
)

cat(
	"  R/results/01_missing_cas_audit.csv\n"
)

cat(
	"  R/results/01_reference_lookup.csv\n"
)

cat(
	"  R/results/01_repeated_profile_audit.csv\n"
)

cat(
	"  R/results/01_missingness_summary.csv\n"
)

cat(
	"  R/results/01_excluded_rows.csv\n"
)