############################################################
# Clean curated SMILES lookup table
# Starting file: ../data/raw/rdkit_descriptors.csv
# Output file: ../data/interim/compound_structure_lookup.csv
############################################################
setwd("R/")

smiles <- read.csv(
	"../data/raw/rdkit_descriptors.csv",
	stringsAsFactors = FALSE
)

############################################################
# Basic column checks
############################################################

required_cols <- c(
	"compound_id",
	"Compound",
	"CAS.No",
	"SMILES_as_reported",
	"SMILES_parent",
	"structure_source",
	"structure_note",
	"structure_status",
	"structure_type"
)

missing_cols <- required_cols[
	!(required_cols %in% names(smiles))
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
# Standardize empty strings and whitespace
############################################################

for (col in names(smiles)) {
	if (is.character(smiles[[col]])) {
		smiles[[col]] <- trimws(smiles[[col]])
		smiles[[col]][smiles[[col]] == ""] <- NA
	}
}

############################################################
# Basic harmonization
############################################################

# Use CAS number as compound_id if compound_id is missing.
smiles$compound_id[
	is.na(smiles$compound_id)
] <- smiles$CAS.No[
	is.na(smiles$compound_id)
]

# Make sure compound_id and CAS.No are character strings.
smiles$compound_id <- as.character(smiles$compound_id)
smiles$CAS.No <- as.character(smiles$CAS.No)

# If SMILES_parent is missing, use SMILES_as_reported.
smiles$SMILES_parent[
	is.na(smiles$SMILES_parent)
] <- smiles$SMILES_as_reported[
	is.na(smiles$SMILES_parent)
]

# Fill default metadata.
smiles$structure_source[
	is.na(smiles$structure_source)
] <- "PubChem"

smiles$structure_note[
	is.na(smiles$structure_note)
] <- "manually checked; duplicate synonym collapsed where applicable"

smiles$structure_status[
	is.na(smiles$structure_status)
] <- "valid"

smiles$structure_type[
	is.na(smiles$structure_type)
] <- "neutral"

############################################################
# Save duplicate audit before collapsing
############################################################

id_counts <- table(smiles$compound_id)

duplicate_ids <- names(id_counts)[
	id_counts > 1
]

duplicate_rows <- smiles[
	smiles$compound_id %in% duplicate_ids,
]

duplicate_rows <- duplicate_rows[
	order(
		duplicate_rows$compound_id,
		duplicate_rows$Compound
	),
]

write.csv(
	duplicate_rows,
	"../data/interim/compound_structure_duplicate_audit_before_collapse.csv",
	row.names = FALSE
)

############################################################
# Check duplicated IDs with different SMILES
############################################################

duplicate_smiles_check <- aggregate(
	SMILES_as_reported ~ compound_id,
	data = duplicate_rows,
	function(x) {
		length(unique(na.omit(x)))
	}
)

names(duplicate_smiles_check)[2] <- "n_unique_SMILES_as_reported"

duplicate_smiles_check <- duplicate_smiles_check[
	order(duplicate_smiles_check$compound_id),
]

write.csv(
	duplicate_smiles_check,
	"../data/interim/compound_structure_duplicate_smiles_check.csv",
	row.names = FALSE
)

conflicting_ids <- duplicate_smiles_check$compound_id[
	duplicate_smiles_check$n_unique_SMILES_as_reported > 1
]

# 71-91-0 is the known intentional conflict:
# tetraethylammonium cation vs tetraethylammonium bromide salt.
allowed_conflicting_ids <- c(
	"71-91-0"
)

unexpected_conflicting_ids <- setdiff(
	conflicting_ids,
	allowed_conflicting_ids
)

if (length(unexpected_conflicting_ids) > 0) {
	stop(
		paste(
			"Unexpected duplicated compound_id with different SMILES_as_reported:",
			paste(unexpected_conflicting_ids, collapse = ", ")
		)
	)
}

############################################################
# Manually resolve 71-91-0
############################################################

# Remove all existing 71-91-0 rows.
smiles_no_tetraethylammonium <- smiles[
	smiles$compound_id != "71-91-0",
]

# Add one curated row with both salt and parent forms.
tetraethylammonium_row <- data.frame(
	compound_id = "71-91-0",
	Compound = "tetraethylammonium bromide",
	CAS.No = "71-91-0",
	SMILES_as_reported = "CC[N+](CC)(CC)CC.[Br-]",
	SMILES_parent = "CC[N+](CC)(CC)CC",
	structure_source = "PubChem",
	structure_note = "reported bromide salt; parent fragment excludes bromide counterion",
	structure_status = "valid_salt",
	structure_type = "salt",
	stringsAsFactors = FALSE
)

smiles_resolved <- rbind(
	smiles_no_tetraethylammonium,
	tetraethylammonium_row
)

############################################################
# Collapse remaining duplicate synonyms
############################################################

# For duplicated compound IDs with identical SMILES, keep one row.
# Because these have been manually checked, different names are treated as synonyms.
smiles_resolved <- smiles_resolved[
	order(
		smiles_resolved$compound_id,
		smiles_resolved$Compound
	),
]

smiles_clean <- smiles_resolved[
	!duplicated(smiles_resolved$compound_id),
]

############################################################
# Improve notes for collapsed synonym IDs
############################################################

collapsed_synonym_ids <- setdiff(
	duplicate_ids,
	"71-91-0"
)

smiles_clean$structure_note[
	smiles_clean$compound_id %in% collapsed_synonym_ids
] <- paste(
	smiles_clean$structure_note[
		smiles_clean$compound_id %in% collapsed_synonym_ids
	],
	"duplicate compound-name synonym collapsed",
	sep = "; "
)

############################################################
# Final checks
############################################################

if (any(duplicated(smiles_clean$compound_id))) {
	stop("Duplicate compound_id values remain after cleaning.")
}

if (any(is.na(smiles_clean$compound_id))) {
	stop("Missing compound_id values remain after cleaning.")
}

if (any(is.na(smiles_clean$SMILES_as_reported))) {
	stop("Missing SMILES_as_reported values remain after cleaning.")
}

if (any(is.na(smiles_clean$SMILES_parent))) {
	stop("Missing SMILES_parent values remain after cleaning.")
}

############################################################
# Reorder columns
############################################################

smiles_clean <- smiles_clean[
	,
	required_cols
]

smiles_clean <- smiles_clean[
	order(smiles_clean$compound_id),
]

############################################################
# Save cleaned structure lookup table
############################################################

write.csv(
	smiles_clean,
	"../data/interim/compound_structure_lookup.csv",
	row.names = FALSE
)

############################################################
# Console summary
############################################################

cat("SMILES curation table cleaned.\n")
cat("Original rows:", nrow(smiles), "\n")
cat("Duplicate compound IDs before collapse:", length(duplicate_ids), "\n")
cat("Final rows:", nrow(smiles_clean), "\n")
cat("Final unique compound IDs:", length(unique(smiles_clean$compound_id)), "\n")
cat("Any duplicated compound IDs remain:", any(duplicated(smiles_clean$compound_id)), "\n")
cat("Output written to: ../data/interim/compound_structure_lookup.csv\n")