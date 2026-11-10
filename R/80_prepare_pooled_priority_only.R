# Prepare the pooled dataset from the project root.
# The runner supplies output directories. Use --read-only-curation=true for fixed crosswalks.
library(readxl)

parse_options <- function(defaults) {
	for (arg in commandArgs(trailingOnly = TRUE)) {
		if (!grepl("^--[^=]+=", arg)) stop("Use --name=value: ", arg)
		key <- sub("^--([^=]+)=.*$", "\\1", arg)
		if (!key %in% names(defaults)) stop("Unknown option: --", key)
		defaults[[key]] <- sub("^--[^=]+=", "", arg)
	}
	defaults
}

options80 <- parse_options(list(
	"output-dir" = "R/results",
	"curation-policy" = "priority",
	"log-tolerance" = "0.01",
	"temperature-tolerance" = "0.2",
	"scenario" = "primary",
	"database-priority" = "SkinPiX,HuSkinDB,Cheruvu",
	"read-only-curation" = "false"
))
output_dir <- options80[["output-dir"]]
curation_policy <- options80[["curation-policy"]]
scenario <- options80[["scenario"]]
database_priority <- trimws(strsplit(options80[["database-priority"]], ",", fixed = TRUE)[[1]])
log_tolerance <- suppressWarnings(as.numeric(options80[["log-tolerance"]]))
temperature_tolerance <- suppressWarnings(as.numeric(options80[["temperature-tolerance"]]))
read_only_option <- options80[["read-only-curation"]]
if (length(database_priority) != 3L ||
	!setequal(database_priority, c("SkinPiX", "HuSkinDB", "Cheruvu"))) {
	stop("--database-priority must list SkinPiX, HuSkinDB and Cheruvu once each.")
}
if (!curation_policy %in% c("priority", "agreement")) stop("Unknown curation policy.")
if (!read_only_option %in% c("true", "false")) stop("--read-only-curation must be true or false.")
if (!nzchar(output_dir) || !nzchar(scenario)) stop("Empty output directory or scenario.")
if (any(!is.finite(c(log_tolerance, temperature_tolerance))) ||
	any(c(log_tolerance, temperature_tolerance) < 0)) {
	stop("Tolerances must be finite nonnegative numbers.")
}
read_only_curation <- read_only_option == "true"

main_output <- identical(
	normalizePath(output_dir, winslash = "/", mustWork = FALSE),
	normalizePath("R/results", winslash = "/", mustWork = FALSE)
)
final_file <- if (main_output) "data/processed/80_pooled.csv" else file.path(output_dir, "80_pooled.csv")
required_dirs <- unique(c(output_dir, dirname(final_file), "data/curation"))
missing_dirs <- required_dirs[!dir.exists(required_dirs)]
if (length(missing_dirs)) stop("Required directories are missing: ", paste(missing_dirs, collapse = ", "))
if (file.exists(final_file) && unlink(final_file) != 0L) stop("Cannot remove stale pooled output: ", final_file)

raw_paths <- c(
	Cheruvu = "data/raw/cheruvu.xlsx",
	HuSkinDB = "data/raw/huskinDB.csv",
	SkinPiX = "data/raw/skinpix.xlsx"
)
source_crosswalk <- "data/curation/80_source_crosswalk.csv"
compound_crosswalk <- "data/curation/80_compound_crosswalk.csv"
structure_file <- "data/raw/rdkit_descriptors.csv"

read_csv <- function(path) {
	read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
}
write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
write_result <- function(x, name) write_csv(x, file.path(output_dir, name))

clean <- function(x) {
	x <- trimws(as.character(x))
	x[tolower(x) %in% c("", "na", "n/a", "unknown", "-", "--")] <- NA
	x
}
number <- function(x) suppressWarnings(as.numeric(as.character(x)))
name_key <- function(x) gsub("[^a-z0-9]", "", iconv(tolower(clean(x)), to = "ASCII//TRANSLIT"))
doi_key <- function(x) {
	x <- sub("^https?://(dx\\.)?doi\\.org/", "", tolower(clean(x)))
	x <- sub("^doi:[[:space:]]*", "", x)
	gsub("[[:space:]]+", "", x)
}
blank <- function(x) {
	x <- as.character(x)
	x[is.na(x)] <- ""
	x
}
author_year <- function(x) {
	vapply(as.character(x), function(z) {
		if (is.na(z)) return(NA_character_)
		year <- regmatches(z, regexpr("(19|20)[0-9]{2}", z))
		author <- regmatches(z, regexpr("[[:alpha:]]+", z))
		if (!length(year) || !length(author)) return(NA_character_)
		paste(tolower(author), year, sep = "_")
	}, character(1))
}

# Keep manual decisions when refreshing crosswalk candidates.
sync_review <- function(template, path, key, manual) {
	if (read_only_curation && !file.exists(path)) stop("Missing curated crosswalk: ", path)
	if (file.exists(path)) {
		old <- read_csv(path)
		if (!all(c(key, manual) %in% names(old))) stop("Missing required columns in ", path)
		old_keys <- trimws(as.character(old[[key]]))
		if (anyDuplicated(old_keys)) stop("Duplicate review keys in ", path)
		j <- match(trimws(as.character(template[[key]])), old_keys)
		found <- !is.na(j)
		if (read_only_curation && any(!found)) stop("New identities need curation in ", path)
		for (field in manual) template[[field]][found] <- old[[field]][j[found]]
	}

	if ("publication_id" %in% names(template)) {
		ids <- clean(template$publication_id)
		is_doi <- !is.na(ids) & grepl(
			"^(doi:|https?://(dx\\.)?doi\\.org/|10\\.[0-9]{4,9}/)", ids, ignore.case = TRUE
		)
		ids[is_doi] <- paste0("DOI:", doi_key(ids[is_doi]))
		template$publication_id <- ids
	}
	template$decision <- tolower(trimws(blank(template$decision)))
	if (any(!template$decision %in% c("keep", "exclude", "review"))) {
		stop("Use keep, exclude, or review in ", path)
	}
	if (!read_only_curation) write_csv(template, path)
	template
}

# Read the raw records and harmonize permeability to log10(cm/h).
ch <- as.data.frame(read_excel(raw_paths[["Cheruvu"]]))
hu <- read_csv(raw_paths[["HuSkinDB"]])
sp <- as.data.frame(read_excel(raw_paths[["SkinPiX"]]))
make_records <- function(dataset, n) {
	data.frame(
		record_id = paste0(dataset, "_", seq_len(n)),
		dataset_origin = dataset,
		original_row = seq_len(n),
		compound_name = NA_character_,
		cas = NA_character_,
		smiles_raw = NA_character_,
		original_reference = NA_character_,
		doi = NA_character_,
		log_kp = NA_real_,
		temperature_c = NA_real_,
		stringsAsFactors = FALSE
	)
}

a <- make_records("Cheruvu", nrow(ch))
a$compound_name <- clean(ch$Compound)
a$cas <- clean(ch[["CAS No"]])
a$original_reference <- clean(ch$Reference)
a$log_kp <- number(ch$logkpl)
temp <- number(ch$Texpi)
a$temperature_c <- ifelse(temp > 100, temp - 273.15, temp)

b <- make_records("HuSkinDB", nrow(hu))
b$compound_name <- clean(hu[["Compound name"]])
b$smiles_raw <- clean(hu$Smiles)
b$original_reference <- clean(hu$reference)
b$doi <- doi_key(hu$DOI)
b$log_kp <- number(hu[["logkp (cm/s)"]]) + log10(3600)
b$temperature_c <- number(hu[["donor/skin surface temperature (°C)"]])

d <- make_records("SkinPiX", nrow(sp))
d$compound_name <- clean(sp[["compound name"]])
d$cas <- clean(sp[["CAS number"]])
d$smiles_raw <- clean(sp$SMILES)
d$original_reference <- paste0(
	blank(sp$author), " ", blank(sp$date), "; publication ID=", blank(sp[["publication ID"]])
)
d$doi <- doi_key(sp$doi)
kp <- number(sp[["Kp (cm/h) (processed)"]])
relation <- clean(sp[["Kp relation"]])
valid_sp <- !is.na(relation) & relation == "=" & is.finite(kp) & kp > 0
d$log_kp[valid_sp] <- log10(kp[valid_sp])
d$temperature_c <- number(sp[["donor/skin surface temperature (°C)"]])

raw <- rbind(a, b, d)
valid <- !is.na(raw$compound_name) & is.finite(raw$log_kp)
if (main_output) write_result(raw[!valid, ], "80_excluded_input_records.csv")
dat <- raw[valid, ]
dat$source_key <- paste(
	dat$dataset_origin,
	ifelse(!is.na(dat$doi), paste0("DOI:", dat$doi), paste0("REF:", blank(dat$original_reference))),
	sep = " || "
)
dat$compound_key <- paste(
	dat$dataset_origin, blank(dat$compound_name), blank(dat$cas), blank(dat$smiles_raw), sep = " || "
)
hirvonen_1998b_salts <- dat$record_id %in% c("Cheruvu_380", "Cheruvu_382")
dat$compound_key[hirvonen_1998b_salts] <- paste0(
	trimws(dat$compound_key[hirvonen_1998b_salts]), " Hirvonen 1998b"
)

# Resolve publication identities. Author/year matches are review candidates.
sources <- unique(dat[c("source_key", "dataset_origin", "original_reference", "doi")])
if (anyDuplicated(sources$source_key)) stop("A source key has conflicting reference metadata.")
sources$candidate_publication_id <- ifelse(
	!is.na(sources$doi), paste0("DOI:", sources$doi), paste0("CHREF:", sources$original_reference)
)
sources$publication_id <- sources$candidate_publication_id
sources$decision <- ifelse(!is.na(sources$doi), "keep", "review")
sources$notes <- ""
known <- sources[!is.na(sources$doi), ]
known$author_year <- author_year(known$original_reference)
sources$author_year_candidates <- ""

for (i in which(sources$dataset_origin == "Cheruvu")) {
	reference <- sources$original_reference[i]
	key <- author_year(reference)
	hits <- unique(known$candidate_publication_id[!is.na(known$author_year) & known$author_year == key])
	hits <- hits[!is.na(hits)]
	sources$author_year_candidates[i] <- paste(hits, collapse = "; ")
	n_years <- length(regmatches(reference, gregexpr("(19|20)[0-9]{2}", reference))[[1]])
	if (n_years > 1) {
		sources$notes[i] <- "Multiple cited publications: resolve attribution."
	} else if (length(hits) == 1) {
		sources$candidate_publication_id[i] <- hits
		sources$publication_id[i] <- hits
		sources$notes[i] <- "Author/year candidate only: verify publication."
	} else {
		sources$notes[i] <- "Verify publication and possible cross-dataset overlap."
	}
}
sources <- sync_review(sources, source_crosswalk, "source_key", c("publication_id", "decision", "notes"))
if (any(sources$decision == "keep" & is.na(clean(sources$publication_id)))) {
	stop("Every kept source needs a publication_id.")
}
if (main_output) write_result(sources, "80_source_crosswalk_snapshot.csv")

# Resolve compound structures from the curated crosswalk.
compounds <- unique(dat[c("compound_key", "dataset_origin", "compound_name", "cas", "smiles_raw")])
compounds$smiles <- compounds$smiles_raw
compounds$decision <- ifelse(!is.na(compounds$smiles_raw), "keep", "review")
compounds$notes <- ""
if (file.exists(structure_file)) {
	structures <- read_csv(structure_file)
	if (!all(c("Compound", "CAS.No", "SMILES_parent") %in% names(structures))) {
		stop("Unexpected structure mapping columns in ", structure_file)
	}
	structure_names <- name_key(structures$Compound)
	structure_cas <- clean(structures$CAS.No)
	for (i in which(is.na(compounds$smiles))) {
		nm <- name_key(compounds$compound_name[i])
		cas <- clean(compounds$cas[i])
		hit_name <- !is.na(structure_names) & structure_names == nm
		hit_cas <- !is.na(structure_cas) & !is.na(cas) & structure_cas == cas
		candidates <- unique(clean(structures$SMILES_parent[hit_name | hit_cas]))
		candidates <- candidates[!is.na(candidates)]
		if (length(candidates) == 1) {
			compounds$smiles[i] <- candidates
			compounds$notes[i] <- "Structure-file candidate: verify name, CAS, chemical form and reconstructed-structure provenance."
		} else {
			compounds$notes[i] <- "Resolve structure manually."
		}
	}
}
compounds <- sync_review(compounds, compound_crosswalk, "compound_key", c("smiles", "decision", "notes"))

# Canonicalize SMILES and compute descriptors while retaining salts and stereochemistry.
structure_input <- unique(data.frame(smiles = clean(compounds$smiles), stringsAsFactors = FALSE))
structure_input <- structure_input[!is.na(structure_input$smiles), , drop = FALSE]
infile <- tempfile(fileext = ".csv")
outfile <- tempfile(fileext = ".csv")
pyfile <- tempfile(fileext = ".py")
write_csv(structure_input, infile)

python_code <- '
import csv
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors

fields = [
	"smiles", "compound_id", "mw", "logp", "hbd", "hba",
	"tpsa", "rotatable_bonds", "structure_error"
]

with open(sys.argv[1], newline="", encoding="utf-8") as src, open(sys.argv[2], "w", newline="", encoding="utf-8") as dst:
	writer = csv.DictWriter(dst, fieldnames=fields)
	writer.writeheader()
	for row in csv.DictReader(src):
		smi = row["smiles"]
		result = dict.fromkeys(fields, "")
		result["smiles"] = smi
		try:
			mol = Chem.MolFromSmiles(smi)
			if mol is None:
				raise ValueError("Invalid SMILES")
			result.update(
				compound_id=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
				mw=Descriptors.MolWt(mol),
				logp=Crippen.MolLogP(mol),
				hbd=Lipinski.NumHDonors(mol),
				hba=Lipinski.NumHAcceptors(mol),
				tpsa=rdMolDescriptors.CalcTPSA(mol),
				rotatable_bonds=Lipinski.NumRotatableBonds(mol)
			)
		except Exception as exc:
			result["structure_error"] = str(exc)
		writer.writerow(result)
'
writeLines(python_code, pyfile)
python <- unname(Sys.which("python"))
if (!nzchar(python)) stop("Activate Python with RDKit installed.")
status <- system2(python, c(shQuote(pyfile), shQuote(infile), shQuote(outfile)))
if (status != 0L) stop("RDKit descriptor calculation failed.")
desc <- read_csv(outfile)
unlink(c(infile, outfile, pyfile))

# Apply manual decisions and protocol exclusions.
si <- match(dat$source_key, sources$source_key)
ci <- match(dat$compound_key, compounds$compound_key)
dat$publication_id <- sources$publication_id[si]
dat$source_decision <- sources$decision[si]
dat$compound_decision <- compounds$decision[ci]
dat$smiles <- compounds$smiles[ci]
di <- match(dat$smiles, desc$smiles)
descriptor_fields <- c("compound_id", "mw", "logp", "hbd", "hba", "tpsa", "rotatable_bonds", "structure_error")
for (field in descriptor_fields) dat[[field]] <- desc[[field]][di]

identity_excluded <- dat$source_decision == "exclude" | dat$compound_decision == "exclude"
# Omit the potentially overlapping Williams and Barry summary.
summary_excluded <- !identity_excluded & dat$record_id %in% "Cheruvu_137"
# Retain Hirvonen 1998 control protocols III and VII.
protocol_excluded <- !identity_excluded &
	dat$record_id %in% paste0("Cheruvu_", c(359:362, 365:370, 373:376))
excluded <- identity_excluded | summary_excluded | protocol_excluded
resolved <- !excluded & dat$source_decision == "keep" & dat$compound_decision == "keep" &
	!is.na(dat$compound_id) & is.finite(dat$mw) & is.finite(dat$logp)
unresolved <- !excluded & !resolved
n_unresolved <- sum(unresolved)
unresolved_file <- file.path(output_dir, "80_unresolved_records.csv")
if (n_unresolved > 0L) {
	write_csv(dat[unresolved, ], unresolved_file)
} else {
	unlink(unresolved_file)
}
pool <- dat[resolved, ]

# Correct nine unit errors verified against the original publications.
unit_corrections <- data.frame(
	record_id = c(paste0("HuSkinDB_", 30:37), "HuSkinDB_501"),
	expected_log_kp = c(-2.839, -4.343, -3.192, -6.979, -3.514, -3.113, -3.830, -2.807, -9.588) + log10(3600),
	original_kp_cm_h = c(1.45e-3, 4.54e-5, 6.43e-4, 1.05e-7, 3.06e-4, 7.71e-4, 1.48e-4, 1.56e-3, 0.93e-3),
	evidence = c(
		rep("Barber 1992; DOI:10.1016/0272-0590(92)90086-w; Table 1, p.495; human stratum corneum", 8),
		"Barber 1995; DOI:10.1016/0378-4274(95)03393-y; Table 2, p.170; human water, pre-exposure"
	),
	stringsAsFactors = FALSE
)
unit_corrections$corrected_log_kp <- log10(unit_corrections$original_kp_cm_h)
correction_rows <- match(unit_corrections$record_id, pool$record_id)
if (anyNA(correction_rows)) {
	stop("Correction records missing from pool: ", paste(unit_corrections$record_id[is.na(correction_rows)], collapse = ", "))
}
unit_corrections$previous_log_kp <- pool$log_kp[correction_rows]
if (any(abs(unit_corrections$previous_log_kp - unit_corrections$expected_log_kp) > 1e-6)) {
	stop("Unexpected values in correction records. Check before proceeding.")
}
pool$log_kp[correction_rows] <- unit_corrections$corrected_log_kp
if (main_output) write_result(unit_corrections, "80_verified_unit_corrections.csv")

# Primary pooling selects one database per publication.
auto_decision <- rep("keep", nrow(pool))
if (curation_policy == "priority") {
	auto_reason <- rep("", nrow(pool))
	for (ii in split(seq_len(nrow(pool)), pool$publication_id)) {
		origins <- unique(pool$dataset_origin[ii])
		preferred <- database_priority[database_priority %in% origins][1]
		auto_decision[ii] <- ifelse(pool$dataset_origin[ii] == preferred, "keep", "exclude")
		auto_reason[ii] <- paste0(
			"Database-priority curation: retained ", preferred, " representation of the publication. ",
			"Selection did not use permeability agreement or temperature."
		)
	}
} else {
	# Agreement scenarios compare records within each publication and compound.
	auto_reason <- rep("Retained: publication-compound group occurs in one database.", nrow(pool))
	group_key <- paste(pool$publication_id, pool$compound_id, sep = "\r")
	groups <- split(seq_len(nrow(pool)), match(group_key, unique(group_key)))
	for (ii in groups) {
		origins <- sort(unique(pool$dataset_origin[ii]))
		if (length(origins) < 2L) next
		conflict_reasons <- character()

		for (datasets in combn(origins, 2, simplify = FALSE)) {
			aa <- ii[pool$dataset_origin[ii] == datasets[1]]
			bb <- ii[pool$dataset_origin[ii] == datasets[2]]
			aa <- aa[order(pool$record_id[aa])]
			bb <- bb[order(pool$record_id[bb])]
			difference <- abs(outer(pool$log_kp[aa], pool$log_kp[bb], "-"))
			temp_a <- pool$temperature_c[aa]
			temp_b <- pool$temperature_c[bb]
			known_temperature <- outer(is.finite(temp_a), is.finite(temp_b), "&")
			temperature_difference <- abs(outer(temp_a, temp_b, "-"))
			temperature_conflict <- known_temperature & temperature_difference > temperature_tolerance + 1e-10

			# Prioritize compatible temperatures, then recorded temperatures, then permeability agreement.
			base <- min(length(aa), length(bb)) + 1
			cost <- temperature_conflict * base^2 + (!known_temperature) * base + difference / (1 + max(difference))
			if (nrow(cost) <= ncol(cost)) {
				rows <- seq_len(nrow(cost))
				columns <- as.integer(clue::solve_LSAP(cost))
			} else {
				columns <- seq_len(ncol(cost))
				rows <- as.integer(clue::solve_LSAP(t(cost)))
			}
			selected <- cbind(rows, columns)

			if (any(difference[selected] > log_tolerance + 1e-10)) {
				conflict_reasons <- c(conflict_reasons, paste0(
					paste(datasets, collapse = " vs "), ": paired log Kp difference > ", log_tolerance
				))
			}
			if (any(temperature_conflict[selected])) {
				conflict_reasons <- c(conflict_reasons, paste0(
					paste(datasets, collapse = " vs "), ": paired temperature difference > ", temperature_tolerance, " C"
				))
			}
		}

		if (length(conflict_reasons)) {
			auto_decision[ii] <- "exclude"
			auto_reason[ii] <- paste(
				"AUTO: excluded entire overlapping group.", paste(unique(conflict_reasons), collapse = "; ")
			)
		} else {
			preferred <- database_priority[database_priority %in% origins][1]
			auto_decision[ii] <- ifelse(pool$dataset_origin[ii] == preferred, "keep", "exclude")
			auto_reason[ii] <- paste0(
				"AUTO: retained complete ", preferred, " representation; other database representations excluded. ",
				"No paired log Kp difference above ", log_tolerance,
				" or paired temperature difference above ", temperature_tolerance, " C. ",
				"Unmatched records retained only in the preferred database."
			)
		}
	}
}
pooled <- pool[auto_decision == "keep", , drop = FALSE]

# Keep manuscript and downstream audits for the primary dataset.
if (main_output) {
	auto_audit <- pool[c(
		"record_id", "dataset_origin", "original_row", "publication_id",
		"compound_id", "compound_name", "log_kp", "temperature_c"
	)]
	auto_audit$decision <- auto_decision
	auto_audit$notes <- auto_reason
	write_result(auto_audit, "80_automatic_measurement_decisions.csv")

	coverage_fields <- c("mw", "logp", "hbd", "hba", "tpsa", "rotatable_bonds", "temperature_c")
	coverage <- data.frame(
		variable = coverage_fields,
		n_available = vapply(coverage_fields, function(v) sum(is.finite(pooled[[v]])), integer(1)),
		n_records = nrow(pooled)
	)
	write_result(coverage, "80_descriptor_coverage.csv")

	input_paths <- c(unname(raw_paths), source_crosswalk, compound_crosswalk)
	if (file.exists(structure_file)) input_paths <- c(input_paths, structure_file)
	write_result(
		data.frame(path = input_paths, md5 = unname(tools::md5sum(input_paths))),
		"80_input_checksums.csv"
	)
}

# Save settings for the manuscript and sensitivity runners.
ready <- n_unresolved == 0L && nrow(pooled) > 0L
settings <- data.frame(
	scenario = scenario,
	curation_policy = curation_policy,
	log_tolerance = if (curation_policy == "agreement") log_tolerance else NA_real_,
	temperature_tolerance = if (curation_policy == "agreement") temperature_tolerance else NA_real_,
	database_priority = paste(database_priority, collapse = " > "),
	read_only_curation = read_only_curation,
	verified_unit_corrections = nrow(unit_corrections),
	analysis_ready = ready,
	n_valid_input = nrow(dat),
	n_excluded_identity = sum(identity_excluded),
	n_excluded_overlapping_summary = sum(summary_excluded),
	n_excluded_protocol = sum(protocol_excluded),
	n_unresolved = n_unresolved,
	n_measurements = nrow(pooled),
	n_compounds = length(unique(pooled$compound_id)),
	n_publications = length(unique(pooled$publication_id)),
	stringsAsFactors = FALSE
)
if (ready) write_csv(pooled, final_file)
write_result(settings, "80_analysis_settings.csv")
if (n_unresolved > 0L) stop("Unresolved identity or structure records. See ", unresolved_file)
if (!ready) stop("No measurements remain after pooling.")
message("Wrote ", nrow(pooled), " measurements to ", final_file, " (", scenario, ").")
