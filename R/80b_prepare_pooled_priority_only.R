# 80_prepare_audit_pooled_data.R
#
# Run from the project root:
# Rscript R/80_prepare_audit_pooled_data.R
#
# Requires readxl and Python/RDKit available as "python".
#
# Review files are preserved in data/curation/.
# Generated audits are written to R/results/80_*.
#
# This script does not fit models or automatically collapse duplicates.

library(readxl)

dir.create("data/curation", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("R/results", recursive = TRUE, showWarnings = FALSE)

# Remove a stale final output: it must pass review again on this run.
final_file <- "data/processed/80_pooled.csv"
if (file.exists(final_file)) unlink(final_file)

read_csv <- function(path) {
  read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA")
  )
}

write_csv <- function(x, path) {
  write.csv(x, path, row.names = FALSE, na = "")
}

clean <- function(x) {
  x <- trimws(as.character(x))
  x[tolower(x) %in% c("", "na", "n/a", "unknown", "-", "--")] <- NA
  x
}

number <- function(x) suppressWarnings(as.numeric(as.character(x)))

name_key <- function(x) {
  x <- iconv(tolower(clean(x)), to = "ASCII//TRANSLIT")
  gsub("[^a-z0-9]", "", x)
}

doi_key <- function(x) {
  x <- tolower(clean(x))
  x <- sub("^https?://(dx\\.)?doi\\.org/", "", x)
  x <- sub("^doi:[[:space:]]*", "", x)
  gsub("[[:space:]]+", "", x)
}

blank <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

# One first-author/year candidate key.
# Used only to propose Cheruvu matches, never to confirm them.
author_year <- function(x) {
  vapply(as.character(x), function(z) {
    if (is.na(z)) return(NA_character_)
    year <- regmatches(z, regexpr("(19|20)[0-9]{2}", z))
    author <- regmatches(z, regexpr("[[:alpha:]]+", z))
    if (!length(year) || !length(author)) return(NA_character_)
    paste(tolower(author), year, sep = "_")
  }, character(1))
}

# Refresh generated fields while preserving specified manual fields.
sync_review <- function(template, path, key, manual) {
  if (file.exists(path)) {
    old <- read_csv(path)

    if (!all(c(key, manual) %in% names(old))) {
      stop("Missing required columns in ", path)
    }
    if (anyDuplicated(old[[key]])) {
      stop("Duplicate review keys in ", path)
    }

    template_keys <- trimws(as.character(template[[key]]))
	old_keys <- trimws(as.character(old[[key]]))

	if (anyDuplicated(old_keys)) {
	  stop("Duplicate review keys after trimming whitespace in ", path)
	}

	j <- match(template_keys, old_keys)
    found <- !is.na(j)

    for (field in manual) {
      template[[field]][found] <- old[[field]][j[found]]
    }
  }

  template$decision <- tolower(trimws(blank(template$decision)))

  if (any(!template$decision %in% c("keep", "exclude", "review"))) {
    stop("Use keep, exclude, or review in ", path)
  }

  write_csv(template, path)
  template
}

# ------------------------------------------------------------
# 1. Read original files and preserve row provenance
# ------------------------------------------------------------

ch <- as.data.frame(read_excel("data/raw/cheruvu.xlsx"))
hu <- read_csv("data/raw/huskinDB.csv")
sp <- as.data.frame(read_excel("data/raw/skinpix.xlsx"))

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
  blank(sp$author), " ", blank(sp$date),
  "; publication ID=", blank(sp[["publication ID"]])
)
d$doi <- doi_key(sp$doi)

kp <- number(sp[["Kp (cm/h) (processed)"]])
relation <- clean(sp[["Kp relation"]])
valid_sp <- !is.na(relation) & relation == "=" &
  is.finite(kp) & kp > 0

d$log_kp[valid_sp] <- log10(kp[valid_sp])
d$temperature_c <- number(
  sp[["donor/skin surface temperature (°C)"]]
)

raw <- rbind(a, b, d)

valid <- !is.na(raw$compound_name) & is.finite(raw$log_kp)
write_csv(raw[!valid, ], "R/results/80_excluded_input_records.csv")
dat <- raw[valid, ]

# Stable keys retain the original reference/name information.
dat$source_key <- paste(
  dat$dataset_origin,
  ifelse(
    !is.na(dat$doi),
    paste0("DOI:", dat$doi),
    paste0("REF:", blank(dat$original_reference))
  ),
  sep = " || "
)

dat$compound_key <- paste(
  dat$dataset_origin, blank(dat$compound_name),
  blank(dat$cas), blank(dat$smiles_raw),
  sep = " || "
)

# ------------------------------------------------------------
# 2. Source crosswalk
# ------------------------------------------------------------

sources <- unique(dat[c(
  "source_key", "dataset_origin", "original_reference", "doi"
)])

if (anyDuplicated(sources$source_key)) {
  stop("A source key has conflicting reference metadata.")
}

sources$candidate_publication_id <- ifelse(
  !is.na(sources$doi),
  paste0("DOI:", sources$doi),
  paste0("CHREF:", sources$original_reference)
)

sources$publication_id <- sources$candidate_publication_id
sources$decision <- ifelse(!is.na(sources$doi), "keep", "review")
sources$notes <- ""

# DOI-bearing references provide author/year candidates.
known <- sources[!is.na(sources$doi), ]
known$author_year <- author_year(known$original_reference)

sources$author_year_candidates <- ""

for (i in which(sources$dataset_origin == "Cheruvu")) {
  key <- author_year(sources$original_reference[i])
  hits <- unique(known$candidate_publication_id[
    !is.na(known$author_year) & known$author_year == key
  ])
  hits <- hits[!is.na(hits)]

  sources$author_year_candidates[i] <- paste(hits, collapse = "; ")

  # Several cited years indicate potentially mixed attribution.
  years <- gregexpr(
    "(19|20)[0-9]{2}", sources$original_reference[i]
  )
  n_years <- length(regmatches(sources$original_reference[i], years)[[1]])

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

sources <- sync_review(
  sources, "data/curation/80_source_crosswalk.csv",
  "source_key", c("publication_id", "decision", "notes")
)

if (any(sources$decision == "keep" &
        is.na(clean(sources$publication_id)))) {
  stop("Every kept source needs a publication_id.")
}

# ------------------------------------------------------------
# 3. Compound crosswalk
# ------------------------------------------------------------

compounds <- unique(dat[c(
  "compound_key", "dataset_origin", "compound_name", "cas", "smiles_raw"
)])

compounds$smiles <- compounds$smiles_raw
compounds$decision <- ifelse(
  !is.na(compounds$smiles_raw), "keep", "review"
)
compounds$notes <- ""

# The existing structure file supplies candidates, not verified mappings.
structure_file <- "data/raw/rdkit_descriptors.csv"

if (file.exists(structure_file)) {
  structures <- read_csv(structure_file)

  required <- c("Compound", "CAS.No", "SMILES_parent")
  if (!all(required %in% names(structures))) {
    stop("Unexpected structure mapping columns in ", structure_file)
  }

  for (i in which(is.na(compounds$smiles))) {
    nm <- name_key(compounds$compound_name[i])
    cas <- clean(compounds$cas[i])

    hit_name <- !is.na(name_key(structures$Compound)) &
      name_key(structures$Compound) == nm
    hit_cas <- !is.na(clean(structures$CAS.No)) &
      !is.na(cas) & clean(structures$CAS.No) == cas

    candidates <- unique(clean(
      structures$SMILES_parent[hit_name | hit_cas]
    ))
    candidates <- candidates[!is.na(candidates)]

    if (length(candidates) == 1) {
      compounds$smiles[i] <- candidates
      compounds$notes[i] <- paste(
        "Structure-file candidate: verify name, CAS, chemical form",
        "and reconstructed-structure provenance."
      )
    } else {
      compounds$notes[i] <- "Resolve structure manually."
    }
  }
}

compounds <- sync_review(
  compounds, "data/curation/80_compound_crosswalk.csv",
  "compound_key", c("smiles", "decision", "notes")
)

# ------------------------------------------------------------
# 4. Standardize structures and calculate descriptors
# ------------------------------------------------------------

structure_input <- unique(data.frame(
  smiles = clean(compounds$smiles),
  stringsAsFactors = FALSE
))
structure_input <- structure_input[
  !is.na(structure_input$smiles), , drop = FALSE
]

infile <- tempfile(fileext = ".csv")
outfile <- tempfile(fileext = ".csv")
pyfile <- tempfile(fileext = ".py")
write_csv(structure_input, infile)

# Preserve molecular fragments, charge, and specified stereochemistry.
# Do not silently convert salts to neutral parent compounds.
python_code <- '
import csv
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors

fields = [
    "smiles", "compound_id", "mw", "logp",
    "hbd", "hba", "tpsa", "rotatable_bonds", "structure_error"
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

            canonical = Chem.MolToSmiles(
                mol, canonical=True, isomericSmiles=True
            )
            result.update(
                compound_id=canonical,
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

status <- system2(
  python,
  c(shQuote(pyfile), shQuote(infile), shQuote(outfile))
)
if (status != 0) stop("RDKit descriptor calculation failed.")

desc <- read_csv(outfile)
unlink(c(infile, outfile, pyfile))

write_csv(desc, "R/results/80_structure_descriptor_audit.csv")

# ------------------------------------------------------------
# 5. Attach decisions and retain only resolved records for preview
# ------------------------------------------------------------

si <- match(dat$source_key, sources$source_key)
ci <- match(dat$compound_key, compounds$compound_key)

dat$publication_id <- sources$publication_id[si]
dat$source_decision <- sources$decision[si]
dat$compound_decision <- compounds$decision[ci]
dat$smiles <- compounds$smiles[ci]

di <- match(dat$smiles, desc$smiles)
descriptor_fields <- c(
  "compound_id", "mw", "logp", "hbd", "hba",
  "tpsa", "rotatable_bonds", "structure_error"
)
for (field in descriptor_fields) dat[[field]] <- desc[[field]][di]

excluded <- dat$source_decision == "exclude" |
  dat$compound_decision == "exclude"

resolved <- !excluded &
  dat$source_decision == "keep" &
  dat$compound_decision == "keep" &
  !is.na(dat$compound_id) &
  is.finite(dat$mw) & is.finite(dat$logp)

unresolved <- !excluded & !resolved
write_csv(dat[unresolved, ], "R/results/80_unresolved_records.csv")

pool <- dat[resolved, ]

# ------------------------------------------------------------
# 6. Database-priority-only sensitivity analysis
# ------------------------------------------------------------

database_priority <- c("SkinPiX", "HuSkinDB", "Cheruvu")

if (any(!pool$dataset_origin %in% database_priority)) {
  stop("Unexpected dataset origin.")
}

auto_decision <- rep("keep", nrow(pool))
auto_reason <- rep("", nrow(pool))
overlap <- rep(FALSE, nrow(pool))

publication_groups <- split(
  seq_len(nrow(pool)),
  pool$publication_id
)

for (ii in publication_groups) {
  origins <- unique(pool$dataset_origin[ii])
  preferred <- database_priority[
    database_priority %in% origins
  ][1]

  overlap[ii] <- length(origins) > 1L

  auto_decision[ii] <- ifelse(
    pool$dataset_origin[ii] == preferred,
    "keep",
    "exclude"
  )

  auto_reason[ii] <- paste0(
    "Database-priority-only sensitivity: retained ",
    preferred,
    " representation of the publication. ",
    "Selection did not use permeability agreement or temperature."
  )
}

auto_audit <- pool[c(
  "record_id", "dataset_origin", "original_row",
  "publication_id", "compound_id", "compound_name",
  "log_kp", "temperature_c"
)]

auto_audit$decision <- auto_decision
auto_audit$notes <- auto_reason

write_csv(
  auto_audit,
  "R/results/80_automatic_measurement_decisions.csv"
)

# Compatibility with Section 8: no manual measurement review.
measurement_review <- auto_audit[overlap, , drop = FALSE]

write_csv(
  pool,
  "R/results/80_resolved_records_before_deduplication.csv"
)

preview <- pool[auto_decision == "keep", , drop = FALSE]

write_csv(
  preview,
  "R/results/80_pooled_preview.csv"
)

cat(
  "\nDatabase-priority-only sensitivity\n",
  "Resolved measurements: ", nrow(pool), "\n",
  "Retained measurements: ", nrow(preview), "\n",
  "Excluded measurements: ",
  sum(auto_decision == "exclude"), "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Source-held-out anchor eligibility
# ------------------------------------------------------------

eligibility <- function(x) {
  ids <- sort(unique(x$publication_id))

  if (!length(ids)) {
    return(data.frame(
      publication_id = character(),
      dataset_origins = character(),
      n_measurements = integer(),
      n_compounds = integer(),
      n_anchors = integer(),
      n_test_compounds = integer(),
      n_test_measurements = integer()
    ))
  }

  do.call(rbind, lapply(ids, function(id) {
    held <- x[x$publication_id == id, ]
    train <- x[x$publication_id != id, ]

    held_compounds <- unique(held$compound_id)
    train_compounds <- unique(train$compound_id)
    anchors <- intersect(held_compounds, train_compounds)
    tests <- setdiff(held_compounds, train_compounds)

    data.frame(
      publication_id = id,
      dataset_origins = paste(
        sort(unique(held$dataset_origin)), collapse = "; "
      ),
      n_measurements = nrow(held),
      n_compounds = length(held_compounds),
      n_anchors = length(anchors),
      n_test_compounds = length(tests),
      n_test_measurements = sum(held$compound_id %in% tests),
      stringsAsFactors = FALSE
    )
  }))
}

eligible <- eligibility(preview)
for (k in 1:5) {
  eligible[[paste0("eligible_", k, "_anchors")]] <-
    eligible$n_anchors >= k & eligible$n_test_compounds >= 1
}

write_csv(eligible, "R/results/80_source_eligibility.csv")

eligible_temp <- eligibility(
  preview[is.finite(preview$temperature_c), ]
)
write_csv(
  eligible_temp,
  "R/results/80_temperature_complete_source_eligibility.csv"
)

panel_summary <- data.frame(
  anchors = 1:5,
  n_sources = vapply(1:5, function(k) {
    sum(eligible$n_anchors >= k & eligible$n_test_compounds >= 1)
  }, integer(1))
)
write_csv(panel_summary, "R/results/80_anchor_panel_eligibility.csv")

# Counts describe resolved, retained records, not an imputed full pool.
coverage <- data.frame(
  variable = c(
    "mw", "logp", "hbd", "hba", "tpsa",
    "rotatable_bonds", "temperature_c"
  ),
  n_available = vapply(
    c("mw", "logp", "hbd", "hba", "tpsa",
      "rotatable_bonds", "temperature_c"),
    function(v) sum(is.finite(preview[[v]])),
    integer(1)
  ),
  n_records = nrow(preview)
)
write_csv(coverage, "R/results/80_descriptor_coverage.csv")

# ------------------------------------------------------------
# 8. Final readiness gate
# ------------------------------------------------------------

n_unresolved <- sum(unresolved)
n_measurement_review <- sum(
  measurement_review$decision == "review"
)

ready <- n_unresolved == 0 &&
  n_measurement_review == 0 &&
  nrow(preview) > 0

if (ready) {
  write_csv(preview, final_file)
}

summary_lines <- c(
  paste("Analysis-ready pooled file:", ready),
  paste("Valid input records:", nrow(dat)),
  paste("Excluded by source/compound decisions:", sum(excluded)),
  paste("Unresolved identity/structure records:", n_unresolved),
  paste("Measurement records awaiting review:", n_measurement_review),
  paste("Retained preview measurements:", nrow(preview)),
  paste("Retained preview compounds:", length(unique(preview$compound_id))),
  paste("Retained preview publications:", length(unique(preview$publication_id))),
  "",
  "Anchor eligibility in the currently resolved preview:",
  capture.output(print(panel_summary, row.names = FALSE)),
  "",
  "Eligibility can change as unresolved records are added.",
  "Publication matches require bibliographic verification.",
  "Shared structures do not establish shared original measurements.",
  "Inspect original experimental conditions before deduplication."
)

writeLines(summary_lines, "R/results/80_merge_status.txt")
cat(paste(summary_lines, collapse = "\n"), "\n")

if (!ready) {
  message(
    "\nReview data/curation/80_* files and rerun.",
    "\nThe preview is provisional; no final pooled dataset was written."
  )
}