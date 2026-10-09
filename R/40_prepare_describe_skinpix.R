# ============================================================
# 40_prepare_describe_skinpix.R
# ============================================================

library(readxl)

dat0 <- as.data.frame(read_excel("data/raw/skinpix.xlsx"))
huskin <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)

n_raw <- nrow(dat0)
dat0$source_row <- seq_len(nrow(dat0))

clean_text <- function(x) {
	x <- trimws(as.character(x))
	x[tolower(x) %in% c("", "na", "n/a", "unknown", "-", "--")] <- NA
	x
}

name_key <- function(x) {
	x <- clean_text(x)
	x <- tolower(iconv(x, to = "ASCII//TRANSLIT"))
	x <- gsub("[^a-z0-9]+", "_", x)
	gsub("^_|_$", "", x)
}

doi_key <- function(x) {
	x <- tolower(clean_text(x))
	x <- gsub("^https?://(dx\\.)?doi\\.org/", "", x)
	gsub("\\s+", "", x)
}

text_cols <- c(
	"publication ID", "SMILES", "CAS number", "compound name", "Kp relation",
	"skin donor type", "skin source site", "skin preparation", "used layer",
	"donor type", "category donor type", "acceptor type", "category acceptor type",
	"cell type", "author", "doi"
)

for (x in text_cols) dat0[[x]] <- clean_text(dat0[[x]])

skinpix_smiles <- data.frame(smiles = unique(dat0$SMILES[!is.na(dat0$SMILES)]))
huskin_smiles <- data.frame(smiles = unique(huskin$smiles[!is.na(huskin$smiles)]))

write.csv(skinpix_smiles, "data/processed/40_skinpix_smiles.csv", row.names = FALSE)
write.csv(huskin_smiles, "data/processed/40_huskindb_smiles.csv", row.names = FALSE)

status <- system2("python", "R/40_make_rdkit_descriptors.py")
if (status != 0) stop("RDKit descriptor script failed.")

skinpix_desc <- read.csv("data/processed/40_skinpix_rdkit_descriptors.csv", stringsAsFactors = FALSE)
huskin_desc <- read.csv("data/processed/40_huskindb_rdkit_descriptors.csv", stringsAsFactors = FALSE)

num_cols <- c(
	"Kp (cm/h) (processed)", "donor pH", "donor/skin surface temperature (°C)",
	"acceptor pH", "acceptor temperature (°C)"
)

for (x in num_cols) dat0[[x]] <- suppressWarnings(as.numeric(dat0[[x]]))

non_equal <- !is.na(dat0[["Kp relation"]]) & dat0[["Kp relation"]] != "="
bad_kp <- is.na(dat0[["Kp (cm/h) (processed)"]]) | dat0[["Kp (cm/h) (processed)"]] <= 0
keep <- !non_equal & !bad_kp

dat <- dat0[keep, ]

dat$smiles_key <- clean_text(dat$SMILES)
dat$name_key <- name_key(dat[["compound name"]])

compound_key <- ifelse(!is.na(dat$smiles_key), paste0("SMILES:", dat$smiles_key), paste0("NAME:", dat$name_key))
if (any(is.na(compound_key) | compound_key == "NAME:")) stop("Missing compound identity.")

dat$compound_id <- paste0("SPX_C", match(compound_key, unique(compound_key)))

source_key <- ifelse(!is.na(dat[["publication ID"]]), paste0("PUB:", dat[["publication ID"]]), paste0("DOI:", doi_key(dat$doi)))
if (any(is.na(source_key) | source_key == "DOI:")) stop("Missing source identity.")

dat$reference_id <- paste0("SPX_R", match(source_key, unique(source_key)))

# HuSkin <-> SkinPiX matching

hcmp <- unique(huskin[c("compound_id", "compound_name", "smiles")])
hcmp$smiles_key <- clean_text(hcmp$smiles)
hcmp$name_key <- name_key(hcmp$compound_name)

smiles_n <- table(hcmp$smiles_key)
name_n <- table(hcmp$name_key)

smiles_map <- hcmp[!is.na(hcmp$smiles_key) & smiles_n[hcmp$smiles_key] == 1, c("smiles_key", "compound_id", "compound_name", "smiles")]
name_map <- hcmp[!is.na(hcmp$name_key) & name_n[hcmp$name_key] == 1, c("name_key", "compound_id", "compound_name", "smiles")]

smiles_lookup <- setNames(smiles_map$compound_id, smiles_map$smiles_key)
name_lookup <- setNames(name_map$compound_id, name_map$name_key)

dat$huskin_compound_id <- unname(smiles_lookup[dat$smiles_key])
dat$match_method <- ifelse(!is.na(dat$huskin_compound_id), "SMILES", NA)

i <- is.na(dat$huskin_compound_id) & !is.na(dat$name_key)
dat$huskin_compound_id[i] <- unname(name_lookup[dat$name_key[i]])
dat$match_method[i & !is.na(dat$huskin_compound_id)] <- "Name"

dat$overlap_huskin <- !is.na(dat$huskin_compound_id)

huskin_doi <- unique(doi_key(huskin$doi))
dat$doi_normalized <- doi_key(dat$doi)
dat$doi_overlaps_huskin <- !is.na(dat$doi_normalized) & dat$doi_normalized %in% huskin_doi

# RDKit descriptors

mw_lookup <- setNames(skinpix_desc$mw_rdkit, skinpix_desc$smiles)
logp_lookup <- setNames(skinpix_desc$logp_rdkit, skinpix_desc$smiles)

dat$mw_rdkit <- unname(mw_lookup[dat$SMILES])
dat$logp_rdkit <- unname(logp_lookup[dat$SMILES])

hmw_lookup <- setNames(huskin_desc$mw_rdkit, huskin_desc$smiles)
hlogp_lookup <- setNames(huskin_desc$logp_rdkit, huskin_desc$smiles)

huskin$mw_rdkit <- unname(hmw_lookup[huskin$smiles])
huskin$logp_rdkit <- unname(hlogp_lookup[huskin$smiles])

write.csv(huskin, "data/processed/40_huskindb_rdkit.csv", row.names = FALSE, na = "")

# SkinPiX processed dataset

out <- data.frame(
	source_row = dat$source_row,
	compound_id = dat$compound_id,
	compound_name = dat[["compound name"]],
	smiles = dat$SMILES,
	cas = dat[["CAS number"]],
	reference_id = dat$reference_id,
	publication_id = dat[["publication ID"]],
	author = dat$author,
	doi = dat$doi,
	doi_normalized = dat$doi_normalized,
	doi_overlaps_huskin = dat$doi_overlaps_huskin,
	huskin_compound_id = dat$huskin_compound_id,
	match_method = dat$match_method,
	overlap_huskin = dat$overlap_huskin,
	kp_cm_h = dat[["Kp (cm/h) (processed)"]],
	log_kp = log10(dat[["Kp (cm/h) (processed)"]]),
	mw_rdkit = dat$mw_rdkit,
	logp_rdkit = dat$logp_rdkit,
	temperature_c = dat[["donor/skin surface temperature (°C)"]],
	skin_donor_type = dat[["skin donor type"]],
	skin_source_site = dat[["skin source site"]],
	skin_preparation = dat[["skin preparation"]],
	used_layer = dat[["used layer"]],
	donor_ph = dat[["donor pH"]],
	donor_type = dat[["donor type"]],
	donor_category = dat[["category donor type"]],
	acceptor_temperature_c = dat[["acceptor temperature (°C)"]],
	acceptor_ph = dat[["acceptor pH"]],
	acceptor_type = dat[["acceptor type"]],
	acceptor_category = dat[["category acceptor type"]],
	cell_type = dat[["cell type"]],
	stringsAsFactors = FALSE
)

if (nrow(out) > n_raw) stop("Preparation increased row count.")

write.csv(out, "data/processed/40_skinpix.csv", row.names = FALSE, na = "")

# Reference-level overlap summary

refs <- unique(out$reference_id)

reference_overlap <- do.call(rbind, lapply(refs, function(ref) {
	x <- out[out$reference_id == ref, ]

	data.frame(
		reference_id = ref,
		publication_id = x$publication_id[1],
		doi = x$doi[1],
		n_rows = nrow(x),
		n_compounds = length(unique(x$compound_id)),
		n_huskin_anchors = length(unique(x$huskin_compound_id[!is.na(x$huskin_compound_id)])),
		n_skinpix_only = length(unique(x$compound_id[!x$overlap_huskin])),
		doi_overlaps_huskin = any(x$doi_overlaps_huskin)
	)
}))

write.csv(reference_overlap, "R/results/40_skinpix_reference_overlap.csv", row.names = FALSE)

# Compound identity audit

audit <- unique(out[out$overlap_huskin, c("compound_id", "compound_name", "smiles", "huskin_compound_id", "match_method")])

huskin_names <- setNames(hcmp$compound_name, hcmp$compound_id)
huskin_smiles_lookup <- setNames(hcmp$smiles, hcmp$compound_id)

audit$huskin_compound_name <- unname(huskin_names[audit$huskin_compound_id])
audit$huskin_smiles <- unname(huskin_smiles_lookup[audit$huskin_compound_id])

audit <- audit[c(
	"compound_id",
	"compound_name",
	"smiles",
	"huskin_compound_id",
	"huskin_compound_name",
	"huskin_smiles",
	"match_method"
)]

write.csv(audit, "R/results/40_skinpix_huskin_overlap_audit.csv", row.names = FALSE)

# Descriptor compatibility audit

sp_desc <- unique(out[out$overlap_huskin, c("compound_id", "compound_name", "huskin_compound_id", "mw_rdkit", "logp_rdkit")])

hu_desc <- aggregate(
	cbind(mw_huskin = huskin$mw, logp_huskin = huskin$log_kow, mw_rdkit_huskin = huskin$mw_rdkit, logp_rdkit_huskin = huskin$logp_rdkit),
	list(huskin_compound_id = huskin$compound_id),
	median,
	na.rm = TRUE
)

descriptor_audit <- merge(sp_desc, hu_desc, by = "huskin_compound_id", all.x = TRUE, sort = FALSE)

descriptor_audit$delta_mw_rdkit <- descriptor_audit$mw_rdkit - descriptor_audit$mw_rdkit_huskin
descriptor_audit$delta_logp_rdkit <- descriptor_audit$logp_rdkit - descriptor_audit$logp_rdkit_huskin
descriptor_audit$delta_huskin_mw_vs_rdkit <- descriptor_audit$mw_huskin - descriptor_audit$mw_rdkit_huskin
descriptor_audit$delta_huskin_logp_vs_rdkit <- descriptor_audit$logp_huskin - descriptor_audit$logp_rdkit_huskin

write.csv(descriptor_audit, "R/results/40_descriptor_compatibility_audit.csv", row.names = FALSE)

# Summary

n_compounds <- length(unique(out$compound_id))
n_refs <- length(unique(out$reference_id))
n_overlap <- length(unique(out$compound_id[out$overlap_huskin]))
n_novel <- length(unique(out$compound_id[!out$overlap_huskin]))

cat(
	"\nSkinPiX dataset prepared\n",
	"Raw rows: ", n_raw, "\n",
	"Non-equality Kp rows: ", sum(non_equal), "\n",
	"Missing/non-positive Kp rows: ", sum(bad_kp), "\n",
	"Total rows removed: ", sum(!keep), "\n",
	"Retained rows: ", nrow(out), "\n",
	"Compounds: ", n_compounds, "\n",
	"References: ", n_refs, "\n",
	"HuSkin-overlapping compounds: ", n_overlap, "\n",
	"SkinPiX-only compounds: ", n_novel, "\n",
	"SMILES matches: ", length(unique(out$compound_id[out$match_method == "SMILES"])), "\n",
	"Name-fallback matches: ", length(unique(out$compound_id[out$match_method == "Name"])), "\n",
	"References with >=1 HuSkin anchor: ", sum(reference_overlap$n_huskin_anchors >= 1), "\n",
	"References with >=3 HuSkin anchors: ", sum(reference_overlap$n_huskin_anchors >= 3), "\n",
	"References with >=5 HuSkin anchors: ", sum(reference_overlap$n_huskin_anchors >= 5), "\n",
	"References with anchors and SkinPiX-only compounds: ", sum(reference_overlap$n_huskin_anchors >= 1 & reference_overlap$n_skinpix_only >= 1), "\n",
	"References with DOI overlap with HuSkin: ", sum(reference_overlap$doi_overlaps_huskin), "\n",
	"SkinPiX missing RDKit MW: ", sum(is.na(out$mw_rdkit)), "\n",
	"SkinPiX missing RDKit logP: ", sum(is.na(out$logp_rdkit)), "\n",
	"HuSkin missing RDKit MW: ", sum(is.na(huskin$mw_rdkit)), "\n",
	"HuSkin missing RDKit logP: ", sum(is.na(huskin$logp_rdkit)), "\n",
	"Median |MW difference| between matched structures: ", round(median(abs(descriptor_audit$delta_mw_rdkit), na.rm = TRUE), 3), "\n",
	"Median |logP difference| between matched structures: ", round(median(abs(descriptor_audit$delta_logp_rdkit), na.rm = TRUE), 3), "\n",
	sep = ""
)