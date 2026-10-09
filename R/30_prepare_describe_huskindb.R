# ============================================================
# 30_prepare_describe_huskindb.R
# ============================================================

dat <- read.csv("data/raw/huskinDB.csv", stringsAsFactors = FALSE, check.names = FALSE)
n_raw <- nrow(dat)
dat$source_row <- seq_len(nrow(dat))

clean_text <- function(x) {
	x <- trimws(as.character(x))
	x[tolower(x) %in% c("", "na", "n/a", "unknown", "-", "--")] <- NA
	x
}

text_cols <- c(
	"Compound name", "Smiles", "skin source type", "skin source site", "used layer",
	"skin preparation", "donor type", "acceptor type", "cell type", "reference", "DOI"
)
for (x in text_cols) dat[[x]] <- clean_text(dat[[x]])

num_cols <- c(
	"logkp (cm/s)", "storage temperature (°C)", "storage duration (days)", "neat",
	"donor/skin surface temperature (°C)", "donor pH", "acceptor temperature (°C)",
	"acceptor pH", "LogP", "molecular weight"
)
for (x in num_cols) dat[[x]] <- suppressWarnings(as.numeric(dat[[x]]))

name_key <- tolower(dat[["Compound name"]])
name_key <- iconv(name_key, to = "ASCII//TRANSLIT")
name_key <- gsub("[^a-z0-9]+", "_", name_key)
name_key <- gsub("^_|_$", "", name_key)

compound_key <- ifelse(!is.na(dat$Smiles), paste0("SMILES:", dat$Smiles), paste0("NAME:", name_key))
keep <- !is.na(dat[["logkp (cm/s)"]]) & !is.na(compound_key) & compound_key != "NAME:"
dat <- dat[keep, ]
compound_key <- compound_key[keep]

dat$compound_id <- paste0("HUC_", match(compound_key, unique(compound_key)))
dat$compound_name <- ave(dat[["Compound name"]], dat$compound_id, FUN = function(x) x[which.max(nchar(x))])

doi_key <- tolower(dat$DOI)
doi_key <- gsub("\\s+", "", doi_key)
source_key <- ifelse(!is.na(doi_key), paste0("DOI:", doi_key), paste0("REF:", dat$reference))
if (any(is.na(source_key) | source_key == "REF:")) stop("Missing source identity.")
dat$reference_id <- paste0("HUR_", match(source_key, unique(source_key)))

out <- data.frame(
	source_row = dat$source_row,
	compound_id = dat$compound_id,
	compound_name = dat$compound_name,
	smiles = dat$Smiles,
	reference_id = dat$reference_id,
	reference_raw = dat$reference,
	doi = dat$DOI,
	log_kp_cm_s = dat[["logkp (cm/s)"]],
	log_kp = dat[["logkp (cm/s)"]] + log10(3600),
	mw = dat[["molecular weight"]],
	log_kow = dat$LogP,
	temperature_c = dat[["donor/skin surface temperature (°C)"]],
	skin_source_type = dat[["skin source type"]],
	skin_source_site = dat[["skin source site"]],
	used_layer = dat[["used layer"]],
	skin_preparation = dat[["skin preparation"]],
	neat = dat$neat,
	donor_ph = dat[["donor pH"]],
	donor_type = dat[["donor type"]],
	acceptor_temperature_c = dat[["acceptor temperature (°C)"]],
	acceptor_ph = dat[["acceptor pH"]],
	acceptor_type = dat[["acceptor type"]],
	cell_type = dat[["cell type"]],
	stringsAsFactors = FALSE
)

if (nrow(out) > n_raw) stop("Preparation increased row count.")

write.csv(out, "data/processed/30_huskindb.csv", row.names = FALSE, na = "")

n_compounds <- length(unique(out$compound_id))
n_references <- length(unique(out$reference_id))
n_repeated <- sum(table(out$compound_id) > 1)
refs_per_compound <- tapply(out$reference_id, out$compound_id, function(x) length(unique(x)))
n_multiref <- sum(refs_per_compound > 1)

compounds_per_reference <- tapply(out$compound_id, out$reference_id, function(x) length(unique(x)))
shared_per_reference <- tapply(out$compound_id, out$reference_id, function(x) sum(unique(x) %in% names(refs_per_compound)[refs_per_compound > 1]))

cat(
	"\nHuSkinDB dataset prepared\n",
	"Raw rows: ", n_raw, "\n",
	"Retained rows: ", nrow(out), "\n",
	"Compounds: ", n_compounds, "\n",
	"References: ", n_references, "\n",
	"Repeated compounds: ", n_repeated, "\n",
	"Multi-reference compounds: ", n_multiref, "\n",
	"References with >=3 shared compounds: ", sum(shared_per_reference >= 3), "\n",
	"References with >=5 shared compounds: ", sum(shared_per_reference >= 5), "\n",
	"Missing MW: ", sum(is.na(out$mw)), "\n",
	"Missing logKow: ", sum(is.na(out$log_kow)), "\n",
	"Missing donor temperature: ", sum(is.na(out$temperature_c)), "\n",
	"logKp converted from log(cm/s) to log(cm/h): +", round(log10(3600), 4), "\n",
	sep = ""
)