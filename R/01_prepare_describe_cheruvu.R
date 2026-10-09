# ============================================================
# 01_prepare_describe_cheruvu.R
# ============================================================

library(readxl)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

dat <- as.data.frame(read_excel("data/raw/Table 1.xlsx"), stringsAsFactors = FALSE)
n_raw <- nrow(dat)
dat$source_row <- seq_len(nrow(dat))

clean_text <- function(x) {
	x <- trimws(as.character(x))
	x[x %in% c("", "NA", "N/A", "N.A.", "-", "--")] <- NA
	x
}

dat$Compound <- clean_text(dat$Compound)
dat$`CAS No` <- clean_text(dat$`CAS No`)
dat$Reference <- clean_text(dat$Reference)

num <- c("MWa", "logKowb", "Mptc", "LogSaqd", "LogSoce", "Hdf", "Hag", "MVh", "Texpi", "Skin thicknessj", "logkpl", "logJmaxm")
for (x in num) dat[[x]] <- suppressWarnings(as.numeric(dat[[x]]))

key <- tolower(dat$Compound)
key <- iconv(key, to = "ASCII//TRANSLIT")
key <- gsub("[^a-z0-9]+", "_", key)
key <- gsub("^_|_$", "", key)

key[key %in% c("beta_estradiol", "estradiol")] <- "estradiol"
key[key %in% c("resorcinol", "resorcinol_1_3_benzenediol")] <- "resorcinol"
key[key %in% c("decanol", "n_decanol")] <- "decanol"
key[key %in% c("hexanol", "n_hexanol")] <- "hexanol"
key[key %in% c("heptanol", "n_heptanol")] <- "heptanol"
key[key %in% c("octanol", "n_octanol")] <- "octanol"
key[key %in% c("nonanol", "n_nonanol")] <- "nonanol"
key[key %in% c("propanol", "n_propanol")] <- "propanol"
key[key %in% c("butanol", "n_butanol")] <- "butanol"
key[key %in% c("pentanol", "n_pentanol")] <- "pentanol"
key[grepl("^5_fluorouracil", key)] <- "5_fluorouracil"

keep <- !is.na(key) & key != ""
dat <- dat[keep, ]
key <- key[keep]

dat$compound_id <- paste0("CMP_", match(key, unique(key)))
dat$compound_name <- ave(dat$Compound, dat$compound_id, FUN = function(x) x[which.max(nchar(x))])
dat$reference_id <- ifelse(is.na(dat$Reference), NA, paste0("REF_", match(dat$Reference, unique(dat$Reference[!is.na(dat$Reference)]))))

temperature_c <- ifelse(dat$Texpi > 100, dat$Texpi - 273.15, dat$Texpi)

out <- data.frame(
	source_row = dat$source_row,
	compound_id = dat$compound_id,
	compound_name = dat$compound_name,
	compound_raw = dat$Compound,
	cas = dat$`CAS No`,
	set = dat$set,
	reference_id = dat$reference_id,
	reference_raw = dat$Reference,
	mw = dat$MWa,
	log_kow = dat$logKowb,
	melting_point_k = dat$Mptc,
	log_saq = dat$LogSaqd,
	log_soc = dat$LogSoce,
	h_donor = dat$Hdf,
	h_acceptor = dat$Hag,
	molecular_volume = dat$MVh,
	temperature_c = temperature_c,
	skin_thickness_mm = dat$`Skin thicknessj`,
	skin_integrity_test = dat$`Skin Integrity testk`,
	skin_selection_criteria = dat$`Skin Selection criteria`,
	log_kp = dat$logkpl,
	log_jmax = dat$logJmaxm,
	stringsAsFactors = FALSE
)

out <- out[!is.na(out$compound_id) & !is.na(out$log_kp), ]

if (nrow(out) > n_raw) stop("Preparation increased row count.")

write.csv(out, "data/processed/01_cheruvu.csv", row.names = FALSE, na = "")

n_compounds <- length(unique(out$compound_id))
n_references <- length(unique(na.omit(out$reference_id)))
n_repeated <- sum(table(out$compound_id) > 1)
n_multiref <- sum(tapply(out$reference_id, out$compound_id, function(x) length(unique(na.omit(x)))) > 1)
within_range <- tapply(out$log_kp, out$compound_id, function(x) diff(range(x, na.rm = TRUE)))
within_range <- within_range[is.finite(within_range) & table(out$compound_id)[names(within_range)] > 1]
n_blank_compound <- n_raw - sum(keep)

cat(
	"\nCheruvu dataset prepared\n",
	"Raw rows: ", n_raw, "\n",
	"Rows removed for missing compound: ", n_blank_compound, "\n",
	"Retained rows: ", nrow(out), "\n",
	"Compounds: ", n_compounds, "\n",
	"References: ", n_references, "\n",
	"Repeated compounds: ", n_repeated, "\n",
	"Multi-reference compounds: ", n_multiref, "\n",
	"Median within-compound logKp range: ", round(median(within_range, na.rm = TRUE), 3), "\n",
	"Missing reference: ", sum(is.na(out$reference_id)), "\n",
	"Missing temperature: ", sum(is.na(out$temperature_c)), "\n",
	"Missing skin thickness: ", sum(is.na(out$skin_thickness_mm)), "\n",
	"Missing S1 descriptors: ", sum(!complete.cases(out[c("mw", "log_kow")])), "\n",
	"Missing S2 descriptors: ", sum(!complete.cases(out[c("mw", "log_kow", "melting_point_k", "log_saq", "h_donor", "h_acceptor")])), "\n",
	"Temperature range (C): ", paste(round(range(out$temperature_c, na.rm = TRUE), 1), collapse = " to "), "\n",
	sep = ""
)