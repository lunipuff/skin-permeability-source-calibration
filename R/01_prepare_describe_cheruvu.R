# 01_prepare_describe_cheruvu.R

library(readxl)

dat <- as.data.frame(read_excel("data/raw/cheruvu.xlsx"), stringsAsFactors = FALSE)
n_raw <- nrow(dat)
dat$source_row <- seq_len(n_raw)

clean_text <- function(x) {
	x <- trimws(as.character(x))
	x[x %in% c("", "NA", "N/A", "N.A.", "-", "--")] <- NA
	x
}

dat$Compound <- clean_text(dat$Compound)
dat$Reference <- clean_text(dat$Reference)

num <- c("MWa", "logKowb", "Mptc", "LogSaqd", "Hdf", "Hag", "Texpi", "Skin thicknessj", "logkpl")
for (x in num) dat[[x]] <- suppressWarnings(as.numeric(dat[[x]]))

key <- iconv(tolower(dat$Compound), to = "ASCII//TRANSLIT")
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
dat$reference_id <- ifelse(is.na(dat$Reference), NA, paste0("REF_", match(dat$Reference, unique(dat$Reference[!is.na(dat$Reference)]))))

out <- data.frame(
	source_row = dat$source_row,
	compound_id = dat$compound_id,
	reference_id = dat$reference_id,
	mw = dat$MWa,
	log_kow = dat$logKowb,
	melting_point_k = dat$Mptc,
	log_saq = dat$LogSaqd,
	h_donor = dat$Hdf,
	h_acceptor = dat$Hag,
	temperature_c = ifelse(dat$Texpi > 100, dat$Texpi - 273.15, dat$Texpi),
	skin_thickness_mm = dat$`Skin thicknessj`,
	skin_integrity_test = dat$`Skin Integrity testk`,
	log_kp = dat$logkpl,
	stringsAsFactors = FALSE
)
out <- out[!is.na(out$compound_id) & !is.na(out$log_kp), ]
if (nrow(out) > n_raw) stop("Preparation increased row count.")

write.csv(out, "data/processed/01_cheruvu.csv", row.names = FALSE, na = "")

cat(
	"\nCheruvu dataset prepared\n",
	"Raw rows: ", n_raw, "\n",
	"Rows removed for missing compound: ", sum(!keep), "\n",
	"Rows removed for missing permeability: ", sum(is.na(dat$logkpl)), "\n",
	"Retained rows: ", nrow(out), "\n",
	"Compounds: ", length(unique(out$compound_id)), "\n",
	"References: ", length(unique(na.omit(out$reference_id))), "\n",
	"Missing reference: ", sum(is.na(out$reference_id)), "\n",
	"Missing temperature: ", sum(is.na(out$temperature_c)), "\n",
	"Missing skin thickness: ", sum(is.na(out$skin_thickness_mm)), "\n",
	"Missing S1 descriptors: ", sum(!complete.cases(out[c("mw", "log_kow")])), "\n",
	"Missing S2 descriptors: ", sum(!complete.cases(out[c("mw", "log_kow", "melting_point_k", "log_saq", "h_donor", "h_acceptor")])), "\n",
	sep = ""
)