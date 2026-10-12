# ============================================================
# 99a_make_tables.R
# ============================================================

fmt <- function(x, digits = 1) {
	round(x, digits)
}

fmt_ci <- function(low, high, digits = 1) {
	paste0(
		"[",
		formatC(low, format = "f", digits = digits),
		", ",
		formatC(high, format = "f", digits = digits),
		"]"
	)
}

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

cheruvu <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)
huskin <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)

m10 <- read.csv("R/results/10_model_comparison.csv", stringsAsFactors = FALSE)
lrt10 <- read.csv("R/results/10_likelihood_ratio_tests.csv", stringsAsFactors = FALSE)
fix10 <- read.csv("R/results/10_fixed_effects.csv", stringsAsFactors = FALSE)
var10 <- read.csv("R/results/10_variance_components.csv", stringsAsFactors = FALSE)

m20 <- read.csv("R/results/20_model_summary.csv", stringsAsFactors = FALSE)
coef20 <- read.csv("R/results/20_coefficients.csv", stringsAsFactors = FALSE)

perf21 <- read.csv("R/results/21_loco_performance.csv", stringsAsFactors = FALSE)
ad21 <- read.csv("R/results/21_applicability_summary.csv", stringsAsFactors = FALSE)
sens21 <- read.csv("R/results/21_sensitivity_performance.csv", stringsAsFactors = FALSE)

rep31 <- read.csv("R/results/31_replication_summary.csv", stringsAsFactors = FALSE)
cal31 <- read.csv("R/results/31_calibration_summary.csv", stringsAsFactors = FALSE)
match31 <- read.csv("R/results/31_calibration_matched_summary.csv", stringsAsFactors = FALSE)
matchref31 <- read.csv("R/results/31_calibration_matched_reference.csv", stringsAsFactors = FALSE)

ref41 <- read.csv("R/results/41_skinpix_calibration_reference.csv", stringsAsFactors = FALSE)
match41 <- read.csv("R/results/41_skinpix_calibration_matched_3anchor.csv", stringsAsFactors = FALSE)

blind50 <- read.csv("R/results/50_source_blind_loco.csv", stringsAsFactors = FALSE)
wrong50 <- read.csv("R/results/50_wrong_source_calibration.csv", stringsAsFactors = FALSE)
wrongref50 <- read.csv("R/results/50_wrong_source_reference.csv", stringsAsFactors = FALSE)

loco60 <- read.csv("R/results/60_loco_bootstrap_ci.csv", stringsAsFactors = FALSE)
var60 <- read.csv("R/results/60_variance_component_ci.csv", stringsAsFactors = FALSE)
cal60 <- read.csv("R/results/60_huskin_calibration_ci.csv", stringsAsFactors = FALSE)
wrong60 <- read.csv("R/results/60_wrong_source_ci.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Table 1. Dataset structure
# ------------------------------------------------------------

dataset_summary <- function(x) {
	counts <- table(x$compound_id)

	c(
		measurements = nrow(x),
		compounds = length(unique(x$compound_id)),
		sources = length(unique(x$reference_id)),
		repeated = sum(counts > 1),
		multisource = sum(tapply(x$reference_id, x$compound_id, function(z) length(unique(z))) > 1)
	)
}

ch <- dataset_summary(cheruvu)
hu <- dataset_summary(huskin)

table1 <- data.frame(
	Dataset = c("Cheruvu", "HuSkinDB", "SkinPiX"),
	Measurements = c(ch["measurements"], hu["measurements"], 199),
	Compounds = c(ch["compounds"], hu["compounds"], 108),
	Sources = c(ch["sources"], hu["sources"], 36),
	Repeated_compounds = c(ch["repeated"], hu["repeated"], 65),
	Multi_source_compounds = c(ch["multisource"], hu["multisource"], 14),
	Role = c(
		"Development",
		"Independent replication and held-out-source calibration",
		"Cross-dataset transfer validation"
	),
	check.names = FALSE
)

write.csv(
	table1,
	"R/tables/Table_1_dataset_structure.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Table 2. Hierarchical decomposition and molecular signal
# ------------------------------------------------------------

ch_ci <- var60[var60$dataset == "Cheruvu", ]
hu_ci <- var60[var60$dataset == "HuSkinDB", ]

table2 <- data.frame(
	Dataset = c("Cheruvu", "Cheruvu", "Cheruvu", "HuSkinDB", "HuSkinDB", "HuSkinDB"),
	Component = c("Compound", "Source", "Residual", "Compound", "Source", "Residual"),
	Estimate_percent = formatC(c(
		ch_ci$estimate_percent[match(c("Compound", "Source", "Residual"), ch_ci$component)],
		hu_ci$estimate_percent[match(c("Compound", "Source", "Residual"), hu_ci$component)]
	), format = "f", digits = 1),
	CI_95 = c(
		fmt_ci(
			ch_ci$ci_low[match(c("Compound", "Source", "Residual"), ch_ci$component)],
			ch_ci$ci_high[match(c("Compound", "Source", "Residual"), ch_ci$component)]
		),
		fmt_ci(
			hu_ci$ci_low[match(c("Compound", "Source", "Residual"), hu_ci$component)],
			hu_ci$ci_high[match(c("Compound", "Source", "Residual"), hu_ci$component)]
		)
	),
	check.names = FALSE
)

names(table2) <- c(
	"Dataset",
	"Component",
	"Component (%)",
	"Component, 95% CI"
)

write.csv(
	table2,
	"R/tables/Table_2_hierarchical_decomposition.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Table 3. Source information and LOCO prediction
# ------------------------------------------------------------

approach_names <- c(
	conventional = "Conventional aggregate QSPR",
	hierarchical_source_blind = "Hierarchical, source blind",
	hierarchical_known_source = "Hierarchical, known source"
)

table3 <- loco60

table3$Dataset <- table3$dataset
table3$Specification <- table3$specification
table3$Approach <- unname(approach_names[table3$approach])
table3$`RMSE (95% CI)` <- paste(
	formatC(table3$RMSE, format = "f", digits = 3),
	fmt_ci(
		table3$RMSE_ci_low,
		table3$RMSE_ci_high,
		3
	)
)

table3$`RMSE improvement vs conventional (%)` <- ifelse(
	table3$approach == "conventional",
	"Reference",
	paste(
		formatC(table3$improvement_percent, format = "f", digits = 1),
		fmt_ci(
			table3$improvement_ci_low,
			table3$improvement_ci_high,
			1
		)
	)
)

table3 <- table3[, c(
	"Dataset",
	"Specification",
	"Approach",
	"RMSE (95% CI)",
	"RMSE improvement vs conventional (%)"
)]

write.csv(
	table3,
	"R/tables/Table_3_source_information_prediction.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Table 4A. Matched HuSkin source calibration
# ------------------------------------------------------------

table4a <- data.frame(
	anchors = cal60$anchors,
	n_references = cal60$n_references,
	RMSE_uncalibrated = fmt(cal60$RMSE_uncalibrated, 3),
	RMSE_calibrated = fmt(cal60$RMSE_calibrated, 3),
	RMSE_improvement_percent = fmt(cal60$improvement_percent, 1),
	RMSE_improvement_95_CI = fmt_ci(
		cal60$improvement_ci_low,
		cal60$improvement_ci_high,
		1
	),
	check.names = FALSE
)

write.csv(
	table4a,
	"R/tables/Table_4A_huskin_matched_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Table 4B. Wrong-source falsification
# ------------------------------------------------------------

table4b <- data.frame(
	anchors = wrong60$anchors,
	n_references = wrong60$n_references,
	correct_RMSE_improvement = fmt(wrong60$correct_improvement_percent, 1),
	correct_RMSE_improvement_95_CI = fmt_ci(
		wrong60$correct_ci_low,
		wrong60$correct_ci_high,
		1
	),
	wrong_RMSE_improvement = fmt(wrong60$wrong_improvement_percent, 1),
	wrong_RMSE_improvement_95_CI = fmt_ci(
		wrong60$wrong_ci_low,
		wrong60$wrong_ci_high,
		1
	),
	RMSE_wrong_minus_correct = fmt(wrong60$RMSE_wrong_minus_correct, 3),
	RMSE_difference_95_CI = fmt_ci(
		wrong60$difference_ci_low,
		wrong60$difference_ci_high,
		3
	),
	check.names = FALSE
)

write.csv(
	table4b,
	"R/tables/Table_4B_wrong_source_falsification.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Table 4C. SkinPiX transfer
# ------------------------------------------------------------

one41 <- ref41[ref41$anchors == 1, ]

all_one <- data.frame(
	cohort = "All eligible sources",
	anchors = 1,
	n_references = length(unique(one41$reference_id)),
	RMSE_uncalibrated = mean(one41$RMSE_uncalibrated),
	RMSE_calibrated = mean(one41$RMSE_calibrated)
)

all_one$RMSE_improvement_percent <- 100 * (
	1 - all_one$RMSE_calibrated / all_one$RMSE_uncalibrated
)

matched41 <- match41[, c(
	"anchors",
	"n_references",
	"RMSE_uncalibrated",
	"RMSE_calibrated",
	"RMSE_improvement_percent"
)]

matched41$cohort <- "Matched three-source cohort"

table4c <- rbind(
	all_one[, c(
		"cohort",
		"anchors",
		"n_references",
		"RMSE_uncalibrated",
		"RMSE_calibrated",
		"RMSE_improvement_percent"
	)],
	matched41[, c(
		"cohort",
		"anchors",
		"n_references",
		"RMSE_uncalibrated",
		"RMSE_calibrated",
		"RMSE_improvement_percent"
	)]
)

table4c$RMSE_uncalibrated <- fmt(table4c$RMSE_uncalibrated, 3)
table4c$RMSE_calibrated <- fmt(table4c$RMSE_calibrated, 3)
table4c$RMSE_improvement_percent <- fmt(table4c$RMSE_improvement_percent, 1)

write.csv(
	table4c,
	"R/tables/Table_4C_skinpix_transfer.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S1. Cheruvu hierarchical models
# ------------------------------------------------------------

write.csv(
	m10,
	"R/tables/Table_S1_cheruvu_hierarchical_models.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S2. Molecular QSPR
# ------------------------------------------------------------

write.csv(
	m20,
	"R/tables/Table_S2A_qspr_model_summary.csv",
	row.names = FALSE
)

write.csv(
	coef20,
	"R/tables/Table_S2B_qspr_coefficients.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S3. Cheruvu validation
# ------------------------------------------------------------

write.csv(
	perf21,
	"R/tables/Table_S3A_cheruvu_loco.csv",
	row.names = FALSE
)

write.csv(
	ad21,
	"R/tables/Table_S3B_cheruvu_applicability.csv",
	row.names = FALSE
)

write.csv(
	sens21,
	"R/tables/Table_S3C_cheruvu_sensitivity.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S4. HuSkin replication
# ------------------------------------------------------------

write.csv(
	rep31,
	"R/tables/Table_S4_huskin_replication.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S5. HuSkin calibration, all eligible sources
# ------------------------------------------------------------

write.csv(
	cal31,
	"R/tables/Table_S5_huskin_all_eligible_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S6. Matched HuSkin calibration by source
# ------------------------------------------------------------

write.csv(
	matchref31,
	"R/tables/Table_S6_huskin_reference_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S7. Applicability-domain summary
# ------------------------------------------------------------

domain31 <- read.csv(
	"R/results/31_applicability_domain.csv",
	stringsAsFactors = FALSE
)

domain41 <- read.csv(
	"R/results/41_skinpix_transfer_domain.csv",
	stringsAsFactors = FALSE
)

s7_ch <- ad21
s7_ch$dataset <- "Cheruvu"

s7_hu <- data.frame(
	dataset = "HuSkinDB",
	specification = "S1",
	n_compounds = nrow(domain31),
	n_in_domain = sum(domain31$in_domain),
	percent_in_domain = 100 * mean(domain31$in_domain)
)

s7_sp <- data.frame(
	dataset = "SkinPiX transfer",
	specification = "HuSkin S1 domain",
	n_compounds = nrow(domain41),
	n_in_domain = sum(domain41$in_domain),
	percent_in_domain = 100 * mean(domain41$in_domain)
)

if (all(c("specification", "n_compounds", "n_in_domain", "percent_in_domain") %in% names(s7_ch))) {
	s7_ch <- s7_ch[, c(
		"dataset",
		"specification",
		"n_compounds",
		"n_in_domain",
		"percent_in_domain"
	)]

	s7 <- rbind(
		s7_ch,
		s7_hu,
		s7_sp
	)

	write.csv(
		s7,
		"R/tables/Table_S7_applicability_domain.csv",
		row.names = FALSE
	)
}

# ------------------------------------------------------------
# Supplementary Tables S8.
# Retain outputs generated by the SkinPiX preparation workflow.
# ------------------------------------------------------------

# Table_S8A_skinpix_transfer_audit.csv
# Table_S8B_skinpix_compound_overlap.csv
#
# These tables depend on the provenance/audit outputs produced by script 40
# and are intentionally not reconstructed from downstream summaries here.

# ------------------------------------------------------------
# Supplementary Table S9. SkinPiX calibration by source
# ------------------------------------------------------------

write.csv(
	ref41,
	"R/tables/Table_S9_skinpix_reference_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S10. Wrong-source falsification by source
# ------------------------------------------------------------

write.csv(
	wrongref50,
	"R/tables/Table_S10_wrong_source_by_reference.csv",
	row.names = FALSE
)

cat(
	"\nManuscript tables complete\n\n",
	"Main tables\n",
	"Table 1: ", nrow(table1), " rows\n",
	"Table 2: ", nrow(table2), " rows\n",
	"Table 3: ", nrow(table3), " rows\n",
	"Table 4A: ", nrow(table4a), " rows\n",
	"Table 4B: ", nrow(table4b), " rows\n",
	"Table 4C: ", nrow(table4c), " rows\n\n",
	"Uncertainty added to Tables 2, 3, 4A, and 4B\n",
	"SkinPiX external transfer remains descriptive without formal CI\n",
	sep = ""
)