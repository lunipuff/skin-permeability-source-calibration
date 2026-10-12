# ============================================================
# 99a_make_tables.R
# ============================================================

fmt <- function(x, digits = 1) {
	formatC(x, format = "f", digits = digits)
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
domain31 <- read.csv("R/results/31_applicability_domain.csv", stringsAsFactors = FALSE)

ref41 <- read.csv("R/results/41_skinpix_calibration_reference.csv", stringsAsFactors = FALSE)
match41 <- read.csv("R/results/41_skinpix_calibration_matched_3anchor.csv", stringsAsFactors = FALSE)
domain41 <- read.csv("R/results/41_skinpix_transfer_domain.csv", stringsAsFactors = FALSE)

wrongref50 <- read.csv("R/results/50_wrong_source_reference.csv", stringsAsFactors = FALSE)

loco60 <- read.csv("R/results/60_loco_bootstrap_ci.csv", stringsAsFactors = FALSE)
var60 <- read.csv("R/results/60_variance_component_ci.csv", stringsAsFactors = FALSE)
cal60 <- read.csv("R/results/60_huskin_calibration_ci.csv", stringsAsFactors = FALSE)
wrong60 <- read.csv("R/results/60_wrong_source_ci.csv", stringsAsFactors = FALSE)

# ============================================================
# Main tables
# ============================================================

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
	`Repeated compounds` = c(ch["repeated"], hu["repeated"], 65),
	`Multi-source compounds` = c(ch["multisource"], hu["multisource"], 14),
	Role = c(
		"Development",
		"Replication and held-out-source calibration",
		"Cross-dataset transfer"
	),
	check.names = FALSE
)

write.csv(
	table1,
	"R/tables/Table_1_dataset_structure.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Table 2. Variance decomposition
# ------------------------------------------------------------

ch_ci <- var60[var60$dataset == "Cheruvu", ]
hu_ci <- var60[var60$dataset == "HuSkinDB", ]

table2 <- data.frame(
	Dataset = c("Cheruvu", "Cheruvu", "Cheruvu", "HuSkinDB", "HuSkinDB", "HuSkinDB"),
	`Variance component` = c("Compound", "Source", "Residual", "Compound", "Source", "Residual"),
	`Estimate (%)` = fmt(c(
		ch_ci$estimate_percent[match(c("Compound", "Source", "Residual"), ch_ci$component)],
		hu_ci$estimate_percent[match(c("Compound", "Source", "Residual"), hu_ci$component)]
	), 1),
	`95% CI` = c(
		fmt_ci(ch_ci$ci_low[match(c("Compound", "Source", "Residual"), ch_ci$component)], ch_ci$ci_high[match(c("Compound", "Source", "Residual"), ch_ci$component)], 1),
		fmt_ci(hu_ci$ci_low[match(c("Compound", "Source", "Residual"), hu_ci$component)], hu_ci$ci_high[match(c("Compound", "Source", "Residual"), hu_ci$component)], 1)
	),
	check.names = FALSE
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
	hierarchical_known_source = "Hierarchical, source aware"
)

table3 <- loco60

table3$Dataset <- table3$dataset
table3$Specification <- table3$specification
table3$Approach <- unname(approach_names[table3$approach])
table3$`RMSE (95% CI)` <- paste(fmt(table3$RMSE, 3), fmt_ci(table3$RMSE_ci_low, table3$RMSE_ci_high, 3))

table3$`RMSE improvement vs conventional (%)` <- ifelse(
	table3$approach == "conventional",
	"Reference",
	paste(fmt(table3$improvement_percent, 1), fmt_ci(table3$improvement_ci_low, table3$improvement_ci_high, 1))
)

table3 <- table3[, c("Dataset", "Specification", "Approach", "RMSE (95% CI)", "RMSE improvement vs conventional (%)")]

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
	RMSE_improvement_95_CI = fmt_ci(cal60$improvement_ci_low, cal60$improvement_ci_high, 1),
	check.names = FALSE
)

write.csv(
	table4a,
	"R/tables/Table_4A_huskin_matched_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Table 4B. Pairwise double-holdout falsification
# ------------------------------------------------------------

table4b <- data.frame(
	anchors = wrong60$anchors,
	n_references = wrong60$n_references,
	matched_RMSE_improvement = fmt(wrong60$correct_improvement_percent, 1),
	matched_RMSE_improvement_95_CI = fmt_ci(wrong60$correct_ci_low, wrong60$correct_ci_high, 1),
	different_source_RMSE_improvement = fmt(wrong60$wrong_improvement_percent, 1),
	different_source_RMSE_improvement_95_CI = fmt_ci(wrong60$wrong_ci_low, wrong60$wrong_ci_high, 1),
	RMSE_wrong_minus_matched = fmt(wrong60$RMSE_wrong_minus_correct, 3),
	RMSE_difference_95_CI = fmt_ci(wrong60$difference_ci_low, wrong60$difference_ci_high, 3),
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

all_one$RMSE_improvement_percent <- 100 * (1 - all_one$RMSE_calibrated / all_one$RMSE_uncalibrated)

matched41 <- match41[, c("anchors", "n_references", "RMSE_uncalibrated", "RMSE_calibrated", "RMSE_improvement_percent")]
matched41$cohort <- "Matched three-source cohort"

table4c <- rbind(
	all_one[, c("cohort", "anchors", "n_references", "RMSE_uncalibrated", "RMSE_calibrated", "RMSE_improvement_percent")],
	matched41[, c("cohort", "anchors", "n_references", "RMSE_uncalibrated", "RMSE_calibrated", "RMSE_improvement_percent")]
)

table4c$RMSE_uncalibrated <- fmt(table4c$RMSE_uncalibrated, 3)
table4c$RMSE_calibrated <- fmt(table4c$RMSE_calibrated, 3)
table4c$RMSE_improvement_percent <- fmt(table4c$RMSE_improvement_percent, 1)

write.csv(
	table4c,
	"R/tables/Table_4C_skinpix_transfer.csv",
	row.names = FALSE
)

# ============================================================
# Supplementary tables
# ============================================================

# ------------------------------------------------------------
# Supplementary Table S1. Cheruvu hierarchical model comparison
# ------------------------------------------------------------

s1 <- m10[, c("model", "AIC", "BIC", "logLik", "df")]

names(s1) <- c(
	"Model",
	"AIC",
	"BIC",
	"Log-likelihood",
	"df"
)

s1$AIC <- fmt(s1$AIC, 2)
s1$BIC <- fmt(s1$BIC, 2)
s1$`Log-likelihood` <- fmt(s1$`Log-likelihood`, 2)

write.csv(
	s1,
	"R/tables/Table_S1_cheruvu_hierarchical_models.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S2A. QSPR model summary
# ------------------------------------------------------------

s2a <- m20[, c(
	"specification",
	"n_compounds",
	"n_observations",
	"compound_variance_baseline",
	"compound_variance_joint",
	"compound_variance_explained",
	"AIC_baseline",
	"AIC_joint"
)]

names(s2a) <- c(
	"Specification",
	"Compounds",
	"Measurements",
	"Baseline compound variance",
	"Joint-model compound variance",
	"Compound variance reduction (%)",
	"Baseline AIC",
	"Joint-model AIC"
)

s2a$`Baseline compound variance` <- fmt(s2a$`Baseline compound variance`, 3)
s2a$`Joint-model compound variance` <- fmt(s2a$`Joint-model compound variance`, 3)
s2a$`Compound variance reduction (%)` <- fmt(s2a$`Compound variance reduction (%)`, 1)
s2a$`Baseline AIC` <- fmt(s2a$`Baseline AIC`, 2)
s2a$`Joint-model AIC` <- fmt(s2a$`Joint-model AIC`, 2)

write.csv(
	s2a,
	"R/tables/Table_S2A_qspr_model_summary.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S2B. QSPR coefficients
# ------------------------------------------------------------

s2b <- coef20

names(s2b)[names(s2b) == "specification"] <- "Specification"
names(s2b)[names(s2b) == "approach"] <- "Approach"
names(s2b)[names(s2b) == "term"] <- "Term"
names(s2b)[names(s2b) == "estimate"] <- "Estimate"
names(s2b)[names(s2b) == "std_error"] <- "SE"
names(s2b)[names(s2b) == "standard_error"] <- "SE"

if ("Estimate" %in% names(s2b)) s2b$Estimate <- fmt(s2b$Estimate, 3)
if ("SE" %in% names(s2b)) s2b$SE <- fmt(s2b$SE, 3)

write.csv(
	s2b,
	"R/tables/Table_S2B_qspr_coefficients.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S3A. Cheruvu LOCO performance
# ------------------------------------------------------------

s3a <- perf21[, c("specification", "approach", "RMSE", "MAE", "R2", "r")]

names(s3a) <- c(
	"Specification",
	"Approach",
	"RMSE",
	"MAE",
	"R2",
	"r"
)

for (v in c("RMSE", "MAE", "R2", "r")) {
	s3a[[v]] <- fmt(s3a[[v]], 3)
}

write.csv(
	s3a,
	"R/tables/Table_S3A_cheruvu_loco.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S3B. Cheruvu applicability domain
# ------------------------------------------------------------

s3b <- ad21[, c("specification", "in_domain", "n_compounds", "MAE_conventional", "MAE_hierarchical")]

names(s3b) <- c(
	"Specification",
	"In domain",
	"Compounds",
	"MAE conventional",
	"MAE hierarchical"
)

s3b$`In domain` <- ifelse(s3b$`In domain`, "Yes", "No")
s3b$`MAE conventional` <- fmt(s3b$`MAE conventional`, 3)
s3b$`MAE hierarchical` <- fmt(s3b$`MAE hierarchical`, 3)

s3b <- s3b[order(s3b$Specification, s3b$`In domain` == "No"), ]

write.csv(
	s3b,
	"R/tables/Table_S3B_cheruvu_applicability.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S3C. Cheruvu sensitivity analyses
# ------------------------------------------------------------

s3c <- sens21[, c("scenario", "specification", "approach", "RMSE", "MAE", "R2", "r")]

names(s3c) <- c(
	"Scenario",
	"Specification",
	"Approach",
	"RMSE",
	"MAE",
	"R2",
	"r"
)

for (v in c("RMSE", "MAE", "R2", "r")) {
	s3c[[v]] <- fmt(s3c[[v]], 3)
}

write.csv(
	s3c,
	"R/tables/Table_S3C_cheruvu_sensitivity.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S4. HuSkinDB replication summary
# ------------------------------------------------------------

s4 <- data.frame(
	Metric = c(
		"Measurements",
		"Compounds",
		"Literature sources",
		"AIC, compound-only hierarchy",
		"AIC, compound + source hierarchy",
		"Compound variance (%)",
		"Source variance (%)",
		"Residual variance (%)",
		"Compound variance reduction by S1 (%)",
		"Delta AIC after adding S1",
		"Temperature complete-case measurements",
		"Temperature coefficient [log10(cm/h) per C]",
		"Temperature p-value",
		"Singular fit"
	),
	Value = c(
		as.character(rep31$n_observations),
		as.character(rep31$n_compounds),
		as.character(rep31$n_references),
		fmt(rep31$AIC_compound_only, 2),
		fmt(rep31$AIC_compound_reference, 2),
		fmt(rep31$compound_variance_percent, 1),
		fmt(rep31$reference_variance_percent, 1),
		fmt(rep31$residual_variance_percent, 1),
		fmt(rep31$qspr_compound_variance_explained, 1),
		fmt(rep31$qspr_delta_AIC, 2),
		as.character(rep31$temperature_n),
		fmt(rep31$temperature_estimate, 3),
		format.pval(rep31$temperature_p, digits = 3, eps = 0.001),
		ifelse(rep31$singular, "Yes", "No")
	),
	check.names = FALSE
)

write.csv(
	s4,
	"R/tables/Table_S4_huskin_replication.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S5. HuSkin calibration, all eligible sources
# ------------------------------------------------------------

s5 <- cal31[, c(
	"task",
	"anchors",
	"n_references",
	"RMSE_uncalibrated",
	"RMSE_calibrated",
	"RMSE_improvement_percent",
	"MAE_uncalibrated",
	"MAE_calibrated",
	"MAE_improvement_percent"
)]

names(s5) <- c(
	"Task",
	"Anchors",
	"Sources",
	"RMSE uncal.",
	"RMSE cal.",
	"RMSE improvement (%)",
	"MAE uncal.",
	"MAE cal.",
	"MAE improvement (%)"
)

s5$`RMSE uncal.` <- fmt(s5$`RMSE uncal.`, 3)
s5$`RMSE cal.` <- fmt(s5$`RMSE cal.`, 3)
s5$`RMSE improvement (%)` <- fmt(s5$`RMSE improvement (%)`, 1)
s5$`MAE uncal.` <- fmt(s5$`MAE uncal.`, 3)
s5$`MAE cal.` <- fmt(s5$`MAE cal.`, 3)
s5$`MAE improvement (%)` <- fmt(s5$`MAE improvement (%)`, 1)

write.csv(
	s5,
	"R/tables/Table_S5_huskin_all_eligible_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S6. Matched HuSkin calibration by source
# ------------------------------------------------------------

s6 <- matchref31[, c(
	"reference_id",
	"anchors",
	"RMSE_uncalibrated",
	"RMSE_calibrated",
	"RMSE_improvement_percent",
	"MAE_improvement_percent"
)]

names(s6) <- c(
	"Source",
	"Anchors",
	"RMSE uncal.",
	"RMSE cal.",
	"RMSE improvement (%)",
	"MAE improvement (%)"
)

s6$`RMSE uncal.` <- fmt(s6$`RMSE uncal.`, 3)
s6$`RMSE cal.` <- fmt(s6$`RMSE cal.`, 3)
s6$`RMSE improvement (%)` <- fmt(s6$`RMSE improvement (%)`, 1)
s6$`MAE improvement (%)` <- fmt(s6$`MAE improvement (%)`, 1)

s6 <- s6[order(s6$Anchors, s6$Source), ]

write.csv(
	s6,
	"R/tables/Table_S6_huskin_reference_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S7. Applicability-domain summary
# ------------------------------------------------------------

ch_domain <- do.call(rbind, lapply(c("S1", "S2"), function(spec) {
	x <- ad21[ad21$specification == spec, ]

	data.frame(
		Dataset = "Cheruvu",
		Specification = spec,
		`In-domain n` = sum(x$n_compounds[x$in_domain]),
		`Out-of-domain n` = sum(x$n_compounds[!x$in_domain]),
		`MAE in domain` = x$MAE_hierarchical[x$in_domain],
		`MAE out of domain` = x$MAE_hierarchical[!x$in_domain],
		check.names = FALSE
	)
}))

hu_domain <- data.frame(
	Dataset = "HuSkinDB",
	Specification = "S1",
	`In-domain n` = sum(domain31$in_domain),
	`Out-of-domain n` = sum(!domain31$in_domain),
	`MAE in domain` = mean(domain31$error_hierarchical[domain31$in_domain]),
	`MAE out of domain` = mean(domain31$error_hierarchical[!domain31$in_domain]),
	check.names = FALSE
)

sp_domain <- data.frame(
	Dataset = "SkinPiX external test",
	Specification = "S1",
	`In-domain n` = sum(domain41$in_domain),
	`Out-of-domain n` = sum(!domain41$in_domain),
	`MAE in domain` = NA_real_,
	`MAE out of domain` = NA_real_,
	check.names = FALSE
)

s7 <- rbind(ch_domain, hu_domain, sp_domain)

s7$`MAE in domain` <- ifelse(is.na(s7$`MAE in domain`), "--", fmt(s7$`MAE in domain`, 3))
s7$`MAE out of domain` <- ifelse(is.na(s7$`MAE out of domain`), "--", fmt(s7$`MAE out of domain`, 3))

write.csv(
	s7,
	"R/tables/Table_S7_applicability_domain.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S8A. SkinPiX transfer audit
# ------------------------------------------------------------

s8a_raw <- read.csv(
	"R/tables/Table_S8A_skinpix_transfer_audit.csv",
	stringsAsFactors = FALSE,
	check.names = FALSE
)

s8a <- s8a_raw
names(s8a)[1:2] <- c("Metric", "Value")

write.csv(
	s8a,
	"R/tables/Table_S8A_skinpix_transfer_audit_display.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S8B. SkinPiX compound-overlap display
# Full structure/SMILES audit remains machine-readable.
# ------------------------------------------------------------

s8b_raw <- read.csv(
	"R/results/40_skinpix_huskin_overlap_audit.csv",
	stringsAsFactors = FALSE,
	check.names = FALSE
)

s8b <- s8b_raw[, c(
	"huskin_compound_id",
	"huskin_compound_name",
	"compound_id",
	"compound_name",
	"match_method"
)]

names(s8b) <- c(
	"HuSkin ID",
	"HuSkin compound",
	"SkinPiX ID",
	"SkinPiX compound",
	"Match method"
)

write.csv(
	s8b,
	"R/tables/Table_S8B_skinpix_compound_overlap_display.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S9. SkinPiX calibration by source
# ------------------------------------------------------------

s9 <- ref41[, c(
	"reference_id",
	"anchors",
	"RMSE_uncalibrated",
	"RMSE_calibrated",
	"RMSE_improvement_percent",
	"MAE_improvement_percent"
)]

names(s9) <- c(
	"Source",
	"Anchors",
	"RMSE uncal.",
	"RMSE cal.",
	"RMSE improvement (%)",
	"MAE improvement (%)"
)

s9$`RMSE uncal.` <- fmt(s9$`RMSE uncal.`, 3)
s9$`RMSE cal.` <- fmt(s9$`RMSE cal.`, 3)
s9$`RMSE improvement (%)` <- fmt(s9$`RMSE improvement (%)`, 1)
s9$`MAE improvement (%)` <- fmt(s9$`MAE improvement (%)`, 1)

s9 <- s9[order(s9$Anchors, s9$Source), ]

write.csv(
	s9,
	"R/tables/Table_S9_skinpix_reference_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S10. Pairwise double-holdout falsification
# ------------------------------------------------------------

s10 <- wrongref50[, c(
	"reference_id",
	"anchors",
	"n_donors",
	"n_test",
	"RMSE_uncalibrated",
	"RMSE_correct",
	"RMSE_wrong",
	"correct_RMSE_improvement",
	"wrong_RMSE_improvement",
	"RMSE_wrong_minus_correct"
)]

names(s10) <- c(
	"Target source",
	"Anchors",
	"Donor sources",
	"Test n",
	"RMSE uncal.",
	"RMSE matched",
	"RMSE different",
	"Matched improvement (%)",
	"Different improvement (%)",
	"Different - matched RMSE"
)

for (v in c("RMSE uncal.", "RMSE matched", "RMSE different", "Different - matched RMSE")) {
	s10[[v]] <- fmt(s10[[v]], 3)
}

s10$`Matched improvement (%)` <- fmt(s10$`Matched improvement (%)`, 1)
s10$`Different improvement (%)` <- fmt(s10$`Different improvement (%)`, 1)

s10 <- s10[order(s10$Anchors, s10$`Target source`), ]

write.csv(
	s10,
	"R/tables/Table_S10_wrong_source_by_reference.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat(
	"\nManuscript tables complete\n\n",
	"Main tables\n",
	"Table 1: ", nrow(table1), " rows\n",
	"Table 2: ", nrow(table2), " rows\n",
	"Table 3: ", nrow(table3), " rows\n",
	"Table 4A: ", nrow(table4a), " rows\n",
	"Table 4B: ", nrow(table4b), " rows\n",
	"Table 4C: ", nrow(table4c), " rows\n\n",
	"Supplementary display tables\n",
	"S1: ", nrow(s1), " rows\n",
	"S2A: ", nrow(s2a), " rows\n",
	"S2B: ", nrow(s2b), " rows\n",
	"S3A: ", nrow(s3a), " rows\n",
	"S3B: ", nrow(s3b), " rows\n",
	"S3C: ", nrow(s3c), " rows\n",
	"S4: ", nrow(s4), " rows\n",
	"S5: ", nrow(s5), " rows\n",
	"S6: ", nrow(s6), " rows\n",
	"S7: ", nrow(s7), " rows\n",
	"S8A: ", nrow(s8a), " rows\n",
	"S8B display: ", nrow(s8b), " rows\n",
	"S9: ", nrow(s9), " rows\n",
	"S10: ", nrow(s10), " rows\n\n",
	"Full SkinPiX overlap audit retained in Table_S8B_skinpix_compound_overlap.csv\n",
	sep = ""
)