# ============================================================
# 99a_make_tables.R
# ============================================================

fmt <- function(x, digits = 3) {
	ifelse(is.na(x), NA, round(x, digits))
}

pct <- function(x, digits = 1) {
	ifelse(is.na(x), NA, round(x, digits))
}

# ------------------------------------------------------------
# Read processed data and existing analysis outputs
# ------------------------------------------------------------

cheruvu <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)
huskin <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)
skinpix <- read.csv("data/processed/40_skinpix.csv", stringsAsFactors = FALSE)

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
perf31 <- read.csv("R/results/31_loco_performance.csv", stringsAsFactors = FALSE)
ad31 <- read.csv("R/results/31_applicability_domain.csv", stringsAsFactors = FALSE)
cal31 <- read.csv("R/results/31_calibration_summary.csv", stringsAsFactors = FALSE)
match31 <- read.csv("R/results/31_calibration_matched_summary.csv", stringsAsFactors = FALSE)
matchref31 <- read.csv("R/results/31_calibration_matched_reference.csv", stringsAsFactors = FALSE)

ref40 <- read.csv("R/results/40_skinpix_reference_overlap.csv", stringsAsFactors = FALSE)
audit40 <- read.csv("R/results/40_skinpix_huskin_overlap_audit.csv", stringsAsFactors = FALSE)
desc40 <- read.csv("R/results/40_descriptor_compatibility_audit.csv", stringsAsFactors = FALSE)

eligible41 <- read.csv("R/results/41_skinpix_eligible_references.csv", stringsAsFactors = FALSE)
cal41 <- read.csv("R/results/41_skinpix_calibration_summary.csv", stringsAsFactors = FALSE)
ref41 <- read.csv("R/results/41_skinpix_calibration_reference.csv", stringsAsFactors = FALSE)
domain41 <- read.csv("R/results/41_skinpix_transfer_domain.csv", stringsAsFactors = FALSE)

if (file.exists("R/results/41_skinpix_calibration_matched_3anchor.csv")) {
	match41 <- read.csv("R/results/41_skinpix_calibration_matched_3anchor.csv", stringsAsFactors = FALSE)
} else {
	match41 <- NULL
}

blind50 <- read.csv("R/results/50_source_blind_loco.csv", stringsAsFactors = FALSE)
wrong50 <- read.csv("R/results/50_wrong_source_calibration.csv", stringsAsFactors = FALSE)
wrongref50 <- read.csv("R/results/50_wrong_source_reference.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Table 1. Dataset structure and study role
# ------------------------------------------------------------

refs_per_compound_cheruvu <- tapply(cheruvu$reference_id, cheruvu$compound_id, function(x) length(unique(x)))
refs_per_compound_huskin <- tapply(huskin$reference_id, huskin$compound_id, function(x) length(unique(x)))
refs_per_compound_skinpix <- tapply(skinpix$reference_id, skinpix$compound_id, function(x) length(unique(x)))

table1 <- data.frame(
	dataset = c("Cheruvu", "HuSkinDB", "SkinPiX"),
	measurements = c(nrow(cheruvu), nrow(huskin), nrow(skinpix)),
	compounds = c(length(unique(cheruvu$compound_id)), length(unique(huskin$compound_id)), length(unique(skinpix$compound_id))),
	sources = c(length(unique(cheruvu$reference_id)), length(unique(huskin$reference_id)), length(unique(skinpix$reference_id))),
	repeated_compounds = c(sum(table(cheruvu$compound_id) > 1), sum(table(huskin$compound_id) > 1), sum(table(skinpix$compound_id) > 1)),
	multi_source_compounds = c(sum(refs_per_compound_cheruvu > 1), sum(refs_per_compound_huskin > 1), sum(refs_per_compound_skinpix > 1)),
	role = c(
		"Development",
		"Independent replication and held-out-source calibration",
		"Cross-dataset transfer validation"
	),
	stringsAsFactors = FALSE
)

write.csv(table1, "R/tables/Table_1_dataset_structure.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table 2. Hierarchical decomposition and molecular signal
# ------------------------------------------------------------

table2 <- data.frame(
	dataset = c("Cheruvu", "HuSkinDB"),
	compound_variance_percent = c(
		var10$percent[var10$component == "Compound"],
		rep31$compound_variance_percent
	),
	source_variance_percent = c(
		var10$percent[var10$component == "Reference"],
		rep31$reference_variance_percent
	),
	residual_variance_percent = c(
		var10$percent[var10$component == "Residual"],
		rep31$residual_variance_percent
	),
	S1_compound_variance_explained_percent = c(
		m20$compound_variance_explained[m20$specification == "S1"],
		rep31$qspr_compound_variance_explained
	),
	S2_compound_variance_explained_percent = c(
		m20$compound_variance_explained[m20$specification == "S2"],
		NA
	),
	temperature_estimate_log10_per_C = c(
		fix10$estimate[fix10$term == "temperature_c"],
		rep31$temperature_estimate
	),
	temperature_p = c(
		lrt10$p[lrt10$comparison == "M1 vs M2"],
		rep31$temperature_p
	)
)

table2[3:6] <- lapply(table2[3:6], pct)
table2$temperature_estimate_log10_per_C <- fmt(table2$temperature_estimate_log10_per_C, 4)
table2$temperature_p <- signif(table2$temperature_p, 3)

write.csv(table2, "R/tables/Table_2_hierarchical_decomposition.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table 3. Source information and unseen-compound prediction
# ------------------------------------------------------------

table3 <- blind50[c("dataset", "specification", "approach", "RMSE", "MAE", "R2", "improvement_vs_conventional")]
table3$approach <- c(
	conventional = "Conventional aggregate QSPR",
	hierarchical_source_blind = "Hierarchical, source-blind",
	hierarchical_known_source = "Hierarchical, known source"
)[table3$approach]

table3$RMSE <- fmt(table3$RMSE)
table3$MAE <- fmt(table3$MAE)
table3$R2 <- fmt(table3$R2)
table3$improvement_vs_conventional <- pct(table3$improvement_vs_conventional)

write.csv(table3, "R/tables/Table_3_source_information_prediction.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table 4A. HuSkin matched calibration
# ------------------------------------------------------------

table4a <- match31[c(
	"anchors",
	"n_references",
	"RMSE_uncalibrated",
	"RMSE_calibrated",
	"RMSE_improvement_percent",
	"MAE_uncalibrated",
	"MAE_calibrated",
	"MAE_improvement_percent"
)]

table4a[c("RMSE_uncalibrated", "RMSE_calibrated", "MAE_uncalibrated", "MAE_calibrated")] <- lapply(
	table4a[c("RMSE_uncalibrated", "RMSE_calibrated", "MAE_uncalibrated", "MAE_calibrated")],
	fmt
)
table4a$RMSE_improvement_percent <- pct(table4a$RMSE_improvement_percent)
table4a$MAE_improvement_percent <- pct(table4a$MAE_improvement_percent)

write.csv(table4a, "R/tables/Table_4A_huskin_matched_calibration.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table 4B. Correct-source versus wrong-source calibration
# ------------------------------------------------------------

table4b <- wrong50[c(
	"anchors",
	"n_references",
	"RMSE_uncalibrated",
	"RMSE_correct",
	"RMSE_wrong",
	"correct_RMSE_improvement",
	"wrong_RMSE_improvement"
)]

table4b[c("RMSE_uncalibrated", "RMSE_correct", "RMSE_wrong")] <- lapply(
	table4b[c("RMSE_uncalibrated", "RMSE_correct", "RMSE_wrong")],
	fmt
)
table4b$correct_RMSE_improvement <- pct(table4b$correct_RMSE_improvement)
table4b$wrong_RMSE_improvement <- pct(table4b$wrong_RMSE_improvement)

write.csv(table4b, "R/tables/Table_4B_wrong_source_falsification.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table 4C. SkinPiX cross-dataset transfer
# ------------------------------------------------------------

table4c_all <- cal41[c(
	"anchors",
	"n_references",
	"RMSE_uncalibrated",
	"RMSE_calibrated",
	"RMSE_improvement_percent",
	"MAE_uncalibrated",
	"MAE_calibrated",
	"MAE_improvement_percent"
)]
table4c_all$cohort <- "All eligible references"

if (!is.null(match41)) {
	table4c_match <- match41[c(
		"anchors",
		"n_references",
		"RMSE_uncalibrated",
		"RMSE_calibrated",
		"RMSE_improvement_percent",
		"MAE_uncalibrated",
		"MAE_calibrated",
		"MAE_improvement_percent"
	)]
	table4c_match$cohort <- "Matched references supporting >=3 anchors"
	table4c <- rbind(table4c_all, table4c_match)
} else {
	table4c <- table4c_all
}

table4c <- table4c[c(
	"cohort",
	"anchors",
	"n_references",
	"RMSE_uncalibrated",
	"RMSE_calibrated",
	"RMSE_improvement_percent",
	"MAE_uncalibrated",
	"MAE_calibrated",
	"MAE_improvement_percent"
)]

table4c[c("RMSE_uncalibrated", "RMSE_calibrated", "MAE_uncalibrated", "MAE_calibrated")] <- lapply(
	table4c[c("RMSE_uncalibrated", "RMSE_calibrated", "MAE_uncalibrated", "MAE_calibrated")],
	fmt
)
table4c$RMSE_improvement_percent <- pct(table4c$RMSE_improvement_percent)
table4c$MAE_improvement_percent <- pct(table4c$MAE_improvement_percent)

write.csv(table4c, "R/tables/Table_4C_skinpix_transfer.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S1. Cheruvu hierarchical model comparison
# ------------------------------------------------------------

table_s1 <- merge(m10, cbind(model = c("M0", "M1", "M2"), lrt10), by = "model", all.x = TRUE)
table_s1$AIC <- fmt(table_s1$AIC, 2)
table_s1$BIC <- fmt(table_s1$BIC, 2)
table_s1$logLik <- fmt(table_s1$logLik, 2)
table_s1$chisq <- fmt(table_s1$chisq, 2)
table_s1$p <- signif(table_s1$p, 3)

write.csv(table_s1, "R/tables/Table_S1_cheruvu_hierarchical_models.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S2. QSPR coefficients and molecular variance explained
# ------------------------------------------------------------

summary_s2 <- m20[c(
	"specification",
	"compound_variance_explained",
	"AIC_baseline",
	"AIC_joint",
	"BIC_baseline",
	"BIC_joint",
	"LRT_chisq",
	"LRT_df",
	"LRT_p"
)]

coef20$estimate <- fmt(coef20$estimate, 4)
coef20$se <- fmt(coef20$se, 4)

write.csv(summary_s2, "R/tables/Table_S2A_qspr_model_summary.csv", row.names = FALSE)
write.csv(coef20, "R/tables/Table_S2B_qspr_coefficients.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S3. Cheruvu validation and sensitivity analyses
# ------------------------------------------------------------

write.csv(perf21, "R/tables/Table_S3A_cheruvu_loco.csv", row.names = FALSE)
write.csv(ad21, "R/tables/Table_S3B_cheruvu_applicability.csv", row.names = FALSE)
write.csv(sens21, "R/tables/Table_S3C_cheruvu_sensitivity.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S4. HuSkin replication
# ------------------------------------------------------------

table_s4 <- data.frame(
	metric = c(
		"Observations",
		"Compounds",
		"References",
		"AIC compound only",
		"AIC compound + reference",
		"Compound variance (%)",
		"Reference variance (%)",
		"Residual variance (%)",
		"MW + logKow compound variance explained (%)",
		"QSPR delta AIC",
		"Temperature subset n",
		"Temperature estimate",
		"Temperature p",
		"Singular fit"
	),
	value = c(
		rep31$n_observations,
		rep31$n_compounds,
		rep31$n_references,
		fmt(rep31$AIC_compound_only, 2),
		fmt(rep31$AIC_compound_reference, 2),
		pct(rep31$compound_variance_percent),
		pct(rep31$reference_variance_percent),
		pct(rep31$residual_variance_percent),
		pct(rep31$qspr_compound_variance_explained),
		fmt(rep31$qspr_delta_AIC, 2),
		rep31$temperature_n,
		fmt(rep31$temperature_estimate, 4),
		signif(rep31$temperature_p, 3),
		rep31$singular
	)
)

write.csv(table_s4, "R/tables/Table_S4_huskin_replication.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S5. HuSkin all-eligible calibration
# ------------------------------------------------------------

write.csv(cal31, "R/tables/Table_S5_huskin_all_eligible_calibration.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S6. HuSkin matched calibration by reference
# ------------------------------------------------------------

write.csv(matchref31, "R/tables/Table_S6_huskin_reference_calibration.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S7. Applicability-domain summary
# ------------------------------------------------------------

cheruvu_ad <- do.call(rbind, lapply(split(ad21, ad21$specification), function(x) {
	inrow <- x[x$in_domain, ]
	outrow <- x[!x$in_domain, ]

	data.frame(
		dataset = "Cheruvu",
		specification = x$specification[1],
		n_in_domain = if (nrow(inrow)) inrow$n_compounds else 0,
		n_out_domain = if (nrow(outrow)) outrow$n_compounds else 0,
		MAE_in_domain = if (nrow(inrow)) inrow$MAE_hierarchical else NA,
		MAE_out_domain = if (nrow(outrow)) outrow$MAE_hierarchical else NA
	)
}))

huskin_ad <- data.frame(
	dataset = "HuSkinDB",
	specification = "S1",
	n_in_domain = sum(ad31$in_domain),
	n_out_domain = sum(!ad31$in_domain),
	MAE_in_domain = mean(ad31$error_hierarchical[ad31$in_domain]),
	MAE_out_domain = mean(ad31$error_hierarchical[!ad31$in_domain])
)

skinpix_ad <- data.frame(
	dataset = "SkinPiX external test",
	specification = "S1",
	n_in_domain = sum(domain41$in_domain),
	n_out_domain = sum(!domain41$in_domain),
	MAE_in_domain = NA,
	MAE_out_domain = NA
)

table_s7 <- rbind(cheruvu_ad, huskin_ad, skinpix_ad)
table_s7$MAE_in_domain <- fmt(table_s7$MAE_in_domain)
table_s7$MAE_out_domain <- fmt(table_s7$MAE_out_domain)

write.csv(table_s7, "R/tables/Table_S7_applicability_domain.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S8. SkinPiX provenance, overlap and descriptor audit
# ------------------------------------------------------------

table_s8a <- data.frame(
	metric = c(
		"Retained measurements",
		"Compounds",
		"References",
		"HuSkin-overlapping compounds",
		"SkinPiX-only compounds",
		"Eligible transfer references",
		"References with DOI overlap",
		"Median absolute RDKit MW difference among matched compounds",
		"Median absolute RDKit logP difference among matched compounds"
	),
	value = c(
		nrow(skinpix),
		length(unique(skinpix$compound_id)),
		length(unique(skinpix$reference_id)),
		length(unique(skinpix$compound_id[skinpix$overlap_huskin])),
		length(unique(skinpix$compound_id[!skinpix$overlap_huskin])),
		sum(eligible41$eligible),
		sum(ref40$doi_overlaps_huskin),
		fmt(median(abs(desc40$delta_mw_rdkit), na.rm = TRUE), 3),
		fmt(median(abs(desc40$delta_logp_rdkit), na.rm = TRUE), 3)
	)
)

write.csv(table_s8a, "R/tables/Table_S8A_skinpix_transfer_audit.csv", row.names = FALSE)
write.csv(audit40, "R/tables/Table_S8B_skinpix_compound_overlap.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S9. SkinPiX calibration by source
# ------------------------------------------------------------

write.csv(ref41, "R/tables/Table_S9_skinpix_reference_calibration.csv", row.names = FALSE)

# ------------------------------------------------------------
# Table S10. Wrong-source falsification by source
# ------------------------------------------------------------

write.csv(wrongref50, "R/tables/Table_S10_wrong_source_by_reference.csv", row.names = FALSE)

cat(
	"\nManuscript tables created\n",
	"Main tables: 4 table groups\n",
	"Supplementary tables: S1-S10\n",
	"Output directory: R/tables\n",
	sep = ""
)