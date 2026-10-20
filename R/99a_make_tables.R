# ============================================================
# 99a_make_tables.R
# ============================================================

fmt <- function(x, digits) {
	formatC(x, format = "f", digits = digits)
}

fmt_ci <- function(low, high, digits) {
	paste0(
		"[",
		fmt(low, digits),
		", ",
		fmt(high, digits),
		"]"
	)
}

fmt_p <- function(x) {
	ifelse(
		is.na(x),
		NA_character_,
		ifelse(
			x < 0.001,
			"<0.001",
			format.pval(x, digits = 3, eps = 0.001)
		)
	)
}

# ------------
# Anchor Check
# ------------

anchor_check <- read.csv("R/results/31_calibration_performance.csv", stringsAsFactors = FALSE)
anchor_check <- anchor_check[anchor_check$anchors == 1, ]
anchor_check <- aggregate(cbind(RMSE_uncalibrated, RMSE_calibrated) ~ reference_id + task, data = anchor_check, FUN = mean)

shared_check <- anchor_check[anchor_check$task == "Known compound", c("reference_id", "RMSE_uncalibrated", "RMSE_calibrated")]
names(shared_check)[2:3] <- c("shared_before", "shared_after")

novel_check <- anchor_check[anchor_check$task == "New compound", c("reference_id", "RMSE_uncalibrated", "RMSE_calibrated")]
names(novel_check)[2:3] <- c("novel_before", "novel_after")

anchor_check <- merge(shared_check, novel_check, by = "reference_id")
anchor_check$shared_improved <- anchor_check$shared_after < anchor_check$shared_before
anchor_check$novel_improved <- anchor_check$novel_after < anchor_check$novel_before
anchor_check$shared_reduction_percent <- 100 * (1 - anchor_check$shared_after / anchor_check$shared_before)
anchor_check$novel_reduction_percent <- 100 * (1 - anchor_check$novel_after / anchor_check$novel_before)

write.csv(anchor_check, "R/tables/Table_anchor_check.csv", row.names = FALSE)
print(table(anchor_check$shared_improved, anchor_check$novel_improved))

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

cheruvu <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)
huskin <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)

m10 <- read.csv("R/results/10_model_comparison.csv", stringsAsFactors = FALSE)

m20 <- read.csv("R/results/20_model_summary.csv", stringsAsFactors = FALSE)
coef20 <- read.csv("R/results/20_coefficients.csv", stringsAsFactors = FALSE)

perf21 <- read.csv("R/results/21_loco_performance.csv", stringsAsFactors = FALSE)
ad21 <- read.csv("R/results/21_applicability_summary.csv", stringsAsFactors = FALSE)
sens21 <- read.csv("R/results/21_sensitivity_performance.csv", stringsAsFactors = FALSE)

rep31 <- read.csv("R/results/31_replication_summary.csv", stringsAsFactors = FALSE)
cal31 <- read.csv("R/results/31_calibration_summary.csv", stringsAsFactors = FALSE)
matchref31 <- read.csv("R/results/31_calibration_matched_reference.csv", stringsAsFactors = FALSE)

comp32 <- read.csv("R/results/32_calibration_method_summary.csv", stringsAsFactors = FALSE)

ref41 <- read.csv("R/results/41_skinpix_calibration_reference.csv", stringsAsFactors = FALSE)
match41 <- read.csv("R/results/41_skinpix_calibration_matched_3anchor.csv", stringsAsFactors = FALSE)

wrongref50 <- read.csv("R/results/50_wrong_source_reference.csv", stringsAsFactors = FALSE)

loco60 <- read.csv("R/results/60_loco_bootstrap_ci.csv", stringsAsFactors = FALSE)
var60 <- read.csv("R/results/60_variance_component_ci.csv", stringsAsFactors = FALSE)
deploy60 <- read.csv("R/results/60_huskin_one_anchor_all_eligible_ci.csv", stringsAsFactors = FALSE)
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
# Table 2. Hierarchical decomposition
# ------------------------------------------------------------

ch_ci <- var60[var60$dataset == "Cheruvu", ]
hu_ci <- var60[var60$dataset == "HuSkinDB", ]

table2 <- data.frame(
	Dataset = c(
		"Cheruvu",
		"Cheruvu",
		"Cheruvu",
		"HuSkinDB",
		"HuSkinDB",
		"HuSkinDB"
	),
	Component = c(
		"Compound",
		"Source",
		"Residual",
		"Compound",
		"Source",
		"Residual"
	),
	`Component (%)` = fmt(
		c(
			ch_ci$estimate_percent[match(c("Compound", "Source", "Residual"), ch_ci$component)],
			hu_ci$estimate_percent[match(c("Compound", "Source", "Residual"), hu_ci$component)]
		),
		1
	),
	`95% CI` = c(
		fmt_ci(
			ch_ci$ci_low[match(c("Compound", "Source", "Residual"), ch_ci$component)],
			ch_ci$ci_high[match(c("Compound", "Source", "Residual"), ch_ci$component)],
			1
		),
		fmt_ci(
			hu_ci$ci_low[match(c("Compound", "Source", "Residual"), hu_ci$component)],
			hu_ci$ci_high[match(c("Compound", "Source", "Residual"), hu_ci$component)],
			1
		)
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

table3 <- data.frame(
	Dataset = loco60$dataset,
	Specification = loco60$specification,
	`Prediction setting` = unname(approach_names[loco60$approach]),
	`RMSE (95% CI)` = paste(
		fmt(loco60$RMSE, 3),
		fmt_ci(loco60$RMSE_ci_low, loco60$RMSE_ci_high, 3)
	),
	`RMSE improvement vs conventional (%)` = ifelse(
		loco60$approach == "conventional",
		"Reference",
		paste(
			fmt(loco60$improvement_percent, 1),
			fmt_ci(loco60$improvement_ci_low, loco60$improvement_ci_high, 1)
		)
	),
	check.names = FALSE
)

write.csv(
	table3,
	"R/tables/Table_3_source_information_prediction.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Table 4A. HuSkin source calibration
# ------------------------------------------------------------

if (nrow(deploy60) != 1 || deploy60$anchors[1] != 1) {
	stop("Expected one all-eligible one-anchor deployment result")
}

matched_improved <- aggregate(
	RMSE_improvement_percent ~ anchors,
	data = matchref31,
	function(x) sum(x > 0, na.rm = TRUE)
)

names(matched_improved)[2] <- "n_references_improved"

primary4a <- data.frame(
	Cohort = "All eligible sources (primary deployment)",
	Anchors = deploy60$anchors,
	Sources = deploy60$n_references,
	`Sources improved` = paste0(deploy60$n_references_improved, "/", deploy60$n_references),
	`RMSE uncalibrated` = fmt(deploy60$RMSE_uncalibrated, 3),
	`RMSE calibrated` = fmt(deploy60$RMSE_calibrated, 3),
	`RMSE improvement (%)` = fmt(deploy60$improvement_percent, 1),
	`95% CI` = fmt_ci(deploy60$improvement_ci_low, deploy60$improvement_ci_high, 1),
	check.names = FALSE
)

matched4a <- data.frame(
	Cohort = "Matched 5-anchor-capable cohort",
	Anchors = cal60$anchors,
	Sources = cal60$n_references,
	`Sources improved` = paste0(
		matched_improved$n_references_improved[match(cal60$anchors, matched_improved$anchors)],
		"/",
		cal60$n_references
	),
	`RMSE uncalibrated` = fmt(cal60$RMSE_uncalibrated, 3),
	`RMSE calibrated` = fmt(cal60$RMSE_calibrated, 3),
	`RMSE improvement (%)` = fmt(cal60$improvement_percent, 1),
	`95% CI` = fmt_ci(cal60$improvement_ci_low, cal60$improvement_ci_high, 1),
	check.names = FALSE
)

table4a <- rbind(
	primary4a,
	matched4a
)

write.csv(
	table4a,
	"R/tables/Table_4A_huskin_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Table 4B. Wrong-source falsification
# ------------------------------------------------------------

fmt_est_ci <- function(x, low, high, digits) {
	paste0(
		fmt(x, digits),
		" ",
		fmt_ci(low, high, digits)
	)
}

table4b <- data.frame(
	Anchors = wrong60$anchors,
	`Target sources` = wrong60$n_references,
	`Target-source RMSE reduction (%)` = fmt_est_ci(
		wrong60$correct_improvement_percent,
		wrong60$correct_ci_low,
		wrong60$correct_ci_high,
		1
	),
	`Different-source RMSE reduction (%)` = fmt_est_ci(
		wrong60$wrong_improvement_percent,
		wrong60$wrong_ci_low,
		wrong60$wrong_ci_high,
		1
	),
	`RMSE difference` = fmt_est_ci(
		wrong60$RMSE_wrong_minus_correct,
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
	Cohort = "All eligible sources",
	Anchors = 1,
	Sources = length(unique(one41$reference_id)),
	`RMSE uncalibrated` = mean(one41$RMSE_uncalibrated),
	`RMSE calibrated` = mean(one41$RMSE_calibrated),
	check.names = FALSE
)

all_one$`RMSE improvement (%)` <- 100 * (
	1 -
	all_one$`RMSE calibrated` /
	all_one$`RMSE uncalibrated`
)

matched41_table <- data.frame(
	Cohort = "Matched three-source cohort",
	Anchors = match41$anchors,
	Sources = match41$n_references,
	`RMSE uncalibrated` = match41$RMSE_uncalibrated,
	`RMSE calibrated` = match41$RMSE_calibrated,
	`RMSE improvement (%)` = match41$RMSE_improvement_percent,
	check.names = FALSE
)

table4c <- rbind(
	all_one,
	matched41_table
)

table4c$`RMSE uncalibrated` <- fmt(table4c$`RMSE uncalibrated`, 3)
table4c$`RMSE calibrated` <- fmt(table4c$`RMSE calibrated`, 3)
table4c$`RMSE improvement (%)` <- fmt(table4c$`RMSE improvement (%)`, 1)

write.csv(
	table4c,
	"R/tables/Table_4C_skinpix_transfer.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S1. Cheruvu hierarchical models
# ------------------------------------------------------------

s1 <- m10

if ("AIC" %in% names(s1)) {
	s1$AIC <- fmt(s1$AIC, 2)
}

if ("BIC" %in% names(s1)) {
	s1$BIC <- fmt(s1$BIC, 2)
}

if ("logLik" %in% names(s1)) {
	s1$logLik <- fmt(s1$logLik, 2)
}

names(s1) <- gsub("_", " ", names(s1), fixed = TRUE)

write.csv(
	s1,
	"R/tables/Table_S1_cheruvu_hierarchical_models.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S2A. Molecular QSPR model summary
# ------------------------------------------------------------

s2a <- m20

if ("compound_variance_explained" %in% names(s2a)) {
	s2a$compound_variance_explained <- fmt(s2a$compound_variance_explained, 1)
}

if ("AIC" %in% names(s2a)) {
	s2a$AIC <- fmt(s2a$AIC, 2)
}

if ("delta_AIC" %in% names(s2a)) {
	s2a$delta_AIC <- fmt(s2a$delta_AIC, 2)
}

names(s2a) <- gsub("_", " ", names(s2a), fixed = TRUE)

write.csv(
	s2a,
	"R/tables/Table_S2A_qspr_model_summary.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S2B. Molecular QSPR coefficients
# ------------------------------------------------------------

s2b <- coef20

if ("estimate" %in% names(s2b)) {
	s2b$estimate <- fmt(s2b$estimate, 4)
}

if ("std_error" %in% names(s2b)) {
	s2b$std_error <- fmt(s2b$std_error, 4)
}

if ("statistic" %in% names(s2b)) {
	s2b$statistic <- fmt(s2b$statistic, 3)
}

if ("p_value" %in% names(s2b)) {
	s2b$p_value <- fmt_p(s2b$p_value)
}

names(s2b) <- gsub("_", " ", names(s2b), fixed = TRUE)

write.csv(
	s2b,
	"R/tables/Table_S2B_qspr_coefficients.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S3A. Cheruvu LOCO
# ------------------------------------------------------------

s3a <- perf21

if ("RMSE" %in% names(s3a)) {
	s3a$RMSE <- fmt(s3a$RMSE, 3)
}

if ("MAE" %in% names(s3a)) {
	s3a$MAE <- fmt(s3a$MAE, 3)
}

if ("R2" %in% names(s3a)) {
	s3a$R2 <- fmt(s3a$R2, 3)
}

if ("improvement_percent" %in% names(s3a)) {
	s3a$improvement_percent <- fmt(s3a$improvement_percent, 1)
}

names(s3a) <- gsub("_", " ", names(s3a), fixed = TRUE)

write.csv(
	s3a,
	"R/tables/Table_S3A_cheruvu_loco.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S3B. Cheruvu applicability
# ------------------------------------------------------------

s3b <- ad21

if ("RMSE" %in% names(s3b)) {
	s3b$RMSE <- fmt(s3b$RMSE, 3)
}

if ("MAE" %in% names(s3b)) {
	s3b$MAE <- fmt(s3b$MAE, 3)
}

if ("percent_in_domain" %in% names(s3b)) {
	s3b$percent_in_domain <- fmt(s3b$percent_in_domain, 1)
}

names(s3b) <- gsub("_", " ", names(s3b), fixed = TRUE)

write.csv(
	s3b,
	"R/tables/Table_S3B_cheruvu_applicability.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S3C. Cheruvu sensitivity
# ------------------------------------------------------------

s3c <- sens21

if ("RMSE" %in% names(s3c)) {
	s3c$RMSE <- fmt(s3c$RMSE, 3)
}

if ("MAE" %in% names(s3c)) {
	s3c$MAE <- fmt(s3c$MAE, 3)
}

if ("improvement_percent" %in% names(s3c)) {
	s3c$improvement_percent <- fmt(s3c$improvement_percent, 1)
}

names(s3c) <- gsub("_", " ", names(s3c), fixed = TRUE)

write.csv(
	s3c,
	"R/tables/Table_S3C_cheruvu_sensitivity.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S4. HuSkin replication
# ------------------------------------------------------------

s4 <- rep31

if ("RMSE" %in% names(s4)) {
	s4$RMSE <- fmt(s4$RMSE, 3)
}

if ("MAE" %in% names(s4)) {
	s4$MAE <- fmt(s4$MAE, 3)
}

if ("improvement_percent" %in% names(s4)) {
	s4$improvement_percent <- fmt(s4$improvement_percent, 1)
}

if ("source_variance_percent" %in% names(s4)) {
	s4$source_variance_percent <- fmt(s4$source_variance_percent, 1)
}

if ("compound_variance_percent" %in% names(s4)) {
	s4$compound_variance_percent <- fmt(s4$compound_variance_percent, 1)
}

if ("residual_variance_percent" %in% names(s4)) {
	s4$residual_variance_percent <- fmt(s4$residual_variance_percent, 1)
}

if ("p_value" %in% names(s4)) {
	s4$p_value <- fmt_p(s4$p_value)
}

names(s4) <- gsub("_", " ", names(s4), fixed = TRUE)

write.csv(
	s4,
	"R/tables/Table_S4_huskin_replication.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S5. HuSkin calibration, all eligible
# ------------------------------------------------------------

s5 <- cal31

if ("RMSE_uncalibrated" %in% names(s5)) {
	s5$RMSE_uncalibrated <- fmt(s5$RMSE_uncalibrated, 3)
}

if ("RMSE_calibrated" %in% names(s5)) {
	s5$RMSE_calibrated <- fmt(s5$RMSE_calibrated, 3)
}

if ("RMSE_improvement_percent" %in% names(s5)) {
	s5$RMSE_improvement_percent <- fmt(s5$RMSE_improvement_percent, 1)
}

if ("MAE_uncalibrated" %in% names(s5)) {
	s5$MAE_uncalibrated <- fmt(s5$MAE_uncalibrated, 3)
}

if ("MAE_calibrated" %in% names(s5)) {
	s5$MAE_calibrated <- fmt(s5$MAE_calibrated, 3)
}

names(s5) <- gsub("_", " ", names(s5), fixed = TRUE)

write.csv(
	s5,
	"R/tables/Table_S5_huskin_all_eligible_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S6. Matched HuSkin calibration by source
# ------------------------------------------------------------

s6 <- matchref31

if ("RMSE_uncalibrated" %in% names(s6)) {
	s6$RMSE_uncalibrated <- fmt(s6$RMSE_uncalibrated, 3)
}

if ("RMSE_calibrated" %in% names(s6)) {
	s6$RMSE_calibrated <- fmt(s6$RMSE_calibrated, 3)
}

if ("RMSE_improvement_percent" %in% names(s6)) {
	s6$RMSE_improvement_percent <- fmt(s6$RMSE_improvement_percent, 1)
}

names(s6) <- gsub("_", " ", names(s6), fixed = TRUE)

write.csv(
	s6,
	"R/tables/Table_S6_huskin_reference_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S7. Applicability-domain summary
# ------------------------------------------------------------

domain31 <- read.csv("R/results/31_applicability_domain.csv", stringsAsFactors = FALSE)
domain41 <- read.csv("R/results/41_skinpix_transfer_domain.csv", stringsAsFactors = FALSE)

s7_ch <- do.call(rbind, lapply(unique(ad21$specification), function(specification) {
	x <- ad21[ad21$specification == specification, ]
	inside <- x$in_domain %in% TRUE
	outside <- x$in_domain %in% FALSE
	data.frame(Dataset = "Cheruvu", Specification = specification, n_in = sum(x$n_compounds[inside]), n_out = sum(x$n_compounds[outside]), mae_in = x$MAE_hierarchical[inside], mae_out = x$MAE_hierarchical[outside], stringsAsFactors = FALSE)
}))

s7_hu <- data.frame(
	Dataset = "HuSkinDB",
	Specification = "S1",
	n_in = sum(domain31$in_domain %in% TRUE),
	n_out = sum(domain31$in_domain %in% FALSE),
	mae_in = mean(domain31$error_hierarchical[domain31$in_domain %in% TRUE]),
	mae_out = mean(domain31$error_hierarchical[domain31$in_domain %in% FALSE]),
	stringsAsFactors = FALSE
)

s7_sp <- data.frame(
	Dataset = "SkinPiX external test",
	Specification = "S1",
	n_in = sum(domain41$in_domain %in% TRUE),
	n_out = sum(domain41$in_domain %in% FALSE),
	mae_in = NA_real_,
	mae_out = NA_real_,
	stringsAsFactors = FALSE
)

s7 <- rbind(s7_ch, s7_hu, s7_sp)
s7$mae_in <- ifelse(is.na(s7$mae_in), "Not reported", fmt(s7$mae_in, 3))
s7$mae_out <- ifelse(is.na(s7$mae_out), "Not reported", fmt(s7$mae_out, 3))
names(s7) <- c("Dataset", "Specification", "In domain n", "Outside domain n", "MAE in domain", "MAE outside domain")

write.csv(s7, "R/tables/Table_S7_applicability_domain.csv", row.names = FALSE)

# ------------------------------------------------------------
# Supplementary Table S8A. SkinPiX transfer audit
# ------------------------------------------------------------

skinpix <- read.csv("data/processed/40_skinpix.csv", stringsAsFactors = FALSE)
eligible41 <- read.csv("R/results/41_skinpix_eligible_references.csv", stringsAsFactors = FALSE)
overlap40 <- read.csv("R/results/40_skinpix_reference_overlap.csv", stringsAsFactors = FALSE)
descriptor40 <- read.csv("R/results/40_descriptor_compatibility_audit.csv", stringsAsFactors = FALSE)

s8a <- data.frame(
	Metric = c(
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
	Value = c(
		nrow(skinpix),
		length(unique(skinpix$compound_id)),
		length(unique(skinpix$reference_id)),
		length(unique(skinpix$compound_id[skinpix$overlap_huskin %in% TRUE])),
		length(unique(skinpix$compound_id[skinpix$overlap_huskin %in% FALSE])),
		length(unique(eligible41$reference_id[eligible41$eligible %in% TRUE])),
		length(unique(overlap40$reference_id[overlap40$doi_overlaps_huskin %in% TRUE])),
		median(abs(descriptor40$delta_mw_rdkit), na.rm = TRUE),
		median(abs(descriptor40$delta_logp_rdkit), na.rm = TRUE)
	),
	stringsAsFactors = FALSE
)

write.csv(s8a, "R/tables/Table_S8A_skinpix_transfer_audit.csv", row.names = FALSE)
write.csv(s8a, "R/tables/Table_S8A_skinpix_transfer_audit_display.csv", row.names = FALSE)

# ------------------------------------------------------------
# Supplementary Table S8B. HuSkin-SkinPiX overlap
# ------------------------------------------------------------

s8b_raw <- read.csv(
	"R/results/40_skinpix_huskin_overlap_audit.csv",
	stringsAsFactors = FALSE,
	check.names = FALSE
)

s8b <- s8b_raw[, c("huskin_compound_id", "huskin_compound_name", "compound_id", "compound_name", "match_method")]

names(s8b) <- c(
	"HuSkin ID",
	"HuSkin compound",
	"SkinPiX ID",
	"SkinPiX compound",
	"Match method"
)

write.csv(
	s8b,
	"R/tables/Table_S8B_skinpix_compound_overlap.csv",
	row.names = FALSE
)

write.csv(
	s8b,
	"R/tables/Table_S8B_skinpix_compound_overlap_display.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S9. SkinPiX calibration by source
# ------------------------------------------------------------

s9 <- ref41

if ("RMSE_uncalibrated" %in% names(s9)) {
	s9$RMSE_uncalibrated <- fmt(s9$RMSE_uncalibrated, 3)
}

if ("RMSE_calibrated" %in% names(s9)) {
	s9$RMSE_calibrated <- fmt(s9$RMSE_calibrated, 3)
}

if ("RMSE_improvement_percent" %in% names(s9)) {
	s9$RMSE_improvement_percent <- fmt(s9$RMSE_improvement_percent, 1)
}

names(s9) <- gsub("_", " ", names(s9), fixed = TRUE)

write.csv(
	s9,
	"R/tables/Table_S9_skinpix_reference_calibration.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S10. Wrong-source falsification by source
# ------------------------------------------------------------

s10 <- wrongref50

if ("RMSE_correct" %in% names(s10)) {
	s10$RMSE_correct <- fmt(s10$RMSE_correct, 3)
}

if ("RMSE_wrong" %in% names(s10)) {
	s10$RMSE_wrong <- fmt(s10$RMSE_wrong, 3)
}

if ("RMSE_wrong_minus_correct" %in% names(s10)) {
	s10$RMSE_wrong_minus_correct <- fmt(s10$RMSE_wrong_minus_correct, 3)
}

if ("correct_improvement_percent" %in% names(s10)) {
	s10$correct_improvement_percent <- fmt(s10$correct_improvement_percent, 1)
}

if ("wrong_improvement_percent" %in% names(s10)) {
	s10$wrong_improvement_percent <- fmt(s10$wrong_improvement_percent, 1)
}

names(s10) <- gsub("_", " ", names(s10), fixed = TRUE)

write.csv(
	s10,
	"R/tables/Table_S10_wrong_source_by_reference.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table S11. Calibration-method comparison
# ------------------------------------------------------------

s11 <- comp32[, c(
	"cohort",
	"anchors",
	"method",
	"n_references",
	"RMSE_uncalibrated",
	"RMSE_calibrated",
	"RMSE_improvement_percent"
)]

s11$RMSE_uncalibrated <- fmt(s11$RMSE_uncalibrated, 3)
s11$RMSE_calibrated <- fmt(s11$RMSE_calibrated, 3)
s11$RMSE_improvement_percent <- fmt(s11$RMSE_improvement_percent, 1)

names(s11) <- c(
	"Cohort",
	"Anchors",
	"Method",
	"Sources",
	"RMSE uncalibrated",
	"RMSE calibrated",
	"RMSE improvement (%)"
)

write.csv(
	s11,
	"R/tables/Table_S11_calibration_method_comparison.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Supplementary Table. Experimental metadata availability
# ------------------------------------------------------------

skinpix <- read.csv("data/processed/40_skinpix.csv", stringsAsFactors = FALSE)

metadata_fields <- c("Skin donor type", "Anatomical site", "Skin layer", "Skin preparation", "Skin thickness", "Diffusion cell type", "Donor medium", "Receptor medium", "Donor or skin temperature", "Receptor temperature", "Donor pH", "Receptor pH", "Skin integrity test")

metadata_columns <- list(
	Cheruvu = c(NA, NA, NA, NA, "skin_thickness_mm", NA, NA, NA, "temperature_c", NA, NA, NA, "skin_integrity_test"),
	HuSkinDB = c("skin_source_type", "skin_source_site", "used_layer", "skin_preparation", NA, "cell_type", "donor_type", "acceptor_type", "temperature_c", "acceptor_temperature_c", "donor_ph", "acceptor_ph", NA),
	SkinPiX = c("skin_donor_type", "skin_source_site", "used_layer", "skin_preparation", NA, "cell_type", "donor_type", "acceptor_type", "temperature_c", "acceptor_temperature_c", "donor_ph", "acceptor_ph", NA)
)

metadata_notes <- c(
	"Donor origin labels require review",
	"Anatomical site labels require review",
	"Layer labels require harmonization",
	"Preparation labels require harmonization",
	"Separate numeric field available only in Cheruvu",
	"Cell labels require harmonization",
	"Medium descriptions require review",
	"Medium descriptions require review",
	"Expressed in degrees Celsius after preprocessing",
	"Expressed in degrees Celsius",
	"Does not describe full medium composition",
	"Does not describe full medium composition",
	"Test descriptions require review"
)

metadata_clean <- function(x) {
	x <- tolower(trimws(as.character(x)))
	x <- gsub("[[:space:]]+", " ", x)
	x[x %in% c("", "na", "n/a", "n.a.", "nan", "unknown", "not reported", "not specified", "-", "--")] <- NA_character_
	x
}

metadata_datasets <- list(Cheruvu = cheruvu, HuSkinDB = huskin, SkinPiX = skinpix)

table_metadata <- data.frame(Metadata = metadata_fields, stringsAsFactors = FALSE)

for (dataset in names(metadata_datasets)) {
	x <- metadata_datasets[[dataset]]
	coverage <- character(length(metadata_fields))
	for (i in seq_along(metadata_fields)) {
		column <- metadata_columns[[dataset]][i]
		if (is.na(column)) {
			coverage[i] <- "Field not provided"
		} else {
			if (!column %in% names(x)) stop(paste("Missing metadata column:", dataset, column))
			values <- metadata_clean(x[[column]])
			n_reported <- sum(!is.na(values))
			coverage[i] <- paste0(n_reported, "/", nrow(x), " (", fmt(100 * n_reported / nrow(x), 1), "%)")
		}
	}
	table_metadata[[dataset]] <- coverage
}

write.csv(table_metadata, "R/tables/Table_metadata_availability.csv", row.names = FALSE, na = "")

# ------------------------------------------------------------
# Supplementary source summary. Primary one-anchor setting.
# ------------------------------------------------------------

summary70 <- read.csv("R/results/70_failure_source_summary.csv", stringsAsFactors = FALSE)
primary70 <- summary70[summary70$anchors == 1, ]
primary70 <- primary70[order(primary70$dataset, primary70$reference_id), ]
table70 <- primary70[c("dataset", "reference_id", "available_anchors", "test_compounds", "median_anchor_replication", "median_test_replication", "RMSE_uncalibrated", "RMSE_calibrated", "test_mean_residual", "offset_mean", "offset_sd", "improved_percent", "anchor_distance_mean", "protocol_fields_varying")]
for (v in names(table70)[7:13]) table70[[v]] <- round(table70[[v]], 3)
names(table70) <- c("Dataset", "Source", "Available anchors", "Test compounds", "Median anchor records", "Median test records", "RMSE before", "Mean RMSE after", "Mean test residual", "Mean offset", "Offset SD", "Selections improved (%)", "Mean anchor distance", "Metadata fields varying")
write.csv(table70, "R/tables/Table_S_calibration_failure.csv", row.names = FALSE, na = "")


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
	"RMSE and MAE shown to 3 decimals\n",
	"AIC and BIC shown to 2 decimals\n",
	"Percentages shown to 1 decimal\n",
	"Coefficients shown to 4 decimals\n",
	"P values shown to 3 significant digits or <0.001\n",
	"Headers cleaned for manuscript display\n",
	sep = ""
)