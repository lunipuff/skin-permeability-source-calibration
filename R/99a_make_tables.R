# 99a_make_tables.R

fmt <- function(x, digits) formatC(x, format = "f", digits = digits)
fmt_ci <- function(low, high, digits) paste0("[", fmt(low, digits), ", ", fmt(high, digits), "]")
fmt_est_ci <- function(x, low, high, digits) paste(fmt(x, digits), fmt_ci(low, high, digits))

write_table <- function(x, filename, clean_names = FALSE) {
	if (clean_names) names(x) <- gsub("_", " ", names(x), fixed = TRUE)
	write.csv(x, file.path("R/tables", filename), row.names = FALSE, na = "")
}

# Inputs

cheruvu <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)
huskin <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)
skinpix <- read.csv("data/processed/40_skinpix.csv", stringsAsFactors = FALSE)
m10 <- read.csv("R/results/10_model_comparison.csv", stringsAsFactors = FALSE)
m20 <- read.csv("R/results/20_model_summary.csv", stringsAsFactors = FALSE)
coef20 <- read.csv("R/results/20_coefficients.csv", stringsAsFactors = FALSE)
ad21 <- read.csv("R/results/21_applicability_summary.csv", stringsAsFactors = FALSE)
sens21 <- read.csv("R/results/21_sensitivity_performance.csv", stringsAsFactors = FALSE)
rep31 <- read.csv("R/results/31_replication_summary.csv", stringsAsFactors = FALSE)
matchref31 <- read.csv("R/results/31_calibration_matched_reference.csv", stringsAsFactors = FALSE)
calibration31 <- read.csv("R/results/31_calibration_performance.csv", stringsAsFactors = FALSE)
domain31 <- read.csv("R/results/31_applicability_domain.csv", stringsAsFactors = FALSE)
comp32 <- read.csv("R/results/32_calibration_method_summary.csv", stringsAsFactors = FALSE)
overlap40 <- read.csv("R/results/40_skinpix_reference_overlap.csv", stringsAsFactors = FALSE)
descriptor40 <- read.csv("R/results/40_descriptor_compatibility_audit.csv", stringsAsFactors = FALSE)
audit40 <- read.csv("R/results/40_skinpix_huskin_overlap_audit.csv", stringsAsFactors = FALSE)
eligible41 <- read.csv("R/results/41_skinpix_eligible_references.csv", stringsAsFactors = FALSE)
ref41 <- read.csv("R/results/41_skinpix_calibration_reference.csv", stringsAsFactors = FALSE)
match41 <- read.csv("R/results/41_skinpix_calibration_matched_3anchor.csv", stringsAsFactors = FALSE)
domain41 <- read.csv("R/results/41_skinpix_transfer_domain.csv", stringsAsFactors = FALSE)
wrongref50 <- read.csv("R/results/50_wrong_source_reference.csv", stringsAsFactors = FALSE)
loco60 <- read.csv("R/results/60_loco_bootstrap_ci.csv", stringsAsFactors = FALSE)
var60 <- read.csv("R/results/60_variance_component_ci.csv", stringsAsFactors = FALSE)
deploy60 <- read.csv("R/results/60_huskin_one_anchor_all_eligible_ci.csv", stringsAsFactors = FALSE)
cal60 <- read.csv("R/results/60_huskin_calibration_ci.csv", stringsAsFactors = FALSE)
wrong60 <- read.csv("R/results/60_wrong_source_ci.csv", stringsAsFactors = FALSE)
summary70 <- read.csv("R/results/70_failure_source_summary.csv", stringsAsFactors = FALSE)
selection70 <- read.csv("R/results/70_failure_anchor_selections.csv", stringsAsFactors = FALSE)

# 01. Hierarchical variance decomposition

components <- c("Compound", "Source", "Residual")

variance_order <- unlist(lapply(c("Cheruvu", "HuSkinDB"), function(dataset) {
	match(paste(dataset, components), paste(var60$dataset, var60$component))
}), use.names = FALSE)

variance <- var60[variance_order, ]

table01 <- data.frame(
	Dataset = variance$dataset,
	Component = variance$component,
	`Component (%)` = fmt(variance$estimate_percent, 1),
	`95% CI` = fmt_ci(variance$ci_low, variance$ci_high, 1),
	check.names = FALSE
)

write_table(table01, "01_hierarchical_decomposition.csv")

# 02. Source information and LOCO prediction

approach_names <- c(
	conventional = "Conventional aggregate QSPR",
	hierarchical_source_blind = "Hierarchical, source blind",
	hierarchical_known_source = "Hierarchical, source aware"
)

table02 <- data.frame(
	Dataset = loco60$dataset,
	Specification = loco60$specification,
	`Prediction setting` = unname(approach_names[loco60$approach]),
	`RMSE (95% CI)` = fmt_est_ci(loco60$RMSE, loco60$RMSE_ci_low, loco60$RMSE_ci_high, 3),
	`RMSE improvement vs conventional (%)` = ifelse(
		loco60$approach == "conventional", "Reference",
		fmt_est_ci(loco60$improvement_percent, loco60$improvement_ci_low, loco60$improvement_ci_high, 1)
	),
	check.names = FALSE
)

write_table(table02, "02_source_information_prediction.csv")

# 03. HuSkinDB source calibration

matched_improved <- aggregate(RMSE_improvement_percent ~ anchors, data = matchref31, function(x) sum(x > 0, na.rm = TRUE))

names(matched_improved)[2] <- "n_improved"

primary03 <- data.frame(
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

matched03 <- data.frame(
	Cohort = paste0("Matched ", max(cal60$anchors), "-anchor-capable cohort"),
	Anchors = cal60$anchors,
	Sources = cal60$n_references,
	`Sources improved` = paste0(matched_improved$n_improved[match(cal60$anchors, matched_improved$anchors)], "/", cal60$n_references),
	`RMSE uncalibrated` = fmt(cal60$RMSE_uncalibrated, 3),
	`RMSE calibrated` = fmt(cal60$RMSE_calibrated, 3),
	`RMSE improvement (%)` = fmt(cal60$improvement_percent, 1),
	`95% CI` = fmt_ci(cal60$improvement_ci_low, cal60$improvement_ci_high, 1),
	check.names = FALSE
)

table03 <- rbind(primary03, matched03)

write_table(table03, "03_huskin_calibration.csv")

# 04. SkinPiX transfer

one41 <- ref41[ref41$anchors == 1, ]

primary04 <- data.frame(
	Cohort = "All eligible sources",
	Anchors = 1,
	Sources = length(unique(one41$reference_id)),
	`RMSE uncalibrated` = mean(one41$RMSE_uncalibrated),
	`RMSE calibrated` = mean(one41$RMSE_calibrated),
	check.names = FALSE
)

primary04$`RMSE improvement (%)` <- 100 * (1 - primary04$`RMSE calibrated` / primary04$`RMSE uncalibrated`)

matched04 <- data.frame(
	Cohort = "Matched cohort",
	Anchors = match41$anchors,
	Sources = match41$n_references,
	`RMSE uncalibrated` = match41$RMSE_uncalibrated,
	`RMSE calibrated` = match41$RMSE_calibrated,
	`RMSE improvement (%)` = match41$RMSE_improvement_percent,
	check.names = FALSE
)

table04 <- rbind(primary04, matched04)

for (v in c("RMSE uncalibrated", "RMSE calibrated")) table04[[v]] <- fmt(table04[[v]], 3)

table04$`RMSE improvement (%)` <- fmt(table04$`RMSE improvement (%)`, 1)

write_table(table04, "04_skinpix_transfer.csv")

# 05. Dataset structure

datasets <- list(Cheruvu = cheruvu, HuSkinDB = huskin, SkinPiX = skinpix)

roles <- c(
	Cheruvu = "Development",
	HuSkinDB = "Replication and held-out-source calibration",
	SkinPiX = "Cross-dataset transfer"
)

table05 <- do.call(rbind, lapply(names(datasets), function(dataset) {
	x <- datasets[[dataset]]
	data.frame(
		Dataset = dataset,
		Measurements = nrow(x),
		Compounds = length(unique(x$compound_id)),
		Sources = length(unique(x$reference_id)),
		`Repeated compounds` = sum(table(x$compound_id) > 1),
		`Multi-source compounds` = sum(tapply(x$reference_id, x$compound_id, function(z) length(unique(z))) > 1),
		Role = unname(roles[dataset]),
		check.names = FALSE
	)
}))

write_table(table05, "05_dataset_structure.csv")

# S1. Cheruvu hierarchical model comparison

s1 <- m10

for (v in c("AIC", "BIC", "logLik")) s1[[v]] <- fmt(s1[[v]], 2)

write_table(s1, "S1_cheruvu_hierarchical_models.csv", clean_names = TRUE)

# S2. Experimental metadata availability

metadata_fields <- c(
	"Skin donor type", "Anatomical site", "Skin layer", "Skin preparation",
	"Skin thickness", "Diffusion cell type", "Donor medium", "Receptor medium",
	"Donor or skin temperature", "Receptor temperature", "Donor pH", "Receptor pH",
	"Skin integrity test"
)

metadata_columns <- list(
	Cheruvu = c(
		NA, NA, NA, NA, "skin_thickness_mm", NA, NA, NA,
		"temperature_c", NA, NA, NA, "skin_integrity_test"
	),
	HuSkinDB = c(
		"skin_source_type", "skin_source_site", "used_layer", "skin_preparation",
		NA, "cell_type", "donor_type", "acceptor_type", "temperature_c",
		"acceptor_temperature_c", "donor_ph", "acceptor_ph", NA
	),
	SkinPiX = c(
		"skin_donor_type", "skin_source_site", "used_layer", "skin_preparation",
		NA, "cell_type", "donor_type", "acceptor_type", "temperature_c",
		"acceptor_temperature_c", "donor_ph", "acceptor_ph", NA
	)
)

metadata_clean <- function(x) {
	x <- tolower(trimws(as.character(x)))
	x <- gsub("[[:space:]]+", " ", x)
	x[x %in% c("", "na", "n/a", "n.a.", "nan", "unknown", "not reported", "not specified", "-", "--")] <- NA_character_
	x
}

s2 <- data.frame(Metadata = metadata_fields, stringsAsFactors = FALSE)

for (dataset in names(datasets)) {
	x <- datasets[[dataset]]
	s2[[dataset]] <- vapply(metadata_columns[[dataset]], function(column) {
		if (is.na(column)) return("Field not provided")
		if (!column %in% names(x)) stop("Missing metadata column: ", dataset, " ", column)
		n_reported <- sum(!is.na(metadata_clean(x[[column]])))
		paste0(n_reported, "/", nrow(x), " (", fmt(100 * n_reported / nrow(x), 1), "%)")
	}, character(1))
}

write_table(s2, "S2_metadata_availability.csv")

# S3. HuSkinDB replication
write_table(rep31, "S3_huskin_replication.csv", clean_names = TRUE)

# S4. Cheruvu sensitivity analyses
s4 <- sens21
for (v in c("RMSE", "MAE", "R2", "r")) s4[[v]] <- fmt(s4[[v]], 3)
write_table(s4, "S4_cheruvu_sensitivity.csv", clean_names = TRUE)

# S5. QSPR model summary
# Retain numeric values because the Rmd formats this transposed table.
write_table(m20, "S5_qspr_model_summary.csv", clean_names = TRUE)

# S6. QSPR coefficients
s6 <- coef20
for (v in c("estimate", "se")) s6[[v]] <- fmt(s6[[v]], 4)
write_table(s6, "S6_qspr_coefficients.csv", clean_names = TRUE)

# S7. Matched HuSkinDB calibration by source
s7 <- matchref31

for (v in c("RMSE_uncalibrated", "RMSE_calibrated", "MAE_uncalibrated", "MAE_calibrated")) s7[[v]] <- fmt(s7[[v]], 3)
for (v in c("RMSE_improvement_percent", "MAE_improvement_percent")) s7[[v]] <- fmt(s7[[v]], 1)

write_table(s7, "S7_huskin_reference_calibration.csv", clean_names = TRUE)

# S8. Conventional versus hierarchical calibration
s8 <- comp32[c("cohort", "anchors", "method", "n_references", "RMSE_uncalibrated", "RMSE_calibrated", "RMSE_improvement_percent")]

for (v in c("RMSE_uncalibrated", "RMSE_calibrated")) s8[[v]] <- fmt(s8[[v]], 3)

s8$RMSE_improvement_percent <- fmt(s8$RMSE_improvement_percent, 1)
names(s8) <- c("Cohort", "Anchors", "Method", "Sources", "RMSE uncalibrated", "RMSE calibrated", "RMSE improvement (%)")

write_table(s8, "S8_calibration_method_comparison.csv")

# S9. Wrong-source control

s9 <- data.frame(
	Anchors = wrong60$anchors,
	`Target sources` = wrong60$n_references,
	`Target-source RMSE reduction (%)` = fmt_est_ci(wrong60$correct_improvement_percent, wrong60$correct_ci_low, wrong60$correct_ci_high, 1),
	`Different-source RMSE reduction (%)` = fmt_est_ci(wrong60$wrong_improvement_percent, wrong60$wrong_ci_low, wrong60$wrong_ci_high, 1),
	`RMSE difference` = fmt_est_ci(wrong60$RMSE_wrong_minus_correct, wrong60$difference_ci_low, wrong60$difference_ci_high, 3),
	check.names = FALSE
)

write_table(s9, "S9_wrong_source_falsification.csv")

# S10. Wrong-source control by target source
# Retain numeric values because the Rmd applies display formatting.
s10 <- wrongref50[c("reference_id", "anchors", "n_donors", "n_test", "RMSE_uncalibrated", "RMSE_correct", "RMSE_wrong", "RMSE_wrong_minus_correct")]
write_table(s10, "S10_wrong_source_by_reference.csv", clean_names = TRUE)

# S11. SkinPiX transfer audit

s11 <- data.frame(
	Metric = c(
		"Retained measurements", "Compounds", "References", "HuSkin-overlapping compounds",
		"SkinPiX-only compounds", "Eligible transfer references", "References with DOI overlap",
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

write_table(s11, "S11_skinpix_transfer_audit.csv")

# S12. HuSkinDB–SkinPiX compound matches
s12 <- audit40[c("huskin_compound_id", "huskin_compound_name", "compound_id", "compound_name", "match_method")]
names(s12) <- c("HuSkin ID", "HuSkin compound", "SkinPiX ID", "SkinPiX compound", "Match method")
write_table(s12, "S12_skinpix_compound_overlap.csv")

# S13. SkinPiX calibration by source
# Retain numeric values because the Rmd applies display formatting.
s13 <- ref41[c(
	"reference_id", "anchors", "RMSE_uncalibrated", "RMSE_calibrated",
	"MAE_uncalibrated", "MAE_calibrated", "RMSE_improvement_percent", "MAE_improvement_percent"
)]
write_table(s13, "S13_skinpix_reference_calibration.csv", clean_names = TRUE)

# S14. Applicability-domain summary

s14_ch <- do.call(rbind, lapply(unique(ad21$specification), function(specification) {
	x <- ad21[ad21$specification == specification, ]
	inside <- x$in_domain %in% TRUE
	outside <- x$in_domain %in% FALSE
	data.frame(
		Dataset = "Cheruvu",
		Specification = specification,
		n_in = sum(x$n_compounds[inside]),
		n_out = sum(x$n_compounds[outside]),
		mae_in = if (any(inside)) x$MAE_hierarchical[inside] else NA_real_,
		mae_out = if (any(outside)) x$MAE_hierarchical[outside] else NA_real_
	)
}))

s14_hu <- data.frame(
	Dataset = "HuSkinDB",
	Specification = "S1",
	n_in = sum(domain31$in_domain %in% TRUE),
	n_out = sum(domain31$in_domain %in% FALSE),
	mae_in = mean(domain31$error_hierarchical[domain31$in_domain %in% TRUE]),
	mae_out = mean(domain31$error_hierarchical[domain31$in_domain %in% FALSE])
)

s14_sp <- data.frame(
	Dataset = "SkinPiX external test",
	Specification = "S1",
	n_in = sum(domain41$in_domain %in% TRUE),
	n_out = sum(domain41$in_domain %in% FALSE),
	mae_in = NA_real_,
	mae_out = NA_real_
)

s14 <- rbind(s14_ch, s14_hu, s14_sp)

for (v in c("mae_in", "mae_out")) s14[[v]] <- ifelse(is.na(s14[[v]]), "Not reported", fmt(s14[[v]], 3))

names(s14) <- c("Dataset", "Specification", "In domain n", "Outside domain n", "MAE in domain", "MAE outside domain")

write_table(s14, "S14_applicability_domain.csv")

# S15a and S15b. Failure summaries and squared-error decomposition

primary70 <- summary70[summary70$anchors == 1, ]
primary70 <- primary70[order(primary70$dataset, primary70$reference_id), ]

primary70$test_bias_squared <- primary70$test_mean_residual^2
primary70$mean_MSE_change <- primary70$offset_variance + primary70$offset_bias_squared - primary70$test_bias_squared

direct_change <- vapply(seq_len(nrow(primary70)), function(i) {
	x <- selection70[
		selection70$dataset == primary70$dataset[i] &
		selection70$reference_id == primary70$reference_id[i] &
		selection70$anchors == 1,
	]
	stopifnot(nrow(x) == primary70$selections[i])
	mean(x$RMSE_calibrated^2) - x$RMSE_uncalibrated[1]^2
}, numeric(1))
stopifnot(max(abs(primary70$mean_MSE_change - direct_change)) < 1e-8)

s15a <- primary70[c(
	"dataset", "reference_id", "available_anchors", "test_compounds",
	"median_anchor_replication", "median_test_replication", "RMSE_uncalibrated",
	"RMSE_calibrated", "test_mean_residual", "offset_mean", "offset_sd",
	"improved_percent", "anchor_distance_mean", "protocol_fields_varying"
)]

for (v in c("RMSE_uncalibrated", "RMSE_calibrated", "test_mean_residual", "offset_mean", "offset_sd", "improved_percent", "anchor_distance_mean")) {
	s15a[[v]] <- round(s15a[[v]], 3)
}

names(s15a) <- c(
	"Dataset", "Source", "Available anchors", "Test compounds",
	"Median anchor records", "Median test records", "RMSE before", "Mean RMSE after",
	"Mean test residual", "Mean offset", "Offset SD", "Selections improved (%)",
	"Mean anchor distance", "Metadata fields varying"
)

write_table(s15a, "S15a_calibration_failure.csv")

s15b <- primary70[c("dataset", "reference_id", "test_bias_squared", "offset_bias_squared", "offset_variance", "mean_MSE_change")]

for (v in c("test_bias_squared", "offset_bias_squared", "offset_variance", "mean_MSE_change")) s15b[[v]] <- round(s15b[[v]], 4)

names(s15b) <- c("Dataset", "Source", "Test bias sq.", "Mismatch sq.", "Offset var.", "Mean MSE change")

write_table(s15b, "S15b_calibration_failure_decomposition.csv")

# Shared-compound check used by manuscript calculations

anchor_check <- aggregate(
	cbind(RMSE_uncalibrated, RMSE_calibrated) ~ reference_id + task,
	data = calibration31[calibration31$anchors == 1, ], mean
)

shared_check <- anchor_check[anchor_check$task == "Known compound", c("reference_id", "RMSE_uncalibrated", "RMSE_calibrated")]
novel_check <- anchor_check[anchor_check$task == "New compound", c("reference_id", "RMSE_uncalibrated", "RMSE_calibrated")]

names(shared_check)[2:3] <- c("shared_before", "shared_after")
names(novel_check)[2:3] <- c("novel_before", "novel_after")
anchor_check <- merge(shared_check, novel_check, by = "reference_id")
anchor_check$shared_improved <- anchor_check$shared_after < anchor_check$shared_before
anchor_check$novel_improved <- anchor_check$novel_after < anchor_check$novel_before

write_table(anchor_check[c("reference_id", "shared_improved", "novel_improved")], "anchor_check.csv")

cat(
	"\nManuscript tables complete\n",
	"Main tables: 5\n",
	"Supplementary table files: 16\n",
	"Shared-compound check sources: ", nrow(anchor_check), "\n",
	"Failure-diagnostic sources: ", nrow(primary70), "\n",
	"Squared-error decomposition verified\n",
	sep = ""
)