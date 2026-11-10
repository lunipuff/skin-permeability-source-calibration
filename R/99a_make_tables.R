# Format the tables and Additional files used by the pooled manuscript.
# Run from the project root after the pooled analyses.

result_dir <- "R/results"
table_dir <- "R/tables"
sensitivity_dir <- file.path(result_dir, "sensitivity")

read_checked <- function(path, required = character()) {
	x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
	if (!all(required %in% names(x))) {
		stop("Missing columns in ", path, ": ", paste(setdiff(required, names(x)), collapse = ", "))
	}
	x
}
read_result <- function(name, required = character()) {
	read_checked(file.path(result_dir, name), required)
}
fmt <- function(x, digits = 3) {
	out <- rep("", length(x))
	ok <- is.finite(x)
	out[ok] <- formatC(x[ok], format = "f", digits = digits)
	out
}
fmt_ci <- function(low, high, digits = 3) {
	out <- rep("", length(low))
	ok <- is.finite(low) & is.finite(high)
	out[ok] <- paste0("[", fmt(low[ok], digits), ", ", fmt(high[ok], digits), "]")
	out
}
fmt_est_ci <- function(x, low, high, digits = 3) {
	trimws(paste(fmt(x, digits), fmt_ci(low, high, digits)))
}
close_enough <- function(x, y) {
	length(x) == length(y) && all(is.finite(x) & is.finite(y) &
		abs(x - y) <= 1e-7 * pmax(1, abs(y)))
}
write_table <- function(x, filename) {
	write.csv(x, file.path(table_dir, filename), row.names = FALSE, na = "")
}

# Require exhaustive calibration for the primary and agreement analyses.
panel_settings_files <- c(
	file.path(result_dir, "81_analysis_settings.csv"),
	list.files(file.path(result_dir, "sensitivity"),
		pattern = "^81_analysis_settings\\.csv$", recursive = TRUE, full.names = TRUE)
)
for (settings_file in panel_settings_files) {
	panel_settings <- read_checked(settings_file, "panel_method")
	if (nrow(panel_settings) != 1L || panel_settings$panel_method != "exhaustive") {
		stop("Rerun exhaustive calibration and sensitivity analyses before reporting: ", settings_file)
	}
}

# Read pooled analysis outputs.
pooled_path <- "data/processed/80_pooled.csv"
pooled <- read_checked(pooled_path, c(
	"record_id", "publication_id", "compound_id", "compound_name", "dataset_origin"
))
structure <- read_result("85_dataset_structure.csv", c(
	"dataset", "raw_measurements", "retained_measurements", "compounds", "publications"
))
metadata <- read_result("85_metadata_availability.csv", c(
	"dataset", "metadata", "field_provided", "n_reported", "n_records", "percent_reported"
))
variance <- read_result("85_variance_components.csv", c(
	"component", "variance", "variance_ci_low", "variance_ci_high",
	"percent", "percent_ci_low", "percent_ci_high"
))
cal_cols <- c("anchors", "n_sources", "n_sources_improved",
	"mean_source_RMSE_uncalibrated", "mean_source_RMSE_calibrated",
	"RMSE_reduction_percent", "reduction_ci_low", "reduction_ci_high")
all_cal <- read_result("81_pooled_calibration_summary_all.csv", c(
	"anchors", "n_sources", "RMSE_reduction_percent", "reduction_ci_low", "reduction_ci_high"
))
matched_cal <- read_result("81_pooled_calibration_summary_matched.csv", cal_cols)
paired <- read_result("81_pooled_one_vs_five_anchors.csv", c(
	"n_sources", "mean_RMSE_decrease", "decrease_ci_low", "decrease_ci_high"
))
by_source <- read_result("81_pooled_calibration_by_source.csv", c(
	"source_id", "anchors", "RMSE_uncalibrated", "RMSE_calibrated",
	"RMSE_reduction_percent", "n_available_anchors", "n_test_compounds", "n_test_measurements", "n_panels"
))
eligibility <- read_result("81_pooled_source_eligibility.csv", c(
	"source_id", "n_measurements", "n_compounds", "n_anchors",
	"n_test_compounds", "n_test_measurements"
))
sensitivity <- read_checked(file.path(sensitivity_dir, "83_sensitivity_comparison.csv"), c(
	"scenario", "curation_policy", "log_tolerance", "n_measurements",
	"one_anchor_publications", "matched_publications", "one_anchor_reduction_percent",
	"one_anchor_ci_low", "one_anchor_ci_high", "one_to_five_RMSE_decrease",
	"one_to_five_ci_low", "one_to_five_ci_high"
))
run_status <- read_checked(file.path(sensitivity_dir, "83_sensitivity_run_status.csv"), c("scenario", "status"))
failure <- read_result("84_failure_source_summary.csv", c(
	"source_id", "anchors", "draws_improved_percent", "RMSE_calibrated",
	"offset_variance", "offset_mismatch_squared",
	"correctable_test_bias_squared", "mean_MSE_change", "dominant_penalty"
))
failure_paired <- read_result("84_one_vs_five_diagnostics.csv", c(
	"source_id", "RMSE_calibrated_one", "RMSE_calibrated_five", "RMSE_decrease",
	"offset_variance_one", "offset_variance_five",
	"offset_mismatch_squared_one", "offset_mismatch_squared_five"
))
corrections <- read_result("80_verified_unit_corrections.csv", c(
	"record_id", "previous_log_kp", "corrected_log_kp", "original_kp_cm_h", "evidence"
))

# Check that upstream results describe the same analysis.
one <- all_cal[all_cal$anchors == 1, , drop = FALSE]
primary_sensitivity <- sensitivity[sensitivity$scenario == "primary", , drop = FALSE]
pooled_counts <- structure[structure$dataset == "Pooled", , drop = FALSE]
if (nrow(one) != 1L || nrow(primary_sensitivity) != 1L ||
	nrow(pooled_counts) != 1L || nrow(paired) != 1L) stop("Unexpected summary rows.")
if (anyDuplicated(pooled$record_id) || anyNA(pooled[c("publication_id", "compound_id")])) {
	stop("Invalid pooled identities.")
}
if (pooled_counts$retained_measurements != nrow(pooled) ||
		pooled_counts$compounds != length(unique(pooled$compound_id)) ||
		pooled_counts$publications != length(unique(pooled$publication_id))) stop("Dataset counts disagree.")
if (primary_sensitivity$n_measurements != nrow(pooled) ||
		primary_sensitivity$one_anchor_publications != one$n_sources ||
		primary_sensitivity$matched_publications != paired$n_sources) {
	stop("Sensitivity sample counts disagree with the main analysis. Rerun 83.")
}
if (!close_enough(
	c(one$RMSE_reduction_percent, one$reduction_ci_low, one$reduction_ci_high,
		paired$mean_RMSE_decrease, paired$decrease_ci_low, paired$decrease_ci_high),
	c(primary_sensitivity$one_anchor_reduction_percent, primary_sensitivity$one_anchor_ci_low,
		primary_sensitivity$one_anchor_ci_high, primary_sensitivity$one_to_five_RMSE_decrease,
		primary_sensitivity$one_to_five_ci_low, primary_sensitivity$one_to_five_ci_high)
)) {
	stop("Sensitivity summary is stale. Rerun 83.")
}
requested <- sensitivity$scenario[sensitivity$scenario != "primary"]
if (!length(requested) || anyDuplicated(run_status$scenario) ||
		!setequal(requested, run_status$scenario) || any(run_status$status != "SUCCESS")) {
	stop("Sensitivity runs are incomplete or do not match the comparison table.")
}
performance_key <- paste(by_source$source_id, by_source$anchors)
diagnostic_key <- paste(failure$source_id, failure$anchors)
if (anyDuplicated(performance_key) ||
	anyDuplicated(diagnostic_key) ||
	!setequal(performance_key, diagnostic_key)) {
	stop("Performance and diagnostic rows do not match. Rerun 84.")
}
diagnostic_index <- match(performance_key, diagnostic_key)
if (!close_enough(by_source$RMSE_calibrated, failure$RMSE_calibrated[diagnostic_index])) {
	stop("Publication-level failure results disagree with 81. Rerun 84.")
}

pooled_md5 <- unname(tools::md5sum(pooled_path))
for (stage in c("81", "85")) {
	settings <- read_result(paste0(stage, "_analysis_settings.csv"))
	hash <- if (stage == "81") settings$input_md5 else settings$md5
	if (length(hash) != 1L || is.na(hash) || hash != pooled_md5) {
		stop("Input fingerprint differs for stage ", stage, ". Rerun that stage.")
	}
}

# Assign the publication labels shared by the tables and figures.
ids <- sort(unique(as.character(pooled$publication_id)))
publication_key <- data.frame(
	publication_label = sprintf("P%03d", seq_along(ids)), source_id = ids,
	stringsAsFactors = FALSE
)
label <- function(x) {
	j <- match(x, publication_key$source_id)
	if (anyNA(j)) stop("Publication missing from primary index: ", paste(x[is.na(j)], collapse = ", "))
	publication_key$publication_label[j]
}

# Dataset structure and agreement-threshold sensitivity.
t1 <- structure
names(t1) <- c("Dataset", "Raw records", "Retained measurements", "Compounds", "Publications")
write_table(t1, "01_dataset_structure.csv")

policy_labels <- ifelse(sensitivity$curation_policy == "priority", "Database priority (primary)",
	paste0("Agreement threshold ", sensitivity$log_tolerance))
t3 <- data.frame(
	Policy = policy_labels, Measurements = sensitivity$n_measurements,
	`Publications: one / five anchors` = paste0(sensitivity$one_anchor_publications, " / ", sensitivity$matched_publications),
	`One-anchor reduction (%) [95% CI]` = fmt_est_ci(sensitivity$one_anchor_reduction_percent, sensitivity$one_anchor_ci_low, sensitivity$one_anchor_ci_high, 1),
	`One-to-five RMSE decrease [95% CI]` = fmt_est_ci(sensitivity$one_to_five_RMSE_decrease, sensitivity$one_to_five_ci_low, sensitivity$one_to_five_ci_high),
	check.names = FALSE
)
write_table(t3, "03a_curation_sensitivity.csv")

# Publication index and descriptive tables.
pub <- merge(publication_key, eligibility, by = "source_id", all.x = TRUE, sort = FALSE)
if (anyNA(pub$n_measurements)) stop("Incomplete publication eligibility index.")
pub$dataset_origins <- vapply(pub$source_id, function(id) {
	paste(sort(unique(pooled$dataset_origin[pooled$publication_id == id])), collapse = "; ")
}, character(1))
pub <- pub[order(pub$publication_label), ]
write_table(pub, "S1_publication_index.csv")

meta_table <- data.frame(Metadata = unique(metadata$metadata), check.names = FALSE)
for (origin in c("Cheruvu", "HuSkinDB", "SkinPiX", "Pooled")) {
	x <- metadata[metadata$dataset == origin, ]
	j <- match(meta_table$Metadata, x$metadata)
	if (anyNA(j) || anyDuplicated(x$metadata)) stop("Incomplete metadata audit: ", origin)
	x <- x[j, ]
	meta_table[[origin]] <- ifelse(x$n_records == 0, "No retained records",
		ifelse(!x$field_provided, "Field not provided",
			paste0(x$n_reported, "/", x$n_records, " (", fmt(x$percent_reported, 1), "%)")))
}
write_table(meta_table, "S2_metadata_availability.csv")

write_table(
	data.frame(
		Component = variance$component,
		`Variance [95% CI]` = fmt_est_ci(variance$variance, variance$variance_ci_low, variance$variance_ci_high),
		`Variance fraction (%) [95% CI]` = fmt_est_ci(
			variance$percent, variance$percent_ci_low, variance$percent_ci_high, 1
		),
		check.names = FALSE
	),
	"S3_variance_components.csv"
)

# Calibration and failure diagnostics.
write_table(data.frame(
	Cohort = rep("Fixed five-anchor-eligible subset", nrow(matched_cal)),
	Anchors = matched_cal$anchors,
	`Publications improved` = paste0(matched_cal$n_sources_improved, "/", matched_cal$n_sources),
	`RMSE before` = fmt(matched_cal$mean_source_RMSE_uncalibrated),
	`RMSE after` = fmt(matched_cal$mean_source_RMSE_calibrated),
	`RMSE reduction (%) [95% CI]` = fmt_est_ci(
		matched_cal$RMSE_reduction_percent, matched_cal$reduction_ci_low, matched_cal$reduction_ci_high, 1
	),
	check.names = FALSE
), "S6_calibration_matched.csv")

write_table(
	data.frame(
		Publication = label(by_source$source_id),
		Anchors = by_source$anchors,
		`Available anchors` = by_source$n_available_anchors,
		`Test compounds` = by_source$n_test_compounds,
		`Test measurements` = by_source$n_test_measurements,
		`Panels evaluated` = by_source$n_panels,
		`RMSE before` = by_source$RMSE_uncalibrated,
		`RMSE after` = by_source$RMSE_calibrated,
		`RMSE reduction (%)` = by_source$RMSE_reduction_percent,
		`Panels improved (%)` = failure$draws_improved_percent[diagnostic_index],
		check.names = FALSE
	),
	"S7_calibration_by_publication.csv"
)

write_table(data.frame(
	Publication = label(failure$source_id), Anchors = failure$anchors,
	`Offset variance` = fmt(failure$offset_variance),
	`Squared mean-offset mismatch` = fmt(failure$offset_mismatch_squared),
	`Squared mean test bias` = fmt(failure$correctable_test_bias_squared),
	`Mean MSE change` = fmt(failure$mean_MSE_change),
	`Larger penalty term` = failure$dominant_penalty, check.names = FALSE
), "S10_failure_decomposition.csv")

write_table(data.frame(
	Publication = label(failure_paired$source_id),
	`One-anchor RMSE` = fmt(failure_paired$RMSE_calibrated_one),
	`Five-anchor RMSE` = fmt(failure_paired$RMSE_calibrated_five),
	`RMSE decrease` = fmt(failure_paired$RMSE_decrease),
	`Offset variance: one` = fmt(failure_paired$offset_variance_one),
	`Offset variance: five` = fmt(failure_paired$offset_variance_five),
	`Squared mismatch: one` = fmt(failure_paired$offset_mismatch_squared_one),
	`Squared mismatch: five` = fmt(failure_paired$offset_mismatch_squared_five), check.names = FALSE
), "S11_one_vs_five_diagnostics.csv")

# Verified unit corrections.
correction_index <- match(corrections$record_id, pooled$record_id)
if (anyNA(correction_index)) {
	stop("Correction records are missing from the pooled dataset.")
}

correction_study <- ifelse(
	grepl("Barber 1992", corrections$evidence, fixed = TRUE),
	"Barber et al. (1992)",
	ifelse(
		grepl("Barber 1995", corrections$evidence, fixed = TRUE),
		"Barber et al. (1995)",
		NA_character_
	)
)
if (anyNA(correction_study)) {
	stop("Missing study citation for a correction record.")
}

write_table(
	data.frame(
		Record = corrections$record_id,
		Compound = pooled$compound_name[correction_index],
		`Previously processed log10 Kp` = fmt(corrections$previous_log_kp, 3),
		`Corrected log10 Kp` = fmt(corrections$corrected_log_kp, 3),
		`Original Kp (cm/h)` = formatC(corrections$original_kp_cm_h, format = "e", digits = 3),
		Study = correction_study,
		check.names = FALSE
	),
	"S12_verified_unit_corrections.csv"
)

# Export the submission files with names that match the manuscript.
additional_dir <- "R/additional_files"
if (!dir.exists(additional_dir) && !dir.create(additional_dir)) {
	stop("Cannot create: ", additional_dir, call. = FALSE)
}
additional_files <- c(
	"Additional_file_1_publication_index.csv" = file.path(table_dir, "S1_publication_index.csv"),
	"Additional_file_2_calibration_by_publication.csv" = file.path(table_dir, "S7_calibration_by_publication.csv"),
	"Additional_file_3_failure_decomposition.csv" = file.path(table_dir, "S10_failure_decomposition.csv"),
	"Additional_file_4_pooled_dataset.csv" = pooled_path
)
copied <- file.copy(
	unname(additional_files),
	file.path(additional_dir, names(additional_files)),
	overwrite = TRUE
)
if (!all(copied)) {
	stop("Could not export: ", paste(names(additional_files)[!copied], collapse = ", "), call. = FALSE)
}

message("Manuscript tables written to ", table_dir, ".")
message("Additional files written to ", additional_dir, ".")
