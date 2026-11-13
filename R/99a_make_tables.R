# Format the tables and Additional files used by the pooled manuscript.
# Run from the project root after the pooled analyses.

result_dir <- "R/results"
table_dir <- "R/tables"
sensitivity_dir <- file.path(result_dir, "sensitivity")
additional_dir <- "R/additional_files"

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
write_additional <- function(x, filename) {
	write.csv(x, file.path(additional_dir, filename), row.names = FALSE, na = "")
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

settings81 <- read_result("81_analysis_settings.csv", c("n_boot", "seed"))
priority <- read_result("priority_sensitivity/87_priority_comparison.csv")
test_sensitivity <- read_result("86_test_set_sensitivity.csv")
sensitivity_influence <- read_result("sensitivity/83_sensitivity_influence.csv")
influence <- read_result("82_source_influence_summary.csv")
influence_gain <- influence[influence$analysis == "RMSE decrease from one to five anchors (log units)", ]
test_one <- test_sensitivity[test_sensitivity$endpoint == "one_anchor" & test_sensitivity$minimum_test_compounds == 2L, ]
test_gain <- test_sensitivity[test_sensitivity$endpoint == "one_to_five" & test_sensitivity$minimum_test_compounds == 2L, ]

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
pub <- merge(publication_key, eligibility, by = "source_id", all.x = TRUE, sort = FALSE)
if (anyNA(pub$n_measurements)) stop("Incomplete publication eligibility index.")
pub$dataset_origins <- vapply(pub$source_id, function(id) {
	paste(sort(unique(pooled$dataset_origin[pooled$publication_id == id])), collapse = "; ")
}, character(1))
pub <- pub[order(pub$publication_label), ]

# Table 1. Dataset provenance.
t1 <- structure
names(t1) <- c("Dataset", "Raw records", "Retained measurements", "Compounds", "Publications")
write_table(t1, "01_dataset_structure.csv")

# Table 2. Sensitivity overview in absolute RMSE units.
# Absolute one-anchor decreases and their source-level bootstrap interval.
t2_one <- by_source[by_source$anchors == 1L, , drop = FALSE]
t2_one <- t2_one[order(t2_one$source_id), ]
t2_one_gain <- t2_one$RMSE_uncalibrated - t2_one$RMSE_calibrated

set.seed(settings81$seed[[1]] + 1L)
t2_boot <- replicate(settings81$n_boot[[1]], {
	i <- sample.int(length(t2_one_gain), replace = TRUE)
	mean(t2_one_gain[i])
})
t2_one_ci <- quantile(t2_boot, c(0.025, 0.975), names = FALSE)

# Each pooling analysis uses its own uncalibrated baseline.
t2_agreement <- sensitivity[
	sensitivity$curation_policy == "agreement",
	, drop = FALSE
]
t2_agreement_one <- with(
	t2_agreement,
	one_anchor_RMSE_before - one_anchor_RMSE_after
)
t2_priority_one <- with(
	priority,
	one_anchor_RMSE_before - one_anchor_RMSE_after
)
t2_restricted_one <- with(
	test_one,
	mean_RMSE_before - mean_RMSE_after
)

# Omit each source from the primary one-anchor evaluation summary.
t2_omission_one <- (
	sum(t2_one_gain) - t2_one_gain
) / (length(t2_one_gain) - 1L)

t2_range <- function(x) {
	paste0(fmt(min(x), 3), " to ", fmt(max(x), 3))
}

sensitivity_overview <- data.frame(
	Analysis = c(
		"Primary analysis",
		"Applied different agreement thresholds to overlapping records (Table~\\ref{tab:table-s5})",
		"Changed which database took priority (Table~\\ref{tab:table-s6})",
		"Kept studies with at least two test compounds (Table~\\ref{tab:table-s7})",
		"Left one study at a time out of the primary evaluation summary (Table~\\ref{tab:table-s8})"
	),
	`One anchor vs none` = c(
		paste0(
			fmt(mean(t2_one_gain), 3), " ",
			fmt_ci(t2_one_ci[1], t2_one_ci[2], 3)
		),
		t2_range(t2_agreement_one),
		t2_range(t2_priority_one),
		fmt(t2_restricted_one, 3),
		t2_range(t2_omission_one)
	),
	`Five anchors vs one` = c(
		paste0(
			fmt(paired$mean_RMSE_decrease, 3), " ",
			fmt_ci(paired$decrease_ci_low, paired$decrease_ci_high, 3)
		),
		t2_range(t2_agreement$one_to_five_RMSE_decrease),
		t2_range(priority$one_to_five_RMSE_decrease),
		fmt(test_gain$effect, 3),
		t2_range(c(
			influence_gain$minimum_omission_effect,
			influence_gain$maximum_omission_effect
		))
	),
	check.names = FALSE
)

write_table(sensitivity_overview, "02_sensitivity_overview.csv")

# Table S1. Variance components.
write_table(
	data.frame(
		Component = variance$component,
		`Variance [95% CI]` = fmt_est_ci(variance$variance, variance$variance_ci_low, variance$variance_ci_high),
		`Variance fraction (%) [95% CI]` = fmt_est_ci(
			variance$percent, variance$percent_ci_low, variance$percent_ci_high, 1
		),
		check.names = FALSE
	),
	"S1_variance_components.csv"
)

# Table S2. Descriptive molecular QSPR fit.
model_data <- pooled[, c("log_kp", "mw", "logp", "compound_id", "publication_id")]
compound_profile <- model_data[!duplicated(model_data$compound_id), c("mw", "logp")]

for (field in c("mw", "logp")) {
	model_data[[paste0(field, "_z")]] <-
		(model_data[[field]] - mean(compound_profile[[field]])) / sd(compound_profile[[field]])
}

pooled_qspr <- lme4::lmer(
	log_kp ~ mw_z + logp_z + (1 | compound_id) + (1 | publication_id),
	data = model_data, REML = TRUE, na.action = stats::na.fail,
	control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
)

if (length(pooled_qspr@optinfo$conv$lme4$messages) ||
	any(pooled_qspr@optinfo$conv$opt != 0)) {
	stop("The pooled QSPR model did not converge. Inspect pooled_qspr before reporting.")
}

coefficient_terms <- c("(Intercept)", "mw_z", "logp_z")
coefficients_qspr <- summary(pooled_qspr)$coefficients[coefficient_terms, , drop = FALSE]

variance_qspr <- as.data.frame(lme4::VarCorr(pooled_qspr))
variance_qspr <- variance_qspr$vcov[
	match(c("compound_id", "publication_id", "Residual"), variance_qspr$grp)
]

# Nakagawa R2 for this Gaussian model with two random intercepts.
fixed_variance <- stats::var(stats::predict(pooled_qspr, re.form = NA))
random_variance <- sum(variance_qspr[1:2])
total_variance <- fixed_variance + sum(variance_qspr)
marginal_r2 <- fixed_variance / total_variance
conditional_r2 <- (fixed_variance + random_variance) / total_variance

write_table(data.frame(
	Measure = c(
		"Measurements", "Compounds", "Publications",
		"Intercept", "Molecular weight (standardized)", "Crippen logP (standardized)",
		"Compound variance", "Publication variance", "Residual variance",
		"Marginal R2", "Conditional R2"
	),
	Estimate = c(
		as.character(c(nrow(model_data), length(unique(model_data$compound_id)),
			length(unique(model_data$publication_id)))),
		fmt(coefficients_qspr[, "Estimate"]), fmt(variance_qspr),
		fmt(c(marginal_r2, conditional_r2))
	),
	SE = c(rep("", 3), fmt(coefficients_qspr[, "Std. Error"]), rep("", 5)),
	check.names = FALSE
), "S2_pooled_qspr_model_fit.csv")

# Table S3. Experimental metadata availability.
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
write_table(meta_table, "S3_metadata_availability.csv")

# Table S4. Anchor-number comparison.
write_table(data.frame(
	Anchors = matched_cal$anchors,
	`Publications improved` = paste0(matched_cal$n_sources_improved, "/", matched_cal$n_sources),
	`RMSE before` = fmt(matched_cal$mean_source_RMSE_uncalibrated),
	`RMSE after` = fmt(matched_cal$mean_source_RMSE_calibrated),
	`RMSE reduction (%) [95% CI]` = fmt_est_ci(
		matched_cal$RMSE_reduction_percent, matched_cal$reduction_ci_low, matched_cal$reduction_ci_high, 1
	),
	check.names = FALSE
), "S4_calibration_matched.csv")

# Table S5. Agreement-threshold sensitivity.
policy_labels <- ifelse(sensitivity$curation_policy == "priority", "Database priority (primary)",
	paste0("Agreement threshold ", sensitivity$log_tolerance))
t3 <- data.frame(
	Policy = policy_labels, Measurements = sensitivity$n_measurements,
	`Publications: one / five anchors` = paste0(sensitivity$one_anchor_publications, " / ", sensitivity$matched_publications),
	`One-anchor reduction (%) [95% CI]` = fmt_est_ci(sensitivity$one_anchor_reduction_percent, sensitivity$one_anchor_ci_low, sensitivity$one_anchor_ci_high, 1),
	`One-to-five RMSE decrease [95% CI]` = fmt_est_ci(sensitivity$one_to_five_RMSE_decrease, sensitivity$one_to_five_ci_low, sensitivity$one_to_five_ci_high),
	check.names = FALSE
)
write_table(t3, "S5_curation_sensitivity.csv")

# Table S6. Database-priority sensitivity.
write_table(data.frame(
	`Database priority` = priority$database_priority,
	Measurements = priority$n_measurements,
	`Sources with one / five anchors` = paste(priority$one_anchor_publications, priority$matched_publications, sep = " / "),
	`One-anchor reduction (%) [95% CI]` = with(priority,
		fmt_est_ci(one_anchor_reduction_percent, one_anchor_ci_low, one_anchor_ci_high, 1)),
	`One-to-five decrease [95% CI]` = with(priority,
		fmt_est_ci(one_to_five_RMSE_decrease, one_to_five_ci_low, one_to_five_ci_high)),
	check.names = FALSE
), "S6_database_priority_sensitivity.csv")

# Table S7. Test-set sensitivity.
s_test <- test_sensitivity

effect_ci <- vapply(seq_len(nrow(s_test)), function(i) {
	digits <- if (s_test$endpoint[i] == "one_anchor") 1L else 3L

	paste(
		fmt(s_test$effect[i], digits),
		fmt_ci(s_test$ci_low[i], s_test$ci_high[i], digits)
	)
}, character(1))

s_test_print <- data.frame(
	Comparison = ifelse(
		s_test$endpoint == "one_anchor",
		"One-anchor reduction (%)",
		"One-to-five decrease (log units)"
	),
	"Minimum test compounds" = s_test$minimum_test_compounds,
	Sources = s_test$n_sources,
	"Test compounds" = s_test$n_test_compounds,
	"Test measurements" = s_test$n_test_measurements,
	"Effect [95% CI]" = effect_ci,
	check.names = FALSE
)

write_table(s_test_print, "S7_test_set_sensitivity.csv")

# Table S8. Evaluation-publication omissions.
# Script 87 records which orders share the same analysis dataset.
primary_group <- unique(priority$analysis_scenario[priority$is_primary_order])
alternative_groups <- setdiff(unique(priority$analysis_scenario), primary_group)
priority_orders <- setNames(
	gsub(" > ", " then ", priority$database_priority, fixed = TRUE),
	priority$scenario
)

priority_omissions <- read_result("priority_sensitivity/87_priority_influence.csv")

# The primary analysis is already included in the script 83 output.
priority_omissions <- priority_omissions[
	priority_omissions$scenario %in% alternative_groups,
	names(sensitivity_influence),
	drop = FALSE
]

endpoint_labels <- c(
	"One-anchor RMSE reduction (%)" = "One anchor (%)",
	"RMSE decrease from one to five anchors (log units)" = "One to five (log units)"
)

stopifnot(
	length(primary_group) == 1L,
	all(table(
		factor(priority_omissions$scenario, levels = alternative_groups),
		factor(priority_omissions$analysis, levels = names(endpoint_labels))
	) == 1L),
	all(priority_omissions$analysis %in% names(endpoint_labels))
)

omissions <- rbind(sensitivity_influence, priority_omissions)
scenario_labels <- c(
	primary = "Primary priority",
	agreement_0p01 = "Agreement 0.01",
	agreement_0p02 = "Agreement 0.02",
	agreement_0p05 = "Agreement 0.05",
	priority_orders[alternative_groups]
)

source_labels <- pub$publication_label[
	match(omissions$most_influential_source, pub$source_id)
]

one_anchor <- omissions$analysis == "One-anchor RMSE reduction (%)"
format_effect <- function(x) ifelse(one_anchor, fmt(x, 1), fmt(x, 3))

s8 <- data.frame(
	`Curation policy` = unname(scenario_labels[omissions$scenario]),
	Comparison = unname(endpoint_labels[omissions$analysis]),
	`Full effect` = format_effect(omissions$full_effect),
	`Omission range` = paste(
		format_effect(omissions$minimum_omission_effect), "to",
		format_effect(omissions$maximum_omission_effect)
	),
	`CIs including zero` = paste0(
		omissions$n_intervals_including_zero, "/", omissions$n_omissions
	),
	`Largest change source` = source_labels,
	check.names = FALSE
)

if (anyNA(s8)) stop("Incomplete labels or results in the omission table.")

shared_orders <- priority$scenario != priority$analysis_scenario
same_dataset_note <- paste0(
	priority_orders[priority$scenario[shared_orders]],
	" used the same dataset as ",
	priority_orders[priority$analysis_scenario[shared_orders]],
	"."
)

write_table(s8, "S8_evaluation_publication_omissions.csv")

# Table S9. One-versus-five diagnostics.
write_table(data.frame(
	Publication = label(failure_paired$source_id),
	`One-anchor RMSE` = fmt(failure_paired$RMSE_calibrated_one),
	`Five-anchor RMSE` = fmt(failure_paired$RMSE_calibrated_five),
	`RMSE decrease` = fmt(failure_paired$RMSE_decrease),
	`Offset variance: one` = fmt(failure_paired$offset_variance_one),
	`Offset variance: five` = fmt(failure_paired$offset_variance_five),
	`Squared mismatch: one` = fmt(failure_paired$offset_mismatch_squared_one),
	`Squared mismatch: five` = fmt(failure_paired$offset_mismatch_squared_five), check.names = FALSE
), "S9_one_vs_five_diagnostics.csv")

# Table S10. Verified permeability corrections.
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
	"S10_verified_unit_corrections.csv"
)

# Additional files. Export directly without duplicate supplementary-table copies.
write_additional(pub, "Additional_file_1_publication_index.csv")
write_additional(
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
	"Additional_file_2_calibration_by_publication.csv"
)

write_additional(data.frame(
	Publication = label(failure$source_id), Anchors = failure$anchors,
	`Offset variance` = fmt(failure$offset_variance),
	`Squared mean-offset mismatch` = fmt(failure$offset_mismatch_squared),
	`Squared mean test bias` = fmt(failure$correctable_test_bias_squared),
	`Mean MSE change` = fmt(failure$mean_MSE_change),
	`Larger penalty term` = failure$dominant_penalty, check.names = FALSE
), "Additional_file_3_failure_decomposition.csv")

if (!file.copy(pooled_path, file.path(additional_dir, "Additional_file_4_pooled_dataset.csv"), overwrite = TRUE)) {
	stop("Could not export the pooled dataset.")
}

# Remove superseded reporting files after the replacement outputs are complete.
obsolete_tables <- c(
	"03a_curation_sensitivity.csv", "87_database_priority_sensitivity.csv",
	"S1_publication_index.csv", "S2_metadata_availability.csv", "S3_variance_components.csv",
	"S6_calibration_matched.csv", "S7_calibration_by_publication.csv", "S10_failure_decomposition.csv",
	"S11_one_vs_five_diagnostics.csv", "S12_verified_unit_corrections.csv", "S13_pooled_qspr_model_fit.csv"
)
unlink(file.path(table_dir, obsolete_tables))
message("Tables 1–2, Supplementary Tables S1–S10 and Additional files 1–4 written.")
