# 99a_make_tables.R
# Run from the project root after 80-85:
# Rscript R/99a_make_tables.R
# Formats saved results only. No models, bootstraps or curation are rerun.

result_dir <- "R/results"
table_dir <- "R/tables"
sensitivity_dir <- file.path(result_dir, "sensitivity")

read_checked <- function(path, required = character()) {
  if (!file.exists(path)) stop("Missing input: ", path)
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
index <- list()
write_table <- function(x, filename, caption, note = "") {
  write.csv(x, file.path(table_dir, filename), row.names = FALSE, na = "")
  index[[length(index) + 1L]] <<- data.frame(
    file = filename, caption = caption, note = note, rows = nrow(x), stringsAsFactors = FALSE
  )
}

# 1. Inputs from the pooled pipeline only.
pooled <- read_checked("data/processed/80_pooled.csv", c(
  "record_id", "publication_id", "compound_id", "dataset_origin"
))
structure <- read_result("85_dataset_structure.csv", c(
  "dataset", "raw_measurements", "retained_measurements", "compounds", "publications"
))
metadata <- read_result("85_metadata_availability.csv", c(
  "dataset", "metadata", "field_provided", "n_reported", "n_records", "percent_reported"
))
variance <- read_result("85_variance_components.csv", c(
  "component", "variance", "variance_ci_low", "variance_ci_high",
  "percent", "percent_ci_low", "percent_ci_high", "n_boot_requested", "n_boot_usable"
))
cal_cols <- c("anchors", "n_sources", "n_sources_improved",
              "mean_source_RMSE_uncalibrated", "mean_source_RMSE_calibrated",
              "RMSE_reduction_percent", "reduction_ci_low", "reduction_ci_high")
all_cal <- read_result("81_pooled_calibration_summary_all.csv", cal_cols)
matched_cal <- read_result("81_pooled_calibration_summary_matched.csv", cal_cols)
paired <- read_result("81_pooled_one_vs_five_anchors.csv", c(
  "n_sources", "mean_RMSE_one_anchor", "mean_RMSE_five_anchors",
  "mean_RMSE_decrease", "decrease_ci_low", "decrease_ci_high"
))
by_source <- read_result("81_pooled_calibration_by_source.csv", c(
  "source_id", "anchors", "RMSE_uncalibrated", "RMSE_calibrated",
  "RMSE_reduction_percent", "n_available_anchors", "n_test_compounds"
))
eligibility <- read_result("81_pooled_source_eligibility.csv", c(
  "source_id", "n_measurements", "n_compounds", "n_anchors",
  "n_test_compounds", "n_test_measurements"
))
fits <- read_result("81_pooled_fit_diagnostics.csv", c("source_id", "success", "singular"))
influence <- read_result("82_source_influence_summary.csv", c(
  "analysis", "full_effect", "minimum_omission_effect", "maximum_omission_effect",
  "most_influential_source", "n_intervals_including_zero", "n_omissions"
))
sensitivity <- read_checked(file.path(sensitivity_dir, "83_sensitivity_comparison.csv"), c(
  "scenario", "curation_policy", "log_tolerance", "temperature_tolerance", "n_measurements",
  "one_anchor_publications", "matched_publications", "one_anchor_reduction_percent",
  "one_anchor_ci_low", "one_anchor_ci_high", "one_to_five_RMSE_decrease",
  "one_to_five_ci_low", "one_to_five_ci_high"
))
sensitivity_influence <- read_checked(file.path(sensitivity_dir, "83_sensitivity_influence.csv"), c("scenario", names(influence)))
run_status <- read_checked(file.path(sensitivity_dir, "83_sensitivity_run_status.csv"), c("scenario", "status"))
failure <- read_result("84_failure_source_summary.csv", c(
  "source_id", "anchors", "available_anchors", "test_compounds", "draws_improved_percent",
  "RMSE_uncalibrated", "RMSE_calibrated", "RMSE_reduction_percent",
  "test_mean_residual", "offset_mean", "offset_sd", "mean_offset_wrong_direction",
  "MSE_uncalibrated", "offset_variance", "offset_mismatch_squared",
  "correctable_test_bias_squared", "mean_MSE_change", "dominant_penalty"
))
failure_overview <- read_result("84_failure_overview.csv", c("anchors", "publications"))
failure_paired <- read_result("84_one_vs_five_diagnostics.csv", c(
  "source_id", "RMSE_calibrated_one", "RMSE_calibrated_five", "RMSE_decrease",
  "offset_variance_one", "offset_variance_five",
  "offset_mismatch_squared_one", "offset_mismatch_squared_five"
))
corrections <- read_result("80_verified_unit_corrections.csv", c(
  "record_id", "previous_log_kp", "corrected_log_kp", "original_kp_cm_h", "evidence"
))
coverage <- read_result("80_descriptor_coverage.csv", c("variable", "n_available", "n_records"))

# 2. Check that the primary outputs agree before writing tables.
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
f1 <- failure[failure$anchors == 1, , drop = FALSE]
j <- match(paste(failure$source_id, failure$anchors), paste(by_source$source_id, by_source$anchors))
if (anyNA(j) || nrow(failure) != nrow(by_source) ||
    !close_enough(failure$RMSE_calibrated, by_source$RMSE_calibrated[j])) {
  stop("Publication-level failure results disagree with 81. Rerun 84.")
}
if (nrow(f1) != one$n_sources ||
    !close_enough(mean(f1$RMSE_calibrated), one$mean_source_RMSE_calibrated)) {
  stop("Failure diagnostics disagree with primary calibration. Rerun 84.")
}
if (!close_enough(sum(variance$percent), 100)) stop("Variance percentages do not sum to 100.")
for (stage in c("81", "85")) {
  settings <- read_result(paste0(stage, "_analysis_settings.csv"))
  hash <- if (stage == "81") settings$input_md5 else settings$md5
  if (length(hash) != 1L || is.na(hash) || hash != unname(tools::md5sum("data/processed/80_pooled.csv"))) {
    stop("Input fingerprint differs for stage ", stage, ". Rerun that stage.")
  }
}

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Stable publication labels for both tables and future figures.
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
calibration_table <- function(x, cohort) {
  data.frame(
    Cohort = rep(cohort, nrow(x)), Anchors = x$anchors,
    Publications = x$n_sources,
    `Publications improved` = paste0(x$n_sources_improved, "/", x$n_sources),
    `RMSE before` = fmt(x$mean_source_RMSE_uncalibrated),
    `RMSE after` = fmt(x$mean_source_RMSE_calibrated),
    `RMSE reduction (%) [95% CI]` = fmt_est_ci(x$RMSE_reduction_percent, x$reduction_ci_low, x$reduction_ci_high, 1),
    check.names = FALSE
  )
}
cal_note <- paste(
  "RMSE is in log10 Kp units, with Kp in cm/h. Metrics are averaged across anchor draws within",
  "publication and then across publications with equal weight. Percentage reduction is the ratio",
  "of mean RMSE values, not the mean of publication-specific percentage reductions.",
  "Confidence intervals resample publications without refitting models."
)

# 3. Main tables.
t1 <- structure
names(t1) <- c("Dataset", "Raw records", "Retained measurements", "Compounds", "Publications")
write_table(t1, "01_dataset_structure.csv", "Dataset provenance and retained pooled sample.",
  paste("Raw records precede eligibility filtering. Publication-priority selection was SkinPiX, HuSkinDB,",
        "then Cheruvu. Compounds can occur in more than one database contribution; compound counts",
        "are not additive. The pooled total counts unique compound and publication identities."))

rows <- rbind(
  calibration_table(one, "All eligible publications"),
  calibration_table(matched_cal[matched_cal$anchors %in% c(1, 5), , drop = FALSE], "Fixed five-anchor-eligible subset")
)
write_table(rows, "02a_primary_calibration.csv", "Primary pooled calibration performance.", cal_note)
write_table(data.frame(
  Publications = paired$n_sources,
  `One-anchor RMSE` = fmt(paired$mean_RMSE_one_anchor),
  `Five-anchor RMSE` = fmt(paired$mean_RMSE_five_anchors),
  `RMSE decrease [95% CI]` = fmt_est_ci(paired$mean_RMSE_decrease, paired$decrease_ci_low, paired$decrease_ci_high),
  check.names = FALSE
), "02b_one_vs_five_anchors.csv", "Paired comparison of one and five anchors.",
  "The same publications and test measurements are compared. Positive differences indicate lower RMSE with five anchors. Units are log10 Kp; this is an absolute difference, not a percentage.")

policy_labels <- ifelse(sensitivity$curation_policy == "priority", "Database priority (primary)",
                        paste0("Agreement threshold ", sensitivity$log_tolerance))
t3 <- data.frame(
  Policy = policy_labels, Measurements = sensitivity$n_measurements,
  `Publications: one / five anchors` = paste0(sensitivity$one_anchor_publications, " / ", sensitivity$matched_publications),
  `One-anchor reduction (%) [95% CI]` = fmt_est_ci(sensitivity$one_anchor_reduction_percent, sensitivity$one_anchor_ci_low, sensitivity$one_anchor_ci_high, 1),
  `One-to-five RMSE decrease [95% CI]` = fmt_est_ci(sensitivity$one_to_five_RMSE_decrease, sensitivity$one_to_five_ci_low, sensitivity$one_to_five_ci_high),
  check.names = FALSE
)
write_table(t3, "03a_curation_sensitivity.csv", "Sensitivity to measurement-curation policy.",
  paste("Agreement thresholds are in log10 units. Recorded temperature tolerances (degrees Celsius):",
        paste(unique(sensitivity$temperature_tolerance[sensitivity$curation_policy == "agreement"]), collapse = ", "),
        "Evaluation publications and test measurements can differ between policies.",
        "The one-to-five comparison uses each policy's fixed five-anchor-eligible subset."))

influence_table <- function(x, include_scenario = FALSE) {
  d <- ifelse(grepl("reduction", x$analysis, fixed = TRUE), 1L, 3L)
  # Endpoints have different units, so each is labelled explicitly.
  effects <- vapply(seq_len(nrow(x)), function(i) fmt(x$full_effect[i], d[i]), character(1))
  ranges <- vapply(seq_len(nrow(x)), function(i) paste0(
    fmt(x$minimum_omission_effect[i], d[i]), " to ", fmt(x$maximum_omission_effect[i], d[i])
  ), character(1))
  out <- data.frame(
    Endpoint = x$analysis, `Full effect` = effects, `Omission estimate range` = ranges,
    `Intervals including zero` = paste0(x$n_intervals_including_zero, "/", x$n_omissions),
    `Most influential publication` = label(x$most_influential_source), check.names = FALSE
  )
  if (include_scenario) out <- cbind(Scenario = x$scenario, out)
  out
}
influence_note <- paste("One evaluation publication is omitted from the summary at a time; models are not refitted.",
                        "Omission ranges contain point estimates, not confidence bounds.",
                        "Publication labels are defined in Supplementary Table S1.")
write_table(influence_table(influence), "03b_publication_influence.csv",
            "Influence of individual evaluation publications in the primary analysis.", influence_note)

# 4. Supplementary tables.
pub <- merge(publication_key, eligibility, by = "source_id", all.x = TRUE, sort = FALSE)
if (anyNA(pub$n_measurements)) stop("Incomplete publication eligibility index.")
pub$dataset_origins <- vapply(pub$source_id, function(id) {
  paste(sort(unique(pooled$dataset_origin[pooled$publication_id == id])), collapse = "; ")
}, character(1))
pub <- pub[order(pub$publication_label), ]
write_table(pub, "S1_publication_index.csv", "Publication identities, provenance and calibration eligibility.",
            "Publication labels are assigned from sorted primary publication identifiers and are shared with the figures.")

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
write_table(meta_table, "S2_metadata_availability.csv", "Experimental metadata availability among retained measurements.",
  "Counts describe reported fields, not harmonised protocol equivalence. Pooled percentages use all retained measurements as the denominator, including contributions without that field.")

write_table(data.frame(
  Component = variance$component,
  `Variance [95% CI]` = fmt_est_ci(variance$variance, variance$variance_ci_low, variance$variance_ci_high),
  `Variance fraction (%) [95% CI]` = fmt_est_ci(variance$percent, variance$percent_ci_low, variance$percent_ci_high, 1),
  `Bootstrap replicates usable / requested` = paste0(variance$n_boot_usable, " / ", variance$n_boot_requested),
  check.names = FALSE
), "S3_variance_components.csv", "Pooled variance decomposition.",
  "REML intercept-only model with compound and publication random intercepts. Molecular descriptors and temperature are omitted. Confidence intervals use a parametric bootstrap with both random effects simulated anew.")

fits <- cbind(Publication = label(fits$source_id), fits)
write_table(fits, "S4_fit_diagnostics.csv", "Diagnostics for publication-held-out calibration fits.")
write_table(calibration_table(all_cal, "All eligible publications at each count"),
            "S5_calibration_all_eligible.csv", "Calibration across all eligible publications at each anchor count.",
            paste(cal_note, "Publication composition changes across anchor counts."))
write_table(calibration_table(matched_cal, "Fixed five-anchor-eligible subset"),
            "S6_calibration_matched.csv", "Anchor-number comparison within a fixed publication subset.", cal_note)

write_table(data.frame(
  Publication = label(by_source$source_id), Anchors = by_source$anchors,
  `Available anchors` = by_source$n_available_anchors, `Test compounds` = by_source$n_test_compounds,
  `RMSE before` = fmt(by_source$RMSE_uncalibrated), `RMSE after` = fmt(by_source$RMSE_calibrated),
  `RMSE reduction (%)` = fmt(by_source$RMSE_reduction_percent, 1), check.names = FALSE
), "S7_calibration_by_publication.csv", "Calibration performance by publication and anchor count.",
  "Values average the sampled panels within each publication. Positive reductions indicate lower RMSE; publication-specific confidence intervals are not calculated.")
write_table(influence_table(sensitivity_influence, TRUE), "S8_sensitivity_influence.csv",
            "Evaluation-publication influence under each curation policy.", influence_note)

write_table(data.frame(
  Publication = label(failure$source_id), Anchors = failure$anchors,
  `Available anchors` = failure$available_anchors, `Test compounds` = failure$test_compounds,
  `RMSE before` = fmt(failure$RMSE_uncalibrated), `RMSE after` = fmt(failure$RMSE_calibrated),
  `Draws improved (%)` = fmt(failure$draws_improved_percent, 1), check.names = FALSE
), "S9_failure_performance.csv", "Prediction performance across sampled anchor panels.",
  "The proportion of draws improving is conditional on the saved panel sampling and is not a prospective success probability.")
write_table(data.frame(
  Publication = label(failure$source_id), Anchors = failure$anchors,
  `Offset variance` = fmt(failure$offset_variance),
  `Squared mean-offset mismatch` = fmt(failure$offset_mismatch_squared),
  `Squared mean test bias` = fmt(failure$correctable_test_bias_squared),
  `Mean MSE change` = fmt(failure$mean_MSE_change),
  `Larger penalty term` = failure$dominant_penalty, check.names = FALSE
), "S10_failure_decomposition.csv", "Decomposition of mean test MSE change across anchor panels.",
  "Mean MSE change equals offset variance plus squared mean-offset mismatch minus squared mean test bias. Offset variance uses denominator n. These are retrospective algebraic terms, not identified experimental causes.")
write_table(data.frame(
  Publication = label(failure_paired$source_id),
  `One-anchor RMSE` = fmt(failure_paired$RMSE_calibrated_one),
  `Five-anchor RMSE` = fmt(failure_paired$RMSE_calibrated_five),
  `RMSE decrease` = fmt(failure_paired$RMSE_decrease),
  `Offset variance: one` = fmt(failure_paired$offset_variance_one),
  `Offset variance: five` = fmt(failure_paired$offset_variance_five),
  `Squared mismatch: one` = fmt(failure_paired$offset_mismatch_squared_one),
  `Squared mismatch: five` = fmt(failure_paired$offset_mismatch_squared_five), check.names = FALSE
), "S11_one_vs_five_diagnostics.csv", "Changes in error, offset variability and mismatch from one to five anchors.",
  "Comparisons use the same publications and test measurements. Increasing panel size changes both anchor averaging and shrinkage.")

write_table(data.frame(
  Record = corrections$record_id,
  `Previous log10 Kp (cm/h)` = fmt(corrections$previous_log_kp, 6),
  `Corrected log10 Kp (cm/h)` = fmt(corrections$corrected_log_kp, 6),
  `Original Kp (cm/h)` = formatC(corrections$original_kp_cm_h, format = "e", digits = 3),
  Evidence = corrections$evidence, check.names = FALSE
), "S12_verified_unit_corrections.csv", "Verified permeability corrections from original publications.",
  "Previous values refer to the processing output before correction, not the original database column. Raw files were preserved.")
write_table(data.frame(
  Variable = coverage$variable, Available = coverage$n_available, Measurements = coverage$n_records,
  `Available (%)` = fmt(ifelse(coverage$n_records > 0, 100*coverage$n_available/coverage$n_records, NA_real_), 1),
  check.names = FALSE
), "S13_descriptor_coverage.csv", "Descriptor and temperature availability in the retained pooled dataset.",
  "Molecular descriptors were calculated consistently from curated structures; availability does not imply inclusion in the primary model.")
write_table(failure_overview, "S14_failure_overview.csv", "Overview of calibration failures by anchor count.",
  "Eligible publication sets can differ across anchor counts. Wrong-direction counts overlap with the mismatch/variability classifications.")

# Sidecar captions/notes for use when writing the revised Rmd.
write.csv(do.call(rbind, index), file.path(table_dir, "table_index.csv"), row.names = FALSE, na = "")
cat("\nPooled manuscript tables written to ", table_dir, ".\n", sep = "")
cat("Captions and footnotes: R/tables/table_index.csv\n")
cat("Update the Rmd to these new filenames; old table imports are not compatible.\n")
