# 83_run_sensitivity.R
# Place in R/ beside 80, 81 and 82. Run from the project root:
# Rscript R/83_run_sensitivity.R
# Optional: --thresholds=0.01,0.02,0.05 --temperature-tolerance=0.2
# No files are moved. Main outputs and manual crosswalks are not modified.

parse_options <- function(defaults) {
  for (arg in commandArgs(trailingOnly = TRUE)) {
    if (!grepl("^--[^=]+=", arg)) stop("Use --name=value: ", arg)
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    if (!key %in% names(defaults)) stop("Unknown option: --", key)
    defaults[[key]] <- sub("^--[^=]+=", "", arg)
  }
  defaults
}

main <- function() {
  options_run <- parse_options(list(
    "thresholds" = "0.01,0.02,0.05", "temperature-tolerance" = "0.2"
  ))
  raw_thresholds <- strsplit(options_run[["thresholds"]], ",", fixed = TRUE)[[1]]
  thresholds <- suppressWarnings(as.numeric(raw_thresholds))
  temperature_tolerance <- suppressWarnings(as.numeric(options_run[["temperature-tolerance"]]))
  if (!length(thresholds) || any(!is.finite(thresholds)) || any(thresholds < 0) ||
      anyDuplicated(thresholds)) stop("Thresholds must be distinct finite nonnegative numbers.")
  if (!is.finite(temperature_tolerance) || temperature_tolerance < 0) stop("Invalid temperature tolerance.")

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(script_arg) != 1L) stop("Run with Rscript R/83_run_sensitivity.R")
  script_path <- sub("^--file=", "", script_arg)
  # Rscript can encode spaces in its own --file argument as ~+~.
  if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
  script_path <- normalizePath(script_path, mustWork = TRUE)
  root <- dirname(dirname(script_path))
  old_wd <- setwd(root)
  on.exit(setwd(old_wd), add = TRUE)

  scripts <- file.path("R", c(
    "80_prepare_pooled_priority_only.R",
    "81_validate_pooled_calibration.R",
    "82_pooled_source_influence.R"
  ))
  if (any(!file.exists(scripts))) stop("Keep the three analysis scripts beside this runner in R/.")
  for (script in scripts) invisible(parse(file = script))
  packages <- c("readxl", "lme4", "clue")
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) stop("Install required packages: ", paste(missing_packages, collapse = ", "))
  python <- unname(Sys.which("python"))
  if (!nzchar(python)) stop("Activate a Python environment with RDKit first.")
  check <- system2(python, c("-c", shQuote("import rdkit")), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(check, "status")) && attr(check, "status") != 0L) stop("Python/RDKit preflight failed.")

  read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required_primary <- file.path("R/results", c(
    "80_analysis_settings.csv", "80_input_checksums.csv",
    "81_pooled_calibration_summary_all.csv", "81_pooled_one_vs_five_anchors.csv",
    "82_source_influence_summary.csv", "81_analysis_settings.csv", "82_analysis_settings.csv"
  ))
  if (any(!file.exists(required_primary))) {
    stop("Run the updated 80, 81 and 82 with no arguments first, to generate the corrected primary analysis.")
  }
  primary_settings <- read_csv(required_primary[1])
  if (nrow(primary_settings) != 1L || primary_settings$curation_policy != "priority" ||
      !isTRUE(primary_settings$analysis_ready) || primary_settings$verified_unit_corrections != 9) {
    stop("Main outputs must be a completed priority-only preparation with the nine unit corrections.")
  }
  primary_inputs <- read_csv(required_primary[2])
  if (any(!file.exists(primary_inputs$path)) ||
      any(unname(tools::md5sum(primary_inputs$path)) != primary_inputs$md5)) {
    stop("Raw data or curation files changed since primary preparation. Rerun main 80, 81 and 82 first.")
  }

  for (settings_file in required_primary[6:7]) {
    provenance <- read_csv(settings_file)
    if (nrow(provenance) != 1L || !file.exists(provenance$input_file) ||
        unname(tools::md5sum(provenance$input_file)) != provenance$input_md5) {
      stop("Main validation inputs changed. Rerun main 80, 81 and 82 first.")
    }
  }

  # Preserve exact bytes of every main output and both manual crosswalks.
  protected <- c(
    list.files("R/results", pattern = "^(80_|81_|82_)", full.names = TRUE),
    "data/processed/80_pooled.csv",
    "data/curation/80_source_crosswalk.csv",
    "data/curation/80_compound_crosswalk.csv"
  )
  protected <- protected[file.exists(protected) & !dir.exists(protected)]
  protected_hash <- tools::md5sum(protected)
  verify_protected <- function() {
    now <- tools::md5sum(protected)
    if (anyNA(now) || !identical(unname(now), unname(protected_hash))) {
      stop("A main output or curated crosswalk changed during sensitivity analysis.")
    }
  }

  base_dir <- "R/results/sensitivity"
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  label <- function(x) gsub(".", "p", format(x, scientific = FALSE, trim = TRUE, digits = 15), fixed = TRUE)
  scenario_names <- paste0("agreement_", vapply(thresholds, label, character(1)))
  if (anyDuplicated(scenario_names)) stop("Threshold labels collide; use fewer decimal places.")
  manifest <- data.frame(
    scenario = scenario_names, curation_policy = "agreement", log_tolerance = thresholds,
    temperature_tolerance = temperature_tolerance,
    output_directory = file.path(base_dir, scenario_names),
    status = "PENDING", started = "", finished = "", stringsAsFactors = FALSE
  )
  manifest_file <- file.path(base_dir, "83_sensitivity_run_status.csv")
  save_manifest <- function() write.csv(manifest, manifest_file, row.names = FALSE)
  save_manifest()
  # Only these old sensitivity summaries are removed, to prevent stale comparisons after a failed rerun.
  comparison_paths <- file.path(base_dir, c("83_sensitivity_comparison.csv", "83_sensitivity_influence.csv"))
  for (path in comparison_paths) if (file.exists(path)) unlink(path)
  write.csv(data.frame(script = scripts, md5 = unname(tools::md5sum(scripts))),
            file.path(base_dir, "83_script_checksums.csv"), row.names = FALSE)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

  run_step <- function(script, settings, log_file) {
    args <- paste0("--", names(settings), "=", unlist(settings, use.names = FALSE))
    message("  ", basename(script))
    status <- system2(rscript, c("--vanilla", shQuote(script), vapply(args, shQuote, character(1))),
                      stdout = log_file, stderr = log_file)
    if (status != 0L) stop("Failed: ", basename(script), ". See ", log_file)
  }

  for (i in seq_along(thresholds)) {
    folder <- manifest$output_directory[i]
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
    # Prevent output symlinks from redirecting scenario writes into primary outputs.
    target <- normalizePath(folder, winslash = "/", mustWork = TRUE)
    boundary <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
    if (!startsWith(paste0(target, "/"), paste0(boundary, "/")) ||
        identical(target, boundary)) stop("Scenario folder is outside the sensitivity directory.")
    for (entry in list.files(folder, all.files = TRUE, full.names = TRUE, no.. = TRUE)) {
      link <- Sys.readlink(entry)
      if (!is.na(link) && nzchar(link)) stop("Remove redirected scenario output: ", entry)
    }
    manifest$status[i] <- "RUNNING"
    manifest$started[i] <- as.character(Sys.time())
    save_manifest()
    message("[", i, "/", length(thresholds), "] ", scenario_names[i])
    tryCatch({
      run_step(scripts[1], list(
        "output-dir" = folder, "curation-policy" = "agreement",
        "log-tolerance" = format(thresholds[i], scientific = FALSE, digits = 15),
        "temperature-tolerance" = temperature_tolerance,
        "scenario" = scenario_names[i], "read-only-curation" = "true"
      ), file.path(folder, "80.log"))
      pooled_file <- file.path(folder, "80_pooled.csv")
      if (!file.exists(pooled_file)) stop("Preparation did not produce an analysis-ready dataset.")
      run_step(scripts[2], list("input-file" = pooled_file, "output-dir" = folder), file.path(folder, "81.log"))
      run_step(scripts[3], list(
        "input-file" = file.path(folder, "81_pooled_calibration_by_source.csv"), "output-dir" = folder
      ), file.path(folder, "82.log"))
      verify_protected()
    }, error = function(e) {
      manifest$status[i] <<- "FAILED"
      manifest$finished[i] <<- as.character(Sys.time())
      save_manifest()
      stop(conditionMessage(e), call. = FALSE)
    })
    manifest$status[i] <- "SUCCESS"
    manifest$finished[i] <- as.character(Sys.time())
    save_manifest()
  }

  summarise_folder <- function(folder) {
    settings <- read_csv(file.path(folder, "80_analysis_settings.csv"))
    all <- read_csv(file.path(folder, "81_pooled_calibration_summary_all.csv"))
    one <- all[all$anchors == 1, , drop = FALSE]
    paired <- read_csv(file.path(folder, "81_pooled_one_vs_five_anchors.csv"))
    if (nrow(settings) != 1L || nrow(one) != 1L || nrow(paired) != 1L) stop("Unexpected summary rows in ", folder)
    data.frame(
      scenario = settings$scenario, curation_policy = settings$curation_policy,
      log_tolerance = settings$log_tolerance, temperature_tolerance = settings$temperature_tolerance,
      n_measurements = settings$n_measurements, n_compounds = settings$n_compounds,
      n_publications = settings$n_publications, one_anchor_publications = one$n_sources,
      one_anchor_improved = one$n_sources_improved,
      one_anchor_RMSE_before = one$mean_source_RMSE_uncalibrated,
      one_anchor_RMSE_after = one$mean_source_RMSE_calibrated,
      one_anchor_reduction_percent = one$RMSE_reduction_percent,
      one_anchor_ci_low = one$reduction_ci_low, one_anchor_ci_high = one$reduction_ci_high,
      matched_publications = paired$n_sources,
      one_to_five_RMSE_decrease = paired$mean_RMSE_decrease,
      one_to_five_ci_low = paired$decrease_ci_low, one_to_five_ci_high = paired$decrease_ci_high
    )
  }
  folders <- c("R/results", manifest$output_directory)
  comparison <- do.call(rbind, lapply(folders, summarise_folder))
  write.csv(comparison, comparison_paths[1], row.names = FALSE, na = "")
  influence <- do.call(rbind, lapply(folders, function(folder) {
    settings <- read_csv(file.path(folder, "80_analysis_settings.csv"))
    cbind(scenario = settings$scenario, read_csv(file.path(folder, "82_source_influence_summary.csv")))
  }))
  write.csv(influence, comparison_paths[2], row.names = FALSE, na = "")
  verify_protected()
  print(comparison, row.names = FALSE)
  message("Complete. Scenario outputs and comparison tables: ", base_dir)
}

main()
