# Rerun pooled calibration at alternative agreement thresholds in existing scenario folders.
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
	options_run <- parse_options(list("thresholds" = "0.01,0.02,0.05", "temperature-tolerance" = "0.2"))
	thresholds <- suppressWarnings(as.numeric(strsplit(options_run[["thresholds"]], ",", fixed = TRUE)[[1]]))
	temperature_tolerance <- suppressWarnings(as.numeric(options_run[["temperature-tolerance"]]))
	if (!length(thresholds) || any(!is.finite(thresholds)) || any(thresholds < 0) || anyDuplicated(thresholds)) {
		stop("Thresholds must be distinct finite nonnegative numbers.")
	}
	if (!is.finite(temperature_tolerance) || temperature_tolerance < 0) stop("Invalid temperature tolerance.")

	script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
	if (length(script_arg) != 1L) stop("Run with Rscript R/83_run_sensitivity.R")
	script_path <- sub("^--file=", "", script_arg)
	if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
	root <- dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
	old_wd <- setwd(root)
	on.exit(setwd(old_wd))

	primary_dir <- "R/results"
	base_dir <- file.path(primary_dir, "sensitivity")
	threshold_label <- function(x) format(x, scientific = FALSE, trim = TRUE, digits = 15)
	scenario_names <- paste0("agreement_", gsub(".", "p", vapply(thresholds, threshold_label, character(1)), fixed = TRUE))
	if (anyDuplicated(scenario_names)) stop("Threshold labels collide; use fewer decimal places.")
	scenario_dirs <- file.path(base_dir, scenario_names)
	missing_dirs <- scenario_dirs[!dir.exists(scenario_dirs)]
	if (length(missing_dirs)) stop("Required scenario directories are missing: ", paste(missing_dirs, collapse = ", "))

	# Keep scenario writes inside their own directories.
	output_dirs <- c(base_dir, scenario_dirs)
	if (any(normalizePath(output_dirs, winslash = "/", mustWork = TRUE) != file.path(root, output_dirs))) {
		stop("Refusing to write through a redirected sensitivity directory.")
	}
	entries <- unlist(lapply(output_dirs, list.files, all.files = TRUE, full.names = TRUE, no.. = TRUE))
	if (length(entries)) {
		links <- Sys.readlink(entries)
		if (any(!is.na(links) & nzchar(links))) stop("Remove redirected sensitivity outputs before running 83.")
	}

	read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
	read_result <- function(folder, name) read_csv(file.path(folder, name))
	required_primary <- c(
		"80_analysis_settings.csv", "80_input_checksums.csv", "81_analysis_settings.csv",
		"81_pooled_calibration_summary_all.csv", "81_pooled_one_vs_five_anchors.csv",
		"82_analysis_settings.csv", "82_source_influence_summary.csv"
	)
	if (any(!file.exists(file.path(primary_dir, required_primary)))) stop("Complete primary scripts 80, 81 and 82 first.")
	primary_settings <- read_result(primary_dir, "80_analysis_settings.csv")
	if (nrow(primary_settings) != 1L || !isTRUE(primary_settings$analysis_ready) ||
		primary_settings$curation_policy != "priority" || primary_settings$verified_unit_corrections != 9) {
		stop("Main outputs must come from the completed priority analysis with the nine unit corrections.")
	}
	primary_inputs <- read_result(primary_dir, "80_input_checksums.csv")
	raw_hashes <- unname(tools::md5sum(primary_inputs$path))
	if (anyNA(raw_hashes) || !identical(raw_hashes, primary_inputs$md5)) {
		stop("Raw data or curation changed. Rerun primary scripts 80, 81 and 82.")
	}
	for (name in c("81_analysis_settings.csv", "82_analysis_settings.csv")) {
		settings <- read_result(primary_dir, name)
		if (nrow(settings) != 1L || !identical(unname(tools::md5sum(settings$input_file)), settings$input_md5)) {
			stop("Primary validation inputs changed. Rerun scripts 81 and 82.")
		}
	}

	protected <- c(
		list.files(primary_dir, pattern = "^(80_|81_|82_)", full.names = TRUE),
		"data/processed/80_pooled.csv",
		"data/curation/80_source_crosswalk.csv",
		"data/curation/80_compound_crosswalk.csv"
	)
	protected_hash <- tools::md5sum(protected)
	if (anyNA(protected_hash)) stop("Primary data or curation files are missing.")
	verify_protected <- function() {
		if (!identical(tools::md5sum(protected), protected_hash)) {
			stop("A primary output or curated crosswalk changed during the sensitivity run.")
		}
	}

	manifest <- data.frame(
		scenario = scenario_names,
		curation_policy = "agreement",
		log_tolerance = thresholds,
		temperature_tolerance = temperature_tolerance,
		output_directory = scenario_dirs,
		status = "PENDING",
		started = "",
		finished = "",
		stringsAsFactors = FALSE
	)
	manifest_file <- file.path(base_dir, "83_sensitivity_run_status.csv")
	save_manifest <- function() write.csv(manifest, manifest_file, row.names = FALSE)
	comparison_file <- file.path(base_dir, "83_sensitivity_comparison.csv")
	influence_file <- file.path(base_dir, "83_sensitivity_influence.csv")
	if (unlink(c(comparison_file, influence_file), expand = FALSE) != 0L) stop("Cannot remove stale sensitivity summaries.")
	save_manifest()

	rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
	run_step <- function(script, settings, log_file) {
		args <- paste0("--", names(settings), "=", unlist(settings, use.names = FALSE))
		message("  ", basename(script))
		status <- system2(rscript, c("--vanilla", shQuote(script), vapply(args, shQuote, character(1))),
			stdout = log_file, stderr = log_file)
		if (status != 0L) stop("Failed: ", basename(script), ". See ", log_file)
	}

	for (i in seq_along(thresholds)) {
		folder <- scenario_dirs[i]
		manifest$status[i] <- "RUNNING"
		manifest$started[i] <- as.character(Sys.time())
		save_manifest()
		message("[", i, "/", length(thresholds), "] ", scenario_names[i])
		tryCatch({
			run_step("R/80_prepare_pooled_priority_only.R", list(
				"output-dir" = folder,
				"curation-policy" = "agreement",
				"log-tolerance" = threshold_label(thresholds[i]),
				"temperature-tolerance" = temperature_tolerance,
				"scenario" = scenario_names[i],
				"read-only-curation" = "true"
			), file.path(folder, "80.log"))
			run_step("R/81_validate_pooled_calibration.R", list(
				"input-file" = file.path(folder, "80_pooled.csv"),
				"output-dir" = folder
			), file.path(folder, "81.log"))
			run_step("R/82_pooled_source_influence.R", list(
				"input-file" = file.path(folder, "81_pooled_calibration_by_source.csv"),
				"output-dir" = folder
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

	# Combine the primary and sensitivity estimates for the manuscript.
	summarise_folder <- function(folder) {
		settings <- read_result(folder, "80_analysis_settings.csv")
		all <- read_result(folder, "81_pooled_calibration_summary_all.csv")
		one <- all[all$anchors == 1L, , drop = FALSE]
		paired <- read_result(folder, "81_pooled_one_vs_five_anchors.csv")
		if (nrow(settings) != 1L || nrow(one) != 1L || nrow(paired) != 1L) stop("Unexpected summary rows in ", folder)
		data.frame(
			scenario = settings$scenario,
			curation_policy = settings$curation_policy,
			log_tolerance = settings$log_tolerance,
			temperature_tolerance = settings$temperature_tolerance,
			n_measurements = settings$n_measurements,
			n_compounds = settings$n_compounds,
			n_publications = settings$n_publications,
			one_anchor_publications = one$n_sources,
			one_anchor_improved = one$n_sources_improved,
			one_anchor_RMSE_before = one$mean_source_RMSE_uncalibrated,
			one_anchor_RMSE_after = one$mean_source_RMSE_calibrated,
			one_anchor_reduction_percent = one$RMSE_reduction_percent,
			one_anchor_ci_low = one$reduction_ci_low,
			one_anchor_ci_high = one$reduction_ci_high,
			matched_publications = paired$n_sources,
			one_to_five_RMSE_decrease = paired$mean_RMSE_decrease,
			one_to_five_ci_low = paired$decrease_ci_low,
			one_to_five_ci_high = paired$decrease_ci_high
		)
	}
	folders <- c(primary_dir, scenario_dirs)
	comparison <- do.call(rbind, lapply(folders, summarise_folder))
	influence <- do.call(rbind, lapply(seq_along(folders), function(i) {
		cbind(scenario = comparison$scenario[i], read_result(folders[i], "82_source_influence_summary.csv"))
	}))
	write.csv(comparison, comparison_file, row.names = FALSE, na = "")
	write.csv(influence, influence_file, row.names = FALSE, na = "")
	message("Completed ", length(thresholds), " agreement scenarios. Results in ", base_dir, ".")
}

main()
