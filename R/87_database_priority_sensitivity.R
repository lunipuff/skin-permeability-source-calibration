# Compare all six database priority orders in existing scenario directories.
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
	rscript_default <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
	options_run <- parse_options(list("run-influence" = "true", "rscript" = rscript_default))
	if (!options_run[["run-influence"]] %in% c("true", "false")) stop("--run-influence must be true or false.")
	run_influence <- options_run[["run-influence"]] == "true"
	rscript <- normalizePath(options_run[["rscript"]], winslash = "/", mustWork = TRUE)

	script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
	if (length(script_arg) != 1L) stop("Run with Rscript R/87_database_priority_sensitivity.R")
	script_path <- sub("^--file=", "", script_arg)
	if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
	root <- dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
	old_wd <- setwd(root)
	on.exit(setwd(old_wd))

	orders <- list(
		c("SkinPiX", "HuSkinDB", "Cheruvu"),
		c("SkinPiX", "Cheruvu", "HuSkinDB"),
		c("HuSkinDB", "SkinPiX", "Cheruvu"),
		c("HuSkinDB", "Cheruvu", "SkinPiX"),
		c("Cheruvu", "SkinPiX", "HuSkinDB"),
		c("Cheruvu", "HuSkinDB", "SkinPiX")
	)
	scenarios <- vapply(orders, function(x) paste0("priority_", paste(tolower(x), collapse = "_")), character(1))
	primary_dir <- "R/results"
	base_dir <- file.path(primary_dir, "priority_sensitivity")
	scenario_dirs <- file.path(base_dir, scenarios)
	table_file <- "R/tables/87_database_priority_sensitivity.csv"

	# Check output locations before starting the runs.
	output_dirs <- c(base_dir, scenario_dirs, dirname(table_file))
	missing_dirs <- output_dirs[!dir.exists(output_dirs)]
	if (length(missing_dirs)) stop("Required output directories are missing: ", paste(missing_dirs, collapse = ", "))
	if (any(normalizePath(output_dirs, winslash = "/", mustWork = TRUE) != file.path(root, output_dirs))) {
		stop("Refusing to write through a redirected priority output directory.")
	}
	entries <- c(table_file, unlist(lapply(
		c(base_dir, scenario_dirs), list.files, all.files = TRUE, full.names = TRUE, no.. = TRUE
	)))
	links <- Sys.readlink(entries)
	if (any(!is.na(links) & nzchar(links))) stop("Remove redirected priority outputs before running 87.")

	read_csv <- function(path) {
		if (!file.exists(path)) stop("Missing input: ", path)
		read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
	}
	read_result <- function(folder, name) read_csv(file.path(folder, name))
	write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")

	# Check that the primary data and calibration use the current inputs.
	primary_input <- "data/processed/80_pooled.csv"
	prep <- read_result(primary_dir, "80_analysis_settings.csv")
	cal_settings <- read_result(primary_dir, "81_analysis_settings.csv")
	if (nrow(prep) != 1L || !isTRUE(prep$analysis_ready) ||
		!isTRUE(prep$curation_policy == "priority") || !isTRUE(prep$verified_unit_corrections == 9)) {
		stop("Complete the corrected primary priority analysis first.")
	}
	if (nrow(cal_settings) != 1L || !isTRUE(cal_settings$panel_method == "exhaustive") ||
		!file.exists(primary_input) ||
		!identical(cal_settings$input_md5, unname(tools::md5sum(primary_input)))) {
		stop("Complete exhaustive 81 on the current primary pooled dataset first.")
	}
	inputs <- read_result(primary_dir, "80_input_checksums.csv")
	if (!nrow(inputs) || !all(c("path", "md5") %in% names(inputs))) stop("Invalid primary input checksums.")
	input_hashes <- unname(tools::md5sum(inputs$path))
	if (anyNA(input_hashes) || !identical(input_hashes, inputs$md5)) {
		stop("Raw inputs or curation changed since 80. Rerun main 80 and 81 first.")
	}
	protected <- unique(c(
		list.files(primary_dir, full.names = TRUE), primary_input, inputs$path,
		list.files("data/curation", recursive = TRUE, full.names = TRUE)
	))
	protected <- protected[!dir.exists(protected)]
	protected_hash <- tools::md5sum(protected)
	verify_protected <- function() {
		if (!identical(tools::md5sum(protected), protected_hash)) {
			stop("A primary output, raw input or curation file changed during 87.")
		}
	}
	on.exit(verify_protected(), add = TRUE, after = FALSE)

	comparison_file <- file.path(base_dir, "87_priority_comparison.csv")
	influence_file <- file.path(base_dir, "87_priority_influence.csv")
	if (unlink(c(comparison_file, influence_file, table_file), expand = FALSE) != 0L) {
		stop("Cannot remove stale priority summaries.")
	}
	manifest <- data.frame(
		scenario = scenarios,
		database_priority = vapply(orders, paste, character(1), collapse = " > "),
		preparation_directory = scenario_dirs,
		dataset_fingerprint = "",
		analysis_scenario = "",
		analysis_directory = "",
		reused_identical_dataset = FALSE,
		status = "PENDING",
		error = "",
		stringsAsFactors = FALSE
	)
	manifest_file <- file.path(base_dir, "87_run_status.csv")
	save_manifest <- function() write_csv(manifest, manifest_file)
	save_manifest()

	run_step <- function(script, arguments, log_file) {
		args <- paste0("--", names(arguments), "=", unlist(arguments, use.names = FALSE))
		status <- system2(rscript, c("--vanilla", shQuote(script), vapply(args, shQuote, character(1))),
			stdout = log_file, stderr = log_file)
		if (status != 0L) stop("Failed ", basename(script), ". See ", log_file)
	}
	canonical_data <- function(path) {
		x <- read_csv(path)
		if (!nrow(x) || !"record_id" %in% names(x) || anyNA(x$record_id) || anyDuplicated(x$record_id)) {
			stop("Invalid pooled dataset: ", path)
		}
		x <- x[order(x$record_id, method = "radix"), sort(names(x)), drop = FALSE]
		rownames(x) <- NULL
		x
	}
	fingerprint <- function(x) {
		temporary <- tempfile(fileext = ".csv")
		on.exit(unlink(temporary))
		write_csv(x, temporary)
		unname(tools::md5sum(temporary))
	}
	primary_data <- canonical_data(primary_input)
	unique_data <- list()

	# Reuse calibration only when pooled datasets are identical.
	for (i in seq_along(orders)) {
		folder <- scenario_dirs[i]
		manifest$status[i] <- "RUNNING"
		save_manifest()
		message("[", i, "/", length(orders), "] ", manifest$database_priority[i])
		tryCatch({
			run_step("R/80_prepare_pooled_priority_only.R", list(
				"output-dir" = folder,
				"curation-policy" = "priority",
				"database-priority" = paste(orders[[i]], collapse = ","),
				"scenario" = scenarios[i],
				"read-only-curation" = "true"
			), file.path(folder, "80.log"))
			pooled_file <- file.path(folder, "80_pooled.csv")
			prepared <- read_result(folder, "80_analysis_settings.csv")
			if (nrow(prepared) != 1L || !isTRUE(prepared$analysis_ready) ||
				!isTRUE(prepared$database_priority == manifest$database_priority[i])) {
				stop("Scenario preparation did not apply the requested priority.")
			}
			data_i <- canonical_data(pooled_file)
			hash <- fingerprint(data_i)
			manifest$dataset_fingerprint[i] <- hash
			if (manifest$database_priority[i] == prep$database_priority && !identical(data_i, primary_data)) {
				stop("The primary-order dataset differs from the main preparation. Regenerate main 80 and 81.")
			}

			owner <- match(hash, manifest$dataset_fingerprint[seq_len(i - 1L)])
			if (!is.na(owner)) {
				analysis_scenario <- manifest$analysis_scenario[owner]
				if (!identical(data_i, unique_data[[analysis_scenario]])) {
					stop("Dataset fingerprint collision. Inspect the prepared data.")
				}
				manifest$analysis_scenario[i] <- analysis_scenario
				manifest$analysis_directory[i] <- manifest$analysis_directory[owner]
				manifest$reused_identical_dataset[i] <- TRUE
				message("  Identical dataset. Using ", analysis_scenario)
			} else {
				manifest$analysis_scenario[i] <- scenarios[i]
				manifest$analysis_directory[i] <- folder
				unique_data[[scenarios[i]]] <- data_i
				run_step("R/81_validate_pooled_calibration.R", list(
					"input-file" = pooled_file, "output-dir" = folder
				), file.path(folder, "81.log"))
				fitted <- read_result(folder, "81_pooled_fit_diagnostics.csv")
				method <- read_result(folder, "81_analysis_settings.csv")
				if (!nrow(fitted) || !"success" %in% names(fitted) || !isTRUE(all(fitted$success)) ||
					nrow(method) != 1L || !isTRUE(method$panel_method == "exhaustive")) {
					stop("Not all eligible folds completed with exhaustive panels. Inspect 81.log.")
				}
				if (run_influence) {
					run_step("R/82_pooled_source_influence.R", list(
						"input-file" = file.path(folder, "81_pooled_calibration_by_source.csv"),
						"output-dir" = folder
					), file.path(folder, "82.log"))
				}
			}
			verify_protected()
			manifest$status[i] <- "SUCCESS"
			save_manifest()
		}, error = function(e) {
			manifest$status[i] <<- "FAILED"
			manifest$error[i] <<- conditionMessage(e)
			save_manifest()
			stop(conditionMessage(e), call. = FALSE)
		})
	}

	# Combine scenario estimates and format the supplementary table.
	comparison <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
		folder <- manifest$analysis_directory[i]
		prepared <- read_result(scenario_dirs[i], "80_analysis_settings.csv")
		all <- read_result(folder, "81_pooled_calibration_summary_all.csv")
		one <- all[all$anchors == 1L, , drop = FALSE]
		pair <- read_result(folder, "81_pooled_one_vs_five_anchors.csv")
		if (nrow(one) != 1L || nrow(pair) != 1L) stop("Unexpected result rows in ", folder)
		data.frame(
			scenario = scenarios[i],
			database_priority = manifest$database_priority[i],
			is_primary_order = manifest$database_priority[i] == prep$database_priority,
			dataset_fingerprint = manifest$dataset_fingerprint[i],
			analysis_scenario = manifest$analysis_scenario[i],
			analysis_directory = folder,
			n_measurements = prepared$n_measurements,
			n_compounds = prepared$n_compounds,
			n_publications = prepared$n_publications,
			one_anchor_publications = one$n_sources,
			one_anchor_improved = one$n_sources_improved,
			one_anchor_RMSE_before = one$mean_source_RMSE_uncalibrated,
			one_anchor_RMSE_after = one$mean_source_RMSE_calibrated,
			one_anchor_reduction_percent = one$RMSE_reduction_percent,
			one_anchor_ci_low = one$reduction_ci_low,
			one_anchor_ci_high = one$reduction_ci_high,
			matched_publications = pair$n_sources,
			one_to_five_RMSE_decrease = pair$mean_RMSE_decrease,
			one_to_five_ci_low = pair$decrease_ci_low,
			one_to_five_ci_high = pair$decrease_ci_high,
			stringsAsFactors = FALSE
		)
	}))
	if (run_influence) {
		influence <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
			cbind(
				scenario = scenarios[i],
				database_priority = manifest$database_priority[i],
				read_result(manifest$analysis_directory[i], "82_source_influence_summary.csv")
			)
		}))
		write_csv(influence, influence_file)
	}

	fmt_ci <- function(x, low, high, digits) {
		paste0(
			formatC(x, format = "f", digits = digits), " [",
			formatC(low, format = "f", digits = digits), ", ",
			formatC(high, format = "f", digits = digits), "]"
		)
	}
	formatted <- data.frame(
		"Database priority" = comparison$database_priority,
		"Measurements" = comparison$n_measurements,
		"Sources with one / five anchors" = paste(comparison$one_anchor_publications, comparison$matched_publications, sep = " / "),
		"One-anchor reduction (%) [95% CI]" = with(comparison,
			fmt_ci(one_anchor_reduction_percent, one_anchor_ci_low, one_anchor_ci_high, 1L)),
		"One-to-five decrease [95% CI]" = with(comparison,
			fmt_ci(one_to_five_RMSE_decrease, one_to_five_ci_low, one_to_five_ci_high, 3L)),
		check.names = FALSE
	)
	write_csv(comparison, comparison_file)
	write_csv(formatted, table_file)
	message("Completed all ", length(orders), " priority orders across ", length(unique_data), " distinct datasets.")
}

main()
