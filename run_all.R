#!/usr/bin/env Rscript

# Run with Rscript run_all.R.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript run_all.R.", call. = FALSE)
runner <- sub("^--file=", "", script_arg)
if (!file.exists(runner)) runner <- gsub("~+~", " ", runner, fixed = TRUE)
root <- dirname(normalizePath(runner, winslash = "/", mustWork = TRUE))
setwd(root)

scripts <- file.path("R", c(
	"80_prepare_pooled_priority_only.R",
	"81_validate_pooled_calibration.R",
	"82_pooled_source_influence.R",
	"83_run_sensitivity.R",
	"84_diagnose_pooled_calibration.R",
	"85_describe_pooled_variation.R",
	"86_test_compound_sensitivity.R",
	"87_database_priority_sensitivity.R",
	"99a_make_tables.R",
	"99b_make_figures.R"
))

# Check inputs and dependencies before clearing generated outputs.
required <- c(
	scripts,
	file.path("data/raw", c("cheruvu.xlsx", "huskinDB.csv", "skinpix.xlsx")),
	file.path("data/curation", c("80_source_crosswalk.csv", "80_compound_crosswalk.csv"))
)
missing <- required[!file.exists(required) | dir.exists(required)]
if (length(missing)) stop("Missing required files:\n", paste(missing, collapse = "\n"), call. = FALSE)
for (script in scripts) invisible(parse(file = script))

packages <- c("readxl", "lme4", "clue", "ggplot2", "patchwork")
installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
if (any(!installed)) stop("Install R packages: ", paste(packages[!installed], collapse = ", "), call. = FALSE)
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) stop("Rscript was not found: ", rscript, call. = FALSE)
python <- unname(Sys.which("python"))
if (!nzchar(python)) stop("Activate the Python environment containing RDKit.", call. = FALSE)
python_info <- suppressWarnings(system2(
	python, c("-c", shQuote("import rdkit; print(rdkit.__version__)")), stdout = TRUE, stderr = TRUE
))
if (!is.null(attr(python_info, "status"))) {
	stop("Python could not load RDKit:\n", paste(python_info, collapse = "\n"), call. = FALSE)
}

# Check every cleanup path before deleting anything.
output_dirs <- c("data/processed", "R/results", "R/tables", "R/figures", "R/additional_files")
for (path in output_dirs) {
	parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
	link <- Sys.readlink(path)
	if (parent != file.path(root, dirname(path)) || (!is.na(link) && nzchar(link))) {
		stop("Refusing to clear a redirected output path: ", path, call. = FALSE)
	}
	if (file.exists(path) && !dir.exists(path)) stop("Expected a directory: ", path, call. = FALSE)
}

# Clear generated outputs and keep the manually prepared graphical abstract.
for (path in output_dirs) {
	if (!dir.exists(path) && !dir.create(path)) stop("Cannot create: ", path, call. = FALSE)
	entries <- list.files(path, full.names = TRUE, all.files = TRUE, no.. = TRUE)
	entries <- setdiff(entries, "R/figures/00_GA.pdf")
	if (length(entries) && unlink(entries, recursive = TRUE, expand = FALSE) != 0L) {
		stop("Could not clear: ", path, call. = FALSE)
	}
}

log_dir <- "R/results/logs"
run_dirs <- c(
	log_dir,
	file.path("R/results/sensitivity", c("agreement_0p01", "agreement_0p02", "agreement_0p05")),
	file.path("R/results/priority_sensitivity", c(
		"priority_skinpix_huskindb_cheruvu",
		"priority_skinpix_cheruvu_huskindb",
		"priority_huskindb_skinpix_cheruvu",
		"priority_huskindb_cheruvu_skinpix",
		"priority_cheruvu_skinpix_huskindb",
		"priority_cheruvu_huskindb_skinpix"
	))
)
for (path in run_dirs) {
	if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) stop("Cannot create: ", path, call. = FALSE)
}

# Run each stage in a separate R session.
status_log <- file.path(log_dir, "run_all.log")
writeLines(c(
	paste("Started:", Sys.time()),
	paste("Python:", python),
	paste("RDKit:", python_info),
	capture.output(sessionInfo())
), status_log)
for (i in seq_along(scripts)) {
	script <- scripts[[i]]
	log_file <- file.path(log_dir, sub("\\.R$", ".log", basename(script)))
	args <- if (i == 1L) "--read-only-curation=true" else character()
	message(sprintf("[%d/%d] %s", i, length(scripts), script))
	status <- system2(rscript, c("--vanilla", shQuote(script), args), stdout = log_file, stderr = log_file)
	if (status != 0L) {
		cat("FAILED ", script, "\n", sep = "", file = status_log, append = TRUE)
		stop("Script failed: ", script, "\nSee: ", log_file, call. = FALSE)
	}
	cat("OK ", script, "\n", sep = "", file = status_log, append = TRUE)
}
cat("COMPLETE ", format(Sys.time()), "\n", sep = "", file = status_log, append = TRUE)
message("Complete. Logs are in ", log_dir, "/.")

thx <- function() {
	cat(paste("Stay hydrated. Have a nice day.", "\n"))
	cat(
		paste0(c(rep(" ", 14), "_", rep(" ", 3), "_", rep(" ", 3), "_"), collapse = ""),
		paste0(c(rep(" ", 13), "/ \\ / \\ / \\"), collapse = ""),
		paste0(c(rep(" ", 13), "\\", rep(" ", 2), "V", rep(" ", 3), "V", rep(" ", 2), "/ __"), collapse = ""),
		paste0(c(rep(" ", 5), "_____.---'", rep(" ", 7), "'-/ .\\____"), collapse = ""),
		paste0(c(rep(" ", 3), ".'", rep(" ", 5), "/", rep(" ", 12), "|", rep(" ", 3), ".\\", rep(" ", 3), "'."), collapse = ""),
		paste0(c(rep(" ", 2), "/", rep(" ", 5), "_'", rep(" ", 14), "\\_|", rep(" ", 2), "|'.", rep(" ", 2), "\\"), collapse = ""),
		paste0(c(" ", "|", rep(" ", 3), "_ '", rep(" ", 2), "/", rep(" ", 16), "- . - . |"), collapse = ""),
		paste0(c(" ", "|_'", rep(" ", 6), "|", rep(" ", 3), "(.", rep(" ", 2), ") (.", rep(" ", 2), ")", rep(" ", 3), "/", rep(" ", 3), "./`|"), collapse = ""),
		paste0(c(rep(" ", 2), "\\", rep(" ", 6), "/", rep(" ", 10), ".", rep(" ", 7), "'--/", rep(" ", 2), "/"), collapse = ""),
		paste0(c(rep(" ", 3), "'.___/", rep(" ", 10), "--", rep(" ", 9), "|_.'"), collapse = ""),
		paste0(c(rep(" ", 7), "\\", rep(" ", 20), "/"), collapse = ""),
		paste0(c(rep(" ", 8), "'---.__________.---'"), collapse = ""),
		paste0(c(rep(" ", 11), "/ |", rep(" ", 10), "| \\"), collapse = ""),
		paste0(c(rep(" ", 10), "/__|", rep(" ", 10), "|__\\"), collapse = ""),
		paste0(c(rep(" ", 11), "\\ /", rep(" ", 10), "\\ /"), collapse = ""),
		paste0(c(rep(" ", 11), "/", rep(" ", 13), "\\"), collapse = ""),
		paste0(c(rep(" ", 11), "'.__", rep(" ", 2), "____", rep(" ", 2), "__.'"), collapse = ""),
		paste0(c(rep(" ", 14), "/__\\", rep(" ", 2), "/__\\"), collapse = ""),
		paste0(c(rep(" ", 14), "'──'", rep(" ", 2), "'──'"), collapse = ""),
		"",
		sep = "\n"
	)
}

thx()
