#!/usr/bin/env Rscript
# Place beside README.md, R/, and data/. Run: Rscript run_all.R
# Rebuilds generated outputs in place; no dated folders are created.
# Clears data/processed/, R/results/, R/tables/, and R/figures/ first.
# Preserves data/raw/, analysis scripts, and 00_GA.pdf in figures/.

main <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) != 1L) {
    stop("Run this file with: Rscript run_all.R", call. = FALSE)
  }
  root <- dirname(normalizePath(sub("^--file=", "", script_arg),
                                winslash = "/", mustWork = TRUE))
  old_wd <- setwd(root)
  on.exit(setwd(old_wd), add = TRUE)

  scripts <- file.path("R", c(
    "01_prepare_describe_cheruvu.R",
    "10_fit_diagnose_cheruvu_hierarchical.R",
    "20_fit_cheruvu_qspr_models.R",
    "21_validate_cheruvu_qspr.R",
    "30_prepare_describe_huskindb.R",
    "31_fit_validate_huskindb.R",
    "32_compare_calibration_methods.R",
    "40_prepare_describe_skinpix.R",
    "41_validate_skinpix_transfer.R",
    "50_falsification_checks.R",
    "60_uncertainty_analysis.R",
    "70_calibration_failure_analysis.R",
    "99a_make_tables.R",
    "99b_make_figures.R"
  ))
  output_dirs <- c("data/processed", "R/results", "R/tables", "R/figures")
  raw_files <- file.path("data/raw", c("cheruvu.xlsx", "huskinDB.csv", "skinpix.xlsx"))

  # Check requirements before deleting any generated files.
  required <- c(scripts, "R/40_make_rdkit_descriptors.py", raw_files)
  missing <- required[!file.exists(required)]
  if (length(missing)) {
    stop("Missing required files:\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
  for (script in scripts) parse(file = script)

  packages <- c("readxl", "lme4", "ggplot2", "patchwork")
  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (any(!installed)) {
    stop("Install required R packages: ", paste(packages[!installed], collapse = ", "),
         call. = FALSE)
  }
  if (!"free" %in% getNamespaceExports("patchwork")) {
    stop("The figure script requires a patchwork version exporting free().", call. = FALSE)
  }

  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  # Script 40 invokes 'python' itself, so check that exact interpreter.
  python <- unname(Sys.which("python"))
  if (!nzchar(python)) {
    stop("'python' is not on PATH. Activate the Python environment containing RDKit.", call. = FALSE)
  }
  python_code <- paste(
    "import sys, ast, rdkit;",
    "ast.parse(open('R/40_make_rdkit_descriptors.py', encoding='utf-8').read());",
    "print('Python executable:', sys.executable);",
    "print('Python:', sys.version.replace(chr(10), ' '));",
    "print('RDKit:', rdkit.__version__)"
  )
  python_info <- system2(python, c("-c", shQuote(python_code)), stdout = TRUE, stderr = TRUE)
  python_status <- attr(python_info, "status")
  if (!is.null(python_status) && python_status != 0L) {
    stop("Python/RDKit preflight failed:\n", paste(python_info, collapse = "\n"), call. = FALSE)
  }

  # Refuse cleanup through redirected output paths outside this project.
  for (path in output_dirs) {
    link <- Sys.readlink(path)
    if (!is.na(link) && nzchar(link)) stop("Output directory is a symbolic link: ", path)
    parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
    if (!startsWith(paste0(parent, "/"), paste0(root, "/"))) {
      stop("Output parent is outside the project: ", path)
    }
    if (file.exists(path) && !dir.exists(path)) stop("Expected a directory: ", path)
  }

  clear_directory <- function(path, preserve_ga = FALSE) {
    entries <- list.files(path, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    for (entry in entries) {
      link <- Sys.readlink(entry)
      is_link <- !is.na(link) && nzchar(link)
      if (dir.exists(entry) && !is_link) {
        clear_directory(entry, preserve_ga)
        if (!length(list.files(entry, all.files = TRUE, no.. = TRUE))) {
          if (unlink(entry, recursive = TRUE) != 0L) stop("Cannot remove: ", entry)
        }
      } else {
        keep <- preserve_ga && !is_link &&
          grepl("^(00_)?GA\\.pdf$", basename(entry), ignore.case = TRUE)
        if (!keep && unlink(entry, recursive = FALSE) != 0L) stop("Cannot remove: ", entry)
      }
    }
  }

  message("Clearing generated outputs; preserving raw data and the graphical abstract.")
  for (path in output_dirs) {
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) stop("Cannot create: ", path)
    clear_directory(path, preserve_ga = identical(path, "R/figures"))
  }

  log_dir <- file.path("R/results", "logs")
  if (!dir.create(log_dir, recursive = TRUE) && !dir.exists(log_dir)) stop("Cannot create log directory.")
  writeLines(c(
    paste("Project:", root),
    paste("Started:", Sys.time()),
    paste("R packages:", paste(paste0(packages, "=", vapply(packages, function(p) {
      as.character(packageVersion(p))
    }, character(1))), collapse = ", ")),
    python_info,
    "", capture.output(sessionInfo())
  ), "R/results/sessionInfo.txt")

  status_file <- file.path(root, "R/results/run_status.txt")
  writeLines(paste("RUNNING", Sys.time()), status_file)
  completed <- FALSE
  on.exit({
    if (!completed) cat("FAILED or INTERRUPTED\n", file = status_file, append = TRUE)
  }, add = TRUE)

  for (i in seq_along(scripts)) {
    script <- scripts[i]
    log_file <- file.path(log_dir, sub("\\.R$", ".log", basename(script)))
    message(sprintf("[%d/%d] %s", i, length(scripts), script))
    # Each analysis script gets its own clean R process and a combined output log.
    status <- system2(rscript, c("--vanilla", shQuote(script)),
                      stdout = log_file, stderr = log_file)
    if (status != 0L) {
      stop("Script failed: ", script, "\nSee: ", log_file, call. = FALSE)
    }
    cat(paste("OK", script), "\n", file = status_file, append = TRUE)
  }
  completed <- TRUE
  cat(paste("SUCCESS", Sys.time()), "\n", file = status_file, append = TRUE)
  message("Complete. Outputs: data/processed/, R/results/, R/tables/, R/figures/.")
  message("Logs and environment versions: R/results/.")
}

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

main()
thx()
