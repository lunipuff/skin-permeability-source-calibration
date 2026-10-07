# ============================================================
# 00_config.R
# Project-wide packages and utility functions
# ============================================================

# ---- packages ----

packages <- c(
  "readxl",
  "readr",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "ggplot2",
  "lme4",
  "lmerTest",
  "broom",
  "broom.mixed",
  "performance"
)

to_install <- packages[
  !packages %in% rownames(installed.packages())
]

if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(
  lapply(
    packages,
    library,
    character.only = TRUE
  )
)


# ---- global options ----

options(
  dplyr.summarise.inform = FALSE,
  readr.show_col_types = FALSE
)

set.seed(2026)


# ---- utility functions ----

# Standardize free-text labels without changing their meaning.
clean_text <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_squish() %>%
    dplyr::na_if("")
}


# Standardize CAS formatting.
# CAS number is used as the primary chemical identifier where available.
clean_cas <- function(x) {
	x <- x %>%
		as.character() %>%
		stringr::str_squish() %>%
		stringr::str_replace_all("\\s+", "")

	missing_tokens <- c(
		"",
		"NA",
		"N/A",
		"N.A.",
		"NONE",
		"-",
		"--"
	)

	x[
		toupper(x) %in%
			missing_tokens
	] <- NA_character_

	x
}


# Stable fallback identifier for compounds without CAS numbers.
normalize_compound_name <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_to_lower() %>%
    stringr::str_squish() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

# Color palette
palette <- function(n = 100) {
	bay_modified <- c(
		"#00496f",
		"#0f85a0",
		"#20b2aa",
		"#f0ead6",
		"#f1d94a",
		"#ed8b00",
		"#dd4124"
	)

	grDevices::colorRampPalette(
		bay_modified
	)(n)
}
