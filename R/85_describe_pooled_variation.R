# 85_describe_pooled_variation.R
#
# Run from the project root after 80:
# Rscript R/85_describe_pooled_variation.R
#
# Descriptive model:
# log_kp ~ 1 + (1 | compound_id) + (1 | publication_id)
#
# Does not modify the pooled dataset or previous analysis outputs.

if (
  !requireNamespace("lme4", quietly = TRUE) ||
  !requireNamespace("readxl", quietly = TRUE)
) {
  stop("Install lme4 and readxl first.")
}

input_file <- "data/processed/80_pooled.csv"
output_dir <- "R/results"

n_boot <- 2000L
seed <- 2026L

if (n_boot < 2L) {
  stop("Use at least two bootstrap replicates.")
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_result <- function(x, filename) {
  write.csv(
    x,
    file.path(output_dir, filename),
    row.names = FALSE,
    na = ""
  )
}

clean <- function(x) {
  x <- trimws(as.character(x))

  x[tolower(x) %in% c(
    "", "na", "n/a", "n.a.", "nan", "unknown",
    "not reported", "not specified", "-", "--"
  )] <- NA_character_

  x
}

# ------------------------------------------------------------
# 1. Read and check pooled data
# ------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Run 80 first.")
}

dat <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c(
  "record_id", "dataset_origin", "original_row",
  "compound_id", "publication_id", "log_kp",
  "temperature_c"
)

if (!all(required %in% names(dat))) {
  stop("Missing pooled input columns.")
}

identity_fields <- setdiff(required, "temperature_c")

if (
  !nrow(dat) ||
  anyNA(dat[identity_fields]) ||
  anyDuplicated(dat$record_id) ||
  any(!is.finite(dat$log_kp))
) {
  stop("Invalid pooled records.")
}

if (
  any(!nzchar(dat$compound_id)) ||
  any(!nzchar(dat$publication_id))
) {
  stop("Missing compound/publication identity.")
}

origins <- c("Cheruvu", "HuSkinDB", "SkinPiX")

if (any(!dat$dataset_origin %in% origins)) {
  stop("Unexpected dataset origin.")
}

raw_paths <- c(
  Cheruvu = "data/raw/cheruvu.xlsx",
  HuSkinDB = "data/raw/huskinDB.csv",
  SkinPiX = "data/raw/skinpix.xlsx"
)

if (any(!file.exists(raw_paths))) {
  stop("Missing original raw files.")
}

# Check that metadata are recovered from the same raw files
# used during pooled preparation.

checksum_file <- file.path(
  output_dir,
  "80_input_checksums.csv"
)

if (file.exists(checksum_file)) {
  checks <- read.csv(
    checksum_file,
    stringsAsFactors = FALSE
  )

  j <- match(unname(raw_paths), checks$path)

  if (
    anyNA(j) ||
    any(
      unname(tools::md5sum(raw_paths)) !=
        checks$md5[j]
    )
  ) {
    stop(
      "Raw files differ from preparation inputs. ",
      "Rerun 80 first."
    )
  }
}

raw <- list(
  Cheruvu = as.data.frame(
    readxl::read_excel(raw_paths["Cheruvu"])
  ),

  HuSkinDB = read.csv(
    raw_paths["HuSkinDB"],
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),

  SkinPiX = as.data.frame(
    readxl::read_excel(raw_paths["SkinPiX"])
  )
)

# ------------------------------------------------------------
# 2. Dataset structure
# ------------------------------------------------------------

structure <- do.call(rbind, lapply(
  origins,
  function(origin) {
    x <- dat[
      dat$dataset_origin == origin,
      ,
      drop = FALSE
    ]

    data.frame(
      dataset = origin,
      raw_measurements = nrow(raw[[origin]]),
      retained_measurements = nrow(x),
      compounds = length(unique(x$compound_id)),
      publications = length(unique(x$publication_id))
    )
  }
))

structure <- rbind(
  structure,
  data.frame(
    dataset = "Pooled",

    raw_measurements = sum(
      vapply(raw, nrow, integer(1))
    ),

    retained_measurements = nrow(dat),
    compounds = length(unique(dat$compound_id)),
    publications = length(unique(dat$publication_id))
  )
)

# Compounds shared between database contributions make the
# database-specific compound counts non-additive.

write_result(
  structure,
  "85_dataset_structure.csv"
)

# ------------------------------------------------------------
# 3. Recover metadata using original row provenance
# ------------------------------------------------------------

fields <- c(
  "Skin donor type",
  "Anatomical site",
  "Skin layer",
  "Skin preparation",
  "Skin thickness",
  "Diffusion cell type",
  "Donor medium",
  "Receptor medium",
  "Donor or skin temperature",
  "Receptor temperature",
  "Donor pH",
  "Receptor pH",
  "Skin integrity test"
)

# NA means that no separate field was supplied.
# These mappings refer to columns in the original raw files.

mapping <- list(
  Cheruvu = c(
    NA, NA, NA, NA,
    "Skin thicknessj",
    NA, NA, NA,
    "Texpi",
    NA, NA, NA,
    "Skin Integrity testk"
  ),

  HuSkinDB = c(
    "skin source type",
    "skin source site",
    "used layer",
    "skin preparation",
    NA,
    "cell type",
    "donor type",
    "acceptor type",
    "donor/skin surface temperature (°C)",
    "acceptor temperature (°C)",
    "donor pH",
    "acceptor pH",
    NA
  ),

  SkinPiX = c(
    "skin donor type",
    "skin source site",
    "used layer",
    "skin preparation",
    NA,
    "cell type",
    "donor type",
    "acceptor type",
    "donor/skin surface temperature (°C)",
    "acceptor temperature (°C)",
    "donor pH",
    "acceptor pH",
    NA
  )
)

metadata <- dat[c(
  "record_id",
  "dataset_origin",
  "publication_id",
  "compound_id"
)]

for (field in fields) {
  metadata[[field]] <- NA_character_
}

availability <- list()

for (origin in origins) {

  ii <- which(dat$dataset_origin == origin)
  rows <- dat$original_row[ii]

  if (any(
    !is.finite(rows) |
      rows != floor(rows) |
      rows < 1 |
      rows > nrow(raw[[origin]])
  )) {
    stop("Invalid raw row indices: ", origin)
  }

  if (
    length(ii) &&
    !identical(
      as.character(dat$record_id[ii]),
      paste0(origin, "_", rows)
    )
  ) {
    stop(
      "Record IDs disagree with original row provenance: ",
      origin
    )
  }

  for (f in seq_along(fields)) {

    column <- mapping[[origin]][f]
    provided <- !is.na(column)

    if (
      provided &&
      !column %in% names(raw[[origin]])
    ) {
      stop(
        "Missing expected raw metadata column: ",
        origin, " / ", column
      )
    }

    values <- if (provided) {
      clean(raw[[origin]][[column]][rows])
    } else {
      rep(NA_character_, length(ii))
    }

    # Use the harmonised Celsius values already in the pool.
    if (fields[f] == "Donor or skin temperature") {
      values <- clean(dat$temperature_c[ii])
    }

    metadata[[fields[f]]][ii] <- values

    availability[[length(availability) + 1L]] <-
      data.frame(
        dataset = origin,
        metadata = fields[f],
        field_provided = provided,
        n_reported = sum(!is.na(values)),
        n_records = length(ii),

        percent_reported = if (length(ii)) {
          100 * mean(!is.na(values))
        } else {
          NA_real_
        }
      )
  }
}

availability <- do.call(rbind, availability)

pool_availability <- do.call(rbind, lapply(
  fields,
  function(field) {

    x <- availability[
      availability$metadata == field,
      ,
      drop = FALSE
    ]

    data.frame(
      dataset = "Pooled",
      metadata = field,

      field_provided = any(
        x$field_provided & x$n_records > 0
      ),

      n_reported = sum(x$n_reported),
      n_records = nrow(dat),

      percent_reported =
        100 * sum(x$n_reported) / nrow(dat)
    )
  }
))

write_result(
  metadata,
  "85_retained_metadata.csv"
)

write_result(
  rbind(availability, pool_availability),
  "85_metadata_availability.csv"
)

# ------------------------------------------------------------
# 4. Within-compound variation and figure data
# ------------------------------------------------------------

# Retain individual measurements rather than averaging them.
# Figure data include compounds with at least two measurements.
# n_publications distinguishes within- from cross-publication
# replication.

compound_summary <- do.call(rbind, lapply(
  split(dat, dat$compound_id),
  function(x) {
    data.frame(
      compound_id = x$compound_id[1],
      n_measurements = nrow(x),
      n_publications = length(unique(x$publication_id)),
      mean_log_kp = mean(x$log_kp),
      minimum = min(x$log_kp),
      maximum = max(x$log_kp),
      range = diff(range(x$log_kp)),
      stringsAsFactors = FALSE
    )
  }
))

rownames(compound_summary) <- NULL

compound_summary <- compound_summary[
  order(
    compound_summary$range,
    compound_summary$compound_id
  ),
]

compound_summary$plot_order <- NA_integer_

repeated <- compound_summary$n_measurements >= 2

compound_summary$plot_order[repeated] <-
  seq_len(sum(repeated))

j <- match(
  dat$compound_id,
  compound_summary$compound_id
)

points <- dat[
  repeated[j],
  ,
  drop = FALSE
]

points$centered_log_kp <-
  points$log_kp -
  compound_summary$mean_log_kp[j[repeated[j]]]

points$compound_order <-
  compound_summary$plot_order[j[repeated[j]]]

write_result(
  compound_summary,
  "85_within_compound_summary.csv"
)

write_result(
  points,
  "85_centered_measurements.csv"
)

# ------------------------------------------------------------
# 5. Pooled variance decomposition
# ------------------------------------------------------------

# Intercept-only descriptive model.
# Molecular descriptors and temperature are not included.

model_data <- dat

model_data$compound_id <-
  factor(model_data$compound_id)

model_data$publication_id <-
  factor(model_data$publication_id)

ctrl <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 100000)
)

fit <- lme4::lmer(
  log_kp ~ 1 +
    (1 | compound_id) +
    (1 | publication_id),
  data = model_data,
  REML = TRUE,
  control = ctrl
)

saveRDS(
  fit,
  file.path(output_dir, "85_pooled_variance_model.rds")
)

extract_variance <- function(model) {

  vc <- as.data.frame(lme4::VarCorr(model))

  v <- setNames(
    vc$vcov[match(
      c("compound_id", "publication_id", "Residual"),
      vc$grp
    )],
    c("Compound", "Publication", "Residual")
  )

  if (
    any(!is.finite(v)) ||
    sum(v) <= 0
  ) {
    stop("Invalid variance components.")
  }

  c(
    setNames(
      v,
      paste0(names(v), "_variance")
    ),

    setNames(
      100 * v / sum(v),
      paste0(names(v), "_percent")
    )
  )
}

estimate <- extract_variance(fit)

# ------------------------------------------------------------
# 6. Parametric-bootstrap confidence intervals
# ------------------------------------------------------------

# Simulate compound effects, publication effects and residuals
# anew in each replicate; refit the same model.

message(
  "Variance bootstrap: ",
  n_boot,
  " replicates."
)

boot <- lme4::bootMer(
  fit,
  FUN = extract_variance,
  nsim = n_boot,
  type = "parametric",
  use.u = FALSE,
  seed = seed
)

usable <- apply(
  boot$t,
  1,
  function(x) all(is.finite(x))
)

bootstrap_notes <- c(
  paste("Requested replicates:", n_boot),
  paste("Usable replicates:", sum(usable)),
  paste("Non-finite replicates:", sum(!usable)),

  "Bootstrap fitting failures:",
  capture.output(print(attr(boot, "bootFail"))),

  "Failure messages:",
  capture.output(print(attr(boot, "boot.fail.msgs"))),

  "All bootstrap messages:",
  capture.output(print(attr(boot, "boot.all.msgs")))
)

writeLines(
  bootstrap_notes,
  file.path(
    output_dir,
    "85_bootstrap_diagnostics.txt"
  )
)

if (sum(usable) < 2L) {
  stop(
    "Too few usable bootstrap replicates; ",
    "inspect diagnostics."
  )
}

ci <- apply(
  boot$t[usable, , drop = FALSE],
  2,
  quantile,
  probs = c(0.025, 0.975),
  names = FALSE
)

variance <- data.frame(
  component = c(
    "Compound",
    "Publication",
    "Residual"
  ),

  variance = unname(estimate[1:3]),
  variance_ci_low = ci[1, 1:3],
  variance_ci_high = ci[2, 1:3],

  percent = unname(estimate[4:6]),
  percent_ci_low = ci[1, 4:6],
  percent_ci_high = ci[2, 4:6],

  n_boot_requested = n_boot,
  n_boot_usable = sum(usable),

  row.names = NULL
)

write_result(
  variance,
  "85_variance_components.csv"
)

# ------------------------------------------------------------
# 7. Model diagnostics and provenance
# ------------------------------------------------------------

write_result(
  data.frame(
    input = normalizePath(input_file),
    md5 = unname(tools::md5sum(input_file)),
    n_boot = n_boot,
    seed = seed,
    REML = TRUE
  ),
  "85_analysis_settings.csv"
)

writeLines(
  c(
    capture.output(summary(fit)),
    paste(
      "Singular fit:",
      lme4::isSingular(fit)
    )
  ),
  file.path(
    output_dir,
    "85_model_diagnostics.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "85_sessionInfo.txt")
)

cat("\nPooled variation analysis complete.\n")

print(structure, row.names = FALSE)
print(variance, row.names = FALSE)

cat(
  "\nPercentages describe an intercept-only model ",
  "without molecular descriptors or temperature.\n",
  "Publication effects do not identify experimental causes ",
  "or stable laboratory effects.\n",
  "Metadata availability does not establish protocol equivalence.\n",
  sep = ""
)