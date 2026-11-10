# Describe dataset composition, metadata availability, and variance in pooled permeability.
input_file <- "data/processed/80_pooled.csv"
output_dir <- "R/results"
n_boot <- 2000L
seed <- 2026L
if (!file.exists(input_file)) stop("Missing pooled dataset. Run script 80 first.")
if (!dir.exists(output_dir)) stop("Required output directory is missing: ", output_dir)

read_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
write_result <- function(x, filename) {
	write.csv(x, file.path(output_dir, filename), row.names = FALSE, na = "")
}
clean <- function(x) {
	x <- trimws(as.character(x))
	x[tolower(x) %in% c(
		"", "na", "n/a", "n.a.", "nan", "unknown", "not reported", "not specified", "-", "--"
	)] <- NA_character_
	x
}

# Read the pooled data and original files used to recover metadata.
dat <- read_csv(input_file)
required <- c("record_id", "dataset_origin", "original_row", "compound_id", "publication_id", "log_kp", "temperature_c")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing pooled input columns: ", paste(missing, collapse = ", "))
if (!nrow(dat) || anyNA(dat[setdiff(required, "temperature_c")]) ||
	anyDuplicated(dat$record_id) || any(!is.finite(dat$log_kp))) stop("Invalid pooled records.")
if (any(!nzchar(dat$compound_id)) || any(!nzchar(dat$publication_id))) {
	stop("Missing compound or publication identity.")
}
raw_paths <- c(
	Cheruvu = "data/raw/cheruvu.xlsx",
	HuSkinDB = "data/raw/huskinDB.csv",
	SkinPiX = "data/raw/skinpix.xlsx"
)
origins <- names(raw_paths)
if (any(!dat$dataset_origin %in% origins)) stop("Unexpected dataset origin.")
checksum_file <- file.path(output_dir, "80_input_checksums.csv")
if (file.exists(checksum_file)) {
	checks <- read_csv(checksum_file)
	expected <- checks$md5[match(unname(raw_paths), checks$path)]
	actual <- unname(tools::md5sum(raw_paths))
	if (anyNA(actual) || !identical(actual, expected)) stop("Raw files differ from preparation inputs. Rerun 80 first.")
}
raw <- list(
	Cheruvu = as.data.frame(readxl::read_excel(raw_paths[["Cheruvu"]])),
	HuSkinDB = read_csv(raw_paths[["HuSkinDB"]]),
	SkinPiX = as.data.frame(readxl::read_excel(raw_paths[["SkinPiX"]]))
)

dataset_structure <- do.call(rbind, lapply(origins, function(origin) {
	x <- dat[dat$dataset_origin == origin, , drop = FALSE]
	data.frame(
		dataset = origin,
		raw_measurements = nrow(raw[[origin]]),
		retained_measurements = nrow(x),
		compounds = length(unique(x$compound_id)),
		publications = length(unique(x$publication_id))
	)
}))
dataset_structure <- rbind(dataset_structure, data.frame(
	dataset = "Pooled",
	raw_measurements = sum(vapply(raw, nrow, integer(1))),
	retained_measurements = nrow(dat),
	compounds = length(unique(dat$compound_id)),
	publications = length(unique(dat$publication_id))
))
write_result(dataset_structure, "85_dataset_structure.csv")

# Map pooled metadata fields to raw columns. NA means no separate field was supplied.
fields <- c(
	"Skin donor type", "Anatomical site", "Skin layer", "Skin preparation",
	"Skin thickness", "Diffusion cell type", "Donor medium", "Receptor medium",
	"Donor or skin temperature", "Receptor temperature", "Donor pH", "Receptor pH",
	"Skin integrity test"
)
mapping <- list(
	Cheruvu = c(NA, NA, NA, NA, "Skin thicknessj", NA, NA, NA, "Texpi", NA, NA, NA, "Skin Integrity testk"),
	HuSkinDB = c(
		"skin source type", "skin source site", "used layer", "skin preparation", NA,
		"cell type", "donor type", "acceptor type", "donor/skin surface temperature (°C)",
		"acceptor temperature (°C)", "donor pH", "acceptor pH", NA
	),
	SkinPiX = c(
		"skin donor type", "skin source site", "used layer", "skin preparation", NA,
		"cell type", "donor type", "acceptor type", "donor/skin surface temperature (°C)",
		"acceptor temperature (°C)", "donor pH", "acceptor pH", NA
	)
)

availability <- list()
for (origin in origins) {
	ii <- which(dat$dataset_origin == origin)
	rows <- dat$original_row[ii]
	if (any(!is.finite(rows) | rows != floor(rows) | rows < 1 | rows > nrow(raw[[origin]]))) {
		stop("Invalid raw row indices: ", origin)
	}
	if (length(ii) && !identical(as.character(dat$record_id[ii]), paste0(origin, "_", rows))) {
		stop("Record IDs disagree with original row provenance: ", origin)
	}
	for (f in seq_along(fields)) {
		column <- mapping[[origin]][f]
		provided <- !is.na(column)
		if (provided && !column %in% names(raw[[origin]])) stop("Missing raw metadata column: ", origin, " / ", column)
		values <- if (fields[f] == "Donor or skin temperature") {
			clean(dat$temperature_c[ii])
		} else if (provided) {
			clean(raw[[origin]][[column]][rows])
		} else {
			rep(NA_character_, length(ii))
		}
		reported <- !is.na(values)
		availability[[length(availability) + 1L]] <- data.frame(
			dataset = origin,
			metadata = fields[f],
			field_provided = provided,
			n_reported = sum(reported),
			n_records = length(ii),
			percent_reported = if (length(ii)) 100 * mean(reported) else NA_real_
		)
	}
}

availability <- do.call(rbind, availability)
pool_availability <- do.call(rbind, lapply(fields, function(field) {
	x <- availability[availability$metadata == field, , drop = FALSE]
	n_reported <- sum(x$n_reported)
	data.frame(
		dataset = "Pooled",
		metadata = field,
		field_provided = any(x$field_provided & x$n_records > 0),
		n_reported = n_reported,
		n_records = nrow(dat),
		percent_reported = 100 * n_reported / nrow(dat)
	)
}))
write_result(rbind(availability, pool_availability), "85_metadata_availability.csv")

# Center measurements within compounds with at least two observations.
compound_summary <- do.call(rbind, lapply(split(dat, dat$compound_id), function(x) {
	data.frame(
		compound_id = x$compound_id[1],
		n_measurements = nrow(x),
		mean_log_kp = mean(x$log_kp),
		range = diff(range(x$log_kp)),
		stringsAsFactors = FALSE
	)
}))
repeated <- compound_summary[compound_summary$n_measurements >= 2L, , drop = FALSE]
repeated <- repeated[order(repeated$range, repeated$compound_id), ]
j <- match(dat$compound_id, repeated$compound_id)
keep <- !is.na(j)
points <- dat[keep, , drop = FALSE]
points$centered_log_kp <- points$log_kp - repeated$mean_log_kp[j[keep]]
points$compound_order <- j[keep]
write_result(points, "85_centered_measurements.csv")

# Fit the intercept-only variance model with compound and publication effects.
model_data <- dat
model_data$compound_id <- factor(model_data$compound_id)
model_data$publication_id <- factor(model_data$publication_id)
ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
fit <- lme4::lmer(
	log_kp ~ 1 + (1 | compound_id) + (1 | publication_id),
	data = model_data,
	REML = TRUE,
	control = ctrl
)

components <- c("Compound", "Publication", "Residual")
extract_variance <- function(model) {
	vc <- as.data.frame(lme4::VarCorr(model))
	v <- vc$vcov[match(c("compound_id", "publication_id", "Residual"), vc$grp)]
	if (any(!is.finite(v)) || sum(v) <= 0) stop("Invalid variance components.")
	c(setNames(v, paste0(components, "_variance")), setNames(100 * v / sum(v), paste0(components, "_percent")))
}
estimate <- extract_variance(fit)

# Simulate new random effects and residuals for each bootstrap refit.
message("Variance bootstrap: ", n_boot, " replicates.")
boot <- lme4::bootMer(fit, FUN = extract_variance, nsim = n_boot, type = "parametric", use.u = FALSE, seed = seed)
usable <- rowSums(!is.finite(boot$t)) == 0L
n_usable <- sum(usable)
if (n_usable < n_boot) {
	message("Usable bootstrap replicates: ", n_usable, "/", n_boot)
	failure_messages <- attr(boot, "boot.fail.msgs")
	if (length(failure_messages)) print(failure_messages)
}
if (n_usable < 2L) stop("Too few usable bootstrap replicates.")
ci <- apply(boot$t[usable, , drop = FALSE], 2, quantile, probs = c(0.025, 0.975), names = FALSE)

variance <- data.frame(
	component = components,
	variance = unname(estimate[1:3]),
	variance_ci_low = ci[1, 1:3],
	variance_ci_high = ci[2, 1:3],
	percent = unname(estimate[4:6]),
	percent_ci_low = ci[1, 4:6],
	percent_ci_high = ci[2, 4:6],
	n_boot_requested = n_boot,
	n_boot_usable = n_usable,
	row.names = NULL
)
write_result(variance, "85_variance_components.csv")
write_result(data.frame(
	input = normalizePath(input_file),
	md5 = unname(tools::md5sum(input_file)),
	n_boot = n_boot,
	seed = seed,
	REML = TRUE
), "85_analysis_settings.csv")
message("Completed pooled variation analysis.")
