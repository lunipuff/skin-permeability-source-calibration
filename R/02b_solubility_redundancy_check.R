############################################################
# 02b_solubility_redundancy_check.R
# Descriptor redundancy screening and red-flag set generation
############################################################

source("R/00_config.R")

############################################################
# Packages
############################################################

required_packages <- c(
	"ggplot2"
)

missing_packages <- required_packages[
	!vapply(
		required_packages,
		requireNamespace,
		logical(1),
		quietly = TRUE
	)
]

if (length(missing_packages) > 0) {
	stop(
		"Missing required packages: ",
		paste(missing_packages, collapse = ", "),
		"\nInstall them with:\ninstall.packages(c(",
		paste(sprintf('"%s"', missing_packages), collapse = ", "),
		"))",
		call. = FALSE
	)
}

############################################################
# Output folders
############################################################

dir.create("results/descriptor_redundancy", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("manuscript/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("manuscript/tables", showWarnings = FALSE, recursive = TRUE)

############################################################
# Output paths
############################################################

path_all_descriptor_correlation_matrix <- "results/descriptor_redundancy/all_descriptor_correlation_matrix.csv"
path_all_descriptor_high_correlations <- "results/descriptor_redundancy/all_descriptor_high_correlations.csv"
path_all_descriptor_vif <- "results/descriptor_redundancy/all_descriptor_vif.csv"
path_all_descriptor_redundancy_summary <- "results/descriptor_redundancy/all_descriptor_redundancy_summary.csv"

path_descriptor_red_flag_sets <- "results/descriptor_redundancy/descriptor_red_flag_sets.csv"
path_descriptor_soft_warning_sets <- "results/descriptor_redundancy/descriptor_soft_warning_sets.csv"

path_solubility_gse_dataset <- "results/descriptor_redundancy/solubility_gse_dataset.csv"
path_solubility_gse_correlation_matrix <- "results/descriptor_redundancy/solubility_gse_correlation_matrix.csv"
path_solubility_gse_regression_summary <- "results/descriptor_redundancy/solubility_gse_regression_summary.csv"
path_solubility_gse_vif <- "results/descriptor_redundancy/solubility_gse_vif.csv"

path_table_descriptor_redundancy_summary <- "tables/tableS1_descriptor_redundancy_summary.csv"
path_manuscript_table_descriptor_redundancy_summary <- "manuscript/tables/tableS1_descriptor_redundancy_summary.csv"

path_fig_descriptor_correlation_png <- "figures/figure2_descriptor_correlation.png"
path_fig_descriptor_correlation_pdf <- "figures/figure2_descriptor_correlation.pdf"
path_manuscript_fig_descriptor_correlation_png <- "manuscript/figures/figure2_descriptor_correlation.png"
path_manuscript_fig_descriptor_correlation_pdf <- "manuscript/figures/figure2_descriptor_correlation.pdf"

path_fig_descriptor_vif_png <- "figures/figureS1_descriptor_vif.png"
path_fig_descriptor_vif_pdf <- "figures/figureS1_descriptor_vif.pdf"
path_manuscript_fig_descriptor_vif_png <- "manuscript/figures/figureS1_descriptor_vif.png"
path_manuscript_fig_descriptor_vif_pdf <- "manuscript/figures/figureS1_descriptor_vif.pdf"

path_fig_solubility_correlation_png <- "figures/figureS2_solubility_descriptor_correlation.png"
path_fig_solubility_correlation_pdf <- "figures/figureS2_solubility_descriptor_correlation.pdf"

path_fig_gse_logS_vs_LogSaqd_png <- "figures/figureS3_gse_logS_vs_LogSaqd.png"
path_fig_gse_logS_vs_LogSaqd_pdf <- "figures/figureS3_gse_logS_vs_LogSaqd.pdf"

############################################################
# Redundancy-screening thresholds
############################################################

hard_pairwise_correlation_threshold <- 0.85
soft_pairwise_correlation_threshold <- 0.70

extreme_vif_threshold <- 100
moderate_vif_threshold <- 5

############################################################
# Helper functions
############################################################

check_required_columns <- function(data, required_cols, data_name = "data") {
	missing_cols <- required_cols[
		!(required_cols %in% names(data))
	]

	if (length(missing_cols) > 0) {
		stop(
			paste0(
				"Missing required columns in ",
				data_name,
				": ",
				paste(missing_cols, collapse = ", ")
			),
			call. = FALSE
		)
	}
}

to_numeric_columns <- function(data, cols) {
	for (col in cols) {
		data[[col]] <- as.numeric(data[[col]])
	}

	data
}

calculate_vif <- function(data, predictor_cols) {
	complete_data <- data[
		stats::complete.cases(data[, predictor_cols, drop = FALSE]),
		predictor_cols,
		drop = FALSE
	]

	if (nrow(complete_data) < length(predictor_cols) + 2) {
		stop(
			"Too few complete observations to calculate VIF.",
			call. = FALSE
		)
	}

	vif_values <- vapply(
		predictor_cols,
		function(response_col) {
			other_cols <- setdiff(predictor_cols, response_col)

			if (length(other_cols) == 0) {
				return(NA_real_)
			}

			formula_text <- paste(
				response_col,
				"~",
				paste(other_cols, collapse = " + ")
			)

			fit <- stats::lm(
				stats::as.formula(formula_text),
				data = complete_data
			)

			r_squared <- summary(fit)$r.squared

			if (is.na(r_squared)) {
				return(NA_real_)
			}

			if (r_squared >= 1) {
				return(Inf)
			}

			1 / (1 - r_squared)
		},
		numeric(1)
	)

	data.frame(
		descriptor = names(vif_values),
		VIF = as.numeric(vif_values),
		stringsAsFactors = FALSE
	)
}

make_correlation_long <- function(correlation_matrix) {
	out <- as.data.frame(
		as.table(correlation_matrix),
		stringsAsFactors = FALSE
	)

	names(out) <- c(
		"descriptor_1",
		"descriptor_2",
		"correlation"
	)

	out
}

make_high_correlation_table <- function(correlation_matrix, threshold) {
	correlation_long <- make_correlation_long(correlation_matrix)

	correlation_long <- correlation_long[
		as.character(correlation_long$descriptor_1) <
			as.character(correlation_long$descriptor_2),
		,
		drop = FALSE
	]

	correlation_long$absolute_correlation <- abs(correlation_long$correlation)

	out <- correlation_long[
		correlation_long$absolute_correlation >= threshold,
		,
		drop = FALSE
	]

	out <- out[
		order(-out$absolute_correlation),
		,
		drop = FALSE
	]

	row.names(out) <- NULL
	out
}

plot_correlation_heatmap <- function(correlation_matrix, title, output_png, output_pdf) {
	correlation_long <- make_correlation_long(correlation_matrix)

	correlation_long$row_index <- match(
		as.character(correlation_long$descriptor_1),
		rownames(correlation_matrix)
	)

	correlation_long$col_index <- match(
		as.character(correlation_long$descriptor_2),
		colnames(correlation_matrix)
	)

	correlation_long <- correlation_long[
		correlation_long$row_index > correlation_long$col_index,
		,
		drop = FALSE
	]

	correlation_long$descriptor_1 <- factor(
		correlation_long$descriptor_1,
		levels = rownames(correlation_matrix)
	)

	correlation_long$descriptor_2 <- factor(
		correlation_long$descriptor_2,
		levels = rev(colnames(correlation_matrix))
	)

	p <- ggplot2::ggplot(
		correlation_long,
		ggplot2::aes(
			x = descriptor_1,
			y = descriptor_2,
			fill = correlation
		)
	) +
		ggplot2::geom_tile(
			color = "white",
			linewidth = 0.25
		) +
		ggplot2::geom_text(
			ggplot2::aes(label = sprintf("%.2f", correlation)),
			size = 3
		) +
		ggplot2::scale_fill_gradient2(
			low = "#2166AC",
			mid = "white",
			high = "#B2182B",
			midpoint = 0,
			limits = c(-1, 1),
			name = "Pearson r"
		) +
		ggplot2::coord_fixed() +
		ggplot2::labs(
			title = title,
			x = NULL,
			y = NULL
		) +
		ggplot2::theme_minimal(base_size = 11) +
		ggplot2::theme(
			axis.text.x = ggplot2::element_text(
				angle = 45,
				hjust = 1,
				vjust = 1
			),
			panel.grid = ggplot2::element_blank(),
			plot.title = ggplot2::element_text(face = "bold")
		)

	ggplot2::ggsave(
		filename = output_png,
		plot = p,
		width = 7.5,
		height = 6.5,
		dpi = 300
	)

	ggplot2::ggsave(
		filename = output_pdf,
		plot = p,
		width = 7.5,
		height = 6.5
	)

	p
}

plot_vif_barplot <- function(vif_table, title, output_png, output_pdf) {
	vif_plot_data <- vif_table

	vif_plot_data <- vif_plot_data[
		order(vif_plot_data$VIF),
		,
		drop = FALSE
	]

	vif_plot_data$descriptor <- factor(
		vif_plot_data$descriptor,
		levels = vif_plot_data$descriptor
	)

	vif_plot_data$label <- ifelse(
		vif_plot_data$VIF >= 1000,
		formatC(vif_plot_data$VIF, format = "e", digits = 2),
		sprintf("%.2f", vif_plot_data$VIF)
	)

	p <- ggplot2::ggplot(
		vif_plot_data,
		ggplot2::aes(
			x = VIF,
			y = descriptor
		)
	) +
		ggplot2::geom_col(width = 0.7) +
		ggplot2::geom_text(
			ggplot2::aes(label = label),
			hjust = -0.1,
			size = 3
		) +
		ggplot2::geom_vline(
			xintercept = moderate_vif_threshold,
			linetype = "dashed",
			linewidth = 0.4
		) +
		ggplot2::geom_vline(
			xintercept = extreme_vif_threshold,
			linetype = "dotted",
			linewidth = 0.4
		) +
		ggplot2::scale_x_log10(
			expand = ggplot2::expansion(mult = c(0.02, 0.25))
		) +
		ggplot2::labs(
			title = title,
			x = "Variance inflation factor, log10 scale",
			y = NULL
		) +
		ggplot2::theme_minimal(base_size = 11) +
		ggplot2::theme(
			panel.grid.minor = ggplot2::element_blank(),
			plot.title = ggplot2::element_text(face = "bold")
		)

	ggplot2::ggsave(
		filename = output_png,
		plot = p,
		width = 7,
		height = 5,
		dpi = 300
	)

	ggplot2::ggsave(
		filename = output_pdf,
		plot = p,
		width = 7,
		height = 5
	)

	p
}

make_predictor_set_table <- function(predictor_sets, reason, source, severity) {
	if (length(predictor_sets) == 0) {
		return(
			data.frame(
				set_id = character(0),
				predictor_set = character(0),
				n_predictors = integer(0),
				reason = character(0),
				source = character(0),
				severity = character(0),
				stringsAsFactors = FALSE
			)
		)
	}

	data.frame(
		set_id = paste0(source, "_", seq_along(predictor_sets)),
		predictor_set = vapply(
			predictor_sets,
			function(x) paste(sort(unique(x)), collapse = ";"),
			character(1)
		),
		n_predictors = vapply(
			predictor_sets,
			function(x) length(unique(x)),
			integer(1)
		),
		reason = reason,
		source = source,
		severity = severity,
		stringsAsFactors = FALSE
	)
}

deduplicate_predictor_set_table <- function(x) {
	if (nrow(x) == 0) {
		return(x)
	}

	x <- x[!duplicated(x$predictor_set), , drop = FALSE]
	x$set_id <- paste0(x$severity, "_", seq_len(nrow(x)))
	row.names(x) <- NULL
	x
}

summarize_lm <- function(model_name, fit) {
	fit_summary <- summary(fit)

	data.frame(
		model = model_name,
		n = stats::nobs(fit),
		r_squared = fit_summary$r.squared,
		adjusted_r_squared = fit_summary$adj.r.squared,
		residual_standard_error = fit_summary$sigma,
		stringsAsFactors = FALSE
	)
}

############################################################
# Load cleaned dataset
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

if (exists("all_predictors")) {
	candidate_descriptors <- all_predictors
} else {
	candidate_descriptors <- c(
		"MWa",
		"logKowb",
		"Mptc",
		"LogSaqd",
		"LogSoce",
		"Hdf",
		"Hag",
		"MVh",
		"Texpi",
		"Skin.thicknessj"
	)
}

check_required_columns(
	df,
	candidate_descriptors,
	data_name = "cleaned dataset"
)

df <- to_numeric_columns(
	df,
	candidate_descriptors
)

descriptor_data <- df[
	stats::complete.cases(df[, candidate_descriptors, drop = FALSE]),
	candidate_descriptors,
	drop = FALSE
]

if (nrow(descriptor_data) == 0) {
	stop(
		"No complete observations available for descriptor redundancy screening.",
		call. = FALSE
	)
}

############################################################
# A. General descriptor correlation screening
############################################################

all_descriptor_correlation <- stats::cor(
	descriptor_data,
	use = "complete.obs",
	method = "pearson"
)

write.csv(
	all_descriptor_correlation,
	path_all_descriptor_correlation_matrix
)

all_descriptor_high_correlations <- make_high_correlation_table(
	all_descriptor_correlation,
	threshold = soft_pairwise_correlation_threshold
)

write.csv(
	all_descriptor_high_correlations,
	path_all_descriptor_high_correlations,
	row.names = FALSE
)

plot_correlation_heatmap(
	correlation_matrix = all_descriptor_correlation,
	title = "Candidate descriptor correlation structure",
	output_png = path_fig_descriptor_correlation_png,
	output_pdf = path_fig_descriptor_correlation_pdf
)

file.copy(
	path_fig_descriptor_correlation_png,
	path_manuscript_fig_descriptor_correlation_png,
	overwrite = TRUE
)

file.copy(
	path_fig_descriptor_correlation_pdf,
	path_manuscript_fig_descriptor_correlation_pdf,
	overwrite = TRUE
)

############################################################
# B. General descriptor VIF screening
############################################################

all_descriptor_vif <- calculate_vif(
	data = df,
	predictor_cols = candidate_descriptors
)

all_descriptor_vif <- all_descriptor_vif[
	order(-all_descriptor_vif$VIF),
	,
	drop = FALSE
]

write.csv(
	all_descriptor_vif,
	path_all_descriptor_vif,
	row.names = FALSE
)

plot_vif_barplot(
	vif_table = all_descriptor_vif,
	title = "Variance inflation factors for candidate descriptors",
	output_png = path_fig_descriptor_vif_png,
	output_pdf = path_fig_descriptor_vif_pdf
)

file.copy(
	path_fig_descriptor_vif_png,
	path_manuscript_fig_descriptor_vif_png,
	overwrite = TRUE
)

file.copy(
	path_fig_descriptor_vif_pdf,
	path_manuscript_fig_descriptor_vif_pdf,
	overwrite = TRUE
)

############################################################
# C. Generate hard red-flag predictor sets
############################################################

hard_pairwise_sets <- all_descriptor_high_correlations[
	all_descriptor_high_correlations$absolute_correlation >=
		hard_pairwise_correlation_threshold,
	,
	drop = FALSE
]

hard_pairwise_predictor_sets <- vector(
	"list",
	nrow(hard_pairwise_sets)
)

if (nrow(hard_pairwise_sets) > 0) {
	for (i in seq_len(nrow(hard_pairwise_sets))) {
		hard_pairwise_predictor_sets[[i]] <- c(
			as.character(hard_pairwise_sets$descriptor_1[i]),
			as.character(hard_pairwise_sets$descriptor_2[i])
		)
	}
}

hard_pairwise_table <- make_predictor_set_table(
	predictor_sets = hard_pairwise_predictor_sets,
	reason = paste0(
		"Pairwise Pearson correlation with |r| >= ",
		hard_pairwise_correlation_threshold,
		". Candidate formulas containing this descriptor combination should be skipped."
	),
	source = "hard_pairwise_correlation",
	severity = "hard"
)

extreme_vif_descriptors <- all_descriptor_vif$descriptor[
	all_descriptor_vif$VIF >= extreme_vif_threshold
]

extreme_vif_predictor_sets <- list()

if (length(extreme_vif_descriptors) >= 2) {
	extreme_vif_predictor_sets <- list(
		extreme_vif_descriptors
	)
}

extreme_vif_table <- make_predictor_set_table(
	predictor_sets = extreme_vif_predictor_sets,
	reason = paste0(
		"Descriptors with VIF >= ",
		extreme_vif_threshold,
		". Candidate formulas containing all descriptors in this multivariable redundancy set should be skipped."
	),
	source = "extreme_vif_group",
	severity = "hard"
)

descriptor_red_flag_sets <- rbind(
	hard_pairwise_table,
	extreme_vif_table
)

descriptor_red_flag_sets <- deduplicate_predictor_set_table(
	descriptor_red_flag_sets
)

write.csv(
	descriptor_red_flag_sets,
	path_descriptor_red_flag_sets,
	row.names = FALSE
)

############################################################
# D. Generate soft-warning predictor sets
############################################################

soft_pairwise_sets <- all_descriptor_high_correlations[
	all_descriptor_high_correlations$absolute_correlation >=
		soft_pairwise_correlation_threshold &
		all_descriptor_high_correlations$absolute_correlation <
			hard_pairwise_correlation_threshold,
	,
	drop = FALSE
]

soft_pairwise_predictor_sets <- vector(
	"list",
	nrow(soft_pairwise_sets)
)

if (nrow(soft_pairwise_sets) > 0) {
	for (i in seq_len(nrow(soft_pairwise_sets))) {
		soft_pairwise_predictor_sets[[i]] <- c(
			as.character(soft_pairwise_sets$descriptor_1[i]),
			as.character(soft_pairwise_sets$descriptor_2[i])
		)
	}
}

soft_pairwise_table <- make_predictor_set_table(
	predictor_sets = soft_pairwise_predictor_sets,
	reason = paste0(
		"Pairwise Pearson correlation with |r| >= ",
		soft_pairwise_correlation_threshold,
		" and < ",
		hard_pairwise_correlation_threshold,
		". Used for reporting and interpretation, not automatic exclusion."
	),
	source = "soft_pairwise_correlation",
	severity = "soft"
)

moderate_vif_descriptors <- all_descriptor_vif$descriptor[
	all_descriptor_vif$VIF >= moderate_vif_threshold &
		all_descriptor_vif$VIF < extreme_vif_threshold
]

moderate_vif_predictor_sets <- as.list(
	moderate_vif_descriptors
)

soft_vif_table <- make_predictor_set_table(
	predictor_sets = moderate_vif_predictor_sets,
	reason = paste0(
		"Descriptor VIF >= ",
		moderate_vif_threshold,
		" and < ",
		extreme_vif_threshold,
		". Used for reporting and interpretation, not automatic exclusion."
	),
	source = "moderate_vif",
	severity = "soft"
)

descriptor_soft_warning_sets <- rbind(
	soft_pairwise_table,
	soft_vif_table
)

descriptor_soft_warning_sets <- deduplicate_predictor_set_table(
	descriptor_soft_warning_sets
)

write.csv(
	descriptor_soft_warning_sets,
	path_descriptor_soft_warning_sets,
	row.names = FALSE
)

############################################################
# E. Targeted solubility/GSE redundancy check
############################################################

solubility_cols <- c(
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce"
)

check_required_columns(
	df,
	solubility_cols,
	data_name = "cleaned dataset"
)

solubility_data <- df[
	stats::complete.cases(df[, solubility_cols, drop = FALSE]),
	solubility_cols,
	drop = FALSE
]

if (nrow(solubility_data) == 0) {
	stop(
		"No complete observations available for solubility/GSE redundancy check.",
		call. = FALSE
	)
}

solubility_data$Mpt_C <- solubility_data$Mptc - 273.15

solubility_data$GSE_logS <- 0.5 -
	0.01 * (solubility_data$Mpt_C - 25) -
	solubility_data$logKowb

write.csv(
	solubility_data,
	path_solubility_gse_dataset,
	row.names = FALSE
)

solubility_correlation_cols <- c(
	"logKowb",
	"Mptc",
	"Mpt_C",
	"LogSaqd",
	"LogSoce",
	"GSE_logS"
)

solubility_gse_correlation <- stats::cor(
	solubility_data[, solubility_correlation_cols, drop = FALSE],
	use = "complete.obs",
	method = "pearson"
)

write.csv(
	solubility_gse_correlation,
	path_solubility_gse_correlation_matrix
)

plot_correlation_heatmap(
	correlation_matrix = solubility_gse_correlation,
	title = "Solubility-related descriptor correlation structure",
	output_png = path_fig_solubility_correlation_png,
	output_pdf = path_fig_solubility_correlation_pdf
)

fit_logsaqd_from_logkow_mpt <- stats::lm(
	LogSaqd ~ logKowb + Mptc,
	data = solubility_data
)

fit_logsoce_from_logkow_mpt <- stats::lm(
	LogSoce ~ logKowb + Mptc,
	data = solubility_data
)

fit_logsaqd_from_gse <- stats::lm(
	LogSaqd ~ GSE_logS,
	data = solubility_data
)

fit_logsoce_from_gse <- stats::lm(
	LogSoce ~ GSE_logS,
	data = solubility_data
)

solubility_gse_regression_summary <- rbind(
	summarize_lm(
		"LogSaqd ~ logKowb + Mptc",
		fit_logsaqd_from_logkow_mpt
	),
	summarize_lm(
		"LogSoce ~ logKowb + Mptc",
		fit_logsoce_from_logkow_mpt
	),
	summarize_lm(
		"LogSaqd ~ GSE_logS",
		fit_logsaqd_from_gse
	),
	summarize_lm(
		"LogSoce ~ GSE_logS",
		fit_logsoce_from_gse
	)
)

write.csv(
	solubility_gse_regression_summary,
	path_solubility_gse_regression_summary,
	row.names = FALSE
)

solubility_gse_vif <- calculate_vif(
	data = solubility_data,
	predictor_cols = solubility_cols
)

solubility_gse_vif <- solubility_gse_vif[
	order(-solubility_gse_vif$VIF),
	,
	drop = FALSE
]

write.csv(
	solubility_gse_vif,
	path_solubility_gse_vif,
	row.names = FALSE
)

gse_plot <- ggplot2::ggplot(
	solubility_data,
	ggplot2::aes(
		x = GSE_logS,
		y = LogSaqd
	)
) +
	ggplot2::geom_point(
		alpha = 0.65,
		size = 1.8
	) +
	ggplot2::geom_smooth(
		method = "lm",
		se = TRUE,
		linewidth = 0.7
	) +
	ggplot2::labs(
		title = "Reported aqueous solubility versus GSE-derived estimate",
		x = "GSE-derived logS",
		y = "Reported LogSaqd"
	) +
	ggplot2::theme_minimal(base_size = 11) +
	ggplot2::theme(
		plot.title = ggplot2::element_text(face = "bold")
	)

ggplot2::ggsave(
	filename = path_fig_gse_logS_vs_LogSaqd_png,
	plot = gse_plot,
	width = 6.5,
	height = 5,
	dpi = 300
)

ggplot2::ggsave(
	filename = path_fig_gse_logS_vs_LogSaqd_pdf,
	plot = gse_plot,
	width = 6.5,
	height = 5
)

############################################################
# F. Manuscript-ready descriptor redundancy summary
############################################################

max_abs_pairwise_correlation <- max(
	abs(all_descriptor_correlation[upper.tri(all_descriptor_correlation)]),
	na.rm = TRUE
)

max_full_descriptor_vif <- max(
	all_descriptor_vif$VIF,
	na.rm = TRUE
)

n_hard_red_flag_sets <- nrow(descriptor_red_flag_sets)
n_soft_warning_sets <- nrow(descriptor_soft_warning_sets)

table_descriptor_redundancy_summary <- data.frame(
	Analysis = c(
		"Candidate descriptors screened",
		"Complete observations used for redundancy screening",
		"Maximum absolute pairwise correlation",
		"Descriptor pairs with |Pearson r| >= 0.85",
		"Descriptor pairs with 0.70 <= |Pearson r| < 0.85",
		"Maximum VIF in full candidate descriptor pool",
		"Hard red-flag predictor sets generated",
		"Soft-warning predictor sets generated",
		"Correlation: LogSaqd vs GSE_logS",
		"Correlation: LogSoce vs GSE_logS",
		"Regression R2: LogSaqd ~ logKowb + Mptc",
		"Regression R2: LogSoce ~ logKowb + Mptc"
	),
	Result = c(
		length(candidate_descriptors),
		nrow(descriptor_data),
		round(max_abs_pairwise_correlation, 3),
		sum(all_descriptor_high_correlations$absolute_correlation >= hard_pairwise_correlation_threshold),
		sum(
			all_descriptor_high_correlations$absolute_correlation >= soft_pairwise_correlation_threshold &
				all_descriptor_high_correlations$absolute_correlation <
					hard_pairwise_correlation_threshold
		),
		round(max_full_descriptor_vif, 3),
		n_hard_red_flag_sets,
		n_soft_warning_sets,
		round(solubility_gse_correlation["LogSaqd", "GSE_logS"], 3),
		round(solubility_gse_correlation["LogSoce", "GSE_logS"], 3),
		round(
			solubility_gse_regression_summary$r_squared[
				solubility_gse_regression_summary$model ==
					"LogSaqd ~ logKowb + Mptc"
			],
			3
		),
		round(
			solubility_gse_regression_summary$r_squared[
				solubility_gse_regression_summary$model ==
					"LogSoce ~ logKowb + Mptc"
			],
			3
		)
	),
	Interpretation = c(
		"Number of candidate molecular and experimental-condition descriptors evaluated.",
		"Rows with complete values for all candidate descriptors.",
		"Strongest pairwise linear descriptor relationship observed in the candidate pool.",
		"Pairs treated as hard red-flag combinations for model-search exclusion.",
		"Pairs retained for model search but flagged for interpretation.",
		"Strongest multivariable redundancy diagnostic in the candidate descriptor pool.",
		"Combinations exported for automatic exclusion from candidate formulas.",
		"Combinations exported for reporting and cautious interpretation only.",
		"Agreement between reported aqueous solubility and GSE-derived solubility estimate.",
		"Agreement between reported octanol solubility and GSE-derived aqueous solubility estimate.",
		"Fraction of variation in aqueous solubility explained by lipophilicity and melting point.",
		"Fraction of variation in octanol solubility explained by lipophilicity and melting point."
	),
	check.names = FALSE,
	stringsAsFactors = FALSE
)

write.csv(
	table_descriptor_redundancy_summary,
	path_table_descriptor_redundancy_summary,
	row.names = FALSE
)

write.csv(
	table_descriptor_redundancy_summary,
	path_manuscript_table_descriptor_redundancy_summary,
	row.names = FALSE
)

write.csv(
	table_descriptor_redundancy_summary,
	path_solubility_redundancy_table,
	row.names = FALSE
)

write.csv(
	table_descriptor_redundancy_summary,
	path_table_solubility_redundancy,
	row.names = FALSE
)

############################################################
# Console summary
############################################################

cat("\nDescriptor redundancy screening complete.\n\n")

cat("Main descriptor screening outputs:\n")
cat("  Correlation matrix: ", path_all_descriptor_correlation_matrix, "\n", sep = "")
cat("  High-correlation table: ", path_all_descriptor_high_correlations, "\n", sep = "")
cat("  VIF table: ", path_all_descriptor_vif, "\n", sep = "")
cat("  Red-flag predictor sets: ", path_descriptor_red_flag_sets, "\n", sep = "")
cat("  Soft-warning predictor sets: ", path_descriptor_soft_warning_sets, "\n", sep = "")
cat("  Descriptor redundancy summary: ", path_all_descriptor_redundancy_summary, "\n", sep = "")
cat("  Manuscript summary table: ", path_manuscript_table_descriptor_redundancy_summary, "\n\n", sep = "")

cat("Main figure outputs:\n")
cat("  Figure 2 PNG: ", path_fig_descriptor_correlation_png, "\n", sep = "")
cat("  Figure 2 PDF: ", path_fig_descriptor_correlation_pdf, "\n", sep = "")
cat("  Manuscript Figure 2 PNG: ", path_manuscript_fig_descriptor_correlation_png, "\n", sep = "")
cat("  Manuscript Figure 2 PDF: ", path_manuscript_fig_descriptor_correlation_pdf, "\n", sep = "")
cat("  VIF figure PNG: ", path_fig_descriptor_vif_png, "\n", sep = "")
cat("  VIF figure PDF: ", path_fig_descriptor_vif_pdf, "\n\n", sep = "")

cat("Solubility/GSE outputs:\n")
cat("  Solubility/GSE dataset: ", path_solubility_gse_dataset, "\n", sep = "")
cat("  Solubility/GSE correlation matrix: ", path_solubility_gse_correlation_matrix, "\n", sep = "")
cat("  Solubility/GSE regression summary: ", path_solubility_gse_regression_summary, "\n", sep = "")
cat("  Solubility/GSE VIF: ", path_solubility_gse_vif, "\n", sep = "")
cat("  Solubility correlation figure PNG: ", path_fig_solubility_correlation_png, "\n", sep = "")
cat("  Solubility correlation figure PDF: ", path_fig_solubility_correlation_pdf, "\n", sep = "")
cat("  GSE plot PNG: ", path_fig_gse_logS_vs_LogSaqd_png, "\n", sep = "")
cat("  GSE plot PDF: ", path_fig_gse_logS_vs_LogSaqd_pdf, "\n\n", sep = "")

cat("Hard red-flag predictor sets:\n")
print(descriptor_red_flag_sets)

cat("\nSoft-warning predictor sets:\n")
print(descriptor_soft_warning_sets)

cat("\nFull descriptor VIF:\n")
print(all_descriptor_vif)

cat("\nDescriptor redundancy summary:\n")
print(table_descriptor_redundancy_summary)