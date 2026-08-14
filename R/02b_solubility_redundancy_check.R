############################################################
# 02b_solubility_redundancy_check.R
# General descriptor redundancy screening and targeted
# solubility/GSE redundancy check
############################################################

source("R/00_config.R")

############################################################
# Required packages
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
# Output paths
############################################################

dir.create("results/descriptor_redundancy", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("manuscript/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("manuscript/tables", showWarnings = FALSE, recursive = TRUE)

if (!exists("path_descriptor_red_flag_sets")) {
	path_descriptor_red_flag_sets <- "results/descriptor_redundancy/descriptor_red_flag_sets.csv"
}

if (!exists("path_descriptor_soft_warning_sets")) {
	path_descriptor_soft_warning_sets <- "results/descriptor_redundancy/descriptor_soft_warning_sets.csv"
}

path_all_descriptor_correlation_matrix <- "results/descriptor_redundancy/all_descriptor_correlation_matrix.csv"
path_all_descriptor_high_correlations <- "results/descriptor_redundancy/all_descriptor_high_correlations.csv"
path_all_descriptor_vif <- "results/descriptor_redundancy/all_descriptor_vif.csv"
path_all_descriptor_summary <- "results/descriptor_redundancy/all_descriptor_redundancy_summary.csv"
path_final_base_predictor_vif <- "results/descriptor_redundancy/final_base_predictor_vif.csv"

path_fig_all_descriptor_correlation_png <- "figures/figure2_descriptor_correlation.png"
path_fig_all_descriptor_correlation_pdf <- "figures/figure2_descriptor_correlation.pdf"
path_manuscript_fig_descriptor_correlation_png <- "manuscript/figures/figure2_descriptor_correlation.png"
path_manuscript_fig_descriptor_correlation_pdf <- "manuscript/figures/figure2_descriptor_correlation.pdf"

path_solubility_dataset <- "results/descriptor_redundancy/solubility_gse_dataset.csv"
path_solubility_correlation_matrix <- "results/descriptor_redundancy/solubility_gse_correlation_matrix.csv"
path_solubility_regression_summary <- "results/descriptor_redundancy/solubility_gse_regression_summary.csv"
path_solubility_vif <- "results/descriptor_redundancy/solubility_gse_vif.csv"

path_table_descriptor_redundancy <- "tables/tableS1_descriptor_redundancy_summary.csv"
path_manuscript_table_descriptor_redundancy <- "manuscript/tables/tableS1_descriptor_redundancy_summary.csv"

path_fig_gse_logS_vs_LogSaqd_png <- "figures/figureS_gse_logS_vs_LogSaqd.png"
path_fig_gse_logS_vs_LogSaqd_pdf <- "figures/figureS_gse_logS_vs_LogSaqd.pdf"
path_fig_solubility_correlation_png <- "figures/figureS_solubility_descriptor_correlation.png"
path_fig_solubility_correlation_pdf <- "figures/figureS_solubility_descriptor_correlation.pdf"

############################################################
# Red-flag thresholds
############################################################

hard_pairwise_correlation_threshold <- 0.85
soft_pairwise_correlation_threshold <- 0.70

extreme_vif_threshold <- 100
moderate_vif_threshold <- 5
extreme_vif_group_correlation_threshold <- 0.50

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
			"Too few complete rows to calculate VIF.",
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

make_high_correlation_table <- function(correlation_matrix, threshold = 0.70) {
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

get_connected_components <- function(nodes, edge_table) {
	if (length(nodes) == 0) {
		return(list())
	}

	remaining <- sort(unique(nodes))
	components <- list()

	while (length(remaining) > 0) {
		current <- remaining[1]
		queue <- current
		component <- character(0)

		while (length(queue) > 0) {
			node <- queue[1]
			queue <- queue[-1]

			if (node %in% component) {
				next
			}

			component <- c(component, node)

			neighbors <- unique(c(
				edge_table$descriptor_2[edge_table$descriptor_1 == node],
				edge_table$descriptor_1[edge_table$descriptor_2 == node]
			))

			neighbors <- neighbors[
				neighbors %in% nodes &
					!(neighbors %in% component)
			]

			queue <- unique(c(queue, neighbors))
		}

		components[[length(components) + 1]] <- sort(unique(component))
		remaining <- setdiff(remaining, component)
	}

	components
}

############################################################
# Load cleaned dataset
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

############################################################
# A. General redundancy screening across all candidate descriptors
############################################################

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
		"No complete rows available for descriptor redundancy screening.",
		call. = FALSE
	)
}

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

all_descriptor_vif <- calculate_vif(
	df,
	candidate_descriptors
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

all_descriptor_summary <- data.frame(
	analysis = c(
		"Number of candidate descriptors screened",
		"Complete observations used for redundancy screening",
		"Number of descriptor pairs with |Pearson r| >= 0.70",
		"Number of descriptor pairs with |Pearson r| >= 0.85",
		"Maximum absolute pairwise correlation",
		"Maximum VIF"
	),
	value = c(
		length(candidate_descriptors),
		nrow(descriptor_data),
		sum(all_descriptor_high_correlations$absolute_correlation >= 0.70),
		sum(all_descriptor_high_correlations$absolute_correlation >= 0.85),
		round(
			max(abs(all_descriptor_correlation[upper.tri(all_descriptor_correlation)])),
			3
		),
		round(max(all_descriptor_vif$VIF, na.rm = TRUE), 3)
	),
	stringsAsFactors = FALSE
)

write.csv(
	all_descriptor_summary,
	path_all_descriptor_summary,
	row.names = FALSE
)

plot_correlation_heatmap(
	all_descriptor_correlation,
	title = "Candidate descriptor correlation structure",
	output_png = path_fig_all_descriptor_correlation_png,
	output_pdf = path_fig_all_descriptor_correlation_pdf
)

file.copy(
	path_fig_all_descriptor_correlation_png,
	path_manuscript_fig_descriptor_correlation_png,
	overwrite = TRUE
)

file.copy(
	path_fig_all_descriptor_correlation_pdf,
	path_manuscript_fig_descriptor_correlation_pdf,
	overwrite = TRUE
)

############################################################
# B. Generate red-flag predictor sets from screening results
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
		"Generated from descriptor screening: empirical pairwise Pearson correlation with |r| >= ",
		hard_pairwise_correlation_threshold,
		"."
	),
	source = "hard_pairwise_correlation",
	severity = "hard"
)

extreme_vif_descriptors <- all_descriptor_vif$descriptor[
	all_descriptor_vif$VIF >= extreme_vif_threshold
]

extreme_vif_edges <- all_descriptor_high_correlations[
	all_descriptor_high_correlations$descriptor_1 %in% extreme_vif_descriptors &
		all_descriptor_high_correlations$descriptor_2 %in% extreme_vif_descriptors &
		all_descriptor_high_correlations$absolute_correlation >=
			extreme_vif_group_correlation_threshold,
	,
	drop = FALSE
]

extreme_vif_components <- get_connected_components(
	nodes = extreme_vif_descriptors,
	edge_table = extreme_vif_edges
)

extreme_vif_components <- extreme_vif_components[
	vapply(
		extreme_vif_components,
		length,
		integer(1)
	) >= 2
]

extreme_vif_group_table <- make_predictor_set_table(
	predictor_sets = extreme_vif_components,
	reason = paste0(
		"Generated from descriptor screening: descriptors with VIF >= ",
		extreme_vif_threshold,
		" connected by pairwise |r| >= ",
		extreme_vif_group_correlation_threshold,
		"."
	),
	source = "extreme_vif_connected_group",
	severity = "hard"
)

descriptor_red_flag_sets <- rbind(
	hard_pairwise_table,
	extreme_vif_group_table
)

descriptor_red_flag_sets <- deduplicate_predictor_set_table(
	descriptor_red_flag_sets
)

write.csv(
	descriptor_red_flag_sets,
	path_descriptor_red_flag_sets,
	row.names = FALSE
)

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
		"Generated from descriptor screening: empirical pairwise Pearson correlation with |r| >= ",
		soft_pairwise_correlation_threshold,
		" and < ",
		hard_pairwise_correlation_threshold,
		"."
	),
	source = "soft_pairwise_correlation",
	severity = "soft"
)

moderate_vif_descriptors <- all_descriptor_vif$descriptor[
	all_descriptor_vif$VIF >= moderate_vif_threshold &
		all_descriptor_vif$VIF < extreme_vif_threshold
]

moderate_vif_sets <- as.list(moderate_vif_descriptors)

soft_vif_table <- make_predictor_set_table(
	predictor_sets = moderate_vif_sets,
	reason = paste0(
		"Generated from descriptor screening: descriptor VIF >= ",
		moderate_vif_threshold,
		" and < ",
		extreme_vif_threshold,
		". Used for reporting only, not automatic model exclusion."
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
# C. Final selected model base-predictor VIF
############################################################

final_base_predictors <- c(
	"MWa",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Texpi"
)

check_required_columns(
	df,
	final_base_predictors,
	data_name = "cleaned dataset"
)

final_base_predictor_vif <- calculate_vif(
	df,
	final_base_predictors
)

final_base_predictor_vif <- final_base_predictor_vif[
	order(-final_base_predictor_vif$VIF),
	,
	drop = FALSE
]

write.csv(
	final_base_predictor_vif,
	path_final_base_predictor_vif,
	row.names = FALSE
)

############################################################
# D. Targeted solubility/GSE redundancy check
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
		"No complete rows available for solubility/GSE redundancy check.",
		call. = FALSE
	)
}

solubility_data$Mpt_C <- solubility_data$Mptc - 273.15

solubility_data$GSE_logS <- 0.5 -
	0.01 * (solubility_data$Mpt_C - 25) -
	solubility_data$logKowb

write.csv(
	solubility_data,
	path_solubility_dataset,
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

solubility_correlation <- stats::cor(
	solubility_data[, solubility_correlation_cols, drop = FALSE],
	use = "complete.obs",
	method = "pearson"
)

write.csv(
	solubility_correlation,
	path_solubility_correlation_matrix
)

plot_correlation_heatmap(
	solubility_correlation,
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

solubility_regression_summary <- rbind(
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
	solubility_regression_summary,
	path_solubility_regression_summary,
	row.names = FALSE
)

solubility_vif <- calculate_vif(
	solubility_data,
	c(
		"logKowb",
		"Mptc",
		"LogSaqd",
		"LogSoce"
	)
)

solubility_vif <- solubility_vif[
	order(-solubility_vif$VIF),
	,
	drop = FALSE
]

write.csv(
	solubility_vif,
	path_solubility_vif,
	row.names = FALSE
)

############################################################
# E. Manuscript-ready descriptor redundancy table
############################################################

table_descriptor_redundancy <- data.frame(
	Analysis = c(
		"Maximum absolute pairwise correlation in full descriptor pool",
		"Number of hard red-flag predictor sets generated",
		"Largest VIF in full descriptor pool",
		"Largest VIF among retained final-model base predictors",
		"Pearson correlation: LogSaqd vs GSE_logS",
		"Pearson correlation: LogSoce vs GSE_logS",
		"Regression: LogSaqd explained by logKowb and Mptc",
		"Regression: LogSoce explained by logKowb and Mptc"
	),
	Result = c(
		round(
			max(abs(all_descriptor_correlation[upper.tri(all_descriptor_correlation)])),
			3
		),
		nrow(descriptor_red_flag_sets),
		round(max(all_descriptor_vif$VIF, na.rm = TRUE), 3),
		round(max(final_base_predictor_vif$VIF, na.rm = TRUE), 3),
		round(solubility_correlation["LogSaqd", "GSE_logS"], 3),
		round(solubility_correlation["LogSoce", "GSE_logS"], 3),
		round(
			solubility_regression_summary$r_squared[
				solubility_regression_summary$model ==
					"LogSaqd ~ logKowb + Mptc"
			],
			3
		),
		round(
			solubility_regression_summary$r_squared[
				solubility_regression_summary$model ==
					"LogSoce ~ logKowb + Mptc"
			],
			3
		)
	),
	Interpretation = c(
		"Strongest pairwise redundancy observed among candidate descriptors.",
		"Number of empirically generated hard exclusion sets used to filter candidate formulas.",
		"Severe multicollinearity diagnostic for the full candidate descriptor pool.",
		"Multicollinearity diagnostic for the retained base predictors in the final selected model.",
		"Agreement between reported aqueous solubility and GSE-derived solubility estimate.",
		"Agreement between reported octanol solubility and GSE-derived aqueous solubility estimate.",
		"Fraction of variation in reported aqueous solubility explained by lipophilicity and melting point.",
		"Fraction of variation in reported octanol solubility explained by lipophilicity and melting point."
	),
	check.names = FALSE,
	stringsAsFactors = FALSE
)

write.csv(
	table_descriptor_redundancy,
	path_table_descriptor_redundancy,
	row.names = FALSE
)

write.csv(
	table_descriptor_redundancy,
	path_manuscript_table_descriptor_redundancy,
	row.names = FALSE
)

############################################################
# F. GSE plot
############################################################

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
# Console summary
############################################################

cat("\nDescriptor redundancy screening complete.\n\n")

cat("General descriptor redundancy outputs:\n")
cat("  Correlation matrix: ", path_all_descriptor_correlation_matrix, "\n", sep = "")
cat("  High-correlation pairs: ", path_all_descriptor_high_correlations, "\n", sep = "")
cat("  Full descriptor VIF: ", path_all_descriptor_vif, "\n", sep = "")
cat("  Final base-predictor VIF: ", path_final_base_predictor_vif, "\n", sep = "")
cat("  Summary: ", path_all_descriptor_summary, "\n", sep = "")
cat("  Red-flag sets: ", path_descriptor_red_flag_sets, "\n", sep = "")
cat("  Soft-warning sets: ", path_descriptor_soft_warning_sets, "\n", sep = "")
cat("  Figure PNG: ", path_fig_all_descriptor_correlation_png, "\n", sep = "")
cat("  Figure PDF: ", path_fig_all_descriptor_correlation_pdf, "\n", sep = "")
cat("  Manuscript figure PNG: ", path_manuscript_fig_descriptor_correlation_png, "\n", sep = "")
cat("  Manuscript figure PDF: ", path_manuscript_fig_descriptor_correlation_pdf, "\n\n", sep = "")

cat("Solubility/GSE redundancy outputs:\n")
cat("  Solubility dataset: ", path_solubility_dataset, "\n", sep = "")
cat("  Correlation matrix: ", path_solubility_correlation_matrix, "\n", sep = "")
cat("  Regression summary: ", path_solubility_regression_summary, "\n", sep = "")
cat("  Solubility VIF: ", path_solubility_vif, "\n", sep = "")
cat("  Descriptor redundancy table: ", path_table_descriptor_redundancy, "\n", sep = "")
cat("  Manuscript descriptor redundancy table: ", path_manuscript_table_descriptor_redundancy, "\n", sep = "")
cat("  GSE plot PNG: ", path_fig_gse_logS_vs_LogSaqd_png, "\n", sep = "")
cat("  GSE plot PDF: ", path_fig_gse_logS_vs_LogSaqd_pdf, "\n", sep = "")
cat("  Solubility correlation PNG: ", path_fig_solubility_correlation_png, "\n", sep = "")
cat("  Solubility correlation PDF: ", path_fig_solubility_correlation_pdf, "\n\n", sep = "")

cat("General redundancy summary:\n")
print(all_descriptor_summary)

cat("\nGenerated hard red-flag sets:\n")
print(descriptor_red_flag_sets)

cat("\nGenerated soft-warning sets:\n")
print(descriptor_soft_warning_sets)

cat("\nFull descriptor VIF:\n")
print(all_descriptor_vif)

cat("\nFinal selected model base-predictor VIF:\n")
print(final_base_predictor_vif)

cat("\nDescriptor redundancy manuscript table:\n")
print(table_descriptor_redundancy)