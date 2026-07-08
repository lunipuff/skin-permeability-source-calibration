############################################################
# 10_partial_effects.R
# Partial-effect interpretation of selected QSPR model
############################################################

source("R/00_config.R")

############################################################
# Load data and selected model
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

selected_model_text <- readLines(
	path_loco_cv_selected_model,
	warn = FALSE
)

selected_model_text <- selected_model_text[1]

selected_formula <- as.formula(
	selected_model_text
)

############################################################
# User settings
############################################################

outcome_col <- "logkpl"

selected_predictors <- get_formula_predictors(
	selected_formula,
	outcome_col = outcome_col
)

predictor_labels <- c(
	MWa = "Molecular weight",
	logKowb = "Octanol-water partition coefficient",
	Mptc = "Melting point",
	LogSaqd = "Aqueous solubility",
	LogSoce = "Organic-solvent solubility",
	Texpi = "Experimental temperature",
	Skin.thicknessj = "Skin thickness",
	Hdf = "Hydrogen-bond donor count",
	Hag = "Hydrogen-bond acceptor count",
	MVh = "Molecular volume"
)

predictor_units <- c(
	MWa = "Da",
	logKowb = "logKow",
	Mptc = "K",
	LogSaqd = "log mol/mL",
	LogSoce = "log mol/mL",
	Texpi = "K",
	Skin.thicknessj = "mm",
	Hdf = "count",
	Hag = "count",
	MVh = "cm3/mol"
)

############################################################
# Check required columns
############################################################

required_cols <- unique(c(
	outcome_col,
	selected_predictors
))

missing_cols <- required_cols[
	!(required_cols %in% names(df))
]

if (length(missing_cols) > 0) {
	stop(
		paste(
			"Missing required columns:",
			paste(missing_cols, collapse = ", ")
		)
	)
}

############################################################
# Prepare model dataset
############################################################

model_df <- df[
	complete.cases(df[, required_cols, drop = FALSE]),
	,
	drop = FALSE
]

############################################################
# Fit selected model on full cleaned dataset
############################################################

selected_fit <- lm(
	selected_formula,
	data = model_df
)

############################################################
# Helper functions
############################################################

get_label <- function(predictor) {
	if (predictor %in% names(predictor_labels)) {
		return(
			unname(
				predictor_labels[predictor]
			)
		)
	}

	predictor
}

get_unit <- function(predictor) {
	if (predictor %in% names(predictor_units)) {
		return(
			unname(
				predictor_units[predictor]
			)
		)
	}

	""
}

make_axis_label <- function(predictor) {
	unit <- get_unit(predictor)

	if (unit == "") {
		return(
			get_label(predictor)
		)
	}

	paste0(
		get_label(predictor),
		" (",
		unit,
		")"
	)
}

make_reference_row <- function(data, predictors) {
	reference_values <- lapply(
		predictors,
		function(p) {
			median(
				data[[p]],
				na.rm = TRUE
			)
		}
	)

	reference_row <- as.data.frame(
		reference_values,
		stringsAsFactors = FALSE
	)

	names(reference_row) <- predictors

	reference_row
}

safe_rbind <- function(x) {
	x <- x[!sapply(x, is.null)]

	if (length(x) == 0) {
		return(data.frame())
	}

	do.call(
		rbind,
		x
	)
}

make_partial_effect_data <- function(fit,
									 data,
									 predictor,
									 predictors,
									 n_grid = 100,
									 lower_quantile = 0.025,
									 upper_quantile = 0.975) {

	reference_row <- make_reference_row(
		data = data,
		predictors = predictors
	)

	x_range <- quantile(
		data[[predictor]],
		probs = c(lower_quantile, upper_quantile),
		na.rm = TRUE
	)

	x_grid <- seq(
		from = x_range[1],
		to = x_range[2],
		length.out = n_grid
	)

	pred_data <- reference_row[
		rep(1, n_grid),
		,
		drop = FALSE
	]

	pred_data[[predictor]] <- x_grid

	pred_mu <- as.numeric(
		predict(
			fit,
			newdata = pred_data
		)
	)

	out <- data.frame(
		predictor = predictor,
		predictor_label = get_label(predictor),
		predictor_unit = get_unit(predictor),
		x = x_grid,
		predicted_logkpl = pred_mu,
		stringsAsFactors = FALSE
	)

	out
}

make_interaction_effect_data <- function(fit,
										 data,
										 x_predictor,
										 modifier_predictor,
										 predictors,
										 n_grid = 100,
										 lower_quantile = 0.025,
										 upper_quantile = 0.975) {

	reference_row <- make_reference_row(
		data = data,
		predictors = predictors
	)

	x_range <- quantile(
		data[[x_predictor]],
		probs = c(lower_quantile, upper_quantile),
		na.rm = TRUE
	)

	x_grid <- seq(
		from = x_range[1],
		to = x_range[2],
		length.out = n_grid
	)

	modifier_values <- quantile(
		data[[modifier_predictor]],
		probs = c(0.10, 0.50, 0.90),
		na.rm = TRUE
	)

	modifier_labels <- c(
		"10th percentile",
		"median",
		"90th percentile"
	)

	out_list <- vector(
		"list",
		length(modifier_values)
	)

	for (i in seq_along(modifier_values)) {

		pred_data <- reference_row[
			rep(1, n_grid),
			,
			drop = FALSE
		]

		pred_data[[x_predictor]] <- x_grid
		pred_data[[modifier_predictor]] <- as.numeric(modifier_values[i])

		pred_mu <- as.numeric(
			predict(
				fit,
				newdata = pred_data
			)
		)

		out_list[[i]] <- data.frame(
			x_predictor = x_predictor,
			modifier_predictor = modifier_predictor,
			x_predictor_label = get_label(x_predictor),
			modifier_predictor_label = get_label(modifier_predictor),
			x = x_grid,
			modifier_level = modifier_labels[i],
			modifier_value = as.numeric(modifier_values[i]),
			predicted_logkpl = pred_mu,
			stringsAsFactors = FALSE
		)
	}

	do.call(
		rbind,
		out_list
	)
}

############################################################
# Generate partial-effect datasets
############################################################

partial_effect_list <- lapply(
	selected_predictors,
	function(p) {
		make_partial_effect_data(
			fit = selected_fit,
			data = model_df,
			predictor = p,
			predictors = selected_predictors,
			n_grid = 100
		)
	}
)

partial_effect_df <- safe_rbind(
	partial_effect_list
)

############################################################
# Generate interaction-effect dataset
############################################################

if (all(c("MWa", "LogSaqd") %in% selected_predictors)) {
	interaction_effect_df <- make_interaction_effect_data(
		fit = selected_fit,
		data = model_df,
		x_predictor = "MWa",
		modifier_predictor = "LogSaqd",
		predictors = selected_predictors,
		n_grid = 100
	)
} else {
	interaction_effect_df <- data.frame()
}

############################################################
# Save effect datasets
############################################################

write.csv(
	partial_effect_df,
	path_partial_effect_dataset,
	row.names = FALSE
)

write.csv(
	interaction_effect_df,
	path_interaction_effect_dataset,
	row.names = FALSE
)

############################################################
# Figure: partial effects
############################################################

make_partial_effect_plot <- function() {

	n_predictors <- length(selected_predictors)

	n_col <- 3
	n_row <- ceiling(n_predictors / n_col)

	op <- par(
		mfrow = c(n_row, n_col),
		mar = c(4.5, 4.8, 3, 1),
		oma = c(0, 0, 2.5, 0)
	)

	on.exit(
		par(op)
	)

	for (p in selected_predictors) {

		plot_df <- partial_effect_df[
			partial_effect_df$predictor == p,
		]

		raw_x <- model_df[[p]]
		raw_y <- model_df[[outcome_col]]

		plot(
			raw_x,
			raw_y,
			pch = 16,
			cex = 0.8,
			xlab = make_axis_label(p),
			ylab = "logKp",
			main = get_label(p)
		)

		lines(
			plot_df$x,
			plot_df$predicted_logkpl,
			lwd = 2
		)

		rug(
			raw_x,
			side = 1
		)
	}

	if (n_predictors < n_row * n_col) {
		plot.new()

		legend(
			"center",
			legend = c(
				"Observed data",
				"Partial-effect prediction"
			),
			pch = c(16, NA),
			lty = c(NA, 1),
			lwd = c(NA, 2),
			bty = "n",
			cex = 1.1
		)

		if ((n_row * n_col - n_predictors) > 1) {
			for (i in seq_len(n_row * n_col - n_predictors - 1)) {
				plot.new()
			}
		}
	}

	mtext(
		"Partial-effect interpretation of selected QSPR model",
		outer = TRUE,
		cex = 1.3,
		font = 2
	)
}

pdf(
	path_fig_partial_effects_pdf,
	width = 11,
	height = 8.5
)

make_partial_effect_plot()

dev.off()

png(
	path_fig_partial_effects_png,
	width = 2750,
	height = 2125,
	res = 250
)

make_partial_effect_plot()

dev.off()

############################################################
# Figure: MWa x LogSaqd interaction
############################################################

make_interaction_plot <- function() {

	if (nrow(interaction_effect_df) == 0) {
		plot.new()
		text(
			0.5,
			0.5,
			"MWa and LogSaqd are not both included in the selected model."
		)
		return(invisible(NULL))
	}

	x_range <- range(
		interaction_effect_df$x,
		na.rm = TRUE
	)

	y_range <- range(
		interaction_effect_df$predicted_logkpl,
		na.rm = TRUE
	)

	par(
		mar = c(5, 5, 3.5, 2)
	)

	plot(
		NA,
		xlim = x_range,
		ylim = y_range,
		xlab = make_axis_label("MWa"),
		ylab = "Predicted logKp",
		main = "Interaction between molecular weight and aqueous solubility"
	)

	line_types <- c(
		"10th percentile" = 2,
		"median" = 1,
		"90th percentile" = 3
	)

	line_order <- c(
		"10th percentile",
		"median",
		"90th percentile"
	)

	for (level in line_order) {

		level_df <- interaction_effect_df[
			interaction_effect_df$modifier_level == level,
		]

		lines(
			level_df$x,
			level_df$predicted_logkpl,
			lwd = 2,
			lty = line_types[level]
		)
	}

	rug_x <- model_df$MWa

	rug_x <- rug_x[
		rug_x >= x_range[1] &
			rug_x <= x_range[2]
	]

	rug(
		rug_x,
		side = 1
	)

	legend_labels <- sapply(
		line_order,
		function(level) {
			level_value <- unique(
				interaction_effect_df$modifier_value[
					interaction_effect_df$modifier_level == level
				]
			)

			paste0(
				level,
				" LogSaqd = ",
				round(level_value[1], 2)
			)
		}
	)

	legend(
		"topright",
		legend = legend_labels,
		lty = line_types[line_order],
		lwd = 2,
		bty = "n"
	)
}

pdf(
	path_fig_interaction_MWa_LogSaqd_pdf,
	width = 7,
	height = 5
)

make_interaction_plot()

dev.off()

png(
	path_fig_interaction_MWa_LogSaqd_png,
	width = 1800,
	height = 1200,
	res = 250
)

make_interaction_plot()

dev.off()

############################################################
# Console summary
############################################################

cat("Partial-effect analysis complete.\n")
cat("Selected model:\n")
cat(selected_model_text, "\n\n")

cat("Selected predictors:\n")
cat(paste(selected_predictors, collapse = ", "), "\n\n")

cat("Partial-effect dataset:", path_partial_effect_dataset, "\n")
cat("Interaction-effect dataset:", path_interaction_effect_dataset, "\n")
cat("Partial-effect figure:", path_fig_partial_effects_png, "\n")
cat("Interaction figure:", path_fig_interaction_MWa_LogSaqd_png, "\n")