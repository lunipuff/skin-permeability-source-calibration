source("R/00_config.R")

############################################################
# Figure: observed vs predicted LOCO-CV predictions
############################################################

if (!requireNamespace("PNWColors", quietly = TRUE)) {
	install.packages("PNWColors")
}

path_loco_cv_predictions <- file.path(
	"results",
	"cross_validation",
	"03_loco_cv_predictions.csv"
)

if (!file.exists(path_loco_cv_predictions)) {
	path_loco_cv_predictions <- file.path(
		"results",
		"cross_validation",
		"03_loco_cv_final_model_predictions.csv"
	)
}

if (!file.exists(path_loco_cv_predictions)) {
	stop(
		"Could not find LOCO-CV prediction file.",
		call. = FALSE
	)
}

est <- read.csv(
	path_loco_cv_predictions,
	stringsAsFactors = FALSE
)

required_cols <- c(
	"logkpl",
	"mu"
)

missing_cols <- required_cols[
	!(required_cols %in% names(est))
]

if (length(missing_cols) > 0) {
	stop(
		"Missing required columns in LOCO-CV prediction file: ",
		paste(missing_cols, collapse = ", "),
		call. = FALSE
	)
}

rmse <- function(observed, predicted) {
	sqrt(
		mean(
			(observed - predicted)^2,
			na.rm = TRUE
		)
	)
}

mae <- function(observed, predicted) {
	mean(
		abs(observed - predicted),
		na.rm = TRUE
	)
}

r2_pred <- function(observed, predicted) {
	1 -
		sum(
			(observed - predicted)^2,
			na.rm = TRUE
		) /
		sum(
			(observed - mean(observed, na.rm = TRUE))^2,
			na.rm = TRUE
		)
}

make_observed_predicted_plot <- function() {
	plot_range <- range(
		c(est$logkpl, est$mu),
		na.rm = TRUE
	)

	residual <- est$logkpl - est$mu

	palette_values <- PNWColors::pnw_palette(
		"Shuksan2",
		n = 101
	)

	max_abs_residual <- max(
		abs(residual),
		na.rm = TRUE
	)

	color_index <- round(
		(residual + max_abs_residual) /
			(2 * max_abs_residual) *
			100
	) + 1

	color_index <- pmax(
		1,
		pmin(
			101,
			color_index
		)
	)

	point_colors <- palette_values[color_index]

	legend_residuals <- c(-1,-0.5,0,0.5,1)

	legend_index <- round(
		(legend_residuals + max_abs_residual) /
			(2 * max_abs_residual) *
			100
	) + 1

	legend_index <- pmax(
		1,
		pmin(
			101,
			legend_index
		)
	)

	legend_colors <- palette_values[legend_index]

	par(mar = c(5, 5, 3, 2))

	plot(
		est$logkpl,
		est$mu,
		pch = 21,
		bg = point_colors,
		col = "black",
		lwd = 0.6,
		xlab = "Observed logKp",
		ylab = "LOCO-CV predicted logKp",
		xlim = plot_range,
		ylim = plot_range,
		main = "Selected model prediction performance"
	)

	abline(
		0,
		1,
		lty = 2
	)

	legend(
		"topleft",
		legend = sprintf(
			"= %.2f",
			legend_residuals
		),
		pt.bg = legend_colors,
		col = "black",
		pch = 21,
		pt.cex = 1.1,
		bty = "n",
		title = "Observed - Predicted"
	)

}

grDevices::pdf(
	file = path_fig_loco_cv_prediction_pdf,
	width = 7,
	height = 5
)

make_observed_predicted_plot()

grDevices::dev.off()

grDevices::png(
	filename = path_fig_loco_cv_prediction_png,
	width = 7,
	height = 5,
	units = "in",
	res = 300
)

make_observed_predicted_plot()

grDevices::dev.off()

############################################################
# 99b_make_manuscript_figures.R
# Generate manuscript and supplementary diagnostic figures
############################################################

source("R/00_config.R")

dir.create(
	"figures",
	showWarnings = FALSE,
	recursive = TRUE
)

############################################################
# Helpers
############################################################

save_base_plot <- function(plot_function, pdf_path, png_path, width, height, png_res = 250) {
	grDevices::pdf(
		file = pdf_path,
		width = width,
		height = height
	)

	plot_function()

	grDevices::dev.off()

	grDevices::png(
		filename = png_path,
		width = width,
		height = height,
		units = "in",
		res = png_res
	)

	plot_function()

	grDevices::dev.off()
}

format_term_label <- function(x) {
	x <- gsub("I\\(Mptc\\^2\\)", "Mptc^2", x)
	x <- gsub("I\\(LogSaqd\\^2\\)", "LogSaqd^2", x)
	x <- gsub("MWa:LogSaqd", "MWa x LogSaqd", x)
	x
}

############################################################
# Supplementary Figure: coefficient stability across folds
############################################################

selected_model <- readLines(
	path_loco_cv_selected_model,
	warn = FALSE
)

selected_model <- trimws(
	gsub(
		"\\s+",
		" ",
		paste(selected_model, collapse = " ")
	)
)

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

if (!("compound_id" %in% names(df))) {
	stop(
		"Cleaned dataset must contain compound_id.",
		call. = FALSE
	)
}

coef_long <- data.frame()

compound_ids <- unique(df$compound_id)

for (compound_id_i in compound_ids) {
	training <- df[
		df$compound_id != compound_id_i,
		,
		drop = FALSE
	]

	fit <- lm(
		as.formula(selected_model),
		data = training
	)

	fold_coef <- coef(fit)

	coef_long <- rbind(
		coef_long,
		data.frame(
			left_out_compound_id = compound_id_i,
			term = names(fold_coef),
			coefficient = as.numeric(fold_coef),
			stringsAsFactors = FALSE
		)
	)
}

coef_long <- coef_long[
	coef_long$term != "(Intercept)",
	,
	drop = FALSE
]

coef_long$term_label <- format_term_label(
	coef_long$term
)

term_summary <- aggregate(
	coefficient ~ term + term_label,
	data = coef_long,
	FUN = mean
)

term_summary <- term_summary[
	order(term_summary$coefficient),
	,
	drop = FALSE
]

coef_long$term_label <- factor(
	coef_long$term_label,
	levels = term_summary$term_label
)

make_coefficient_stability_plot <- function() {
	op <- par(
		mar = c(8, 5, 3, 2)
	)

	on.exit(
		par(op),
		add = TRUE
	)

	boxplot(
		coefficient ~ term_label,
		data = coef_long,
		las = 2,
		xlab = "",
		ylab = "Coefficient estimate",
		main = "Selected model coefficient estimates across folds",
		col = "white",
		border = "black"
	)

	abline(
		h = 0,
		lty = 2,
		col = "black"
	)
}

save_base_plot(
	plot_function = make_coefficient_stability_plot,
	pdf_path = path_fig_selected_model_coefficient_boxplot_pdf,
	png_path = path_fig_selected_model_coefficient_boxplot_png,
	width = 7,
	height = 5
)

############################################################
# Supplementary Figures: LOCO-CV residual diagnostics
############################################################

if (!file.exists(path_loco_cv_error_diagnostic_dataset)) {
	stop(
		"Could not find LOCO-CV error diagnostic dataset. Run R/05_loco_cv_error_diagnostics.R first.",
		call. = FALSE
	)
}

diagnostic_df <- read.csv(
	path_loco_cv_error_diagnostic_dataset,
	stringsAsFactors = FALSE
)

required_error_cols <- c(
	"logkpl",
	"loco_mu",
	"loco_residual",
	"loco_abs_error"
)

missing_error_cols <- required_error_cols[
	!(required_error_cols %in% names(diagnostic_df))
]

if (length(missing_error_cols) > 0) {
	stop(
		"Missing required columns in LOCO-CV error diagnostic dataset: ",
		paste(missing_error_cols, collapse = ", "),
		call. = FALSE
	)
}

selected_model_predictors <- c(
	"MWa",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Texpi"
)

diagnostic_predictors <- selected_model_predictors[
	selected_model_predictors %in% names(diagnostic_df)
]

missing_diagnostic_predictors <- selected_model_predictors[
	!(selected_model_predictors %in% names(diagnostic_df))
]

if (length(missing_diagnostic_predictors) > 0) {
	warning(
		"Selected-model predictors missing from LOCO-CV error diagnostic dataset: ",
		paste(missing_diagnostic_predictors, collapse = ", "),
		call. = FALSE
	)
}

if (length(diagnostic_predictors) == 0) {
	stop(
		"No selected-model predictors were found in the LOCO-CV error diagnostic dataset.",
		call. = FALSE
	)
}

make_residual_by_descriptor_plot <- function() {
	n_plots <- length(diagnostic_predictors)
	n_row <- ceiling(n_plots / 3)

	op <- par(
		mfrow = c(n_row, 3),
		mar = c(4, 4.5, 2, 1),
		oma = c(0, 0, 2, 0)
	)

	on.exit(
		par(op),
		add = TRUE
	)

	for (dp in diagnostic_predictors) {
		x <- diagnostic_df[[dp]]
		y <- diagnostic_df$loco_residual

		plot(
			x,
			y,
			pch = 21,
			bg = "white",
			col = "black",
			cex = 0.7,
			lwd = 0.4,
			xlab = dp,
			ylab = "Prediction error"
		)

		abline(
			h = 0,
			lty = 2,
			col = "black"
		)

		if (sum(complete.cases(x, y)) > 5) {
			lines(
				lowess(x, y),
				lwd = 2,
				col = "black"
			)
		}
	}

	mtext(
		"Prediction error by descriptor",
		outer = TRUE,
		font = 2,
		cex = 1.1
	)
}

make_absolute_error_by_descriptor_plot <- function() {
	n_plots <- length(diagnostic_predictors)
	n_row <- ceiling(n_plots / 3)

	op <- par(
		mfrow = c(n_row, 3),
		mar = c(4, 4.5, 2, 1),
		oma = c(0, 0, 2, 0)
	)

	on.exit(
		par(op),
		add = TRUE
	)

	for (dp in diagnostic_predictors) {
		x <- diagnostic_df[[dp]]
		y <- diagnostic_df$loco_abs_error

		plot(
			x,
			y,
			pch = 21,
			bg = "white",
			col = "black",
			cex = 0.7,
			lwd = 0.4,
			xlab = dp,
			ylab = "Absolute error"
		)

		abline(
			h = 1,
			lty = 2,
			col = "black"
		)

		if (sum(complete.cases(x, y)) > 5) {
			lines(
				lowess(x, y),
				lwd = 2,
				col = "black"
			)
		}
	}

	mtext(
		"Absolute prediction error by descriptor",
		outer = TRUE,
		font = 2,
		cex = 1.1
	)
}

make_error_density_plot <- function() {
	error <- diagnostic_df$loco_residual
	error <- error[!is.na(error)]

	error_sd <- sd(
		error,
		na.rm = TRUE
	)

	d_obs <- density(error)

	x <- seq(
		min(
			error,
			-4 * error_sd,
			na.rm = TRUE
		),
		max(
			error,
			4 * error_sd,
			na.rm = TRUE
		),
		length.out = 1000
	)

	y_expected <- dnorm(
		x,
		mean = 0,
		sd = error_sd
	)

	op <- par(
		mar = c(5, 5, 3, 2)
	)

	on.exit(
		par(op),
		add = TRUE
	)

	plot(
		x,
		y_expected,
		type = "l",
		lwd = 2,
		col = "black",
		ylim = range(
			c(y_expected, d_obs$y),
			na.rm = TRUE
		),
		xlab = "Prediction error",
		ylab = "Density",
		main = "Prediction-error density"
	)

	lines(
		d_obs$x,
		d_obs$y,
		lwd = 2,
		lty = 2,
		col = "black"
	)

	abline(
		v = 0,
		lty = 3,
		col = "black"
	)

	legend(
		"topright",
		legend = c(
			"Normal density using observed error SD",
			"Observed error density"
		),
		lwd = c(2, 2),
		lty = c(1, 2),
		col = "black",
		bty = "n"
	)
}

save_base_plot(
	plot_function = make_residual_by_descriptor_plot,
	pdf_path = path_fig_loco_cv_residual_by_descriptor_pdf,
	png_path = path_fig_loco_cv_residual_by_descriptor_png,
	width = 12,
	height = 9
)

save_base_plot(
	plot_function = make_absolute_error_by_descriptor_plot,
	pdf_path = path_fig_loco_cv_absolute_error_by_descriptor_pdf,
	png_path = path_fig_loco_cv_absolute_error_by_descriptor_png,
	width = 12,
	height = 9
)

save_base_plot(
	plot_function = make_error_density_plot,
	pdf_path = path_fig_loco_cv_error_density_pdf,
	png_path = path_fig_loco_cv_error_density_png,
	width = 7,
	height = 5
)

############################################################
# Supplementary Figure: heteroscedastic uncertainty diagnostics
############################################################

if (!file.exists(path_hetero_cv_predictions_best)) {
	stop(
		"Could not find heteroscedastic prediction file. Run R/06_heteroscedastic_uncertainty_model.R first.",
		call. = FALSE
	)
}

best_pred <- read.csv(
	path_hetero_cv_predictions_best,
	stringsAsFactors = FALSE
)

required_hetero_cols <- c(
	"observed",
	"mu",
	"sigma",
	"lower95",
	"upper95",
	"diff",
	"abs_diff",
	"covered95"
)

missing_hetero_cols <- required_hetero_cols[
	!(required_hetero_cols %in% names(best_pred))
]

if (length(missing_hetero_cols) > 0) {
	stop(
		"Missing required columns in heteroscedastic prediction file: ",
		paste(missing_hetero_cols, collapse = ", "),
		call. = FALSE
	)
}

best_hetero_name <- "Selected variance model"

if (file.exists(path_hetero_selected_model)) {
	best_hetero_name <- readLines(
		path_hetero_selected_model,
		warn = FALSE
	)

	best_hetero_name <- paste(
		best_hetero_name,
		collapse = " "
	)
}

best_pred$high_error <- abs(best_pred$diff) > 1

make_hetero_diagnostic_plot <- function() {
	op <- par(
		mfrow = c(2, 2),
		mar = c(5, 5, 3, 2)
	)

	on.exit(
		par(op),
		add = TRUE
	)

	point_col <- ifelse(
		best_pred$high_error,
		"black",
		"grey70"
	)

	plot_range <- range(
		c(best_pred$observed, best_pred$mu),
		na.rm = TRUE
	)

	plot(
		best_pred$observed,
		best_pred$mu,
		pch = 21,
		bg = point_col,
		col = "black",
		cex = 0.8,
		lwd = 0.4,
		xlab = "Observed logKp",
		ylab = "Predicted logKp",
		xlim = plot_range,
		ylim = plot_range,
		main = paste("Observed vs predicted:", best_hetero_name)
	)

	abline(
		0,
		1,
		lty = 2,
		lwd = 2
	)

	legend(
		"topleft",
		legend = c(
			"|error| <= 1",
			"|error| > 1"
		),
		pt.bg = c(
			"grey70",
			"black"
		),
		col = "black",
		pch = 21,
		bty = "n"
	)

	plot(
		best_pred$sigma,
		best_pred$abs_diff,
		pch = 21,
		bg = point_col,
		col = "black",
		cex = 0.8,
		lwd = 0.4,
		xlab = "Predicted uncertainty",
		ylab = "Absolute prediction error",
		main = "Uncertainty calibration"
	)

	if (sum(complete.cases(best_pred$sigma, best_pred$abs_diff)) > 5) {
		lines(
			lowess(best_pred$sigma, best_pred$abs_diff),
			lwd = 2,
			col = "black"
		)
	}

	plot(
		best_pred$mu,
		best_pred$sigma,
		pch = 21,
		bg = point_col,
		col = "black",
		cex = 0.8,
		lwd = 0.4,
		xlab = "Predicted logKp",
		ylab = "Predicted uncertainty",
		main = "Predicted uncertainty"
	)

	if (sum(complete.cases(best_pred$mu, best_pred$sigma)) > 5) {
		lines(
			lowess(best_pred$mu, best_pred$sigma),
			lwd = 2,
			col = "black"
		)
	}

	ord <- order(best_pred$mu)

	plot(
		seq_along(ord),
		best_pred$observed[ord],
		pch = 21,
		bg = ifelse(
			best_pred$covered95[ord],
			"white",
			"black"
		),
		col = "black",
		cex = 0.7,
		lwd = 0.4,
		ylim = range(
			c(
				best_pred$lower95,
				best_pred$upper95,
				best_pred$observed
			),
			na.rm = TRUE
		),
		xlab = "Observations ordered by predicted logKp",
		ylab = "logKp",
		main = "95% prediction intervals"
	)

	segments(
		x0 = seq_along(ord),
		y0 = best_pred$lower95[ord],
		x1 = seq_along(ord),
		y1 = best_pred$upper95[ord],
		col = "grey70"
	)

	points(
		seq_along(ord),
		best_pred$mu[ord],
		pch = 16,
		col = "black",
		cex = 0.5
	)

	points(
		seq_along(ord),
		best_pred$observed[ord],
		pch = 21,
		bg = ifelse(
			best_pred$covered95[ord],
			"white",
			"black"
		),
		col = "black",
		cex = 0.7,
		lwd = 0.4
	)

	legend(
		"topleft",
		legend = c(
			"Predicted mean",
			"Observed covered",
			"Observed not covered"
		),
		col = "black",
		pt.bg = c(
			"black",
			"white",
			"black"
		),
		pch = c(
			16,
			21,
			21
		),
		bty = "n"
	)
}

save_base_plot(
	plot_function = make_hetero_diagnostic_plot,
	pdf_path = path_fig_hetero_diagnostics_pdf,
	png_path = path_fig_hetero_diagnostics_png,
	width = 10,
	height = 10
)

cat("\nSupplementary diagnostic figures generated:\n")
cat("  ", path_fig_selected_model_coefficient_boxplot_png, "\n", sep = "")
cat("  ", path_fig_loco_cv_residual_by_descriptor_png, "\n", sep = "")
cat("  ", path_fig_loco_cv_absolute_error_by_descriptor_png, "\n", sep = "")
cat("  ", path_fig_loco_cv_error_density_png, "\n", sep = "")
cat("  ", path_fig_hetero_diagnostics_png, "\n", sep = "")

############################################################
# Interaction-effect AB panel: MWa x LogSaqd
############################################################

interaction_ab_png <- file.path(
	"figures",
	"figure_interaction_MWa_LogSaqd_AB.png"
)

interaction_ab_pdf <- file.path(
	"figures",
	"figure_interaction_MWa_LogSaqd_AB.pdf"
)

modeling_data <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE,
	check.names = FALSE
)

selected_formula <- logkpl ~ MWa +
	log(Mptc) +
	LogSaqd +
	LogSoce +
	log(Texpi) +
	I(Mptc^2) +
	I(LogSaqd^2) +
	MWa:LogSaqd

model_variables <- all.vars(selected_formula)

interaction_data <- modeling_data[
	complete.cases(modeling_data[, model_variables]),
	,
	drop = FALSE
]

selected_fit <- lm(
	selected_formula,
	data = interaction_data
)

predictors_for_reference <- c(
	"MWa",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Texpi"
)

reference_values <- vapply(
	interaction_data[, predictors_for_reference, drop = FALSE],
	median,
	numeric(1),
	na.rm = TRUE
)

mwa_grid <- seq(
	quantile(
		interaction_data$MWa,
		0.025,
		na.rm = TRUE
	),
	quantile(
		interaction_data$MWa,
		0.975,
		na.rm = TRUE
	),
	length.out = 160
)

logsaqd_grid <- seq(
	quantile(
		interaction_data$LogSaqd,
		0.025,
		na.rm = TRUE
	),
	quantile(
		interaction_data$LogSaqd,
		0.975,
		na.rm = TRUE
	),
	length.out = 160
)

interaction_grid <- expand.grid(
	MWa = mwa_grid,
	LogSaqd = logsaqd_grid
)

for (current_predictor in names(reference_values)) {
	if (!(current_predictor %in% names(interaction_grid))) {
		interaction_grid[[current_predictor]] <- reference_values[[current_predictor]]
	}
}

interaction_grid$predicted_logkpl <- predict(
	selected_fit,
	newdata = interaction_grid
)

prediction_matrix <- matrix(
	interaction_grid$predicted_logkpl,
	nrow = length(mwa_grid),
	ncol = length(logsaqd_grid)
)

interaction_points <- interaction_data[
	interaction_data$MWa >= min(mwa_grid) &
	interaction_data$MWa <= max(mwa_grid) &
	interaction_data$LogSaqd >= min(logsaqd_grid) &
	interaction_data$LogSaqd <= max(logsaqd_grid),
	,
	drop = FALSE
]

rug_mwa <- interaction_points$MWa
rug_logsaqd <- interaction_points$LogSaqd

logsaqd_slices <- as.numeric(
	quantile(
		interaction_data$LogSaqd,
		probs = c(
			0.10,
			0.50,
			0.90
		),
		na.rm = TRUE
	)
)

names(logsaqd_slices) <- c(
	"Low aqueous solubility",
	"Median aqueous solubility",
	"High aqueous solubility"
)

slice_prediction_data <- do.call(
	rbind,
	lapply(
		names(logsaqd_slices),
		function(slice_name) {
			new_data <- data.frame(
				MWa = mwa_grid,
				Mptc = reference_values[["Mptc"]],
				LogSaqd = logsaqd_slices[[slice_name]],
				LogSoce = reference_values[["LogSoce"]],
				Texpi = reference_values[["Texpi"]]
			)

			data.frame(
				MWa = mwa_grid,
				LogSaqd_level = slice_name,
				LogSaqd_value = logsaqd_slices[[slice_name]],
				predicted_logkpl = predict(
					selected_fit,
					newdata = new_data
				),
				stringsAsFactors = FALSE
			)
		}
	)
)

interaction_palette <- PNWColors::pnw_palette(
	"Bay",
	n = 100
)

zlim_shared <- range(
	c(
		prediction_matrix,
		interaction_points$logkpl,
		slice_prediction_data$predicted_logkpl
	),
	finite = TRUE
)

map_to_palette <- function(values, zlim, palette_values) {
	color_index <- round(
		(values - zlim[1]) /
			diff(zlim) *
			(length(palette_values) - 1)
	) + 1

	color_index <- pmax(
		1,
		pmin(
			length(palette_values),
			color_index
		)
	)

	palette_values[color_index]
}

interaction_point_fill <- map_to_palette(
	values = interaction_points$logkpl,
	zlim = zlim_shared,
	palette_values = interaction_palette
)

draw_horizontal_colorbar <- function(zlim, palette_values) {
	par(
		mar = c(2.6, 12, 0.8, 8)
	)

	plot(
		NA,
		xlim = zlim,
		ylim = c(0, 1),
		type = "n",
		xaxs = "i",
		yaxs = "i",
		axes = FALSE,
		xlab = "",
		ylab = ""
	)

	breaks <- seq(
		zlim[1],
		zlim[2],
		length.out = length(palette_values) + 1
	)

	for (i in seq_along(palette_values)) {
		rect(
			xleft = breaks[i],
			ybottom = 0,
			xright = breaks[i + 1],
			ytop = 1,
			col = palette_values[i],
			border = NA
		)
	}

	box()

	colorbar_ticks <- pretty(
		zlim,
		n = 6
	)

	axis(
		side = 1,
		at = colorbar_ticks,
		labels = sprintf("%.0f", colorbar_ticks),
		cex.axis = 0.85
	)

	mtext(
		"logKp",
		side = 3,
		line = 0.2,
		cex = 0.95
	)
}

make_panel_b_background <- function(x_range, y_range, zlim, palette_values) {
	y_breaks <- seq(
		y_range[1],
		y_range[2],
		length.out = 101
	)

	y_centers <- (
		y_breaks[-1] +
		y_breaks[-length(y_breaks)]
	) / 2

	bg_colors <- map_to_palette(
		values = y_centers,
		zlim = zlim,
		palette_values = palette_values
	)

	for (i in seq_along(y_centers)) {
		rect(
			xleft = x_range[1],
			ybottom = y_breaks[i],
			xright = x_range[2],
			ytop = y_breaks[i + 1],
			col = bg_colors[i],
			border = NA
		)
	}
}

make_interaction_ab_plot <- function() {
	layout(
		mat = matrix(
			c(
				1, 2,
				3, 3
			),
			nrow = 2,
			byrow = TRUE
		),
		widths = c(1, 1),
		heights = c(1, 0.15)
	)

	line_types <- c(
		"Low aqueous solubility" = 2,
		"Median aqueous solubility" = 1,
		"High aqueous solubility" = 3
	)

	old_par <- par(
		oma = c(0.5, 0.2, 2.6, 0.2)
	)

	on.exit({
		par(old_par)
		layout(1)
	})

	########################################################
	# Panel A
	########################################################

	par(
		mar = c(5.0, 5.8, 3.2, 1.2)
	)

	image(
		x = mwa_grid,
		y = logsaqd_grid,
		z = prediction_matrix,
		col = interaction_palette,
		zlim = zlim_shared,
		xlab = "Molecular weight (Da)",
		ylab = "Aqueous solubility, LogSaqd (log mol/mL)",
		main = "A. Predicted interaction surface",
		cex.main = 1.15,
		cex.lab = 1.0,
		cex.axis = 0.9,
		useRaster = TRUE
	)

	contour(
		x = mwa_grid,
		y = logsaqd_grid,
		z = prediction_matrix,
		add = TRUE,
		drawlabels = TRUE,
		col = grDevices::adjustcolor(
			"black",
			alpha.f = 0.55
		),
		lwd = 0.7
	)

	points(
		x = interaction_points$MWa,
		y = interaction_points$LogSaqd,
		pch = 21,
		cex = 0.95,
		lwd = 0.8,
		col = "black",
		bg = interaction_point_fill
	)

	rug(
		rug_mwa,
		side = 1,
		col = grDevices::adjustcolor(
			"black",
			alpha.f = 0.35
		)
	)

	rug(
		rug_logsaqd,
		side = 2,
		col = grDevices::adjustcolor(
			"black",
			alpha.f = 0.35
		)
	)

	box()

	########################################################
	# Panel B
	########################################################

	par(
		mar = c(5.0, 5.8, 3.2, 1.2)
	)

	x_range_b <- range(mwa_grid)
	y_range_b <- range(
		slice_prediction_data$predicted_logkpl,
		finite = TRUE
	) + c(-0.5, 0.5)

	plot(
		NA,
		xlim = x_range_b,
		ylim = y_range_b,
		xlab = "Molecular weight (Da)",
		ylab = "Predicted logKp",
		main = "B. Conditional molecular-weight effects",
		cex.main = 1.15,
		cex.lab = 1.0,
		cex.axis = 0.9,
		xaxs = "i",
		yaxs = "i"
	)

	make_panel_b_background(
		x_range = x_range_b,
		y_range = y_range_b,
		zlim = zlim_shared,
		palette_values = interaction_palette
	)

	box()

	for (slice_name in names(logsaqd_slices)) {
		current_slice <- slice_prediction_data[
			slice_prediction_data$LogSaqd_level == slice_name,
			,
			drop = FALSE
		]

		lines(
			current_slice$MWa,
			current_slice$predicted_logkpl,
			lwd = 2.0,
			lty = line_types[[slice_name]],
			col = "black"
		)
	}

	rug(
		rug_mwa,
		side = 1,
		col = grDevices::adjustcolor(
			"black",
			alpha.f = 0.35
		)
	)

	legend(
		"topright",
		inset = 0.01,
		legend = c(
			paste0(
				"Low aqueous solubility, LogSaqd = ",
				sprintf("%.2f", logsaqd_slices[1])
			),
			paste0(
				"Median aqueous solubility, LogSaqd = ",
				sprintf("%.2f", logsaqd_slices[2])
			),
			paste0(
				"High aqueous solubility, LogSaqd = ",
				sprintf("%.2f", logsaqd_slices[3])
			)
		),
		lty = line_types[names(logsaqd_slices)],
		lwd = 2.0,
		col = "black",
		bty = "n",
		cex = 0.88
	)

	########################################################
	# Horizontal color scale
	########################################################

	draw_horizontal_colorbar(
		zlim = zlim_shared,
		palette_values = interaction_palette
	)
	
	mtext(
		"Interaction between molecular weight and aqueous solubility",
		outer = TRUE,
		side = 3,
		line = 0.6,
		font = 2,
		cex = 1.25
	)
}

grDevices::png(
	filename = interaction_ab_png,
	width = 10,
	height = 5.8,
	units = "in",
	res = 300
)

make_interaction_ab_plot()

grDevices::dev.off()

grDevices::pdf(
	file = interaction_ab_pdf,
	width = 10,
	height = 5.8
)

make_interaction_ab_plot()

grDevices::dev.off()

message(
	"Generated interaction-effect AB panel: ",
	interaction_ab_png
)

############################################################
# Supplementary heteroscedastic uncertainty diagnostics
############################################################

hetero_predictions <- read.csv(
	path_hetero_cv_predictions_best,
	stringsAsFactors = FALSE,
	check.names = FALSE
)

required_hetero_cols <- c(
	"observed",
	"mu",
	"sigma"
)

missing_hetero_cols <- required_hetero_cols[
	!(required_hetero_cols %in% names(hetero_predictions))
]

if (length(missing_hetero_cols) > 0) {
	stop(
		"Missing required columns in heteroscedastic prediction file: ",
		paste(missing_hetero_cols, collapse = ", "),
		call. = FALSE
	)
}

hetero_predictions$error <- hetero_predictions$observed -
	hetero_predictions$mu

hetero_predictions$absolute_error <- abs(
	hetero_predictions$error
)

hetero_predictions$lower_95 <- hetero_predictions$mu -
	1.96 * hetero_predictions$sigma

hetero_predictions$upper_95 <- hetero_predictions$mu +
	1.96 * hetero_predictions$sigma

hetero_predictions$covered_95 <- hetero_predictions$observed >=
	hetero_predictions$lower_95 &
	hetero_predictions$observed <=
	hetero_predictions$upper_95

hetero_predictions$high_error <- hetero_predictions$absolute_error > 1

hetero_predictions_ordered <- hetero_predictions[
	order(
		hetero_predictions$mu
	),
	,
	drop = FALSE
]

hetero_predictions_ordered$order_id <- seq_len(
	nrow(hetero_predictions_ordered)
)

make_heteroscedastic_diagnostic_plot <- function() {
	old_par <- par(
		mfrow = c(2, 2),
		mar = c(4.8, 4.8, 3.2, 1.2),
		oma = c(0, 0, 2.2, 0)
	)

	on.exit(par(old_par))

	########################################################
	# Panel A: observed versus predicted
	########################################################

	plot_range <- range(
		c(
			hetero_predictions$observed,
			hetero_predictions$mu
		),
		finite = TRUE
	)

	plot(
		hetero_predictions$observed,
		hetero_predictions$mu,
		xlim = plot_range,
		ylim = plot_range,
		xlab = "Observed logKp",
		ylab = "Predicted logKp",
		main = "Observed versus predicted",
		pch = 21,
		bg = ifelse(
			hetero_predictions$high_error,
			"black",
			"white"
		),
		col = "black",
		cex = 0.8,
		lwd = 0.7
	)

	abline(
		0,
		1,
		lty = 2,
		lwd = 1
	)

	legend(
		"topleft",
		legend = c(
			"|error| <= 1",
			"|error| > 1"
		),
		pch = 21,
		pt.bg = c(
			"white",
			"black"
		),
		col = "black",
		bty = "n",
		cex = 0.8
	)

	########################################################
	# Panel B: uncertainty calibration
	########################################################

	plot(
		hetero_predictions$sigma,
		hetero_predictions$absolute_error,
		xlab = "Predicted sigma",
		ylab = "Absolute prediction error",
		main = "Uncertainty calibration",
		pch = 21,
		bg = "white",
		col = "black",
		cex = 0.8,
		lwd = 0.7
	)

	if (sum(complete.cases(
		hetero_predictions$sigma,
		hetero_predictions$absolute_error
	)) > 5) {
		lines(
			lowess(
				hetero_predictions$sigma,
				hetero_predictions$absolute_error
			),
			lwd = 2
		)
	}

	########################################################
	# Panel C: predicted uncertainty versus predicted mean
	########################################################

	plot(
		hetero_predictions$mu,
		hetero_predictions$sigma,
		xlab = "Predicted logKp",
		ylab = "Predicted sigma",
		main = "Predicted uncertainty",
		pch = 21,
		bg = "white",
		col = "black",
		cex = 0.8,
		lwd = 0.7
	)

	if (sum(complete.cases(
		hetero_predictions$mu,
		hetero_predictions$sigma
	)) > 5) {
		lines(
			lowess(
				hetero_predictions$mu,
				hetero_predictions$sigma
			),
			lwd = 2
		)
	}

	########################################################
	# Panel D: 95% prediction intervals
	########################################################

	plot(
		hetero_predictions_ordered$order_id,
		hetero_predictions_ordered$observed,
		type = "n",
		xlab = "Observations ordered by predicted logKp",
		ylab = "logKp",
		main = "95% prediction intervals",
		ylim = range(
			c(
				hetero_predictions_ordered$lower_95,
				hetero_predictions_ordered$upper_95,
				hetero_predictions_ordered$observed
			),
			finite = TRUE
		)
	)

	segments(
		x0 = hetero_predictions_ordered$order_id,
		y0 = hetero_predictions_ordered$lower_95,
		x1 = hetero_predictions_ordered$order_id,
		y1 = hetero_predictions_ordered$upper_95,
		col = grDevices::adjustcolor(
			"grey50",
			alpha.f = 0.45
		),
		lwd = 0.8
	)

	points(
		hetero_predictions_ordered$order_id,
		hetero_predictions_ordered$mu,
		pch = 16,
		col = "black",
		cex = 0.45
	)

	points(
		hetero_predictions_ordered$order_id,
		hetero_predictions_ordered$observed,
		pch = 21,
		bg = ifelse(
			hetero_predictions_ordered$covered_95,
			"white",
			"black"
		),
		col = "black",
		cex = 0.65,
		lwd = 0.7
	)

	legend(
		"topleft",
		legend = c(
			"Predicted mean",
			"Observed covered",
			"Observed not covered"
		),
		pch = c(
			16,
			21,
			21
		),
		pt.bg = c(
			"black",
			"white",
			"black"
		),
		col = "black",
		bty = "n",
		cex = 0.75
	)

	mtext(
		"Heteroscedastic uncertainty diagnostics",
		outer = TRUE,
		font = 2,
		cex = 1.2
	)
}

grDevices::png(
	filename = path_fig_hetero_diagnostics_png,
	width = 8.5,
	height = 7.0,
	units = "in",
	res = 300
)

make_heteroscedastic_diagnostic_plot()

grDevices::dev.off()

grDevices::pdf(
	file = path_fig_hetero_diagnostics_pdf,
	width = 8.5,
	height = 7.0
)

make_heteroscedastic_diagnostic_plot()

grDevices::dev.off()

message(
	"Generated heteroscedastic uncertainty diagnostics: ",
	path_fig_hetero_diagnostics_png
)

############################################################
# Supplementary LOCO-CV error-density figure
############################################################

error_diagnostic_data <- read.csv(
	path_loco_cv_error_diagnostic_dataset,
	stringsAsFactors = FALSE,
	check.names = FALSE
)

message(
	"LOCO-CV error diagnostic columns: ",
	paste(
		names(error_diagnostic_data),
		collapse = ", "
	)
)

prediction_error <- error_diagnostic_data[["loco_residual"]]
prediction_error <- prediction_error[
	is.finite(prediction_error)
]

make_loco_cv_error_density_plot <- function() {
	error_density <- density(
		prediction_error,
		na.rm = TRUE
	)

	error_mean <- mean(
		prediction_error,
		na.rm = TRUE
	)

	error_sd <- stats::sd(
		prediction_error,
		na.rm = TRUE
	)

	x_range <- range(
		c(
			error_density$x,
			qnorm(
				c(
					0.001,
					0.999
				),
				mean = error_mean,
				sd = error_sd
			)
		),
		finite = TRUE
	)

	x_grid <- seq(
		x_range[1],
		x_range[2],
		length.out = 500
	)

	normal_density <- dnorm(
		x_grid,
		mean = error_mean,
		sd = error_sd
	)

	y_range <- range(
		c(
			error_density$y,
			normal_density
		),
		finite = TRUE
	)

	old_par <- par(
		mar = c(4.8, 4.8, 3.2, 1.2)
	)

	on.exit(par(old_par))

	plot(
		error_density$x,
		error_density$y,
		type = "l",
		lwd = 2.2,
		col = "black",
		xlim = x_range,
		ylim = c(0,1),
		xlab = "Prediction error",
		ylab = "Density",
		main = "Prediction-error density"
	)

	lines(
		x_grid,
		normal_density,
		lwd = 2,
		lty = 2,
		col = "black"
	)

	abline(
		v = 0,
		lty = 3,
		col = "black"
	)

	rug(
		prediction_error,
		col = grDevices::adjustcolor(
			"black",
			alpha.f = 0.35
		)
	)

	legend(
		"topright",
		legend = c(
			"Observed error density",
			"Normal density using observed error SD"
		),
		lty = c(
			1,
			2
		),
		lwd = c(
			2.2,
			2
		),
		col = "black",
		bty = "n",
		cex = 0.85
	)

	box()
}

grDevices::png(
	filename = path_fig_loco_cv_error_density_png,
	width = 6.5,
	height = 4.8,
	units = "in",
	res = 300
)

make_loco_cv_error_density_plot()

grDevices::dev.off()

grDevices::pdf(
	file = path_fig_loco_cv_error_density_pdf,
	width = 6.5,
	height = 4.8
)

make_loco_cv_error_density_plot()

grDevices::dev.off()

message(
	"Generated LOCO-CV error-density figure: ",
	path_fig_loco_cv_error_density_png
)