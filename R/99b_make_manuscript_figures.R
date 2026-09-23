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