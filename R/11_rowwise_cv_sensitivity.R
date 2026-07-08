############################################################
# 11_rowwise_cv_sensitivity.R
# Compare row-wise CV with leave-one-compound-out CV
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
id_col <- "compound_id"

k_folds <- 5
n_repeats <- 50

set.seed(202604)

############################################################
# Check required columns
############################################################

selected_model_vars <- all.vars(
	selected_formula
)

required_cols <- unique(c(
	id_col,
	"CAS.No",
	"Compound",
	selected_model_vars
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
# Prepare common complete-case dataset
############################################################

model_df <- df[
	complete.cases(df[, required_cols, drop = FALSE]),
	,
	drop = FALSE
]

model_df$row_id <- seq_len(
	nrow(model_df)
)

############################################################
# Helper functions
############################################################

safe_rbind <- function(x) {
	x <- x[!sapply(x, is.null)]

	if (length(x) == 0) {
		return(data.frame())
	}

	all_cols <- unique(
		unlist(
			lapply(x, names)
		)
	)

	x2 <- lapply(
		x,
		function(d) {
			missing_cols <- setdiff(
				all_cols,
				names(d)
			)

			if (length(missing_cols) > 0) {
				for (m in missing_cols) {
					d[[m]] <- NA
				}
			}

			d[
				,
				all_cols,
				drop = FALSE
			]
		}
	)

	do.call(
		rbind,
		x2
	)
}

summarize_predictions <- function(pred_df,
								  validation_scheme,
								  repeat_id = NA) {

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	data.frame(
		validation_scheme = validation_scheme,
		repeat_id = repeat_id,
		n_observations = nrow(pred_df),
		n_compounds = length(unique(pred_df$compound_id)),
		RMSE = rmse(pred_df$observed, pred_df$mu),
		MAE = mae(pred_df$observed, pred_df$mu),
		R2_pred = r2_pred(pred_df$observed, pred_df$mu),
		R = cor(
			pred_df$observed,
			pred_df$mu,
			use = "complete.obs"
		),
		mean_error = mean(pred_df$residual, na.rm = TRUE),
		median_abs_error = median(pred_df$abs_error, na.rm = TRUE),
		n_abs_error_gt_1 = sum(pred_df$abs_error > 1, na.rm = TRUE),
		proportion_abs_error_gt_1 = mean(pred_df$abs_error > 1, na.rm = TRUE),
		stringsAsFactors = FALSE
	)
}

############################################################
# Row-wise k-fold CV for one repeat
############################################################

run_rowwise_cv_one_repeat <- function(data,
									  model_formula,
									  repeat_id,
									  k_folds = 5,
									  id_col = "compound_id",
									  outcome_col = "logkpl") {

	n <- nrow(data)

	fold_id <- sample(
		rep(
			seq_len(k_folds),
			length.out = n
		)
	)

	pred_list <- vector(
		"list",
		k_folds
	)

	for (fold in seq_len(k_folds)) {

		training <- data[
			fold_id != fold,
			,
			drop = FALSE
		]

		test <- data[
			fold_id == fold,
			,
			drop = FALSE
		]

		fit <- lm(
			model_formula,
			data = training
		)

		mu <- as.numeric(
			predict(
				fit,
				newdata = test
			)
		)

		compound_seen_in_training <- test[[id_col]] %in% training[[id_col]]

		pred_list[[fold]] <- data.frame(
			validation_scheme = "rowwise_cv",
			repeat_id = repeat_id,
			fold_id = fold,
			row_id = test$row_id,
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			observed = test[[outcome_col]],
			mu = mu,
			compound_seen_in_training = compound_seen_in_training,
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(
		pred_list
	)

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	pred_df
}

############################################################
# Leave-one-compound-out CV on the same complete-case dataset
############################################################

run_loco_cv <- function(data,
						model_formula,
						id_col = "compound_id",
						outcome_col = "logkpl") {

	ids <- unique(
		data[[id_col]]
	)

	pred_list <- vector(
		"list",
		length(ids)
	)

	for (i in seq_along(ids)) {

		id <- ids[i]

		training <- data[
			data[[id_col]] != id,
			,
			drop = FALSE
		]

		test <- data[
			data[[id_col]] == id,
			,
			drop = FALSE
		]

		fit <- lm(
			model_formula,
			data = training
		)

		mu <- as.numeric(
			predict(
				fit,
				newdata = test
			)
		)

		pred_list[[i]] <- data.frame(
			validation_scheme = "leave_one_compound_out",
			repeat_id = NA,
			fold_id = i,
			row_id = test$row_id,
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			observed = test[[outcome_col]],
			mu = mu,
			compound_seen_in_training = FALSE,
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(
		pred_list
	)

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	pred_df
}

############################################################
# Run LOCO-CV
############################################################

cat("Running leave-one-compound-out CV on common complete-case dataset.\n")

loco_predictions <- run_loco_cv(
	data = model_df,
	model_formula = selected_formula,
	id_col = id_col,
	outcome_col = outcome_col
)

loco_summary <- summarize_predictions(
	pred_df = loco_predictions,
	validation_scheme = "leave_one_compound_out",
	repeat_id = NA
)

############################################################
# Run repeated row-wise CV
############################################################

rowwise_prediction_list <- vector(
	"list",
	n_repeats
)

rowwise_summary_list <- vector(
	"list",
	n_repeats
)

for (repeat_id in seq_len(n_repeats)) {

	cat(
		"Running row-wise CV repeat",
		repeat_id,
		"/",
		n_repeats,
		"\n"
	)

	this_pred <- run_rowwise_cv_one_repeat(
		data = model_df,
		model_formula = selected_formula,
		repeat_id = repeat_id,
		k_folds = k_folds,
		id_col = id_col,
		outcome_col = outcome_col
	)

	this_summary <- summarize_predictions(
		pred_df = this_pred,
		validation_scheme = "rowwise_cv",
		repeat_id = repeat_id
	)

	rowwise_prediction_list[[repeat_id]] <- this_pred
	rowwise_summary_list[[repeat_id]] <- this_summary
}

rowwise_predictions <- safe_rbind(
	rowwise_prediction_list
)

rowwise_summary_by_repeat <- safe_rbind(
	rowwise_summary_list
)

############################################################
# Average row-wise predictions across repeats
############################################################

rowwise_prediction_mean <- aggregate(
	cbind(mu, observed) ~ row_id + compound_id + CAS.No + Compound,
	data = rowwise_predictions,
	FUN = mean
)

rowwise_prediction_mean$validation_scheme <- "rowwise_cv_mean"
rowwise_prediction_mean$repeat_id <- NA
rowwise_prediction_mean$compound_seen_in_training <- NA

rowwise_prediction_mean$residual <- rowwise_prediction_mean$observed -
	rowwise_prediction_mean$mu

rowwise_prediction_mean$abs_error <- abs(
	rowwise_prediction_mean$residual
)

rowwise_mean_summary <- summarize_predictions(
	pred_df = rowwise_prediction_mean,
	validation_scheme = "rowwise_cv_mean",
	repeat_id = NA
)

############################################################
# Leakage summary for row-wise CV
############################################################

rowwise_leakage_summary <- aggregate(
	compound_seen_in_training ~ repeat_id,
	data = rowwise_predictions,
	FUN = mean
)

rowwise_leakage_mean <- mean(
	rowwise_leakage_summary$compound_seen_in_training,
	na.rm = TRUE
)

rowwise_leakage_sd <- sd(
	rowwise_leakage_summary$compound_seen_in_training,
	na.rm = TRUE
)

############################################################
# Compile validation-scheme summary
############################################################

rowwise_repeat_summary <- data.frame(
	validation_scheme = "rowwise_cv_repeats",
	repeat_id = NA,
	n_observations = nrow(model_df),
	n_compounds = length(unique(model_df$compound_id)),
	RMSE = mean(rowwise_summary_by_repeat$RMSE, na.rm = TRUE),
	RMSE_sd = sd(rowwise_summary_by_repeat$RMSE, na.rm = TRUE),
	MAE = mean(rowwise_summary_by_repeat$MAE, na.rm = TRUE),
	MAE_sd = sd(rowwise_summary_by_repeat$MAE, na.rm = TRUE),
	R2_pred = mean(rowwise_summary_by_repeat$R2_pred, na.rm = TRUE),
	R2_pred_sd = sd(rowwise_summary_by_repeat$R2_pred, na.rm = TRUE),
	R = mean(rowwise_summary_by_repeat$R, na.rm = TRUE),
	R_sd = sd(rowwise_summary_by_repeat$R, na.rm = TRUE),
	mean_error = mean(rowwise_summary_by_repeat$mean_error, na.rm = TRUE),
	median_abs_error = mean(rowwise_summary_by_repeat$median_abs_error, na.rm = TRUE),
	n_abs_error_gt_1 = mean(rowwise_summary_by_repeat$n_abs_error_gt_1, na.rm = TRUE),
	proportion_abs_error_gt_1 = mean(rowwise_summary_by_repeat$proportion_abs_error_gt_1, na.rm = TRUE),
	proportion_compound_seen_in_training = rowwise_leakage_mean,
	proportion_compound_seen_in_training_sd = rowwise_leakage_sd,
	stringsAsFactors = FALSE
)

loco_summary$RMSE_sd <- NA
loco_summary$MAE_sd <- NA
loco_summary$R2_pred_sd <- NA
loco_summary$R_sd <- NA
loco_summary$proportion_compound_seen_in_training <- 0
loco_summary$proportion_compound_seen_in_training_sd <- NA

rowwise_mean_summary$RMSE_sd <- NA
rowwise_mean_summary$MAE_sd <- NA
rowwise_mean_summary$R2_pred_sd <- NA
rowwise_mean_summary$R_sd <- NA
rowwise_mean_summary$proportion_compound_seen_in_training <- rowwise_leakage_mean
rowwise_mean_summary$proportion_compound_seen_in_training_sd <- rowwise_leakage_sd

validation_scheme_summary <- safe_rbind(
	list(
		loco_summary,
		rowwise_repeat_summary,
		rowwise_mean_summary
	)
)

validation_scheme_summary <- validation_scheme_summary[
	,
	c(
		"validation_scheme",
		"n_observations",
		"n_compounds",
		"RMSE",
		"RMSE_sd",
		"MAE",
		"MAE_sd",
		"R2_pred",
		"R2_pred_sd",
		"R",
		"R_sd",
		"mean_error",
		"median_abs_error",
		"n_abs_error_gt_1",
		"proportion_abs_error_gt_1",
		"proportion_compound_seen_in_training",
		"proportion_compound_seen_in_training_sd"
	)
]

############################################################
# Combine prediction outputs
############################################################

validation_predictions <- safe_rbind(
	list(
		loco_predictions,
		rowwise_predictions
	)
)

############################################################
# Save outputs
############################################################

write.csv(
	validation_predictions,
	path_rowwise_cv_predictions,
	row.names = FALSE
)

write.csv(
	rowwise_summary_by_repeat,
	path_rowwise_cv_summary,
	row.names = FALSE
)

write.csv(
	validation_scheme_summary,
	path_validation_scheme_summary,
	row.names = FALSE
)

############################################################
# Figure: validation-scheme RMSE
############################################################

make_validation_scheme_rmse_plot <- function() {

	plot_df <- validation_scheme_summary[
		validation_scheme_summary$validation_scheme %in% c(
			"leave_one_compound_out",
			"rowwise_cv_repeats"
		),
	]

	plot_df$plot_label <- c(
		"LOCO-CV",
		"Row-wise CV"
	)

	par(
		mar = c(5, 5, 3, 2)
	)

	bar_centers <- barplot(
		plot_df$RMSE,
		names.arg = plot_df$plot_label,
		ylab = "RMSE",
		main = "Validation-scheme sensitivity"
	)

	rowwise_index <- which(
		plot_df$validation_scheme == "rowwise_cv_repeats"
	)

	if (length(rowwise_index) == 1 &&
		!is.na(plot_df$RMSE_sd[rowwise_index])) {

		arrows(
			x0 = bar_centers[rowwise_index],
			y0 = plot_df$RMSE[rowwise_index] -
				plot_df$RMSE_sd[rowwise_index],
			x1 = bar_centers[rowwise_index],
			y1 = plot_df$RMSE[rowwise_index] +
				plot_df$RMSE_sd[rowwise_index],
			angle = 90,
			code = 3,
			length = 0.05
		)
	}

	legend(
		"topright",
		legend = "Row-wise error bar = SD across repeats",
		bty = "n"
	)
}

pdf(
	path_fig_validation_scheme_rmse_pdf,
	width = 6,
	height = 5
)

make_validation_scheme_rmse_plot()

dev.off()

png(
	path_fig_validation_scheme_rmse_png,
	width = 1500,
	height = 1250,
	res = 250
)

make_validation_scheme_rmse_plot()

dev.off()

############################################################
# Figure: row-wise observed vs predicted
############################################################

make_rowwise_observed_predicted_plot <- function() {

	plot_range <- range(
		c(
			rowwise_prediction_mean$observed,
			rowwise_prediction_mean$mu
		),
		na.rm = TRUE
	)

	par(
		mar = c(5, 5, 3, 2)
	)

	plot(
		rowwise_prediction_mean$observed,
		rowwise_prediction_mean$mu,
		pch = 16,
		xlab = "Observed logKp",
		ylab = "Row-wise CV predicted logKp",
		xlim = plot_range,
		ylim = plot_range,
		main = "Row-wise cross-validation"
	)

	abline(
		0,
		1,
		lty = 2,
		lwd = 2
	)

	legend(
		"topleft",
		legend = paste0(
			"RMSE = ",
			round(rowwise_mean_summary$RMSE, 3)
		),
		bty = "n"
	)
}

pdf(
	path_fig_rowwise_observed_predicted_pdf,
	width = 6,
	height = 6
)

make_rowwise_observed_predicted_plot()

dev.off()

png(
	path_fig_rowwise_observed_predicted_png,
	width = 1500,
	height = 1500,
	res = 250
)

make_rowwise_observed_predicted_plot()

dev.off()

############################################################
# Console summary
############################################################

cat("\nRow-wise validation sensitivity analysis complete.\n")
cat("Selected model:\n")
cat(selected_model_text, "\n\n")

cat("Common complete-case observations:", nrow(model_df), "\n")
cat("Common complete-case compounds:", length(unique(model_df$compound_id)), "\n")
cat("Row-wise CV repeats:", n_repeats, "\n")
cat("Row-wise CV folds:", k_folds, "\n\n")

cat("Mean proportion of row-wise test observations whose compound also appeared in training:\n")
cat(round(rowwise_leakage_mean, 3), "\n\n")

cat("Validation-scheme summary:\n")
print(validation_scheme_summary)

cat("\nOutput table:", path_validation_scheme_summary, "\n")