############################################################
# 08_benchmark_models.R
# Compare selected QSPR model with benchmark models
############################################################

source("R/00_config.R")

############################################################
# Load cleaned dataset and selected model
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

selected_model <- readLines(
	path_loco_cv_selected_model,
	warn = FALSE
)

selected_model <- selected_model[1]

############################################################
# User settings
############################################################

outcome_col <- "logkpl"
id_col <- "compound_id"

############################################################
# Extract predictors from selected model
############################################################

selected_model_predictors <- get_formula_predictors(
	selected_model,
	outcome_col = outcome_col
)

linear_selected_formula <- as.formula(
	paste(
		outcome_col,
		"~",
		paste(selected_model_predictors, collapse = " + ")
	)
)

############################################################
# Benchmark model formulas
############################################################

benchmark_formulas <- list(
	null_mean = formula_null,
	potts_guy = formula_potts_guy,
	linear_selected_predictors = linear_selected_formula,
	selected_model = as.formula(selected_model),
	extended_descriptor = formula_extended_descriptor
)

############################################################
# Check required columns
############################################################

required_cols <- unique(c(
	id_col,
	"CAS.No",
	"Compound",
	outcome_col,
	all_predictors,
	unlist(
		lapply(
			benchmark_formulas,
			all.vars
		)
	)
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
# Helper function: safe rbind
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

############################################################
# Helper function: summarize predictions
############################################################

summarize_predictions <- function(pred_df, model_name, formula_text) {
	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	data.frame(
		model = model_name,
		formula = formula_text,
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
# LOCO-CV for one formula-based lm model
############################################################

run_loco_cv_formula_model <- function(data,
									  model_name,
									  model_formula,
									  id_col = "compound_id",
									  outcome_col = "logkpl") {

	model_vars <- all.vars(model_formula)

	required_model_cols <- unique(c(
		id_col,
		"CAS.No",
		"Compound",
		outcome_col,
		model_vars
	))

	model_data <- data[
		complete.cases(data[, required_model_cols, drop = FALSE]),
		,
		drop = FALSE
	]

	ids <- unique(model_data[[id_col]])

	pred_list <- vector(
		"list",
		length(ids)
	)

	for (i in seq_along(ids)) {
		id <- ids[i]

		training <- model_data[
			model_data[[id_col]] != id,
			,
			drop = FALSE
		]

		test <- model_data[
			model_data[[id_col]] == id,
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
			model = model_name,
			df_row = as.integer(rownames(test)),
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			observed = test[[outcome_col]],
			mu = mu,
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(pred_list)

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	formula_text <- paste(
		deparse(model_formula),
		collapse = " "
	)

	summary_df <- summarize_predictions(
		pred_df = pred_df,
		model_name = model_name,
		formula_text = formula_text
	)

	list(
		prediction = pred_df,
		summary = summary_df
	)
}

############################################################
# LOCO-CV random forest benchmark
############################################################

run_loco_cv_random_forest <- function(data,
									  predictors,
									  model_name = "random_forest_core",
									  id_col = "compound_id",
									  outcome_col = "logkpl",
									  ntree = 1000,
									  nodesize = 5) {

	if (!requireNamespace("randomForest", quietly = TRUE)) {
		warning(
			"Package 'randomForest' is not installed. ",
			"Skipping random forest benchmark. Install with install.packages('randomForest')."
		)

		return(NULL)
	}

	required_model_cols <- unique(c(
		id_col,
		"CAS.No",
		"Compound",
		outcome_col,
		predictors
	))

	model_data <- data[
		complete.cases(data[, required_model_cols, drop = FALSE]),
		,
		drop = FALSE
	]

	ids <- unique(model_data[[id_col]])

	pred_list <- vector(
		"list",
		length(ids)
	)

	formula_rf <- as.formula(
		paste(
			outcome_col,
			"~",
			paste(predictors, collapse = " + ")
		)
	)

	mtry_value <- floor(
		sqrt(length(predictors))
	)

	if (mtry_value < 1) {
		mtry_value <- 1
	}

	for (i in seq_along(ids)) {
		id <- ids[i]

		cat(
			"Random forest LOCO-CV fold",
			i,
			"/",
			length(ids),
			"\n"
		)

		training <- model_data[
			model_data[[id_col]] != id,
			,
			drop = FALSE
		]

		test <- model_data[
			model_data[[id_col]] == id,
			,
			drop = FALSE
		]

		set.seed(202604 + i)

		fit <- randomForest::randomForest(
			formula_rf,
			data = training,
			ntree = ntree,
			mtry = mtry_value,
			nodesize = nodesize,
			importance = TRUE
		)

		mu <- as.numeric(
			predict(
				fit,
				newdata = test
			)
		)

		pred_list[[i]] <- data.frame(
			model = model_name,
			df_row = as.integer(rownames(test)),
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			observed = test[[outcome_col]],
			mu = mu,
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(pred_list)

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	formula_text <- paste(
		outcome_col,
		"~",
		paste(predictors, collapse = " + ")
	)

	summary_df <- summarize_predictions(
		pred_df = pred_df,
		model_name = model_name,
		formula_text = formula_text
	)

	summary_df$ntree <- ntree
	summary_df$mtry <- mtry_value
	summary_df$nodesize <- nodesize

	list(
		prediction = pred_df,
		summary = summary_df
	)
}

############################################################
# Run formula-based benchmark models
############################################################

benchmark_prediction_list <- list()
benchmark_summary_list <- list()

for (model_name in names(benchmark_formulas)) {

	cat(
		"Running benchmark model:",
		model_name,
		"\n"
	)

	result <- run_loco_cv_formula_model(
		data = df,
		model_name = model_name,
		model_formula = benchmark_formulas[[model_name]],
		id_col = id_col,
		outcome_col = outcome_col
	)

	benchmark_prediction_list[[model_name]] <- result$prediction
	benchmark_summary_list[[model_name]] <- result$summary
}

############################################################
# Run random forest benchmark
############################################################

rf_result <- run_loco_cv_random_forest(
	data = df,
	predictors = selected_model_predictors,
	model_name = "random_forest_selected_predictors",
	id_col = id_col,
	outcome_col = outcome_col,
	ntree = 1000,
	nodesize = 5
)

if (!is.null(rf_result)) {
	benchmark_prediction_list[["random_forest_selected_predictors"]] <- rf_result$prediction
	benchmark_summary_list[["random_forest_selected_predictors"]] <- rf_result$summary
}

############################################################
# Combine benchmark results
############################################################

benchmark_predictions <- safe_rbind(
	benchmark_prediction_list
)

benchmark_summary <- safe_rbind(
	benchmark_summary_list
)

benchmark_summary <- benchmark_summary[
	order(
		benchmark_summary$RMSE,
		benchmark_summary$MAE
	),
]

benchmark_summary$model_rank <- seq_len(
	nrow(benchmark_summary)
)

############################################################
# Save benchmark results
############################################################

write.csv(
	benchmark_predictions,
	path_benchmark_model_predictions,
	row.names = FALSE
)

write.csv(
	benchmark_summary,
	path_benchmark_model_summary,
	row.names = FALSE
)

write.csv(
	benchmark_summary,
	path_table_benchmark_model_summary,
	row.names = FALSE
)

############################################################
# Figure: benchmark RMSE
############################################################

make_benchmark_rmse_plot <- function() {

	plot_df <- benchmark_summary[
		order(
			benchmark_summary$RMSE,
			decreasing = TRUE
		),
	]

	par(
		mar = c(5, 12, 3, 2)
	)

	barplot(
		plot_df$RMSE,
		names.arg = plot_df$model,
		horiz = TRUE,
		las = 1,
		xlab = "LOCO-CV RMSE",
		main = "Benchmark model comparison"
	)
}

pdf(
	path_fig_benchmark_rmse_pdf,
	width = 8,
	height = 5
)

make_benchmark_rmse_plot()

dev.off()

png(
	path_fig_benchmark_rmse_png,
	width = 2000,
	height = 1250,
	res = 250
)

make_benchmark_rmse_plot()

dev.off()

############################################################
# Figure: observed vs predicted for top benchmark models
############################################################

make_benchmark_observed_predicted_plot <- function() {

	top_models <- benchmark_summary$model[
		seq_len(
			min(4, nrow(benchmark_summary))
		)
	]

	plot_df <- benchmark_predictions[
		benchmark_predictions$model %in% top_models,
	]

	plot_range <- range(
		c(plot_df$observed, plot_df$mu),
		na.rm = TRUE
	)

	op <- par(
		mfrow = c(2, 2),
		mar = c(5, 5, 3, 2)
	)

	on.exit(
		par(op)
	)

	for (model_name in top_models) {

		model_df <- plot_df[
			plot_df$model == model_name,
		]

		plot(
			model_df$observed,
			model_df$mu,
			pch = 16,
			xlab = "Observed logKp",
			ylab = "LOCO-CV predicted logKp",
			xlim = plot_range,
			ylim = plot_range,
			main = model_name
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
				round(
					rmse(model_df$observed, model_df$mu),
					3
				)
			),
			bty = "n"
		)
	}
}

pdf(
	path_fig_benchmark_observed_predicted_pdf,
	width = 9,
	height = 9
)

make_benchmark_observed_predicted_plot()

dev.off()

png(
	path_fig_benchmark_observed_predicted_png,
	width = 2250,
	height = 2250,
	res = 250
)

make_benchmark_observed_predicted_plot()

dev.off()

############################################################
# Console summary
############################################################

cat("Selected model predictors used for linear and random forest benchmarks:\n")
cat(paste(selected_model_predictors, collapse = ", "), "\n")
cat("Benchmark model comparison complete.\n")
cat("Models evaluated:\n")
print(benchmark_summary[, c(
	"model",
	"n_observations",
	"n_compounds",
	"RMSE",
	"MAE",
	"R2_pred",
	"model_rank"
)])

