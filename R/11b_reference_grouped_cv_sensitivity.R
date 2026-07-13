############################################################
# 11b_reference_grouped_cv_sensitivity.R
# Re-evaluate selected model using leave-one-reference-out CV
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
selected_formula <- as.formula(selected_model_text)

############################################################
# User settings
############################################################

outcome_col <- "logkpl"
id_col <- "compound_id"
reference_col <- "Reference"

############################################################
# Check required columns
############################################################

selected_model_vars <- all.vars(selected_formula)

required_cols <- unique(c(
	id_col,
	reference_col,
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

model_df$row_id <- seq_len(nrow(model_df))

model_df[[reference_col]] <- trimws(
	as.character(model_df[[reference_col]])
)

model_df[[reference_col]][
	is.na(model_df[[reference_col]]) |
	model_df[[reference_col]] == ""
] <- "Unknown reference"

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

			d[, all_cols, drop = FALSE]
		}
	)

	do.call(rbind, x2)
}

summarize_predictions <- function(pred_df, validation_scheme) {
	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	data.frame(
		validation_scheme = validation_scheme,
		n_observations = nrow(pred_df),
		n_compounds = length(unique(pred_df$compound_id)),
		n_references = length(unique(pred_df$Reference)),
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

summarize_by_reference <- function(pred_df) {
	split_pred <- split(
		pred_df,
		pred_df$Reference
	)

	summary_list <- lapply(
		names(split_pred),
		function(ref_name) {
			this_df <- split_pred[[ref_name]]
			this_df$residual <- this_df$observed - this_df$mu
			this_df$abs_error <- abs(this_df$residual)

			data.frame(
				Reference = ref_name,
				n_observations = nrow(this_df),
				n_compounds = length(unique(this_df$compound_id)),
				RMSE = rmse(this_df$observed, this_df$mu),
				MAE = mae(this_df$observed, this_df$mu),
				R2_pred = ifelse(
					nrow(this_df) > 1 &&
						sum((this_df$observed - mean(this_df$observed, na.rm = TRUE))^2, na.rm = TRUE) > 0,
					r2_pred(this_df$observed, this_df$mu),
					NA
				),
				R = ifelse(
					nrow(this_df) > 2,
					cor(this_df$observed, this_df$mu, use = "complete.obs"),
					NA
				),
				mean_error = mean(this_df$residual, na.rm = TRUE),
				median_abs_error = median(this_df$abs_error, na.rm = TRUE),
				n_abs_error_gt_1 = sum(this_df$abs_error > 1, na.rm = TRUE),
				proportion_abs_error_gt_1 = mean(this_df$abs_error > 1, na.rm = TRUE),
				stringsAsFactors = FALSE
			)
		}
	)

	out <- safe_rbind(summary_list)

	out <- out[
		order(out$RMSE, decreasing = TRUE),
		,
		drop = FALSE
	]

	row.names(out) <- NULL
	out
}

############################################################
# Leave-one-reference-out CV
############################################################

run_leave_one_reference_out_cv <- function(
	data,
	model_formula,
	reference_col = "Reference",
	id_col = "compound_id",
	outcome_col = "logkpl"
) {
	references <- unique(data[[reference_col]])

	pred_list <- vector(
		"list",
		length(references)
	)

	failed_folds <- list()

	for (i in seq_along(references)) {
		ref <- references[i]

		cat(
			"Running leave-one-reference-out CV fold",
			i,
			"/",
			length(references),
			":",
			ref,
			"\n"
		)

		training <- data[
			data[[reference_col]] != ref,
			,
			drop = FALSE
		]

		test <- data[
			data[[reference_col]] == ref,
			,
			drop = FALSE
		]

		fit <- tryCatch(
			lm(
				model_formula,
				data = training
			),
			error = function(e) {
				failed_folds[[length(failed_folds) + 1]] <<- data.frame(
					Reference = ref,
					n_test_observations = nrow(test),
					n_test_compounds = length(unique(test[[id_col]])),
					error_message = conditionMessage(e),
					stringsAsFactors = FALSE
				)

				return(NULL)
			}
		)

		if (is.null(fit)) {
			next
		}

		mu <- tryCatch(
			as.numeric(
				predict(
					fit,
					newdata = test
				)
			),
			error = function(e) {
				failed_folds[[length(failed_folds) + 1]] <<- data.frame(
					Reference = ref,
					n_test_observations = nrow(test),
					n_test_compounds = length(unique(test[[id_col]])),
					error_message = conditionMessage(e),
					stringsAsFactors = FALSE
				)

				return(rep(NA_real_, nrow(test)))
			}
		)

		pred_list[[i]] <- data.frame(
			validation_scheme = "leave_one_reference_out",
			fold_id = i,
			Reference = test[[reference_col]],
			row_id = test$row_id,
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			observed = test[[outcome_col]],
			mu = mu,
			reference_seen_in_training = FALSE,
			compound_seen_in_training = test[[id_col]] %in% training[[id_col]],
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(pred_list)

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	failed_df <- safe_rbind(failed_folds)

	list(
		predictions = pred_df,
		failed_folds = failed_df
	)
}

############################################################
# Run leave-one-reference-out CV
############################################################

cat("Running leave-one-reference-out CV on common complete-case dataset.\n")

reference_cv_result <- run_leave_one_reference_out_cv(
	data = model_df,
	model_formula = selected_formula,
	reference_col = reference_col,
	id_col = id_col,
	outcome_col = outcome_col
)

reference_cv_predictions <- reference_cv_result$predictions
reference_cv_failed_folds <- reference_cv_result$failed_folds

reference_cv_predictions <- reference_cv_predictions[
	complete.cases(reference_cv_predictions[, c("observed", "mu")]),
	,
	drop = FALSE
]

reference_cv_overall_summary <- summarize_predictions(
	pred_df = reference_cv_predictions,
	validation_scheme = "leave_one_reference_out"
)

reference_cv_by_reference <- summarize_by_reference(
	reference_cv_predictions
)

############################################################
# Add reference-level training/test overlap diagnostics
############################################################

reference_overlap_list <- lapply(
	unique(model_df[[reference_col]]),
	function(ref) {
		test <- model_df[
			model_df[[reference_col]] == ref,
			,
			drop = FALSE
		]

		training <- model_df[
			model_df[[reference_col]] != ref,
			,
			drop = FALSE
		]

		data.frame(
			Reference = ref,
			n_test_observations = nrow(test),
			n_test_compounds = length(unique(test[[id_col]])),
			proportion_test_compounds_seen_in_training = mean(
				test[[id_col]] %in% training[[id_col]]
			),
			stringsAsFactors = FALSE
		)
	}
)

reference_overlap <- safe_rbind(reference_overlap_list)

reference_cv_by_reference <- merge(
	reference_cv_by_reference,
	reference_overlap,
	by = "Reference",
	all.x = TRUE
)

reference_cv_by_reference <- reference_cv_by_reference[
	order(reference_cv_by_reference$RMSE, decreasing = TRUE),
	,
	drop = FALSE
]

############################################################
# Save outputs
############################################################

write.csv(
	reference_cv_predictions,
	path_reference_cv_predictions,
	row.names = FALSE
)

write.csv(
	reference_cv_overall_summary,
	path_reference_cv_overall_summary,
	row.names = FALSE
)

write.csv(
	reference_cv_by_reference,
	path_reference_cv_by_reference,
	row.names = FALSE
)

if (nrow(reference_cv_failed_folds) > 0) {
	write.csv(
		reference_cv_failed_folds,
		"results/cross_validation/11b_reference_grouped_cv_failed_folds.csv",
		row.names = FALSE
	)
}

############################################################
# Figure: top reference-level RMSE values
############################################################

make_reference_cv_rmse_plot <- function(top_n = 20) {
	plot_df <- reference_cv_by_reference[
		order(reference_cv_by_reference$RMSE, decreasing = TRUE),
		,
		drop = FALSE
	]

	plot_df <- plot_df[
		seq_len(min(top_n, nrow(plot_df))),
		,
		drop = FALSE
	]

	plot_df$Reference <- factor(
		plot_df$Reference,
		levels = rev(plot_df$Reference)
	)

	par(mar = c(5, 11, 3, 2))

	barplot(
		rev(plot_df$RMSE),
		names.arg = rev(as.character(plot_df$Reference)),
		horiz = TRUE,
		las = 1,
		xlab = "RMSE",
		main = "Leave-one-reference-out CV: highest-error references"
	)
}

pdf(
	path_fig_reference_cv_rmse_pdf,
	width = 8,
	height = 7
)
make_reference_cv_rmse_plot()
dev.off()

png(
	path_fig_reference_cv_rmse_png,
	width = 2000,
	height = 1750,
	res = 250
)
make_reference_cv_rmse_plot()
dev.off()

############################################################
# Figure: observed vs predicted
############################################################

make_reference_cv_observed_predicted_plot <- function() {
	plot_range <- range(
		c(
			reference_cv_predictions$observed,
			reference_cv_predictions$mu
		),
		na.rm = TRUE
	)

	par(mar = c(5, 5, 3, 2))

	plot(
		reference_cv_predictions$observed,
		reference_cv_predictions$mu,
		pch = 16,
		xlab = "Observed logKp",
		ylab = "Leave-one-reference-out CV predicted logKp",
		xlim = plot_range,
		ylim = plot_range,
		main = "Leave-one-reference-out cross-validation"
	)

	abline(0, 1, lty = 2, lwd = 2)

	legend(
		"topleft",
		legend = paste0(
			"RMSE = ",
			round(reference_cv_overall_summary$RMSE, 3)
		),
		bty = "n"
	)
}

pdf(
	path_fig_reference_cv_observed_predicted_pdf,
	width = 6,
	height = 6
)
make_reference_cv_observed_predicted_plot()
dev.off()

png(
	path_fig_reference_cv_observed_predicted_png,
	width = 1500,
	height = 1500,
	res = 250
)
make_reference_cv_observed_predicted_plot()
dev.off()

############################################################
# Console summary
############################################################

cat("\nLeave-one-reference-out validation sensitivity analysis complete.\n")
cat("Selected model:\n")
cat(selected_model_text, "\n\n")

cat("Common complete-case observations:", nrow(model_df), "\n")
cat("Common complete-case compounds:", length(unique(model_df[[id_col]])), "\n")
cat("References:", length(unique(model_df[[reference_col]])), "\n\n")

cat("Overall leave-one-reference-out summary:\n")
print(reference_cv_overall_summary)

cat("\nHighest-error references:\n")
print(
	head(
		reference_cv_by_reference[
			,
			c(
				"Reference",
				"n_observations",
				"n_compounds",
				"RMSE",
				"MAE",
				"R2_pred",
				"proportion_test_compounds_seen_in_training"
			)
		],
		10
	)
)

if (nrow(reference_cv_failed_folds) > 0) {
	cat("\nFailed folds were recorded in:\n")
	cat("results/cross_validation/11b_reference_grouped_cv_failed_folds.csv\n")
}

cat("\nOutput tables:\n")
cat("  ", path_reference_cv_predictions, "\n")
cat("  ", path_reference_cv_overall_summary, "\n")
cat("  ", path_reference_cv_by_reference, "\n")

############################################################
# Compile main validation-design sensitivity table
# LOCO-CV, row-wise CV, leave-one-reference-out CV
############################################################

############################################################
# Check required output from script 11
############################################################

if (!file.exists(path_validation_scheme_summary)) {
	stop(
		paste0(
			"Missing validation summary from R/11_rowwise_cv_sensitivity.R: ",
			path_validation_scheme_summary,
			"\nRun R/11_rowwise_cv_sensitivity.R first."
		),
		call. = FALSE
	)
}

if (!file.exists(path_rowwise_cv_predictions)) {
	stop(
		paste0(
			"Missing validation predictions from R/11_rowwise_cv_sensitivity.R: ",
			path_rowwise_cv_predictions,
			"\nRun R/11_rowwise_cv_sensitivity.R first."
		),
		call. = FALSE
	)
}

validation_scheme_summary_11 <- read.csv(
	path_validation_scheme_summary,
	stringsAsFactors = FALSE
)

validation_predictions_11 <- read.csv(
	path_rowwise_cv_predictions,
	stringsAsFactors = FALSE
)

############################################################
# Helper functions
############################################################

pull_validation_row <- function(summary_df, scheme_name) {
	this_row <- summary_df[
		summary_df$validation_scheme == scheme_name,
		,
		drop = FALSE
	]

	if (nrow(this_row) != 1) {
		stop(
			paste0(
				"Expected exactly one row for validation_scheme = '",
				scheme_name,
				"', found ",
				nrow(this_row),
				"."
			),
			call. = FALSE
		)
	}

	this_row
}

safe_column <- function(df, col, default = NA) {
	if (col %in% names(df)) {
		return(df[[col]])
	}

	rep(default, nrow(df))
}

count_unique_fold_groups <- function(df) {
	if (!("fold_id" %in% names(df))) {
		return(NA_integer_)
	}

	if ("repeat_id" %in% names(df)) {
		fold_key <- paste(
			df$repeat_id,
			df$fold_id,
			sep = "_"
		)

		fold_key <- fold_key[
			!is.na(df$fold_id)
		]

		return(length(unique(fold_key)))
	}

	length(unique(df$fold_id))
}

############################################################
# Pull summary rows from script 11
############################################################

loco_row <- pull_validation_row(
	validation_scheme_summary_11,
	"leave_one_compound_out"
)

rowwise_row <- pull_validation_row(
	validation_scheme_summary_11,
	"rowwise_cv_repeats"
)

############################################################
# Pull prediction rows from script 11
############################################################

loco_predictions_11 <- validation_predictions_11[
	validation_predictions_11$validation_scheme == "leave_one_compound_out",
	,
	drop = FALSE
]

rowwise_predictions_11 <- validation_predictions_11[
	validation_predictions_11$validation_scheme == "rowwise_cv",
	,
	drop = FALSE
]

if (nrow(loco_predictions_11) == 0) {
	stop("No LOCO-CV predictions found in script 11 output.", call. = FALSE)
}

if (nrow(rowwise_predictions_11) == 0) {
	stop("No row-wise CV predictions found in script 11 output.", call. = FALSE)
}

############################################################
# Methodological checks for script 11
############################################################

loco_n_folds <- length(unique(loco_predictions_11$fold_id))
loco_n_compounds_from_predictions <- length(unique(loco_predictions_11$compound_id))
loco_n_rows <- nrow(loco_predictions_11)

if (loco_n_folds != loco_n_compounds_from_predictions) {
	stop(
		paste0(
			"LOCO-CV check failed: number of folds = ",
			loco_n_folds,
			", but number of unique compounds = ",
			loco_n_compounds_from_predictions,
			"."
		),
		call. = FALSE
	)
}

if (loco_n_folds == loco_n_rows) {
	message(
		"LOCO-CV check: number of folds equals number of prediction rows. ",
		"This means each held-out compound has one observation in this complete-case dataset."
	)
} else {
	message(
		"LOCO-CV check: number of folds does not equal number of prediction rows. ",
		"This is expected when some compounds have repeated observations. ",
		"Folds equal unique compounds, not rows."
	)
}

rowwise_n_repeats <- length(unique(rowwise_predictions_11$repeat_id))
rowwise_n_fold_groups <- count_unique_fold_groups(rowwise_predictions_11)
rowwise_n_rows <- nrow(rowwise_predictions_11)
rowwise_n_obs_per_repeat <- length(unique(rowwise_predictions_11$row_id))

rowwise_folds_per_repeat <- tapply(
	rowwise_predictions_11$fold_id,
	rowwise_predictions_11$repeat_id,
	function(x) length(unique(x))
)

if (length(unique(rowwise_folds_per_repeat)) != 1) {
	stop(
		"Row-wise CV check failed: not all repeats have the same number of folds.",
		call. = FALSE
	)
}

rowwise_k_folds_detected <- unique(rowwise_folds_per_repeat)

if (rowwise_n_rows != rowwise_n_obs_per_repeat * rowwise_n_repeats) {
	stop(
		paste0(
			"Row-wise CV check failed: prediction rows = ",
			rowwise_n_rows,
			", but n_observations × n_repeats = ",
			rowwise_n_obs_per_repeat * rowwise_n_repeats,
			"."
		),
		call. = FALSE
	)
}

message(
	"Row-wise CV check: detected ",
	rowwise_n_repeats,
	" repeats × ",
	rowwise_k_folds_detected,
	" folds = ",
	rowwise_n_fold_groups,
	" repeat-fold validation groups."
)

############################################################
# Pull leave-one-reference-out summary and predictions
############################################################

if (!exists("reference_cv_overall_summary")) {
	if (!file.exists(path_reference_cv_overall_summary)) {
		stop(
			paste0(
				"Missing leave-one-reference-out summary: ",
				path_reference_cv_overall_summary
			),
			call. = FALSE
		)
	}

	reference_cv_overall_summary <- read.csv(
		path_reference_cv_overall_summary,
		stringsAsFactors = FALSE
	)
}

if (!exists("reference_cv_predictions")) {
	if (!file.exists(path_reference_cv_predictions)) {
		stop(
			paste0(
				"Missing leave-one-reference-out predictions: ",
				path_reference_cv_predictions
			),
			call. = FALSE
		)
	}

	reference_cv_predictions <- read.csv(
		path_reference_cv_predictions,
		stringsAsFactors = FALSE
	)
}

loro_row <- reference_cv_overall_summary

if (nrow(loro_row) != 1) {
	stop(
		paste0(
			"Expected one leave-one-reference-out summary row, found ",
			nrow(loro_row),
			"."
		),
		call. = FALSE
	)
}

############################################################
# Methodological checks for LORO-CV
############################################################

loro_n_folds <- length(unique(reference_cv_predictions$fold_id))
loro_n_references <- length(unique(reference_cv_predictions$Reference))
loro_n_rows <- nrow(reference_cv_predictions)

if (loro_n_folds != loro_n_references) {
	stop(
		paste0(
			"LORO-CV check failed: number of folds = ",
			loro_n_folds,
			", but number of unique references = ",
			loro_n_references,
			"."
		),
		call. = FALSE
	)
}

if (loro_n_folds == loro_n_rows) {
	message(
		"LORO-CV check: number of folds equals number of prediction rows. ",
		"This means each held-out reference has one observation in this complete-case dataset."
	)
} else {
	message(
		"LORO-CV check: number of folds does not equal number of prediction rows. ",
		"This is expected when references contain multiple observations. ",
		"Folds equal unique references, not rows."
	)
}

############################################################
# Build main manuscript table
############################################################

validation_design_sensitivity <- data.frame(
	validation_design = c(
		"LOCO-CV",
		"5-fold row-wise CV",
		"Leave-one-reference-out CV"
	),
	held_out_unit = c(
		"Compound",
		"Observation",
		"Literature reference"
	),
	role = c(
		"Primary compound-level validation",
		"Validation-design sensitivity",
		"Reference-level sensitivity"
	),
	n_observations = c(
		loco_row$n_observations,
		rowwise_row$n_observations,
		loro_row$n_observations
	),
	n_compounds = c(
		loco_row$n_compounds,
		rowwise_row$n_compounds,
		loro_row$n_compounds
	),
	n_validation_units = c(
		loco_n_folds,
		rowwise_n_fold_groups,
		loro_n_folds
	),
	n_repeats = c(
		NA,
		rowwise_n_repeats,
		NA
	),
	RMSE = c(
		loco_row$RMSE,
		rowwise_row$RMSE,
		loro_row$RMSE
	),
	MAE = c(
		loco_row$MAE,
		rowwise_row$MAE,
		loro_row$MAE
	),
	R2_pred = c(
		loco_row$R2_pred,
		rowwise_row$R2_pred,
		loro_row$R2_pred
	),
	median_abs_error = c(
		loco_row$median_abs_error,
		rowwise_row$median_abs_error,
		loro_row$median_abs_error
	),
	proportion_abs_error_gt_1 = c(
		loco_row$proportion_abs_error_gt_1,
		rowwise_row$proportion_abs_error_gt_1,
		loro_row$proportion_abs_error_gt_1
	),
	proportion_compound_seen_in_training = c(
		safe_column(loco_row, "proportion_compound_seen_in_training", 0),
		safe_column(rowwise_row, "proportion_compound_seen_in_training"),
		mean(reference_cv_predictions$compound_seen_in_training, na.rm = TRUE)
	),
	stringsAsFactors = FALSE
)

############################################################
# Round numeric columns
############################################################

numeric_cols <- vapply(
	validation_design_sensitivity,
	is.numeric,
	logical(1)
)

validation_design_sensitivity[numeric_cols] <- lapply(
	validation_design_sensitivity[numeric_cols],
	function(x) {
		round(x, 3)
	}
)

############################################################
# Save output
############################################################

write.csv(
	validation_design_sensitivity,
	path_validation_design_sensitivity_table,
	row.names = FALSE
)

cat("\nMain validation-design sensitivity table written to:\n")
cat("  ", path_validation_design_sensitivity_table, "\n\n")

cat("Main validation-design sensitivity table:\n")
print(validation_design_sensitivity)