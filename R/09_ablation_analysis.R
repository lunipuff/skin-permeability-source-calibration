############################################################
# 09_ablation_analysis.R
# Descriptor-group ablation analysis for selected model
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
# Use common complete-case dataset for fair ablation
############################################################

model_df <- df[
	complete.cases(df[, required_cols, drop = FALSE]),
	,
	drop = FALSE
]

############################################################
# Helper functions
############################################################

get_formula_terms <- function(model_formula) {
	attr(
		terms(model_formula),
		"term.labels"
	)
}

get_term_variables <- function(term_label) {
	all.vars(
		as.formula(
			paste(
				"~",
				term_label
			)
		)
	)
}

make_formula_from_terms <- function(term_labels,
									outcome_col = "logkpl") {

	if (length(term_labels) == 0) {
		return(
			as.formula(
				paste(
					outcome_col,
					"~ 1"
				)
			)
		)
	}

	reformulate(
		termlabels = term_labels,
		response = outcome_col
	)
}

remove_terms_containing_variables <- function(model_formula,
											  variables,
											  outcome_col = "logkpl") {

	term_labels <- get_formula_terms(
		model_formula
	)

	keep_terms <- sapply(
		term_labels,
		function(term_label) {
			term_vars <- get_term_variables(
				term_label
			)

			!any(term_vars %in% variables)
		}
	)

	make_formula_from_terms(
		term_labels = term_labels[keep_terms],
		outcome_col = outcome_col
	)
}

remove_interaction_terms <- function(model_formula,
									 outcome_col = "logkpl") {

	term_labels <- get_formula_terms(
		model_formula
	)

	keep_terms <- !grepl(
		":",
		term_labels,
		fixed = TRUE
	)

	make_formula_from_terms(
		term_labels = term_labels[keep_terms],
		outcome_col = outcome_col
	)
}

make_linear_formula_from_selected_model <- function(model_formula,
													outcome_col = "logkpl") {

	predictors <- all.vars(
		model_formula
	)

	predictors <- predictors[
		predictors != outcome_col
	]

	predictors <- unique(
		predictors
	)

	make_formula_from_terms(
		term_labels = predictors,
		outcome_col = outcome_col
	)
}

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
								  model_name,
								  model_formula,
								  ablation_group,
								  ablation_description) {

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	data.frame(
		model = model_name,
		ablation_group = ablation_group,
		ablation_description = ablation_description,
		formula = paste(
			deparse(model_formula),
			collapse = " "
		),
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
# LOCO-CV for one ablation model
############################################################

run_loco_cv_formula_model <- function(data,
									  model_name,
									  model_formula,
									  ablation_group,
									  ablation_description,
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
			model = model_name,
			ablation_group = ablation_group,
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			observed = test[[outcome_col]],
			mu = mu,
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(
		pred_list
	)

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(pred_df$residual)

	summary_df <- summarize_predictions(
		pred_df = pred_df,
		model_name = model_name,
		model_formula = model_formula,
		ablation_group = ablation_group,
		ablation_description = ablation_description
	)

	list(
		prediction = pred_df,
		summary = summary_df
	)
}

############################################################
# Define ablation models
############################################################

ablation_models <- list()

ablation_models[["selected_model"]] <- list(
	formula = selected_formula,
	ablation_group = "none",
	description = "Full selected interpretable QSPR model"
)

ablation_models[["remove_MWa_group"]] <- list(
	formula = remove_terms_containing_variables(
		selected_formula,
		variables = c("MWa"),
		outcome_col = outcome_col
	),
	ablation_group = "MWa",
	description = "Remove molecular-weight terms and interactions involving MWa"
)

ablation_models[["remove_Mptc_group"]] <- list(
	formula = remove_terms_containing_variables(
		selected_formula,
		variables = c("Mptc"),
		outcome_col = outcome_col
	),
	ablation_group = "Mptc",
	description = "Remove melting-point terms"
)

ablation_models[["remove_LogSaqd_group"]] <- list(
	formula = remove_terms_containing_variables(
		selected_formula,
		variables = c("LogSaqd"),
		outcome_col = outcome_col
	),
	ablation_group = "LogSaqd",
	description = "Remove aqueous-solubility terms and interactions involving LogSaqd"
)

ablation_models[["remove_LogSoce_group"]] <- list(
	formula = remove_terms_containing_variables(
		selected_formula,
		variables = c("LogSoce"),
		outcome_col = outcome_col
	),
	ablation_group = "LogSoce",
	description = "Remove organic-solvent solubility term"
)

ablation_models[["remove_Texpi_group"]] <- list(
	formula = remove_terms_containing_variables(
		selected_formula,
		variables = c("Texpi"),
		outcome_col = outcome_col
	),
	ablation_group = "Texpi",
	description = "Remove experimental-temperature term"
)

ablation_models[["remove_interaction_terms"]] <- list(
	formula = remove_interaction_terms(
		selected_formula,
		outcome_col = outcome_col
	),
	ablation_group = "interaction",
	description = "Remove interaction terms only"
)

ablation_models[["linear_selected_predictors"]] <- list(
	formula = make_linear_formula_from_selected_model(
		selected_formula,
		outcome_col = outcome_col
	),
	ablation_group = "nonlinear_terms",
	description = "Use only linear terms for the base predictors in the selected model"
)

############################################################
# Run ablation analysis
############################################################

prediction_list <- list()
summary_list <- list()

for (model_name in names(ablation_models)) {

	cat(
		"Running ablation model:",
		model_name,
		"\n"
	)

	this_model <- ablation_models[[model_name]]

	result <- run_loco_cv_formula_model(
		data = model_df,
		model_name = model_name,
		model_formula = this_model$formula,
		ablation_group = this_model$ablation_group,
		ablation_description = this_model$description,
		id_col = id_col,
		outcome_col = outcome_col
	)

	prediction_list[[model_name]] <- result$prediction
	summary_list[[model_name]] <- result$summary
}

ablation_predictions <- safe_rbind(
	prediction_list
)

ablation_summary <- safe_rbind(
	summary_list
)

############################################################
# Add performance differences relative to selected model
############################################################

selected_rmse <- ablation_summary$RMSE[
	ablation_summary$model == "selected_model"
]

selected_mae <- ablation_summary$MAE[
	ablation_summary$model == "selected_model"
]

selected_r2 <- ablation_summary$R2_pred[
	ablation_summary$model == "selected_model"
]

ablation_summary$delta_RMSE <- ablation_summary$RMSE - selected_rmse
ablation_summary$delta_MAE <- ablation_summary$MAE - selected_mae
ablation_summary$delta_R2_pred <- ablation_summary$R2_pred - selected_r2

ablation_summary <- ablation_summary[
	order(
		ablation_summary$RMSE,
		ablation_summary$MAE
	),
]

ablation_summary$model_rank <- seq_len(
	nrow(ablation_summary)
)

############################################################
# Save outputs
############################################################

write.csv(
	ablation_predictions,
	path_ablation_predictions,
	row.names = FALSE
)

write.csv(
	ablation_summary,
	path_ablation_summary,
	row.names = FALSE
)

write.csv(
	ablation_summary,
	path_table_ablation_summary,
	row.names = FALSE
)

############################################################
# Figure: ablation RMSE
############################################################

make_ablation_rmse_plot <- function() {

	plot_df <- ablation_summary[
		order(
			ablation_summary$RMSE,
			decreasing = TRUE
		),
	]

	par(
		mar = c(5, 14, 3, 2)
	)

	barplot(
		plot_df$RMSE,
		names.arg = plot_df$model,
		horiz = TRUE,
		las = 1,
		xlab = "LOCO-CV RMSE",
		main = "Ablation analysis of selected interpretable model"
	)

	abline(
		v = selected_rmse,
		lty = 2,
		lwd = 2
	)

	legend(
		"bottomright",
		legend = "Selected model RMSE",
		lty = 2,
		lwd = 2,
		bty = "n"
	)
}

pdf(
	path_fig_ablation_rmse_pdf,
	width = 9,
	height = 5
)

make_ablation_rmse_plot()

dev.off()

png(
	path_fig_ablation_rmse_png,
	width = 2250,
	height = 1250,
	res = 250
)

make_ablation_rmse_plot()

dev.off()

############################################################
# Console summary
############################################################

cat("\nAblation analysis complete.\n")
cat("Selected model:\n")
cat(selected_model_text, "\n\n")

cat("Common complete-case observations:", nrow(model_df), "\n")
cat("Common complete-case compounds:", length(unique(model_df$compound_id)), "\n\n")

print(
	ablation_summary[
		,
		c(
			"model",
			"ablation_group",
			"n_observations",
			"n_compounds",
			"RMSE",
			"MAE",
			"R2_pred",
			"delta_RMSE",
			"delta_MAE",
			"delta_R2_pred"
		)
	]
)