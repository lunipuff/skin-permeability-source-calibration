############################################################
# 02_initial_predictor_screening.R
# Screen simple interpretable predictor combinations
############################################################

source("R/00_config.R")

############################################################
# Load cleaned dataset
############################################################

df <- read.csv(path_cleaned_dataset)

n_observations <- nrow(df)

############################################################
# Check required columns
############################################################

required_cols <- c("logkpl", all_predictors)

missing_cols <- required_cols[!(required_cols %in% names(df))]

if (length(missing_cols) > 0) {
	stop(
		paste(
			"Missing required columns:",
			paste(missing_cols, collapse = ", ")
		)
	)
}

############################################################
# Prepare predictor combinations
############################################################

variable_predictors <- all_predictors[
	!(all_predictors %in% fixed_predictors)
]

############################################################
# Helper function: check redundant predictor combinations
############################################################

has_redundant_combination <- function(predictors, redundant_sets) {
	for (redundant_set in redundant_sets) {
		if (all(redundant_set %in% predictors)) {
			return(TRUE)
		}
	}

	FALSE
}

############################################################
# Helper function: fit and summarize one model
############################################################

fit_screening_model <- function(predictors, data) {

	formula_text <- paste(
		"logkpl ~",
		paste(predictors, collapse = " + ")
	)

	fit <- lm(
		as.formula(formula_text),
		data = data
	)

	predicted <- predict(
		fit,
		newdata = data
	)

	observed <- data$logkpl

	coef_values <- coef(fit)

	model_r <- cor(
		observed,
		predicted,
		use = "complete.obs"
	)

	result <- data.frame(
		formula = formula_text,
		logLik = as.numeric(logLik(fit)),
		RMSE = rmse(observed, predicted),
		MAE = mae(observed, predicted),
		R = model_r,
		R2 = model_r^2,
		AIC = AIC(fit),
		n_predictors = length(predictors),
		intercept = coef_values["(Intercept)"],
		sigma = sigma(fit),
		sigma2 = sigma(fit)^2,
		stringsAsFactors = FALSE
	)

	for (p in all_predictors) {
		result[[p]] <- ifelse(
			p %in% names(coef_values),
			coef_values[p],
			NA_real_
		)
	}

	result
}

############################################################
# Screen all predictor combinations
############################################################

screening_results <- data.frame()

for (n_extra in 0:length(variable_predictors)) {

	predictor_combinations <- combn(
		variable_predictors,
		n_extra,
		simplify = FALSE
	)

	for (candidate_extra_predictors in predictor_combinations) {

		predictors <- c(
			fixed_predictors,
			candidate_extra_predictors
		)

#		if (has_redundant_combination(predictors, redundant_predictor_sets)) {
#			next
#		}

		model_result <- fit_screening_model(
			predictors = predictors,
			data = df
		)

		screening_results <- rbind(
			screening_results,
			model_result
		)
	}
}

############################################################
# Rank models
############################################################

screening_results <- screening_results[
	order(
		screening_results$RMSE,
		screening_results$AIC
	),
]

screening_results$model_rank <- seq_len(
	nrow(screening_results)
)

############################################################
# Select best model for each predictor count
############################################################

selected_models <- data.frame()

for (n_pred in sort(unique(screening_results$n_predictors))) {

	subset_results <- screening_results[
		screening_results$n_predictors == n_pred,
	]

	subset_results <- subset_results[
		order(
			subset_results$RMSE,
			subset_results$AIC
		),
	]

	selected_models <- rbind(
		selected_models,
		subset_results[1, ]
	)
}

############################################################
# Save results
############################################################

write.csv(
	screening_results,
	path_predictor_screening_full,
	row.names = FALSE
)

write.csv(
	selected_models,
	path_predictor_screening_selected,
	row.names = FALSE
)

############################################################
# Make RMSE screening figure
############################################################

make_predictor_screening_rmse_plot <- function() {

	par(mar = c(5, 5, 3, 2))

	boxplot(
		RMSE ~ n_predictors,
		data = screening_results,
		xlab = "Number of predictors",
		ylab = "Apparent RMSE",
		main = "Initial predictor screening"
	)

	points(
		selected_models$n_predictors,
		selected_models$RMSE,
		pch = 16
	)

	lines(
		selected_models$n_predictors,
		selected_models$RMSE,
		lwd = 2
	)
}

pdf(
	path_fig_predictor_screening_rmse_by_n_pdf,
	width = 7,
	height = 5
)

make_predictor_screening_rmse_plot()

dev.off()

png(
	path_fig_predictor_screening_rmse_by_n_png,
	width = 1800,
	height = 1200,
	res = 250
)

make_predictor_screening_rmse_plot()

dev.off()

############################################################
# Console summary
############################################################

cat("Initial predictor screening complete.\n")
cat("Total observations:", n_observations, "\n")
cat("Total models screened:", nrow(screening_results), "\n")
cat("Best model formula:\n")
cat(screening_results$formula[1], "\n")
cat("Best model apparent RMSE:", screening_results$RMSE[1], "\n")
cat("Best model apparent MAE:", screening_results$MAE[1], "\n")
cat("Best model apparent R:", screening_results$R[1], "\n")
cat("Best model AIC:", screening_results$AIC[1], "\n")