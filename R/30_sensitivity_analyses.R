# ============================================================
# 30_sensitivity_analyses.R
#
# Purpose:
# 	Test robustness of the extended physicochemical
# 	hierarchical QSPR analysis.
#
# Primary molecular specification tested here:
#
# 	S2:
# 		MW + logKow + melting point +
# 		logSaq + HBD + HBA
#
# Sensitivity analyses:
#
# 	S1. Multi-reference compounds only
#
# 	S2. Exclude observations with
# 		|standardized residual| >3
#
# 	S3. Mean rather than median descriptor aggregation
#
# 	S4. Expanded nonlinear molecular specification:
#
# 		MW + logKow + melting point +
# 		logSaq + logSaq^2 +
# 		HBA + logKow:HBA
#
# For each sensitivity, compare:
#
# 	Conventional compound-mean QSPR
# 	vs
# 	Joint hierarchical QSPR
#
# using leave-one-compound-out cross-validation.
#
# This script DOES NOT:
# 	- select predictors
# 	- optimize formulas
# 	- redefine compound identity
# 	- replace the primary S1/S2 analysis
#
# Run from project root:
# 	source("R/30_sensitivity_analyses.R")
# ============================================================


source("R/00_config.R")


# ============================================================
# 1. Import data
# ============================================================

dat <- readr::read_csv(
	"R/results/10_model_data.csv",
	show_col_types = FALSE
)


diagnostics <- readr::read_csv(
	"R/results/11_observation_diagnostics.csv",
	show_col_types = FALSE
)


primary_performance <- readr::read_csv(
	"R/results/20_loco_performance.csv",
	show_col_types = FALSE
)


primary_variance <- readr::read_csv(
	"R/results/20_compound_variance_reduction.csv",
	show_col_types = FALSE
)


# ============================================================
# 2. Define molecular specifications
# ============================================================

all_descriptor_variables <- c(
	"mw",
	"log_kow",
	"melting_point_K",
	"log_saq",
	"h_donor",
	"h_acceptor"
)


primary_descriptors <- c(
	"mw",
	"log_kow",
	"melting_point_K",
	"log_saq",
	"h_donor",
	"h_acceptor"
)


nonlinear_descriptors <- c(
	"mw",
	"log_kow",
	"melting_point_K",
	"log_saq",
	"h_acceptor"
)


# ============================================================
# 3. Model control
# ============================================================

model_control <- lme4::lmerControl(
	optimizer = "bobyqa",
	optCtrl = list(
		maxfun = 100000
	)
)


# ============================================================
# 4. Performance functions
# ============================================================

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
		abs(
			observed - predicted
		),
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
			(
				observed -
					mean(
						observed,
						na.rm = TRUE
					)
			)^2,
			na.rm = TRUE
		)
}


# ============================================================
# 5. Summarize descriptors to one row per compound
# ============================================================

summarize_descriptors <- function(
	data,
	descriptors,
	method = "median"
) {

	if (!method %in% c("median", "mean")) {
		stop(
			"Descriptor method must be 'median' or 'mean'."
		)
	}


	if (method == "median") {

		out <- data %>%
			group_by(
				compound_id
			) %>%
			summarise(
				compound_name =
					first(
						compound_name
					),

				canonical_name =
					first(
						canonical_name
					),

				across(
					all_of(
						descriptors
					),
					~ if (
						all(is.na(.))
					) {
						NA_real_
					} else {
						median(
							.,
							na.rm = TRUE
						)
					}
				),

				.groups = "drop"
			)

	} else {

		out <- data %>%
			group_by(
				compound_id
			) %>%
			summarise(
				compound_name =
					first(
						compound_name
					),

				canonical_name =
					first(
						canonical_name
					),

				across(
					all_of(
						descriptors
					),
					~ if (
						all(is.na(.))
					) {
						NA_real_
					} else {
						mean(
							.,
							na.rm = TRUE
						)
					}
				),

				.groups = "drop"
			)
	}


	out
}


# ============================================================
# 6. Construct compound endpoint and row-level QSPR data
# ============================================================

prepare_qspr_data <- function(
	data,
	descriptors,
	descriptor_method = "median"
) {

	compound_descriptors <- summarize_descriptors(
		data,
		descriptors,
		method = descriptor_method
	)


	compound_endpoint <- data %>%
		group_by(
			compound_id
		) %>%
		summarise(
			compound_name =
				first(
					compound_name
				),

			canonical_name =
				first(
					canonical_name
				),

			n_observations =
				n(),

			n_references =
				n_distinct(
					reference_id
				),

			raw_mean_logKp =
				mean(
					logKp,
					na.rm = TRUE
				),

			.groups = "drop"
		) %>%
		left_join(
			compound_descriptors %>%
				select(
					compound_id,
					all_of(
						descriptors
					)
				),
			by = "compound_id"
		) %>%
		filter(
			if_all(
				all_of(
					descriptors
				),
				~ !is.na(.)
			)
		)


	eligible_ids <-
		compound_endpoint$compound_id


	row_data <- data %>%
		filter(
			compound_id %in%
				eligible_ids
		) %>%
		select(
			-all_of(
				all_descriptor_variables
			)
		) %>%
		left_join(
			compound_descriptors %>%
				select(
					compound_id,
					all_of(
						descriptors
					)
				),
			by = "compound_id"
		)


	if (
		anyDuplicated(
			compound_endpoint$compound_id
		) > 0
	) {
		stop(
			"compound_endpoint contains duplicate compound_id values."
		)
	}


	list(
		compound =
			compound_endpoint,

		rows =
			row_data
	)
}


# ============================================================
# 7. Scale descriptors
#
# Scaling is calculated across compounds rather than rows so
# highly replicated compounds do not receive extra weight.
# ============================================================

scale_qspr_data <- function(
	compound_data,
	row_data,
	descriptors
) {

	scaling <- compound_data %>%
		summarise(
			across(
				all_of(
					descriptors
				),
				list(
					mean =
						~ mean(
							.,
							na.rm = TRUE
						),

					sd =
						~ sd(
							.,
							na.rm = TRUE
						)
				)
			)
		)


	for (variable in descriptors) {

		variable_mean <-
			scaling[[paste0(variable, "_mean")]]


		variable_sd <-
			scaling[[paste0(variable, "_sd")]]


		if (
			is.na(variable_sd) ||
			variable_sd == 0
		) {
			stop(
				"Cannot scale descriptor ",
				variable,
				"."
			)
		}


		z_name <-
			paste0(
				variable,
				"_z"
			)


		compound_data[[z_name]] <-
			(
				compound_data[[variable]] -
					variable_mean
			) /
				variable_sd


		row_data[[z_name]] <-
			(
				row_data[[variable]] -
					variable_mean
			) /
				variable_sd
	}


	list(
		compound =
			compound_data,

		rows =
			row_data
	)
}


# ============================================================
# 8. Construct model formulas
# ============================================================

make_formula <- function(
	descriptors,
	nonlinear = FALSE
) {

	z_terms <- paste0(
		descriptors,
		"_z"
	)


	if (!nonlinear) {

		molecular_terms <- paste(
			z_terms,
			collapse = " + "
		)

	} else {

		molecular_terms <- paste(
			"mw_z",
			"log_kow_z",
			"melting_point_K_z",
			"log_saq_z",
			"I(log_saq_z^2)",
			"h_acceptor_z",
			"log_kow_z:h_acceptor_z",
			sep = " + "
		)
	}


	raw_formula <- as.formula(
		paste(
			"raw_mean_logKp ~",
			molecular_terms
		)
	)


	joint_formula <- as.formula(
		paste(
			"logKp ~",
			molecular_terms,
			"+ temperature_c_qspr",
			"+ (1 | compound_id)",
			"+ (1 | reference_id)"
		)
	)


	list(
		raw =
			raw_formula,

		joint =
			joint_formula
	)
}


# ============================================================
# 9. Fit full-data sensitivity models
# ============================================================

fit_sensitivity <- function(
	data,
	descriptors,
	descriptor_method,
	nonlinear,
	scenario
) {

	prepared <- prepare_qspr_data(
		data,
		descriptors,
		descriptor_method
	)


	scaled <- scale_qspr_data(
		prepared$compound,
		prepared$rows,
		descriptors
	)


	compound_data <-
		scaled$compound


	row_data <-
		scaled$rows


	temperature_reference <- mean(
		row_data$temperature_K,
		na.rm = TRUE
	)


	row_data <- row_data %>%
		mutate(
			temperature_c_qspr =
				temperature_K -
					temperature_reference,

			compound_id =
				factor(
					compound_id
				),

			reference_id =
				factor(
					reference_id
				)
		)


	formulas <- make_formula(
		descriptors,
		nonlinear =
			nonlinear
	)


	raw_model <- lm(
		formulas$raw,
		data =
			compound_data
	)


	baseline_model <- lmerTest::lmer(
		logKp ~
			temperature_c_qspr +
			(1 | compound_id) +
			(1 | reference_id),
		data =
			row_data,
		REML =
			TRUE,
		control =
			model_control
	)


	joint_model <- lmerTest::lmer(
		formulas$joint,
		data =
			row_data,
		REML =
			TRUE,
		control =
			model_control
	)


	get_compound_variance <- function(model) {

		as.data.frame(
			VarCorr(model)
		) %>%
			filter(
				grp ==
					"compound_id"
			) %>%
			pull(
				vcov
			)
	}


	baseline_variance <-
		get_compound_variance(
			baseline_model
		)


	joint_variance <-
		get_compound_variance(
			joint_model
		)


	summary <- tibble(
		scenario =
			scenario,

		n_observations =
			nrow(
				row_data
			),

		n_compounds =
			n_distinct(
				row_data$compound_id
			),

		n_references =
			n_distinct(
				row_data$reference_id
			),

		baseline_compound_variance =
			baseline_variance,

		joint_compound_variance =
			joint_variance,

		percent_compound_variance_explained =
			100 *
				(
					1 -
						joint_variance /
							baseline_variance
				)
	)


	list(
		summary =
			summary,

		raw_model =
			raw_model,

		baseline_model =
			baseline_model,

		joint_model =
			joint_model,

		compound_data =
			compound_data,

		row_data =
			row_data
	)
}


# ============================================================
# 10. LOCO cross-validation
# ============================================================

run_loco <- function(
	data,
	descriptors,
	descriptor_method,
	nonlinear,
	scenario
) {

	prepared <- prepare_qspr_data(
		data,
		descriptors,
		descriptor_method
	)


	compound_data <-
		prepared$compound


	compound_ids <-
		compound_data$compound_id


	formulas <- make_formula(
		descriptors,
		nonlinear =
			nonlinear
	)


	predictions <- list()

	failures <- list()

	prediction_counter <- 1

	failure_counter <- 1


	for (
		i in seq_along(
			compound_ids
		)
	) {

		held_out <-
			compound_ids[i]


		message(
			scenario,
			" LOCO ",
			i,
			"/",
			length(
				compound_ids
			)
		)


		train_compounds <- compound_data %>%
			filter(
				compound_id !=
					held_out
			)


		test_compound <- compound_data %>%
			filter(
				compound_id ==
					held_out
			)


		train_rows <- data %>%
			filter(
				compound_id %in%
					train_compounds$compound_id
			)


		test_rows <- data %>%
			filter(
				compound_id ==
					held_out
			)


		# ----------------------------------------------------
		# Training-fold descriptor scaling
		# ----------------------------------------------------

		fold_means <- train_compounds %>%
			summarise(
				across(
					all_of(
						descriptors
					),
					~ mean(
						.,
						na.rm = TRUE
					)
				)
			)


		fold_sds <- train_compounds %>%
			summarise(
				across(
					all_of(
						descriptors
					),
					~ sd(
						.,
						na.rm = TRUE
					)
				)
			)


		for (variable in descriptors) {

			z_name <-
				paste0(
					variable,
					"_z"
				)


			train_compounds[[z_name]] <-
				(
					train_compounds[[variable]] -
						fold_means[[variable]]
				) /
					fold_sds[[variable]]


			test_compound[[z_name]] <-
				(
					test_compound[[variable]] -
						fold_means[[variable]]
				) /
					fold_sds[[variable]]
		}


		# ----------------------------------------------------
		# Conventional compound-mean model
		# ----------------------------------------------------

		raw_model <- tryCatch(
			lm(
				formulas$raw,
				data =
					train_compounds
			),
			error = function(e) {
				NULL
			}
		)


		if (is.null(raw_model)) {

			failures[[failure_counter]] <- tibble(
				scenario =
					scenario,

				compound_id =
					held_out,

				approach =
					"Compound-mean QSPR"
			)

			failure_counter <-
				failure_counter + 1

			next
		}


		raw_prediction <- tryCatch(
			as.numeric(
				predict(
					raw_model,
					newdata =
						test_compound
				)
			),
			error = function(e) {
				NA_real_
			}
		)


		# ----------------------------------------------------
		# Attach scaled descriptors to row-level data
		# ----------------------------------------------------

		train_descriptor_data <- train_compounds %>%
			select(
				compound_id,
				all_of(
					paste0(
						descriptors,
						"_z"
					)
				)
			)


		test_descriptor_data <- test_compound %>%
			select(
				compound_id,
				all_of(
					paste0(
						descriptors,
						"_z"
					)
				)
			)


		train_rows <- train_rows %>%
			select(
				-all_of(
					all_descriptor_variables
				)
			) %>%
			left_join(
				train_descriptor_data,
				by =
					"compound_id"
			)


		test_rows <- test_rows %>%
			select(
				-all_of(
					all_descriptor_variables
				)
			) %>%
			left_join(
				test_descriptor_data,
				by =
					"compound_id"
			)


		fold_temperature_reference <- mean(
			train_rows$temperature_K,
			na.rm = TRUE
		)


		train_rows <- train_rows %>%
			mutate(
				temperature_c_qspr =
					temperature_K -
						fold_temperature_reference,

				compound_id =
					factor(
						compound_id
					),

				reference_id =
					factor(
						reference_id
					)
			)


		test_rows <- test_rows %>%
			mutate(
				temperature_c_qspr =
					temperature_K -
						fold_temperature_reference
			)


		# ----------------------------------------------------
		# Joint hierarchical QSPR
		# ----------------------------------------------------

		joint_model <- tryCatch(
			lme4::lmer(
				formulas$joint,
				data =
					train_rows,
				REML =
					TRUE,
				control =
					model_control
			),
			error = function(e) {
				NULL
			}
		)


		if (is.null(joint_model)) {

			failures[[failure_counter]] <- tibble(
				scenario =
					scenario,

				compound_id =
					held_out,

				approach =
					"Joint hierarchical QSPR"
			)

			failure_counter <-
				failure_counter + 1

			next
		}


		joint_prediction <- tryCatch(
			as.numeric(
				predict(
					joint_model,
					newdata =
						test_rows,
					re.form =
						~ (1 | reference_id),
					allow.new.levels =
						TRUE
				)
			),
			error = function(e) {
				rep(
					NA_real_,
					nrow(
						test_rows
					)
				)
			}
		)


		predictions[[prediction_counter]] <- tibble(
			scenario =
				scenario,

			compound_id =
				held_out,

			source_row =
				test_rows$source_row,

			observed_logKp =
				test_rows$logKp,

			compound_mean_prediction =
				raw_prediction,

			joint_hierarchical_prediction =
				joint_prediction
		)


		prediction_counter <-
			prediction_counter + 1
	}


	predictions <- bind_rows(
		predictions
	)


	failures <- bind_rows(
		failures
	)


	performance <- bind_rows(

		tibble(
			scenario =
				scenario,

			approach =
				"Compound-mean QSPR",

			n =
				sum(
					complete.cases(
						predictions$observed_logKp,
						predictions$compound_mean_prediction
					)
				),

			RMSE =
				rmse(
					predictions$observed_logKp,
					predictions$compound_mean_prediction
				),

			MAE =
				mae(
					predictions$observed_logKp,
					predictions$compound_mean_prediction
				),

			R2_pred =
				r2_pred(
					predictions$observed_logKp,
					predictions$compound_mean_prediction
				),

			Pearson_r =
				cor(
					predictions$observed_logKp,
					predictions$compound_mean_prediction,
					use = "complete.obs"
				)
		),

		tibble(
			scenario =
				scenario,

			approach =
				"Joint hierarchical QSPR",

			n =
				sum(
					complete.cases(
						predictions$observed_logKp,
						predictions$joint_hierarchical_prediction
					)
				),

			RMSE =
				rmse(
					predictions$observed_logKp,
					predictions$joint_hierarchical_prediction
				),

			MAE =
				mae(
					predictions$observed_logKp,
					predictions$joint_hierarchical_prediction
				),

			R2_pred =
				r2_pred(
					predictions$observed_logKp,
					predictions$joint_hierarchical_prediction
				),

			Pearson_r =
				cor(
					predictions$observed_logKp,
					predictions$joint_hierarchical_prediction,
					use = "complete.obs"
				)
		)
	)


	compound_errors <- predictions %>%
		group_by(
			compound_id
		) %>%
		summarise(
			compound_mean_RMSE =
				rmse(
					observed_logKp,
					compound_mean_prediction
				),

			compound_mean_MAE =
				mae(
					observed_logKp,
					compound_mean_prediction
				),

			joint_RMSE =
				rmse(
					observed_logKp,
					joint_hierarchical_prediction
				),

			joint_MAE =
				mae(
					observed_logKp,
					joint_hierarchical_prediction
				),

			.groups = "drop"
		)


	compound_balanced <- tibble(
		scenario =
			scenario,

		compound_mean_mean_RMSE =
			mean(
				compound_errors$compound_mean_RMSE,
				na.rm = TRUE
			),

		joint_mean_RMSE =
			mean(
				compound_errors$joint_RMSE,
				na.rm = TRUE
			),

		compound_mean_mean_MAE =
			mean(
				compound_errors$compound_mean_MAE,
				na.rm = TRUE
			),

		joint_mean_MAE =
			mean(
				compound_errors$joint_MAE,
				na.rm = TRUE
			)
	) %>%
		mutate(
			RMSE_improvement_percent =
				100 *
					(
						compound_mean_mean_RMSE -
							joint_mean_RMSE
					) /
						compound_mean_mean_RMSE,

			MAE_improvement_percent =
				100 *
					(
						compound_mean_mean_MAE -
							joint_mean_MAE
					) /
						compound_mean_mean_MAE
		)


	list(
		predictions =
			predictions,

		failures =
			failures,

		performance =
			performance,

		compound_balanced =
			compound_balanced
	)
}


# ============================================================
# 11. Sensitivity 1: multi-reference compounds only
# ============================================================

multi_reference_ids <- dat %>%
	group_by(
		compound_id
	) %>%
	summarise(
		n_references =
			n_distinct(
				reference_id
			),
		.groups = "drop"
	) %>%
	filter(
		n_references >= 2
	) %>%
	pull(
		compound_id
	)


dat_multireference <- dat %>%
	filter(
		compound_id %in%
			multi_reference_ids
	)


s1_fit <- fit_sensitivity(
	dat_multireference,
	descriptors =
		primary_descriptors,
	descriptor_method =
		"median",
	nonlinear =
		FALSE,
	scenario =
		"Multi-reference compounds"
)


s1_cv <- run_loco(
	dat_multireference,
	descriptors =
		primary_descriptors,
	descriptor_method =
		"median",
	nonlinear =
		FALSE,
	scenario =
		"Multi-reference compounds"
)


# ============================================================
# 12. Sensitivity 2: exclude extreme residuals
# ============================================================

extreme_rows <- diagnostics %>%
	filter(
		absolute_standardized_residual > 3
	) %>%
	pull(
		source_row
	)


dat_no_extreme <- dat %>%
	filter(
		!source_row %in%
			extreme_rows
	)


s2_fit <- fit_sensitivity(
	dat_no_extreme,
	descriptors =
		primary_descriptors,
	descriptor_method =
		"median",
	nonlinear =
		FALSE,
	scenario =
		"Exclude |standardized residual| >3"
)


s2_cv <- run_loco(
	dat_no_extreme,
	descriptors =
		primary_descriptors,
	descriptor_method =
		"median",
	nonlinear =
		FALSE,
	scenario =
		"Exclude |standardized residual| >3"
)


# ============================================================
# 13. Sensitivity 3: mean descriptor aggregation
# ============================================================

s3_fit <- fit_sensitivity(
	dat,
	descriptors =
		primary_descriptors,
	descriptor_method =
		"mean",
	nonlinear =
		FALSE,
	scenario =
		"Mean descriptor values"
)


s3_cv <- run_loco(
	dat,
	descriptors =
		primary_descriptors,
	descriptor_method =
		"mean",
	nonlinear =
		FALSE,
	scenario =
		"Mean descriptor values"
)


# ============================================================
# 14. Sensitivity 4: expanded nonlinear specification
#
# This is deliberately treated as a sensitivity analysis
# rather than the primary molecular model.
# ============================================================

s4_fit <- fit_sensitivity(
	dat,
	descriptors =
		nonlinear_descriptors,
	descriptor_method =
		"median",
	nonlinear =
		TRUE,
	scenario =
		"Expanded nonlinear specification"
)


s4_cv <- run_loco(
	dat,
	descriptors =
		nonlinear_descriptors,
	descriptor_method =
		"median",
	nonlinear =
		TRUE,
	scenario =
		"Expanded nonlinear specification"
)


# ============================================================
# 15. Combine model summaries
# ============================================================

sensitivity_model_summary <- bind_rows(
	s1_fit$summary,
	s2_fit$summary,
	s3_fit$summary,
	s4_fit$summary
)


readr::write_csv(
	sensitivity_model_summary,
	"R/results/30_sensitivity_model_summary.csv"
)


# ============================================================
# 16. Combine LOCO results
# ============================================================

sensitivity_loco_performance <- bind_rows(
	s1_cv$performance,
	s2_cv$performance,
	s3_cv$performance,
	s4_cv$performance
)


readr::write_csv(
	sensitivity_loco_performance,
	"R/results/30_sensitivity_loco_performance.csv"
)


sensitivity_compound_balanced <- bind_rows(
	s1_cv$compound_balanced,
	s2_cv$compound_balanced,
	s3_cv$compound_balanced,
	s4_cv$compound_balanced
)


readr::write_csv(
	sensitivity_compound_balanced,
	"R/results/30_sensitivity_compound_balanced.csv"
)


# ============================================================
# 17. Combine failures
# ============================================================

sensitivity_failures <- bind_rows(
	s1_cv$failures,
	s2_cv$failures,
	s3_cv$failures,
	s4_cv$failures
)


readr::write_csv(
	sensitivity_failures,
	"R/results/30_sensitivity_failures.csv"
)


# ============================================================
# 18. Combine sensitivity predictions
# ============================================================

sensitivity_predictions <- bind_rows(
	s1_cv$predictions,
	s2_cv$predictions,
	s3_cv$predictions,
	s4_cv$predictions
)


readr::write_csv(
	sensitivity_predictions,
	"R/results/30_sensitivity_predictions.csv"
)


# ============================================================
# 19. Extract primary S2 results for comparison
# ============================================================

primary_s2_performance <- primary_performance %>%
	filter(
		specification ==
			"S2"
	)


primary_raw_rmse <- primary_s2_performance %>%
	filter(
		approach ==
			"Compound-mean QSPR"
	) %>%
	pull(
		RMSE
	)


primary_joint_rmse <- primary_s2_performance %>%
	filter(
		approach ==
			"Joint hierarchical QSPR"
	) %>%
	pull(
		RMSE
	)


primary_s2_variance <- primary_variance %>%
	filter(
		specification ==
			"S2"
	) %>%
	pull(
		percent_compound_variance_explained
	)


primary_summary <- tibble(
	scenario =
		"Primary S2 analysis",

	percent_compound_variance_explained =
		primary_s2_variance,

	raw_RMSE =
		primary_raw_rmse,

	joint_RMSE =
		primary_joint_rmse,

	RMSE_improvement_percent =
		100 *
			(
				primary_raw_rmse -
					primary_joint_rmse
			) /
				primary_raw_rmse
)


# ============================================================
# 20. Create compact overall sensitivity summary
# ============================================================

sensitivity_performance_wide <- sensitivity_loco_performance %>%
	select(
		scenario,
		approach,
		RMSE
	) %>%
	pivot_wider(
		names_from =
			approach,
		values_from =
			RMSE
	) %>%
	rename(
		raw_RMSE =
			`Compound-mean QSPR`,

		joint_RMSE =
			`Joint hierarchical QSPR`
	)


sensitivity_summary <- sensitivity_model_summary %>%
	select(
		scenario,
		percent_compound_variance_explained
	) %>%
	left_join(
		sensitivity_performance_wide,
		by =
			"scenario"
	) %>%
	mutate(
		RMSE_improvement_percent =
			100 *
				(
					raw_RMSE -
						joint_RMSE
				) /
					raw_RMSE
	)


overall_summary <- bind_rows(
	primary_summary,
	sensitivity_summary
)


readr::write_csv(
	overall_summary,
	"R/results/30_sensitivity_summary.csv"
)


# ============================================================
# 21. Save sensitivity models
# ============================================================

saveRDS(
	list(
		multireference =
			s1_fit,

		no_extreme =
			s2_fit,

		mean_descriptors =
			s3_fit,

		nonlinear =
			s4_fit
	),
	"R/results/30_sensitivity_models.rds"
)


# ============================================================
# 22. Console summary
# ============================================================

cat("\n")
cat("========================================\n")
cat("30_sensitivity_analyses.R complete\n")
cat("========================================\n\n")

cat("Overall sensitivity summary:\n\n")

print(
	overall_summary
)

cat("\nCompound-balanced performance:\n\n")

print(
	sensitivity_compound_balanced
)

cat("\nLOCO failures: ")

if (nrow(sensitivity_failures) == 0) {
	cat("none\n")
} else {
	cat(
		nrow(
			sensitivity_failures
		),
		"\n"
	)
}

cat("\nOutputs written to:\n")

cat(
	"  R/results/30_sensitivity_summary.csv\n"
)

cat(
	"  R/results/30_sensitivity_model_summary.csv\n"
)

cat(
	"  R/results/30_sensitivity_loco_performance.csv\n"
)

cat(
	"  R/results/30_sensitivity_compound_balanced.csv\n"
)

cat(
	"  R/results/30_sensitivity_failures.csv\n"
)

cat(
	"  R/results/30_sensitivity_predictions.csv\n"
)

cat(
	"  R/results/30_sensitivity_models.rds\n"
)