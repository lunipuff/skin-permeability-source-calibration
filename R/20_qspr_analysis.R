# ============================================================
# 20_qspr_analysis.R
#
# Purpose:
# 	Test whether explicitly modeling the measurement hierarchy
# 	improves QSPR analysis relative to conventional
# 	compound-level aggregation.
#
# Molecular specifications are fixed a priori.
# No predictor selection is performed.
#
# Specification S1:
# 	Size-lipophilicity model
#
# 		MW + logKow
#
# Specification S2:
# 	Extended physicochemical model
#
# 		MW + logKow + melting point +
# 		log aqueous solubility + HBD + HBA
#
# For each molecular specification, two approaches are fitted:
#
# 	Conventional QSPR:
#
# 		mean(logKp) ~ molecular descriptors
#
# 	Joint hierarchical QSPR:
#
# 		logKp ~ molecular descriptors +
# 			temperature +
# 			(1 | compound) +
# 			(1 | reference)
#
# Main questions:
#
# 	1. How much between-compound variance is explained by
# 	   each fixed molecular specification?
#
# 	2. Does hierarchical modeling improve prediction relative
# 	   to compound-mean QSPR when molecular information is
# 	   held constant?
#
# 	3. Does the hierarchical advantage persist across a
# 	   simple and a broader molecular specification?
#
# This script DOES NOT:
# 	- select predictors
# 	- drop predictors based on p-values
# 	- optimize molecular formulas
# 	- add nonlinear terms or interactions
# 	- remove outliers
#
# Run from the project root:
# 	source("R/20_qspr_analysis.R")
# ============================================================


source("R/00_config.R")


# ============================================================
# 1. Import prepared model data
# ============================================================

dat <- readr::read_csv(
	"R/results/10_model_data.csv",
	show_col_types = FALSE
)


# ============================================================
# 2. Define fixed molecular specifications
# ============================================================

descriptor_variables <- c(
	"mw",
	"log_kow",
	"melting_point_K",
	"log_saq",
	"h_donor",
	"h_acceptor"
)


specifications <- list(

	S1 = c(
		"mw",
		"log_kow"
	),

	S2 = c(
		"mw",
		"log_kow",
		"melting_point_K",
		"log_saq",
		"h_donor",
		"h_acceptor"
	)
)


specification_labels <- c(
	S1 = "Size-lipophilicity",
	S2 = "Extended physicochemical"
)


required_variables <- c(
	"compound_id",
	"compound_name",
	"canonical_name",
	"reference_id",
	"logKp",
	"temperature_K",
	descriptor_variables
)


missing_variables <- setdiff(
	required_variables,
	names(dat)
)


if (length(missing_variables) > 0) {
	stop(
		"Required variables missing: ",
		paste(
			missing_variables,
			collapse = ", "
		)
	)
}


# ============================================================
# 3. Audit descriptor consistency within canonical compounds
# ============================================================

descriptor_consistency <- dat %>%
	group_by(
		compound_id,
		compound_name,
		canonical_name
	) %>%
	summarise(
		across(
			all_of(
				descriptor_variables
			),
			~ n_distinct(
				.,
				na.rm = TRUE
			),
			.names = "n_{.col}"
		),
		.groups = "drop"
	) %>%
	mutate(
		descriptor_inconsistency =
			if_any(
				starts_with("n_"),
				~ . > 1
			)
	)


readr::write_csv(
	descriptor_consistency,
	"R/results/20_descriptor_consistency.csv"
)


# ============================================================
# 4. Export actual inconsistent descriptor values
# ============================================================

descriptor_value_audit <- dat %>%
	select(
		compound_id,
		compound_name,
		canonical_name,
		all_of(
			descriptor_variables
		)
	) %>%
	pivot_longer(
		cols =
			all_of(
				descriptor_variables
			),
		names_to =
			"descriptor",
		values_to =
			"value"
	) %>%
	group_by(
		compound_id,
		compound_name,
		canonical_name,
		descriptor
	) %>%
	summarise(
		n_values =
			n_distinct(
				value,
				na.rm = TRUE
			),

		values =
			paste(
				sort(
					unique(
						value[!is.na(value)]
					)
				),
				collapse = " | "
			),

		.groups = "drop"
	) %>%
	filter(
		n_values > 1
	) %>%
	arrange(
		compound_name,
		descriptor
	)


readr::write_csv(
	descriptor_value_audit,
	"R/results/20_descriptor_value_audit.csv"
)


# ============================================================
# 5. Create one molecular profile per canonical compound
#
# Median values are used when source rows contain slightly
# different physicochemical values for the same compound.
#
# Descriptor-mean sensitivity is handled later in script 30.
# ============================================================

compound_descriptors <- dat %>%
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
				descriptor_variables
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


if (
	anyDuplicated(
		compound_descriptors$compound_id
	) > 0
) {
	stop(
		"compound_descriptors contains duplicate compound_id values."
	)
}


# ============================================================
# 6. Construct conventional compound-level endpoint
# ============================================================

compound_endpoint <- dat %>%
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

		raw_sd_logKp =
			if_else(
				n() >= 2,
				sd(
					logKp,
					na.rm = TRUE
				),
				NA_real_
			),

		.groups = "drop"
	) %>%
	left_join(
		compound_descriptors %>%
			select(
				compound_id,
				all_of(
					descriptor_variables
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


# ============================================================
# 7. QSPR sample summary by specification
# ============================================================

qspr_sample_summary <- purrr::imap_dfr(
	specifications,
	function(spec_descriptors, spec_id) {

		eligible_compounds <- compound_endpoint %>%
			filter(
				if_all(
					all_of(
						spec_descriptors
					),
					~ !is.na(.)
				)
			)


		eligible_ids <-
			eligible_compounds$compound_id


		eligible_rows <- dat %>%
			filter(
				compound_id %in%
					eligible_ids
			)


		tibble(
			specification =
				spec_id,

			specification_label =
				specification_labels[[spec_id]],

			n_descriptors =
				length(
					spec_descriptors
				),

			descriptors =
				paste(
					spec_descriptors,
					collapse = " + "
				),

			n_compounds =
				nrow(
					eligible_compounds
				),

			n_observations =
				nrow(
					eligible_rows
				),

			n_references =
				n_distinct(
					eligible_rows$reference_id
				),

			n_compounds_excluded =
				nrow(
					compound_endpoint
				) -
					nrow(
						eligible_compounds
					)
		)
	}
)


readr::write_csv(
	qspr_sample_summary,
	"R/results/20_qspr_sample_summary.csv"
)


# ============================================================
# 8. Utility functions
# ============================================================

rmse <- function(observed, predicted) {
	sqrt(
		mean(
			(
				observed -
					predicted
			)^2,
			na.rm = TRUE
		)
	)
}


mae <- function(observed, predicted) {
	mean(
		abs(
			observed -
				predicted
		),
		na.rm = TRUE
	)
}


r2_pred <- function(observed, predicted) {
	1 -
		sum(
			(
				observed -
					predicted
			)^2,
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


model_control <- lme4::lmerControl(
	optimizer = "bobyqa",
	optCtrl = list(
		maxfun = 100000
	)
)


# ============================================================
# 9. Helper: add standardized descriptors
#
# Scaling parameters are estimated at the compound level so
# replicated compounds do not receive additional weight.
# ============================================================

scale_qspr_data <- function(
	compound_data,
	row_data,
	spec_descriptors
) {

	scaling <- compound_data %>%
		summarise(
			across(
				all_of(
					spec_descriptors
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


	for (
		variable in spec_descriptors
	) {

		variable_mean <-
			scaling[[paste0(variable,"_mean")]]


		variable_sd <-
			scaling[[paste0(variable,"_sd")]]


		if (
			is.na(variable_sd) ||
			variable_sd == 0
		) {
			stop(
				"Cannot scale descriptor ",
				variable,
				": zero or missing SD."
			)
		}


		z_name <-
			paste0(
				variable,
				"_z"
			)


		compound_data[[z_name]] <- (compound_data[[variable]] - variable_mean) / variable_sd


		row_data[[z_name]] <- (row_data[[variable]] - variable_mean) / variable_sd
	}


	list(
		compound =
			compound_data,

		rows =
			row_data
	)
}


# ============================================================
# 10. Helper: construct formulas
# ============================================================

make_formulas <- function(
	spec_descriptors
) {

	z_terms <- paste0(
		spec_descriptors,
		"_z"
	)


	molecular_terms <- paste(
		z_terms,
		collapse = " + "
	)


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
# 11. Fit full-data models for each specification
# ============================================================

full_model_results <- list()


for (
	spec_id in names(
		specifications
	)
) {

	message(
		"Fitting full-data specification ",
		spec_id,
		": ",
		specification_labels[[spec_id]]
	)


	spec_descriptors <-
		specifications[[spec_id]]


	eligible_compounds <- compound_endpoint %>%
		filter(
			if_all(
				all_of(
					spec_descriptors
				),
				~ !is.na(.)
			)
		)


	eligible_ids <-
		eligible_compounds$compound_id


	row_data <- dat %>%
		filter(
			compound_id %in%
				eligible_ids
		) %>%
		select(
			-all_of(
				descriptor_variables
			)
		) %>%
		left_join(
			compound_descriptors %>%
				select(
					compound_id,
					all_of(
						spec_descriptors
					)
				),
			by = "compound_id"
		)


	scaled <- scale_qspr_data(
		eligible_compounds,
		row_data,
		spec_descriptors
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


	formulas <- make_formulas(
		spec_descriptors
	)


	raw_model <- lm(
		formulas$raw,
		data = compound_data
	)


	baseline_model <- lmerTest::lmer(
		logKp ~
			temperature_c_qspr +
			(1 | compound_id) +
			(1 | reference_id),
		data = row_data,
		REML = TRUE,
		control = model_control
	)


	joint_model <- lmerTest::lmer(
		formulas$joint,
		data = row_data,
		REML = TRUE,
		control = model_control
	)


	baseline_ml <- update(
		baseline_model,
		REML = FALSE
	)


	joint_ml <- update(
		joint_model,
		REML = FALSE
	)


	full_model_results[[spec_id]] <- list(
		specification =
			spec_id,

		specification_label =
			specification_labels[[spec_id]],

		descriptors =
			spec_descriptors,

		temperature_reference =
			temperature_reference,

		compound_data =
			compound_data,

		row_data =
			row_data,

		raw_model =
			raw_model,

		baseline_model =
			baseline_model,

		joint_model =
			joint_model,

		baseline_ml =
			baseline_ml,

		joint_ml =
			joint_ml
	)
}


# ============================================================
# 12. Full-model likelihood comparisons
#
# The comparison asks whether adding the fixed molecular
# descriptor panel improves the hierarchical model.
#
# Likelihood-ratio statistics are calculated directly from
# the ML log-likelihoods to avoid dependence on package-
# specific anova() column names.
# ============================================================

model_comparison <- purrr::imap_dfr(
	full_model_results,
	function(x, spec_id) {

		baseline_loglik <-
			as.numeric(
				logLik(
					x$baseline_ml
				)
			)


		joint_loglik <-
			as.numeric(
				logLik(
					x$joint_ml
				)
			)


		baseline_df <-
			attr(
				logLik(
					x$baseline_ml
				),
				"df"
			)


		joint_df <-
			attr(
				logLik(
					x$joint_ml
				),
				"df"
			)


		lrt_chisq <-
			2 *
				(
					joint_loglik -
						baseline_loglik
				)


		lrt_df <-
			joint_df -
				baseline_df


		lrt_p <-
			pchisq(
				lrt_chisq,
				df = lrt_df,
				lower.tail = FALSE
			)


		tibble(
			specification =
				spec_id,

			specification_label =
				x$specification_label,

			n_compounds =
				nrow(
					x$compound_data
				),

			n_observations =
				nrow(
					x$row_data
				),

			baseline_logLik =
				baseline_loglik,

			joint_logLik =
				joint_loglik,

			baseline_AIC =
				AIC(
					x$baseline_ml
				),

			joint_AIC =
				AIC(
					x$joint_ml
				),

			delta_AIC =
				AIC(
					x$joint_ml
				) -
					AIC(
						x$baseline_ml
					),

			baseline_BIC =
				BIC(
					x$baseline_ml
				),

			joint_BIC =
				BIC(
					x$joint_ml
				),

			delta_BIC =
				BIC(
					x$joint_ml
				) -
					BIC(
						x$baseline_ml
					),

			LRT_chisq =
				lrt_chisq,

			LRT_df =
				lrt_df,

			LRT_p =
				lrt_p
		)
	}
)


readr::write_csv(
	model_comparison,
	"R/results/20_model_comparison.csv"
)


# ============================================================
# 13. Compound variance explained
#
# Molecular explanatory power is quantified as the reduction
# in compound random-effect variance after adding the fixed
# molecular descriptor panel.
# ============================================================

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


compound_variance_reduction <- purrr::imap_dfr(
	full_model_results,
	function(x, spec_id) {

		baseline_variance <-
			get_compound_variance(
				x$baseline_model
			)


		joint_variance <-
			get_compound_variance(
				x$joint_model
			)


		tibble(
			specification =
				spec_id,

			specification_label =
				x$specification_label,

			baseline_compound_variance =
				baseline_variance,

			joint_compound_variance =
				joint_variance,

			proportion_compound_variance_explained =
				1 -
					joint_variance /
						baseline_variance,

			percent_compound_variance_explained =
				100 *
					(
						1 -
							joint_variance /
								baseline_variance
					)
		)
	}
)


readr::write_csv(
	compound_variance_reduction,
	"R/results/20_compound_variance_reduction.csv"
)


# ============================================================
# 14. Export model coefficients
#
# Coefficients are reported for interpretation.
# They are NOT used for variable selection.
# ============================================================

qspr_coefficients <- purrr::imap_dfr(
	full_model_results,
	function(x, spec_id) {

		raw_coefficients <- broom::tidy(
			x$raw_model,
			conf.int = TRUE
		) %>%
			transmute(
				specification =
					spec_id,

				specification_label =
					x$specification_label,

				approach =
					"Compound-mean QSPR",

				term,
				estimate,
				std.error,
				statistic,
				p.value,
				conf.low,
				conf.high
			)


		joint_coefficients <- broom.mixed::tidy(
			x$joint_model,
			effects = "fixed",
			conf.int = TRUE
		) %>%
			transmute(
				specification =
					spec_id,

				specification_label =
					x$specification_label,

				approach =
					"Joint hierarchical QSPR",

				term,
				estimate,
				std.error,
				statistic,
				p.value,
				conf.low,
				conf.high
			)


		bind_rows(
			raw_coefficients,
			joint_coefficients
		)
	}
)


readr::write_csv(
	qspr_coefficients,
	"R/results/20_qspr_coefficients.csv"
)


# ============================================================
# 15. Leave-one-compound-out cross-validation
#
# All observations from the held-out compound are removed.
#
# Descriptor scaling is recalculated within each training
# fold.
#
# The molecular specification itself never changes.
#
# For joint-model predictions:
# 	- the held-out compound random effect is zero because the
# 	  compound is unseen;
# 	- reference effects are retained for references observed
# 	  in the training fold;
# 	- unseen references receive random effect zero.
# ============================================================

cv_predictions <- list()

cv_failures <- list()

prediction_counter <- 1

failure_counter <- 1


for (
	spec_id in names(
		specifications
	)
) {

	spec_descriptors <-
		specifications[[spec_id]]


	spec_label <-
		specification_labels[[spec_id]]


	eligible_compounds <- compound_endpoint %>%
		filter(
			if_all(
				all_of(
					spec_descriptors
				),
				~ !is.na(.)
			)
		)


	compound_ids <-
		eligible_compounds$compound_id


	formulas <- make_formulas(
		spec_descriptors
	)


	for (
		i in seq_along(
			compound_ids
		)
	) {

		held_out <-
			compound_ids[i]


		message(
			spec_id,
			" LOCO ",
			i,
			"/",
			length(
				compound_ids
			),
			": ",
			held_out
		)


		train_compounds <- eligible_compounds %>%
			filter(
				compound_id !=
					held_out
			)


		test_compound <- eligible_compounds %>%
			filter(
				compound_id ==
					held_out
			)


		train_rows <- dat %>%
			filter(
				compound_id %in%
					train_compounds$compound_id
			)


		test_rows <- dat %>%
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
						spec_descriptors
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
						spec_descriptors
					),
					~ sd(
						.,
						na.rm = TRUE
					)
				)
			)


		for (
			variable in spec_descriptors
		) {

			z_name <-
				paste0(
					variable,
					"_z"
				)


			train_compounds[[z_name]] <- (
					train_compounds[[variable]] - fold_means[[variable]]
				) /
					fold_sds[[variable]]


			test_compound[[z_name]] <- (
					test_compound[[variable]] - fold_means[[variable]]
				) /
					fold_sds[[variable]]
		}


		# ----------------------------------------------------
		# Conventional compound-mean QSPR
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

			cv_failures[[failure_counter]] <- tibble(
				specification =
					spec_id,

				specification_label =
					spec_label,

				compound_id =
					held_out,

				approach =
					"Compound-mean QSPR",

				reason =
					"Model fit failed"
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
						spec_descriptors,
						"_z"
					)
				)
			)


		test_descriptor_data <- test_compound %>%
			select(
				compound_id,
				all_of(
					paste0(
						spec_descriptors,
						"_z"
					)
				)
			)


		train_rows <- train_rows %>%
			select(
				-all_of(
					descriptor_variables
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
					descriptor_variables
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

			cv_failures[[failure_counter]] <- tibble(
				specification =
					spec_id,

				specification_label =
					spec_label,

				compound_id =
					held_out,

				approach =
					"Joint hierarchical QSPR",

				reason =
					"Model fit failed"
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


		# ----------------------------------------------------
		# Store held-out predictions
		# ----------------------------------------------------

		cv_predictions[[prediction_counter]] <- tibble(
			specification =
				spec_id,

			specification_label =
				spec_label,

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
}


cv_predictions <- bind_rows(
	cv_predictions
)


cv_failures <- bind_rows(
	cv_failures
)


readr::write_csv(
	cv_predictions,
	"R/results/20_loco_predictions.csv"
)


readr::write_csv(
	cv_failures,
	"R/results/20_loco_failures.csv"
)


# ============================================================
# 16. Row-level LOCO performance
#
# These metrics evaluate prediction of the individual
# literature measurements for an unseen compound.
# ============================================================

loco_performance <- purrr::map_dfr(
	names(
		specifications
	),
	function(spec_id) {

		x <- cv_predictions %>%
			filter(
				specification ==
					spec_id
			)


		bind_rows(

			tibble(
				specification =
					spec_id,

				specification_label =
					specification_labels[[spec_id]],

				approach =
					"Compound-mean QSPR",

				n =
					sum(
						complete.cases(
							x$observed_logKp,
							x$compound_mean_prediction
						)
					),

				RMSE =
					rmse(
						x$observed_logKp,
						x$compound_mean_prediction
					),

				MAE =
					mae(
						x$observed_logKp,
						x$compound_mean_prediction
					),

				R2_pred =
					r2_pred(
						x$observed_logKp,
						x$compound_mean_prediction
					),

				Pearson_r =
					cor(
						x$observed_logKp,
						x$compound_mean_prediction,
						use =
							"complete.obs"
					)
			),

			tibble(
				specification =
					spec_id,

				specification_label =
					specification_labels[[spec_id]],

				approach =
					"Joint hierarchical QSPR",

				n =
					sum(
						complete.cases(
							x$observed_logKp,
							x$joint_hierarchical_prediction
						)
					),

				RMSE =
					rmse(
						x$observed_logKp,
						x$joint_hierarchical_prediction
					),

				MAE =
					mae(
						x$observed_logKp,
						x$joint_hierarchical_prediction
					),

				R2_pred =
					r2_pred(
						x$observed_logKp,
						x$joint_hierarchical_prediction
					),

				Pearson_r =
					cor(
						x$observed_logKp,
						x$joint_hierarchical_prediction,
						use =
							"complete.obs"
					)
			)
		)
	}
)


readr::write_csv(
	loco_performance,
	"R/results/20_loco_performance.csv"
)


# ============================================================
# 17. Compound-balanced LOCO performance
#
# Errors are first summarized within each held-out compound,
# then averaged across compounds. Each chemical therefore
# contributes equally regardless of replication count.
# ============================================================

loco_compound_errors <- cv_predictions %>%
	group_by(
		specification,
		specification_label,
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


readr::write_csv(
	loco_compound_errors,
	"R/results/20_loco_compound_errors.csv"
)


loco_compound_balanced <- loco_compound_errors %>%
	group_by(
		specification,
		specification_label
	) %>%
	summarise(
		compound_mean_mean_RMSE =
			mean(
				compound_mean_RMSE,
				na.rm = TRUE
			),

		joint_mean_RMSE =
			mean(
				joint_RMSE,
				na.rm = TRUE
			),

		compound_mean_mean_MAE =
			mean(
				compound_mean_MAE,
				na.rm = TRUE
			),

		joint_mean_MAE =
			mean(
				joint_MAE,
				na.rm = TRUE
			),

		compound_mean_median_MAE =
			median(
				compound_mean_MAE,
				na.rm = TRUE
			),

		joint_median_MAE =
			median(
				joint_MAE,
				na.rm = TRUE
			),

		.groups = "drop"
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


readr::write_csv(
	loco_compound_balanced,
	"R/results/20_loco_compound_balanced_performance.csv"
)


# ============================================================
# 18. Direct hierarchy advantage summary
# ============================================================

hierarchy_advantage <- loco_performance %>%
	select(
		specification,
		specification_label,
		approach,
		RMSE,
		MAE,
		R2_pred,
		Pearson_r
	) %>%
	pivot_wider(
		names_from =
			approach,
		values_from =
			c(
				RMSE,
				MAE,
				R2_pred,
				Pearson_r
			)
	) %>%
	mutate(
		RMSE_improvement_percent =
			100 *
				(
					`RMSE_Compound-mean QSPR` -
						`RMSE_Joint hierarchical QSPR`
				) /
					`RMSE_Compound-mean QSPR`,

		MAE_improvement_percent =
			100 *
				(
					`MAE_Compound-mean QSPR` -
						`MAE_Joint hierarchical QSPR`
				) /
					`MAE_Compound-mean QSPR`,

		R2_improvement =
			`R2_pred_Joint hierarchical QSPR` -
				`R2_pred_Compound-mean QSPR`
	)


readr::write_csv(
	hierarchy_advantage,
	"R/results/20_hierarchy_advantage.csv"
)


# ============================================================
# 19. Save full-data model objects
# ============================================================

saveRDS(
	full_model_results,
	"R/results/20_qspr_models.rds"
)


# ============================================================
# 20. Export molecular specification definitions
# ============================================================

specification_definition <- tibble(
	specification = c(
		"S1",
		"S2"
	),

	specification_label = c(
		"Size-lipophilicity",
		"Extended physicochemical"
	),

	descriptors = c(
		"MW + logKow",

		"MW + logKow + melting point + logSaq + HBD + HBA"
	),

	selection = c(
		"None",
		"None"
	),

	rationale = c(
		"Minimal molecular size and lipophilicity representation",

		"Prespecified interpretable physicochemical panel representing size, lipophilicity, solid-state constraint, aqueous solubility, and hydrogen-bonding capacity"
	)
)


readr::write_csv(
	specification_definition,
	"R/results/20_specification_definition.csv"
)


# ============================================================
# 21. Console summary
# ============================================================

cat("\n")
cat("========================================\n")
cat("20_qspr_analysis.R complete\n")
cat("========================================\n\n")

cat("QSPR samples:\n\n")

print(
	qspr_sample_summary
)

cat("\nCompound variance explained:\n\n")

print(
	compound_variance_reduction
)

cat("\nLOCO performance:\n\n")

print(
	loco_performance
)

cat("\nHierarchy advantage:\n\n")

print(
	hierarchy_advantage
)

cat("\nCompound-balanced LOCO performance:\n\n")

print(
	loco_compound_balanced
)

cat("\nLOCO failures: ")

if (nrow(cv_failures) == 0) {
	cat("none\n")
} else {
	cat(
		nrow(
			cv_failures
		),
		"\n"
	)
}

cat("\nOutputs written to:\n")

cat(
	"  R/results/20_specification_definition.csv\n"
)

cat(
	"  R/results/20_qspr_sample_summary.csv\n"
)

cat(
	"  R/results/20_descriptor_consistency.csv\n"
)

cat(
	"  R/results/20_descriptor_value_audit.csv\n"
)

cat(
	"  R/results/20_model_comparison.csv\n"
)

cat(
	"  R/results/20_compound_variance_reduction.csv\n"
)

cat(
	"  R/results/20_qspr_coefficients.csv\n"
)

cat(
	"  R/results/20_loco_predictions.csv\n"
)

cat(
	"  R/results/20_loco_failures.csv\n"
)

cat(
	"  R/results/20_loco_performance.csv\n"
)

cat(
	"  R/results/20_loco_compound_errors.csv\n"
)

cat(
	"  R/results/20_loco_compound_balanced_performance.csv\n"
)

cat(
	"  R/results/20_hierarchy_advantage.csv\n"
)

cat(
	"  R/results/20_qspr_models.rds\n"
)