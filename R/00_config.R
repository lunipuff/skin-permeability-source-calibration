library("dplyr")
library("readxl")
library("knitr")

############################################################
# 00_config.R
# Project-wide settings for skin permeability QSPR analysis
############################################################


############################################################
# File paths
############################################################

path_cleaned_dataset <- "data/processed/skin_permeability_modeling_dataset.csv"

# Paths for workflow
path_workflow_svg <- "figures/figure1_workflow.svg"
path_workflow_pdf <- "figures/figure1_workflow.pdf"
path_workflow_png <- "manuscript/figures/figure1_workflow.png"

# Paths for cleaning
path_raw_excel <- "data/raw/Table 1.xlsx"
path_interim_raw_csv <- "data/interim/raw_dataset.csv"
path_removed_rows <- "data/interim/removed_rows_log.csv"
path_cleaning_summary <- "data/interim/cleaning_summary.csv"
path_descriptor_inconsistency_log <- "data/interim/descriptor_inconsistency_within_merge_groups.csv"
path_required_missingness_summary <- "results/cleaning/required_missingness_summary.csv"
path_cleaning_flow <- "results/cleaning/cleaning_flow.csv"

# Cleaning flow table
path_cleaning_flow_table <- "tables/table_cleaning_flow.csv"
path_cleaning_flow_results <- "results/model_comparison/01b_cleaning_flow_table.csv"

# Solubility redundancy check paths
path_solubility_redundancy_table <- "results/predictor_screening/solubility_redundancy_check.csv"
path_table_solubility_redundancy <- "tables/tableS_solubility_redundancy_check.csv"
path_fig_gse_logs_vs_logsaqd_pdf <- "figures/figureS_gse_logS_vs_LogSaqd.pdf"
path_fig_gse_logs_vs_logsaqd_png <- "figures/figureS_gse_logS_vs_LogSaqd.png"
path_fig_solubility_correlation <- "figures/figureS_solubility_descriptor_correlation.pdf"

# Paths for predictor screening
path_predictor_screening_full <- "results/predictor_screening/initial_predictor_screening_full.csv"
path_predictor_screening_selected <- "results/predictor_screening/initial_predictor_screening_selected.csv"
path_fig_predictor_screening_rmse_by_n_pdf <- "figures/figureS_predictor_screening_RMSE_by_n_predictors.pdf"
path_fig_predictor_screening_rmse_by_n_png <- "figures/figureS_predictor_screening_RMSE_by_n_predictors.png"

# LOCO-CV candidate model search paths
path_loco_cv_linear <- "results/cross_validation/03_loco_cv_linear_models.csv"
path_loco_cv_quadratic <- "results/cross_validation/03_loco_cv_quadratic_models.csv"
path_loco_cv_log <- "results/cross_validation/03_loco_cv_log_models.csv"
path_loco_cv_additive <- "results/cross_validation/03_loco_cv_additive_models.csv"
path_loco_cv_interaction <- "results/cross_validation/03_loco_cv_interaction_models.csv"
path_loco_cv_final_results <- "results/cross_validation/03_loco_cv_candidate_model_search.csv"
path_loco_cv_final_predictions <- "results/cross_validation/03_loco_cv_final_model_predictions.csv"
path_loco_cv_selected_model <- "results/cross_validation/03_selected_model.txt"
path_fig_loco_cv_prediction_pdf <- "figures/figure_loco_cv_observed_vs_predicted.pdf"
path_fig_loco_cv_prediction_png <- "figures/figure_loco_cv_observed_vs_predicted.png"

# LOCO-CV coefficient stability paths
path_loco_cv_coefficient_mean <- "results/cross_validation/04_loco_cv_coefficient_mean.csv"
path_loco_cv_coefficient_sd <- "results/cross_validation/04_loco_cv_coefficient_sd.csv"
path_loco_cv_coefficient_rsd <- "results/cross_validation/04_loco_cv_coefficient_rsd.csv"
path_loco_cv_coefficient_long <- "results/cross_validation/04_loco_cv_coefficient_long.csv"
path_loco_cv_selected_model_coefficient_summary <- "tables/tableS_selected_model_coefficient_stability.csv"
path_fig_selected_model_coefficient_rsd_pdf <- "figures/figureS_selected_model_coefficient_rsd.pdf"
path_fig_selected_model_coefficient_rsd_png <- "figures/figureS_selected_model_coefficient_rsd.png"
path_fig_selected_model_coefficient_boxplot_pdf <- "figures/figureS_selected_model_coefficient_boxplot.pdf"
path_fig_selected_model_coefficient_boxplot_png <- "figures/figureS_selected_model_coefficient_boxplot.png"

# LOCO-CV error diagnostic paths
path_loco_cv_error_diagnostic_dataset <- "results/cross_validation/05_loco_cv_error_diagnostic_dataset.csv"
path_loco_cv_error_summary <- "results/cross_validation/05_loco_cv_error_summary.csv"
path_loco_cv_high_error_compounds <- "tables/tableS_high_error_compounds.csv"
path_fig_loco_cv_error_observed_predicted_pdf <- "figures/figureS_loco_cv_error_observed_predicted.pdf"
path_fig_loco_cv_error_observed_predicted_png <- "figures/figureS_loco_cv_error_observed_predicted.png"
path_fig_loco_cv_residual_by_descriptor_pdf <- "figures/figureS_loco_cv_residual_by_descriptor.pdf"
path_fig_loco_cv_residual_by_descriptor_png <- "figures/figureS_loco_cv_residual_by_descriptor.png"
path_fig_loco_cv_absolute_error_by_descriptor_pdf <- "figures/figureS_loco_cv_absolute_error_by_descriptor.pdf"
path_fig_loco_cv_absolute_error_by_descriptor_png <- "figures/figureS_loco_cv_absolute_error_by_descriptor.png"
path_fig_loco_cv_error_density_pdf <- "figures/figureS_loco_cv_error_density.pdf"
path_fig_loco_cv_error_density_png <- "figures/figureS_loco_cv_error_density.png"
path_fig_loco_cv_high_error_compounds_pdf <- "figures/figureS_loco_cv_high_error_compounds.pdf"
path_fig_loco_cv_high_error_compounds_png <- "figures/figureS_loco_cv_high_error_compounds.png"

# Heteroscedastic uncertainty model paths
path_hetero_full_fit_summary <- "results/uncertainty/06_heteroscedastic_full_fit_summary.csv"
path_hetero_cv_summary <- "results/uncertainty/06_heteroscedastic_cv_summary.csv"
path_hetero_cv_predictions_all <- "results/uncertainty/06_heteroscedastic_cv_predictions_all_models.csv"
path_hetero_cv_predictions_best <- "results/uncertainty/06_heteroscedastic_cv_predictions_best_model.csv"
path_hetero_fold_info <- "results/uncertainty/06_heteroscedastic_fold_info.csv"
path_hetero_selected_model <- "results/uncertainty/06_selected_heteroscedastic_model.txt"
path_table_hetero_cv_summary <- "tables/tableS_heteroscedastic_cv_summary.csv"
path_fig_hetero_diagnostics_pdf <- "figures/figureS_heteroscedastic_uncertainty_diagnostics.pdf"
path_fig_hetero_diagnostics_png <- "figures/figureS_heteroscedastic_uncertainty_diagnostics.png"

# Descriptor relationship and domain overview paths
path_descriptor_relationship_dataset <- "results/applicability_domain/07_descriptor_relationship_dataset.csv"
path_descriptor_central_domain_dataset <- "results/applicability_domain/07_central_descriptor_domain_dataset.csv"
path_descriptor_central_domain_summary <- "results/applicability_domain/07_central_descriptor_domain_summary.csv"
path_fig_descriptor_logkp_relationships_pdf <- "figures/figureS_descriptor_logKp_relationships.pdf"
path_fig_descriptor_logkp_relationships_png <- "figures/figureS_descriptor_logKp_relationships.png"

# Benchmark model paths
path_benchmark_model_summary <- "results/benchmarks/08_benchmark_model_summary.csv"
path_benchmark_model_predictions <- "results/benchmarks/08_benchmark_model_predictions.csv"
path_table_benchmark_model_summary <- "tables/table_benchmark_model_summary.csv"
path_fig_benchmark_rmse_pdf <- "figures/figure_benchmark_model_RMSE.pdf"
path_fig_benchmark_rmse_png <- "figures/figure_benchmark_model_RMSE.png"
path_fig_benchmark_observed_predicted_pdf <- "figures/figure_benchmark_observed_vs_predicted.pdf"
path_fig_benchmark_observed_predicted_png <- "figures/figure_benchmark_observed_vs_predicted.png"

# RDKit descriptor benchmark paths
path_rdkit_descriptors_as_reported <- "data/processed/rdkit_descriptors_as_reported.csv"
path_rdkit_descriptors_parent <- "data/processed/rdkit_descriptors_parent.csv"
path_rdkit_descriptor_log <- "data/interim/rdkit_descriptor_calculation_log.csv"
path_rdkit_benchmark_summary <- "results/benchmarks/08b_rdkit_benchmark_summary.csv"
path_rdkit_benchmark_predictions <- "results/benchmarks/08b_rdkit_benchmark_predictions.csv"
path_table_rdkit_benchmark_summary <- "tables/tableS_rdkit_benchmark_summary.csv"
path_fig_rdkit_benchmark_rmse_pdf <- "figures/figureS_rdkit_benchmark_RMSE.pdf"
path_fig_rdkit_benchmark_rmse_png <- "figures/figureS_rdkit_benchmark_RMSE.png"
path_fig_rdkit_observed_predicted_pdf <- "figures/figureS_rdkit_observed_vs_predicted.pdf"
path_fig_rdkit_observed_predicted_png <- "figures/figureS_rdkit_observed_vs_predicted.png"

# Compiled benchmark table paths
path_table_benchmark_all <- "tables/tableS_all_benchmark_models.csv"
path_table_benchmark_main <- "tables/table_main_benchmark_models.csv"

# Ablation analysis paths
path_ablation_summary <- "results/ablation/09_ablation_summary.csv"
path_ablation_predictions <- "results/ablation/09_ablation_predictions.csv"
path_table_ablation_summary <- "tables/tableS_ablation_summary.csv"
path_fig_ablation_rmse_pdf <- "figures/figureS_ablation_RMSE.pdf"
path_fig_ablation_rmse_png <- "figures/figureS_ablation_RMSE.png"

# Partial-effect analysis paths
path_partial_effect_dataset <- "results/interpretability/10_partial_effect_dataset.csv"
path_interaction_effect_dataset <- "results/interpretability/10_interaction_effect_dataset.csv"
path_fig_partial_effects_pdf <- "figures/figure_partial_effects_selected_model.pdf"
path_fig_partial_effects_png <- "figures/figure_partial_effects_selected_model.png"
path_fig_interaction_MWa_LogSaqd_pdf <- "figures/figure_interaction_MWa_LogSaqd.pdf"
path_fig_interaction_MWa_LogSaqd_png <- "figures/figure_interaction_MWa_LogSaqd.png"

# Row-wise validation sensitivity paths
path_rowwise_cv_predictions <- "results/cross_validation/11_rowwise_cv_predictions.csv"
path_rowwise_cv_summary <- "results/cross_validation/11_rowwise_cv_summary.csv"
path_validation_scheme_summary <- "tables/table_validation_scheme_sensitivity.csv"
path_fig_validation_scheme_rmse_pdf <- "figures/figure_validation_scheme_RMSE.pdf"
path_fig_validation_scheme_rmse_png <- "figures/figure_validation_scheme_RMSE.png"
path_fig_rowwise_observed_predicted_pdf <- "figures/figureS_rowwise_cv_observed_vs_predicted.pdf"
path_fig_rowwise_observed_predicted_png <- "figures/figureS_rowwise_cv_observed_vs_predicted.png"

# Ref group validation sensitivity analysis
path_reference_cv_predictions <- "results/cross_validation/11b_reference_grouped_cv_predictions.csv"
path_reference_cv_overall_summary <- "results/cross_validation/11b_reference_grouped_cv_overall_summary.csv"
path_reference_cv_by_reference <- "results/cross_validation/11b_reference_grouped_cv_by_reference.csv"
path_fig_reference_cv_rmse_pdf <- "figures/figure_reference_grouped_cv_rmse.pdf"
path_fig_reference_cv_rmse_png <- "figures/figure_reference_grouped_cv_rmse.png"
path_fig_reference_cv_observed_predicted_pdf <- "figures/figure_reference_grouped_cv_observed_predicted.pdf"
path_fig_reference_cv_observed_predicted_png <- "figures/figure_reference_grouped_cv_observed_predicted.png"

path_validation_design_sensitivity_table <- "tables/table_validation_design_sensitivity.csv"

# Applicability-domain summary paths
path_ad_prediction_dataset <- "results/applicability_domain/12_applicability_domain_prediction_dataset.csv"
path_ad_descriptor_summary <- "results/applicability_domain/12_descriptor_domain_summary.csv"
path_ad_error_by_domain <- "results/applicability_domain/12_error_by_domain_class.csv"
path_ad_high_error_summary <- "results/applicability_domain/12_high_error_domain_summary.csv"
path_table_ad_descriptor_summary <- "tables/tableS_descriptor_domain_summary.csv"
path_table_ad_error_by_domain <- "tables/table_applicability_domain_error_summary.csv"
path_table_ad_high_error_summary <- "tables/tableS_high_error_domain_summary.csv"
path_fig_ad_error_by_domain_pdf <- "figures/figure_applicability_domain_error_by_domain.pdf"
path_fig_ad_error_by_domain_png <- "figures/figure_applicability_domain_error_by_domain.png"
path_fig_ad_abs_error_by_predictor_pdf <- "figures/figureS_applicability_domain_abs_error_by_predictor.pdf"
path_fig_ad_abs_error_by_predictor_png <- "figures/figureS_applicability_domain_abs_error_by_predictor.png"

# Final model export paths
path_final_model_coefficients <- "tables/table_final_model_coefficients.csv"
path_final_model_equation <- "results/model_comparison/13_final_model_equation.txt"
path_final_model_summary <- "results/model_comparison/13_final_model_summary.csv"
path_final_model_manuscript_values <- "tables/table_manuscript_key_results.csv"
path_table_final_model_performance <- "tables/table_final_model_performance_summary.csv"
path_table_final_model_selected_outputs <- "tables/table_final_model_selected_outputs.csv"
path_fig_final_model_coefficients_pdf <- "figures/figureS_final_model_coefficients.pdf"
path_fig_final_model_coefficients_png <- "figures/figureS_final_model_coefficients.png"

############################################################
# Predictor screening settings
############################################################

all_predictors <- c(
	"MWa",
	"logKowb",
	"Mptc",
	"LogSaqd",
	"LogSoce",
	"Hdf",
	"Hag",
	"MVh",
	"Texpi",
	"Skin.thicknessj"
)

fixed_predictors <- c("MWa")

# Solubility, lipophilicity, and melting point are related through
# established solubility models, including the General Solubility Equation.
# To reduce redundancy during simple-model screening, models containing
# logKowb, Mptc, and a solubility descriptor simultaneously are excluded.

redundant_predictor_sets <- list(
  c("logKowb", "Mptc", "LogSaqd"),
  c("logKowb", "Mptc", "LogSoce")
)

############################################################
# Functions
############################################################

rmse <- function(obs, pred) {
	sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

mae <- function(obs, pred) {
	mean(abs(obs - pred), na.rm = TRUE)
}

r2_pred <- function(obs, pred) {
	1 - sum((obs - pred)^2, na.rm = TRUE) /
		sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
}

bias <- function(obs, pred) {
	mean(pred - obs, na.rm = TRUE)
}

############################################################
# Benchmark model formulas and helper functions
############################################################

formula_null <- logkpl ~ 1

formula_potts_guy <- logkpl ~ MWa + logKowb

formula_linear_core <- logkpl ~ MWa + logKowb + Mptc + LogSaqd + LogSoce + Texpi

formula_quadratic_core <- logkpl ~ MWa + logKowb + I(logKowb^2) +
	Mptc + I(Mptc^2) +
	LogSaqd + I(LogSaqd^2) +
	Texpi + I(Texpi^2)

formula_extended_descriptor <- logkpl ~ MWa + logKowb + I(logKowb^2) +
	Mptc + I(Mptc^2) +
	LogSaqd + I(LogSaqd^2) +
	LogSoce +
	Hdf +
	Hag +
	MVh +
	Texpi + I(Texpi^2) +
	Skin.thicknessj

get_formula_predictors <- function(model_formula, outcome_col = "logkpl") {
	predictors <- all.vars(
		as.formula(model_formula)
	)

	predictors <- predictors[
		predictors != outcome_col
	]

	unique(predictors)
}