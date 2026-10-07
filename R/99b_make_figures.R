# ============================================================
# 99b_make_figures.R
#
# Purpose:
# 	Generate manuscript-ready figures from finalized analysis
# 	outputs.
#
# Main figures:
#
# 	Figure 1
# 		Variance decomposition
#
# 	Figure 2
# 		Compound-disjoint observed vs predicted logKp
#
# 	Figure 3
# 		LOCO RMSE across fixed molecular specifications
#
# Supplementary figures:
#
# 	Figure S1
# 		Reference-level random effects
#
# 	Figure S2
# 		Sensitivity analyses
#
# 	Figure S3
# 		Residual Q-Q plot
#
# 	Figure S4
# 		Residuals versus fitted values
#
# This script DOES NOT:
# 	- fit models
# 	- recalculate analyses
#
# Outputs:
# 	manuscript/figures/
#
# Run from project root:
# 	source("R/99b_make_figures.R")
# ============================================================


source("R/00_config.R")


# ============================================================
# 1. Create output directory
# ============================================================

dir.create(
	"manuscript/figures",
	recursive = TRUE,
	showWarnings = FALSE
)


# ============================================================
# 2. Color palette
# ============================================================

bay <- palette(7)

col_darkblue <- bay[1]
col_blue <- bay[2]
col_teal <- bay[3]
col_beige <- bay[4]
col_yellow <- bay[5]
col_orange <- bay[6]
col_red <- bay[7]


approach_colors <- c(
	"Compound-mean QSPR" = col_orange,
	"Joint hierarchical QSPR" = col_teal
)


variance_colors <- c(
	"Compound" = col_darkblue,
	"Reference" = col_teal,
	"Residual" = col_orange
)


# ============================================================
# 3. Shared figure theme
# ============================================================

manuscript_theme <- theme_bw(
	base_size = 11
) +
	theme(
		panel.grid.minor =
			element_blank(),

		panel.grid.major =
			element_line(
				colour = "grey90",
				linewidth = 0.4
			),

		plot.title =
			element_text(
				face = "bold",
				size = 13,
				margin = margin(
					b = 10
				)
			),

		plot.subtitle =
			element_text(
				size = 10
			),

		axis.title =
			element_text(
				size = 11
			),

		axis.text =
			element_text(
				size = 10
			),

		strip.text =
			element_text(
				face = "bold",
				size = 10
			),

		legend.position =
			"bottom",

		legend.title =
			element_blank(),

		plot.margin =
			margin(
				10,
				12,
				10,
				10
			)
	)


# ============================================================
# 4. Figure 1: Variance decomposition
#
# Primary hierarchical model M2.
# ============================================================

variance_data <- readr::read_csv(
	"R/results/10_variance_proportions.csv",
	show_col_types = FALSE
) %>%
	filter(
		model == "M2"
	) %>%
	mutate(
		component =
			factor(
				component,
				levels = c(
					"Compound",
					"Reference",
					"Residual"
				)
			)
	)


figure_1 <- ggplot(
	variance_data,
	aes(
		x = component,
		y = percent,
		fill = component
	)
) +
	geom_col(
		width = 0.68
	) +
	geom_text(
		aes(
			label =
				paste0(
					round(
						percent,
						1
					),
					"%"
				)
		),
		vjust = -0.45,
		size = 4
	) +
	scale_fill_manual(
		values =
			variance_colors
	) +
	labs(
		x = NULL,
		y = "Variance (%)",
		title = "Variance decomposition of logKp"
	) +
	coord_cartesian(
		ylim = c(
			0,
			max(
				variance_data$percent
			) *
				1.14
		)
	) +
	guides(
		fill = "none"
	) +
	manuscript_theme


ggsave(
	"manuscript/figures/figure_1_variance_decomposition.png",
	figure_1,
	width = 6.2,
	height = 4.5,
	dpi = 400
)


ggsave(
	"manuscript/figures/figure_1_variance_decomposition.pdf",
	figure_1,
	width = 6.2,
	height = 4.5
)


# ============================================================
# 5. Figure 2: LOCO observed versus predicted
#
# Extended physicochemical specification S2.
# ============================================================

loco_predictions <- readr::read_csv(
	"R/results/20_loco_predictions.csv",
	show_col_types = FALSE
)


prediction_plot_data <- loco_predictions %>%
	filter(
		specification == "S2"
	) %>%
	select(
		compound_id,
		source_row,
		observed_logKp,
		compound_mean_prediction,
		joint_hierarchical_prediction
	) %>%
	pivot_longer(
		cols = c(
			compound_mean_prediction,
			joint_hierarchical_prediction
		),
		names_to = "approach",
		values_to = "predicted_logKp"
	) %>%
	mutate(
		approach = recode(
			approach,
			compound_mean_prediction =
				"Compound-mean QSPR",
			joint_hierarchical_prediction =
				"Joint hierarchical QSPR"
		),

		approach =
			factor(
				approach,
				levels = c(
					"Compound-mean QSPR",
					"Joint hierarchical QSPR"
				)
			)
	)


prediction_limits <- range(
	c(
		prediction_plot_data$observed_logKp,
		prediction_plot_data$predicted_logKp
	),
	na.rm = TRUE
)


figure_2 <- ggplot(
	prediction_plot_data,
	aes(
		x = observed_logKp,
		y = predicted_logKp,
		colour = approach
	)
) +
	geom_abline(
		intercept = 0,
		slope = 1,
		linetype = 2,
		colour = "grey30"
	) +
	geom_point(
		alpha = 0.55,
		size = 1.6
	) +
	facet_wrap(
		~ approach,
		nrow = 1
	) +
	scale_colour_manual(
		values =
			approach_colors
	) +
	coord_equal(
		xlim =
			prediction_limits,
		ylim =
			prediction_limits
	) +
	labs(
		x = "Observed logKp",
		y = "LOCO predicted logKp",
		title = "Compound-disjoint prediction: extended physicochemical model"
	) +
	guides(
		colour = "none"
	) +
	manuscript_theme


ggsave(
	"manuscript/figures/figure_2_loco_predictions_s2.png",
	figure_2,
	width = 8.5,
	height = 4.5,
	dpi = 400
)


ggsave(
	"manuscript/figures/figure_2_loco_predictions_s2.pdf",
	figure_2,
	width = 8.5,
	height = 4.5
)


# ============================================================
# 6. Figure 3: LOCO RMSE across molecular specifications
# ============================================================

loco_performance <- readr::read_csv(
	"R/results/20_loco_performance.csv",
	show_col_types = FALSE
) %>%
	mutate(
		specification_label =
			factor(
				specification_label,
				levels = c(
					"Size-lipophilicity",
					"Extended physicochemical"
				)
			),

		approach =
			factor(
				approach,
				levels = c(
					"Compound-mean QSPR",
					"Joint hierarchical QSPR"
				)
			)
	)


figure_3 <- ggplot(
	loco_performance,
	aes(
		x = specification_label,
		y = RMSE,
		fill = approach
	)
) +
	geom_col(
		position =
			position_dodge(
				width = 0.8
			),
		width = 0.7
	) +
	geom_text(
		aes(
			label =
				sprintf(
					"%.3f",
					RMSE
				)
		),
		position =
			position_dodge(
				width = 0.8
			),
		vjust = -0.45,
		size = 3.4
	) +
	scale_fill_manual(
		values =
			approach_colors
	) +
	labs(
		x = NULL,
		y = "LOCO RMSE",
		title = "Compound-disjoint prediction across QSPR approaches"
	) +
	coord_cartesian(
		ylim = c(
			0,
			max(
				loco_performance$RMSE
			) *
				1.15
		)
	) +
	manuscript_theme


ggsave(
	"manuscript/figures/figure_3_qspr_rmse.png",
	figure_3,
	width = 7,
	height = 4.5,
	dpi = 400
)


ggsave(
	"manuscript/figures/figure_3_qspr_rmse.pdf",
	figure_3,
	width = 7,
	height = 4.5
)


# ============================================================
# 7. Supplementary Figure S1:
# Reference-level random effects
# ============================================================

reference_effects <- readr::read_csv(
	"R/results/11_reference_random_effects.csv",
	show_col_types = FALSE
) %>%
	arrange(
		random_effect
	) %>%
	mutate(
		reference_label =
			if_else(
				is.na(reference_raw),
				reference_id,
				reference_raw
			),

		reference_label =
			factor(
				reference_label,
				levels =
					reference_label
			)
	)


figure_s1 <- ggplot(
	reference_effects,
	aes(
		x = random_effect,
		y = reference_label,
		colour = random_effect
	)
) +
	geom_vline(
		xintercept = 0,
		linetype = 2,
		colour = "grey30"
	) +
	geom_errorbarh(
		aes(
			xmin = conf_low,
			xmax = conf_high
		),
		height = 0,
		colour = "grey45",
		linewidth = 0.45
	) +
	geom_point(
		size = 1.7
	) +
	scale_colour_gradientn(
		colours = palette(100)
	) +
	labs(
		x = "Reference-level effect on logKp",
		y = NULL,
		title = "Estimated reference-level heterogeneity"
	) +
	guides(
		colour = "none"
	) +
	manuscript_theme +
	theme(
		axis.text.y =
			element_text(
				size = 5.5
			)
	)


ggsave(
	"manuscript/figures/figure_s1_reference_effects.png",
	figure_s1,
	width = 8,
	height = 10,
	dpi = 400
)


ggsave(
	"manuscript/figures/figure_s1_reference_effects.pdf",
	figure_s1,
	width = 8,
	height = 10
)


# ============================================================
# 8. Supplementary Figure S2:
# Sensitivity analyses
# ============================================================

sensitivity_summary <- readr::read_csv(
	"R/results/30_sensitivity_summary.csv",
	show_col_types = FALSE
) %>%
	mutate(
		scenario =
			factor(
				scenario,
				levels =
					rev(
						scenario
					)
			)
	)


figure_s2 <- ggplot(
	sensitivity_summary,
	aes(
		x = RMSE_improvement_percent,
		y = scenario,
		fill = scenario
	)
) +
	geom_vline(
		xintercept = 0,
		linetype = 2,
		colour = "grey30"
	) +
	geom_col(
		width = 0.65
	) +
	geom_text(
		aes(
			label =
				paste0(
					round(
						RMSE_improvement_percent,
						1
					),
					"%"
				)
		),
		hjust = -0.15,
		size = 3.3
	) +
	scale_fill_manual(
		values = setNames(
			palette(
				nrow(
					sensitivity_summary
				)
			),
			levels(
				sensitivity_summary$scenario
			)
		)
	) +
	labs(
		x = "Reduction in LOCO RMSE (%)",
		y = NULL,
		title = "Hierarchical predictive advantage across sensitivity analyses"
	) +
	coord_cartesian(
		xlim = c(
			min(
				0,
				min(
					sensitivity_summary$RMSE_improvement_percent,
					na.rm = TRUE
				)
			),

			max(
				sensitivity_summary$RMSE_improvement_percent,
				na.rm = TRUE
			) *
				1.18
		)
	) +
	guides(
		fill = "none"
	) +
	manuscript_theme


ggsave(
	"manuscript/figures/figure_s2_sensitivity_rmse.png",
	figure_s2,
	width = 7.5,
	height = 4.8,
	dpi = 400
)


ggsave(
	"manuscript/figures/figure_s2_sensitivity_rmse.pdf",
	figure_s2,
	width = 7.5,
	height = 4.8
)


# ============================================================
# 9. Supplementary Figure S3:
# Residual Q-Q plot
# ============================================================

residual_qq <- readr::read_csv(
	"R/results/11_residual_qq.csv",
	show_col_types = FALSE
)


figure_s3 <- ggplot(
	residual_qq,
	aes(
		x = expected,
		y = observed
	)
) +
	geom_abline(
		intercept = 0,
		slope = 1,
		linetype = 2,
		colour = "grey30"
	) +
	geom_point(
		size = 1.5,
		colour = col_blue,
		alpha = 0.75
	) +
	labs(
		x = "Expected normal quantile",
		y = "Observed standardized residual quantile",
		title = "Residual Q-Q plot"
	) +
	manuscript_theme


ggsave(
	"manuscript/figures/figure_s3_residual_qq.png",
	figure_s3,
	width = 5,
	height = 5,
	dpi = 400
)


ggsave(
	"manuscript/figures/figure_s3_residual_qq.pdf",
	figure_s3,
	width = 5,
	height = 5
)


# ============================================================
# 10. Supplementary Figure S4:
# Residuals versus fitted values
# ============================================================

observation_diagnostics <- readr::read_csv(
	"R/results/11_observation_diagnostics.csv",
	show_col_types = FALSE
)


figure_s4 <- ggplot(
	observation_diagnostics,
	aes(
		x = fitted,
		y = residual
	)
) +
	geom_hline(
		yintercept = 0,
		linetype = 2,
		colour = "grey30"
	) +
	geom_point(
		alpha = 0.5,
		size = 1.5,
		colour = col_blue
	) +
	geom_smooth(
		method = "loess",
		se = FALSE,
		colour = col_orange,
		linewidth = 0.9
	) +
	labs(
		x = "Fitted logKp",
		y = "Residual",
		title = "Residuals versus fitted values"
	) +
	manuscript_theme


ggsave(
	"manuscript/figures/figure_s4_residual_vs_fitted.png",
	figure_s4,
	width = 6,
	height = 4.5,
	dpi = 400
)


ggsave(
	"manuscript/figures/figure_s4_residual_vs_fitted.pdf",
	figure_s4,
	width = 6,
	height = 4.5
)


# ============================================================
# 11. Figure manifest
# ============================================================

figure_manifest <- tibble(
	figure = c(
		"Figure 1",
		"Figure 2",
		"Figure 3",
		"Figure S1",
		"Figure S2",
		"Figure S3",
		"Figure S4"
	),

	file = c(
		"figure_1_variance_decomposition",
		"figure_2_loco_predictions_s2",
		"figure_3_qspr_rmse",
		"figure_s1_reference_effects",
		"figure_s2_sensitivity_rmse",
		"figure_s3_residual_qq",
		"figure_s4_residual_vs_fitted"
	),

	title = c(
		"Variance decomposition of logKp",
		"Compound-disjoint prediction using the extended physicochemical model",
		"Compound-disjoint prediction across QSPR approaches",
		"Estimated reference-level heterogeneity",
		"Hierarchical predictive advantage across sensitivity analyses",
		"Residual Q-Q plot",
		"Residuals versus fitted values"
	)
)


readr::write_csv(
	figure_manifest,
	"manuscript/figures/figure_manifest.csv"
)


# ============================================================
# 12. Console summary
# ============================================================

cat("\n")
cat("========================================\n")
cat("99b_make_figures.R complete\n")
cat("========================================\n\n")

cat(
	"Figures written to manuscript/figures/\n\n"
)

print(
	figure_manifest
)