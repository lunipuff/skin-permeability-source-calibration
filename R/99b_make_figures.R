# ============================================================
# 99b_make_figures.R
# ============================================================

library(ggplot2)
library(patchwork)

source("R/00_config.R")

pal <- palette(7)

c_blue <- pal[1]
c_cyan <- pal[2]
c_teal <- pal[3]
c_cream <- pal[4]
c_yellow <- pal[5]
c_orange <- pal[6]
c_red <- pal[7]
c_grey <- "grey65"
c_grey_dark <- "grey40"
c_grey_light <- "grey82"

theme_manuscript <- function(base_size = 11) {
	theme_classic(base_size = base_size) +
		theme(
			plot.title = element_text(face = "bold", size = base_size + 1.5, margin = margin(b = 5)),
			plot.subtitle = element_text(size = base_size - 0.5, margin = margin(b = 9)),
			axis.title = element_text(size = base_size),
			axis.text = element_text(size = base_size - 0.5),
			strip.background = element_blank(),
			strip.text = element_text(face = "bold", size = base_size),
			legend.title = element_blank(),
			legend.key = element_blank(),
			legend.text = element_text(size = base_size - 1),
			plot.tag = element_text(face = "bold", size = base_size + 2),
			plot.margin = margin(14, 20, 14, 20)
		)
}

save_fig <- function(plot, filename, width, height) {
	ggsave(
		paste0("R/figures/", filename, ".pdf"),
		plot = plot,
		width = width,
		height = height,
		units = "in"
	)
}

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

matchref31 <- read.csv("R/results/31_calibration_matched_reference.csv", stringsAsFactors = FALSE)

eligible41 <- read.csv("R/results/41_skinpix_eligible_references.csv", stringsAsFactors = FALSE)
match41 <- read.csv("R/results/41_skinpix_calibration_matched_3anchor.csv", stringsAsFactors = FALSE)
ref41 <- read.csv("R/results/41_skinpix_calibration_reference.csv", stringsAsFactors = FALSE)
pred41 <- read.csv("R/results/41_skinpix_one_anchor_predictions.csv", stringsAsFactors = FALSE)
domain41 <- read.csv("R/results/41_skinpix_transfer_domain.csv", stringsAsFactors = FALSE)

ad21 <- read.csv("R/results/21_applicability_domain.csv", stringsAsFactors = FALSE)
ad31 <- read.csv("R/results/31_applicability_domain.csv", stringsAsFactors = FALSE)

loco60 <- read.csv("R/results/60_loco_bootstrap_ci.csv", stringsAsFactors = FALSE)
var60 <- read.csv("R/results/60_variance_component_ci.csv", stringsAsFactors = FALSE)
cal60 <- read.csv("R/results/60_huskin_calibration_ci.csv", stringsAsFactors = FALSE)
wrong60 <- read.csv("R/results/60_wrong_source_ci.csv", stringsAsFactors = FALSE)

# ============================================================
# Figure 1
# Variance decomposition
# ============================================================

fig1_dat <- var60[, c("dataset", "component", "estimate_percent")]
names(fig1_dat)[3] <- "percent"

fig1_dat$dataset[fig1_dat$dataset == "Cheruvu"] <- "Cheruvu et al., 2022"
fig1_dat$dataset <- factor(fig1_dat$dataset, levels = c("HuSkinDB", "Cheruvu et al., 2022"))
fig1_dat$component <- factor(fig1_dat$component, levels = c("Residual", "Source", "Compound"))
fig1_dat$label <- paste0(round(fig1_dat$percent, 1), "%")
fig1_dat$text_color <- ifelse(fig1_dat$component == "Compound", "white", "black")

fig1 <- ggplot(fig1_dat, aes(x = percent, y = dataset, fill = component)) +
	geom_col(width = 0.6) +
	geom_text(
		aes(label = label, color = text_color),
		position = position_stack(vjust = 0.5),
		fontface = "bold",
		size = 4
	) +
	scale_color_identity() +
	scale_fill_manual(
		values = c(
			"Compound" = c_blue,
			"Source" = c_teal,
			"Residual" = c_yellow
		),
		breaks = c("Compound", "Source", "Residual")
	) +
	scale_x_continuous(
		breaks = seq(0, 100, 20),
		expand = c(0, 0)
	) +
	coord_cartesian(xlim = c(0, 100), clip = "off") +
	labs(
		title = "Source-associated variation is substantial in two independent datasets",
		x = "Modeled random/residual variance (%)",
		y = NULL
	) +
	theme_manuscript(12) +
	theme(
		legend.position = "top",
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(18, 26, 14, 26)
	)

save_fig(fig1, "Figure_1_variance_decomposition", 8.6, 3.5)

# ============================================================
# Figure 2
# Source-blind versus known-source prediction
# ============================================================

fig2_dat <- loco60

fig2_dat$group <- paste(fig2_dat$dataset, fig2_dat$specification)
fig2_dat$group[fig2_dat$group == "Cheruvu S1"] <- "Cheruvu et al., 2022 — S1"
fig2_dat$group[fig2_dat$group == "Cheruvu S2"] <- "Cheruvu et al., 2022 — S2"
fig2_dat$group[fig2_dat$group == "HuSkinDB S1"] <- "HuSkinDB — S1"

fig2_dat$stage <- factor(
	fig2_dat$approach,
	levels = c(
		"conventional",
		"hierarchical_source_blind",
		"hierarchical_known_source"
	),
	labels = c(
		"Conventional\naggregate",
		"Hierarchical\nsource-blind",
		"Hierarchical\nknown-source"
	)
)

fig2_dat$rmse_label <- sprintf("%.3f", fig2_dat$RMSE)

fig2_dat$gain_label <- ifelse(
	fig2_dat$approach == "conventional",
	"",
	paste0(
		ifelse(fig2_dat$improvement_percent >= 0, "−", "+"),
		abs(round(fig2_dat$improvement_percent, 1)),
		"% RMSE"
	)
)

ylim2 <- range(c(fig2_dat$RMSE_ci_low, fig2_dat$RMSE_ci_high))
ylim2 <- c(max(0, ylim2[1] - 0.08), ylim2[2] + 0.11)

fig2 <- ggplot(fig2_dat, aes(x = stage, y = RMSE, group = group)) +
	geom_line(linewidth = 0.8, color = c_grey) +
	geom_errorbar(
		aes(ymin = RMSE_ci_low, ymax = RMSE_ci_high),
		width = 0.12,
		linewidth = 0.7,
		color = c_grey_dark
	) +
	geom_point(
		aes(fill = stage),
		shape = 21,
		size = 4.2,
		stroke = 0.6,
		color = "white"
	) +
	geom_text(
		aes(label = rmse_label),
		vjust = -1.45,
		size = 3.2
	) +
	geom_text(
		aes(label = gain_label),
		vjust = 1.9,
		size = 3,
		fontface = "bold"
	) +
	facet_wrap(~group, nrow = 1) +
	scale_fill_manual(
		values = c(
			"Conventional\naggregate" = c_cyan,
			"Hierarchical\nsource-blind" = c_teal,
			"Hierarchical\nknown-source" = c_yellow
		)
	) +
	coord_cartesian(ylim = ylim2, clip = "off") +
	labs(
		title = "Source context accounts for the reproducible predictive improvement",
		subtitle = "Points show LOCO RMSE; bars show 95% compound-bootstrap intervals",
		x = NULL,
		y = "Leave-one-compound-out RMSE"
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "none",
		panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.5),
		axis.text.x = element_text(size = 8.2, lineheight = 0.95),
		plot.margin = margin(16, 24, 18, 24)
	)

save_fig(fig2, "Figure_2_source_information_prediction", 10.5, 4.6)

# ============================================================
# Figure 3
# SkinPiX external transfer
# ============================================================

# ------------------------------------------------------------
# A. One-anchor performance in each unseen SkinPiX source
# ------------------------------------------------------------

p3a_dat <- ref41[ref41$anchors == 1, ]
p3a_dat$improvement <- 100 * (1 - p3a_dat$RMSE_calibrated / p3a_dat$RMSE_uncalibrated)
p3a_dat <- p3a_dat[order(p3a_dat$RMSE_uncalibrated, decreasing = TRUE), ]
p3a_dat$source <- paste0("SkinPiX source ", seq_len(nrow(p3a_dat)), "\n(", p3a_dat$reference_id, ")")
p3a_dat$source <- factor(p3a_dat$source, levels = rev(p3a_dat$source))

p3a <- ggplot(p3a_dat, aes(y = source)) +
	geom_segment(
		aes(x = RMSE_uncalibrated, xend = RMSE_calibrated, yend = source),
		color = c_grey_light,
		linewidth = 2
	) +
	geom_point(
		aes(x = RMSE_uncalibrated),
		shape = 21,
		size = 4.1,
		fill = "white",
		color = c_grey_dark,
		stroke = 1
	) +
	geom_point(
		aes(x = RMSE_calibrated),
		shape = 21,
		size = 4.3,
		fill = c_teal,
		color = c_teal
	) +
	geom_text(
		aes(
			x = pmax(RMSE_uncalibrated, RMSE_calibrated),
			label = paste0(
				ifelse(improvement >= 0, "↓ ", "↑ "),
				abs(round(improvement, 1)),
				"%"
			)
		),
		hjust = -0.15,
		size = 3.2,
		fontface = "bold",
		color = c_blue
	) +
	scale_x_continuous(expand = expansion(mult = c(0.03, 0.20))) +
	labs(
		title = "One shared compound calibrates an unseen SkinPiX source",
		subtitle = "Open circle: uncalibrated HuSkin prediction; filled circle: one-anchor calibration",
		x = "RMSE on SkinPiX-only compounds",
		y = NULL
	) +
	theme_manuscript(10) +
	theme(
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(12, 30, 10, 18)
	)

# ------------------------------------------------------------
# B. Compound-level absolute prediction error
# ------------------------------------------------------------

p3b_dat <- aggregate(
	cbind(abs_error_uncalibrated, abs_error_calibrated) ~ compound_id,
	data = pred41,
	mean
)

p3b_dat$improved <- p3b_dat$abs_error_calibrated < p3b_dat$abs_error_uncalibrated
n_improved <- sum(p3b_dat$improved)
n_total <- nrow(p3b_dat)

lim3b <- max(c(p3b_dat$abs_error_uncalibrated, p3b_dat$abs_error_calibrated))
lim3b <- lim3b * 1.03

p3b <- ggplot(
	p3b_dat,
	aes(x = abs_error_uncalibrated, y = abs_error_calibrated)
) +
	geom_abline(
		intercept = 0,
		slope = 1,
		linetype = 2,
		color = c_grey
	) +
	geom_point(
		aes(fill = improved),
		shape = 21,
		size = 2.8,
		color = "white",
		stroke = 0.35,
		alpha = 0.85
	) +
	scale_fill_manual(
		values = c(
			"TRUE" = c_teal,
			"FALSE" = c_red
		),
		labels = c(
			"TRUE" = "Lower error",
			"FALSE" = "Higher error"
		)
	) +
	coord_equal(xlim = c(0, lim3b), ylim = c(0, lim3b), clip = "off") +
	labs(
		title = "Most external compounds have lower error after calibration",
		subtitle = paste0(n_improved, "/", n_total, " SkinPiX-only compounds fall below the identity line"),
		x = "Absolute error before calibration",
		y = "Absolute error after calibration"
	) +
	theme_manuscript(9.8) +
	theme(
		legend.position = "bottom",
		legend.direction = "horizontal",
		plot.margin = margin(12, 16, 12, 16)
	)

# ------------------------------------------------------------
# C. Same three sources with 1-3 anchors
# ------------------------------------------------------------

matched3_refs <- eligible41$reference_id[eligible41$eligible & eligible41$n_anchors >= 3]
p3c_ind <- ref41[ref41$reference_id %in% matched3_refs & ref41$anchors <= 3, ]

source_cols <- setNames(c(c_cyan, c_teal, c_orange), matched3_refs)

p3c <- ggplot() +
	geom_hline(yintercept = 0, linetype = 3, color = c_grey) +
	geom_line(
		data = p3c_ind,
		aes(x = anchors, y = RMSE_improvement_percent, group = reference_id, color = reference_id),
		linewidth = 0.85
	) +
	geom_point(
		data = p3c_ind,
		aes(x = anchors, y = RMSE_improvement_percent, color = reference_id),
		size = 2.6
	) +
	geom_line(
		data = match41,
		aes(x = anchors, y = RMSE_improvement_percent),
		color = c_blue,
		linewidth = 1.6
	) +
	geom_point(
		data = match41,
		aes(x = anchors, y = RMSE_improvement_percent),
		color = c_blue,
		size = 3.8
	) +
	geom_text(
		data = match41,
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			label = paste0(round(RMSE_improvement_percent, 1), "%")
		),
		vjust = -1,
		color = c_blue,
		fontface = "bold",
		size = 3
	) +
	scale_color_manual(values = source_cols) +
	scale_x_continuous(breaks = 1:3) +
	labs(
		title = "Matched external sources suggest additional benefit from more anchors",
		subtitle = "Colored lines: individual sources; dark line: same three-source cohort",
		x = "Number of shared anchor compounds",
		y = "RMSE improvement (%)"
	) +
	theme_manuscript(9.8) +
	theme(
		legend.position = "bottom",
		legend.text = element_text(size = 7.8),
		plot.margin = margin(12, 16, 12, 16)
	)

fig3 <- p3a / (p3b | p3c) +
	plot_layout(heights = c(0.85, 1.15)) +
	plot_annotation(
		title = "Calibration transfers a HuSkin-trained model to independent SkinPiX measurements",
		subtitle = "Model fitting used HuSkinDB only; external evaluation used 54 SkinPiX-only compounds from four unseen sources",
		tag_levels = "A",
		theme = theme(
			plot.title = element_text(face = "bold", size = 13.5),
			plot.subtitle = element_text(size = 10),
			plot.margin = margin(14, 18, 10, 18)
		)
	)

save_fig(fig3, "Figure_3_skinpix_transfer", 10.5, 8.6)

# ============================================================
# Figure 4
# Source-specific falsification
# ============================================================

# ------------------------------------------------------------
# A. Matched-source versus different-source calibration
# ------------------------------------------------------------

p4a_same <- data.frame(
	anchors = wrong60$anchors,
	calibration = "Matched-source calibration",
	improvement = wrong60$correct_improvement_percent,
	ci_low = wrong60$correct_ci_low,
	ci_high = wrong60$correct_ci_high
)

p4a_wrong <- data.frame(
	anchors = wrong60$anchors,
	calibration = "Different-source calibration",
	improvement = wrong60$wrong_improvement_percent,
	ci_low = wrong60$wrong_ci_low,
	ci_high = wrong60$wrong_ci_high
)

p4a_dat <- rbind(p4a_same, p4a_wrong)

p4a_dat$calibration <- factor(
	p4a_dat$calibration,
	levels = c(
		"Matched-source calibration",
		"Different-source calibration"
	)
)

plot4a <- ggplot(
	p4a_dat,
	aes(x = anchors, y = improvement, color = calibration)
) +
	geom_hline(yintercept = 0, linetype = 3, color = c_grey) +
	geom_errorbar(
		aes(ymin = ci_low, ymax = ci_high),
		width = 0.12,
		linewidth = 0.7
	) +
	geom_line(linewidth = 1.25) +
	geom_point(size = 3.7) +
	geom_text(
		aes(label = paste0(round(improvement, 1), "%")),
		vjust = ifelse(p4a_dat$improvement >= 0, -1.25, 1.55),
		fontface = "bold",
		size = 3,
		show.legend = FALSE
	) +
	scale_color_manual(
		values = c(
			"Matched-source calibration" = c_blue,
			"Different-source calibration" = c_red
		)
	) +
	scale_x_continuous(breaks = 1:5) +
	labs(
		title = "Matched-source and different-source calibration diverge",
		subtitle = "Bars show 95% source-bootstrap intervals",
		x = "Number of anchor compounds",
		y = "RMSE improvement vs no calibration (%)"
	) +
	theme_manuscript(10) +
	theme(
		legend.position = "bottom",
		legend.text = element_text(size = 8.2),
		plot.margin = margin(12, 20, 12, 20)
	)

# ------------------------------------------------------------
# B. Direct paired falsification contrast
# ------------------------------------------------------------

plot4b <- ggplot(
	wrong60,
	aes(x = anchors, y = RMSE_wrong_minus_correct)
) +
	geom_hline(yintercept = 0, linetype = 3, color = c_grey) +
	geom_errorbar(
		aes(ymin = difference_ci_low, ymax = difference_ci_high),
		width = 0.13,
		linewidth = 0.75,
		color = c_grey_dark
	) +
	geom_line(
		color = c_teal,
		linewidth = 1.3
	) +
	geom_point(
		color = c_teal,
		size = 3.8
	) +
	geom_text(
		aes(label = sprintf("%.3f", RMSE_wrong_minus_correct)),
		vjust = -1.2,
		size = 3,
		fontface = "bold",
		color = c_blue
	) +
	scale_x_continuous(breaks = 1:5) +
	labs(
		title = "Correct source information outperforms a borrowed offset",
		subtitle = "Positive values favor matched-source calibration; all 95% intervals remain above zero",
		x = "Number of anchor compounds",
		y = "RMSE: different-source minus matched-source"
	) +
	theme_manuscript(10) +
	theme(
		plot.margin = margin(12, 20, 12, 20)
	)

fig4 <- (plot4a | plot4b) +
	plot_annotation(
		title = "The calibration signal is specific to the experimental source",
		subtitle = "Replacing the held-out source offset with information from another source removes or reverses the predictive benefit",
		tag_levels = "A",
		theme = theme(
			plot.title = element_text(face = "bold", size = 13.5),
			plot.subtitle = element_text(size = 10),
			plot.margin = margin(14, 18, 10, 18)
		)
	)

save_fig(fig4, "Figure_4_source_specific_falsification", 10.5, 5.5)

# ============================================================
# Figure S1
# Applicability-domain diagnostics
# ============================================================

ad_ch_s1 <- ad21[ad21$specification == "S1", ]
ad_ch_s1$group <- "Cheruvu et al., 2022 — S1"
ad_ch_s1$ratio <- ad_ch_s1$mahalanobis / qchisq(0.95, 2)

ad_ch_s2 <- ad21[ad21$specification == "S2", ]
ad_ch_s2$group <- "Cheruvu et al., 2022 — S2"
ad_ch_s2$ratio <- ad_ch_s2$mahalanobis / qchisq(0.95, 6)

ad_hu <- ad31
ad_hu$group <- "HuSkinDB"
ad_hu$ratio <- ad_hu$mahalanobis / qchisq(0.95, 2)

ad_sp <- domain41
ad_sp$group <- "SkinPiX external transfer"
ad_sp$ratio <- ad_sp$mahalanobis / qchisq(0.95, 2)

ad_plot <- rbind(
	ad_ch_s1[c("group", "compound_id", "ratio", "in_domain")],
	ad_ch_s2[c("group", "compound_id", "ratio", "in_domain")],
	ad_hu[c("group", "compound_id", "ratio", "in_domain")],
	ad_sp[c("group", "compound_id", "ratio", "in_domain")]
)

ad_plot$ratio <- pmax(ad_plot$ratio, 0.001)

ad_plot$group <- factor(
	ad_plot$group,
	levels = c(
		"SkinPiX external transfer",
		"HuSkinDB",
		"Cheruvu et al., 2022 — S2",
		"Cheruvu et al., 2022 — S1"
	)
)

counts <- tapply(
	ad_plot$in_domain,
	ad_plot$group,
	function(x) paste0(sum(x), "/", length(x), " in-domain")
)

group_labels <- setNames(
	paste0(names(counts), "    ", counts),
	names(counts)
)

figs1 <- ggplot(ad_plot, aes(x = ratio, y = group)) +
	geom_vline(
		xintercept = 1,
		linetype = 2,
		color = c_red,
		linewidth = 0.9
	) +
	geom_jitter(
		aes(fill = in_domain),
		shape = 21,
		height = 0.13,
		width = 0,
		size = 2.2,
		alpha = 0.75,
		color = "white",
		stroke = 0.2
	) +
	scale_fill_manual(
		values = c(
			"TRUE" = c_teal,
			"FALSE" = c_red
		),
		labels = c(
			"TRUE" = "In-domain",
			"FALSE" = "Out-of-domain"
		)
	) +
	scale_x_log10(
		breaks = c(0.01, 0.1, 1, 10, 100),
		labels = c("0.01", "0.1", "1", "10", "100")
	) +
	scale_y_discrete(labels = group_labels) +
	labs(
		title = "Chemical applicability-domain diagnostics",
		subtitle = "Distance is normalized to the fold-specific 95% Mahalanobis threshold; values above 1 exceed the distance criterion",
		x = "Mahalanobis distance / 95% threshold",
		y = NULL
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "top",
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(18, 34, 16, 34)
	)

save_fig(figs1, "Figure_S1_applicability_domain", 9.2, 5.1)

# ============================================================
# Figure S2
# HuSkin matched-source calibration
# ============================================================

s2_ind <- matchref31
s2_end <- s2_ind[s2_ind$anchors == max(s2_ind$anchors), ]

figs2 <- ggplot() +
	geom_hline(yintercept = 0, linetype = 3, color = c_grey) +
	geom_ribbon(
		data = cal60,
		aes(
			x = anchors,
			ymin = improvement_ci_low,
			ymax = improvement_ci_high
		),
		fill = c_blue,
		alpha = 0.12
	) +
	geom_line(
		data = s2_ind[s2_ind$reference_id != "HUR_48", ],
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			group = reference_id,
			color = "Individual held-out source"
		),
		linewidth = 0.65,
		alpha = 0.6
	) +
	geom_point(
		data = s2_ind[s2_ind$reference_id != "HUR_48", ],
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			color = "Individual held-out source"
		),
		size = 1.8,
		alpha = 0.65
	) +
	geom_line(
		data = s2_ind[s2_ind$reference_id == "HUR_48", ],
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			color = "HUR_48"
		),
		linewidth = 1
	) +
	geom_point(
		data = s2_ind[s2_ind$reference_id == "HUR_48", ],
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			color = "HUR_48"
		),
		size = 2.4,
		show.legend = FALSE
	) +
	geom_line(
		data = cal60,
		aes(
			x = anchors,
			y = improvement_percent,
			color = "Mean source-level RMSE"
		),
		linewidth = 1.7
	) +
	geom_point(
		data = cal60,
		aes(
			x = anchors,
			y = improvement_percent,
			color = "Mean source-level RMSE"
		),
		size = 3.8
	) +
	geom_text(
		data = cal60,
		aes(
			x = anchors,
			y = improvement_percent,
			label = paste0(round(improvement_percent, 1), "%")
		),
		vjust = -1,
		color = c_blue,
		fontface = "bold",
		size = 3
	) +
	geom_text(
		data = s2_end,
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			label = reference_id,
			color = ifelse(reference_id == "HUR_48", "HUR_48", "Individual held-out source")
		),
		hjust = -0.12,
		size = 2.7,
		show.legend = FALSE
	) +
	scale_color_manual(
		values = c(
			"Individual held-out source" = c_grey_dark,
			"Mean source-level RMSE" = c_blue,
			"HUR_48" = c_red
		),
		breaks = c(
			"Mean source-level RMSE",
			"Individual held-out source",
			"HUR_48"
		)
	) +
	scale_x_continuous(
		breaks = 1:5,
		limits = c(0.9, 5.8)
	) +
	labs(
		title = "Calibration benefit across the nine matched HuSkin sources",
		subtitle = "Dark line shows reduction in mean source-level RMSE; shading is the 95% source-bootstrap interval",
		x = "Number of anchor compounds",
		y = "RMSE improvement (%)"
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "top",
		plot.margin = margin(18, 48, 14, 20)
	)

save_fig(figs2, "Figure_S2_huskin_source_calibration", 9.2, 5.8)

cat(
	"\nManuscript figures created\n",
	"Figure 1: variance decomposition\n",
	"Figure 2: source-information prediction with compound-bootstrap intervals\n",
	"Figure 3: SkinPiX external transfer\n",
	"Figure 4: source-specific falsification with source-bootstrap intervals\n",
	"Figure S1: applicability-domain diagnostics\n",
	"Figure S2: matched HuSkin calibration with source-bootstrap interval\n",
	sep = ""
)