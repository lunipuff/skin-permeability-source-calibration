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
		plot.subtitle = element_text(size = base_size - 0.5, margin = margin(b = 10)),
		axis.title = element_text(size = base_size),
		axis.text = element_text(size = base_size - 0.5),
		strip.background = element_blank(),
		strip.text = element_text(face = "bold", size = base_size),
		legend.title = element_blank(),
		legend.key = element_blank(),
		legend.text = element_text(size = base_size - 1),
		plot.tag = element_text(face = "bold", size = base_size + 2),
		plot.margin = margin(12, 18, 12, 18)
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

var10 <- read.csv("R/results/10_variance_components.csv", stringsAsFactors = FALSE)
rep31 <- read.csv("R/results/31_replication_summary.csv", stringsAsFactors = FALSE)

blind50 <- read.csv("R/results/50_source_blind_loco.csv", stringsAsFactors = FALSE)
wrong50 <- read.csv("R/results/50_wrong_source_calibration.csv", stringsAsFactors = FALSE)
wrongref50 <- read.csv("R/results/50_wrong_source_reference.csv", stringsAsFactors = FALSE)

match31 <- read.csv("R/results/31_calibration_matched_summary.csv", stringsAsFactors = FALSE)
matchref31 <- read.csv("R/results/31_calibration_matched_reference.csv", stringsAsFactors = FALSE)

eligible41 <- read.csv("R/results/41_skinpix_eligible_references.csv", stringsAsFactors = FALSE)
match41 <- read.csv("R/results/41_skinpix_calibration_matched_3anchor.csv", stringsAsFactors = FALSE)
ref41 <- read.csv("R/results/41_skinpix_calibration_reference.csv", stringsAsFactors = FALSE)
pred41 <- read.csv("R/results/41_skinpix_one_anchor_predictions.csv", stringsAsFactors = FALSE)
domain41 <- read.csv("R/results/41_skinpix_transfer_domain.csv", stringsAsFactors = FALSE)

ad21 <- read.csv("R/results/21_applicability_domain.csv", stringsAsFactors = FALSE)
ad31 <- read.csv("R/results/31_applicability_domain.csv", stringsAsFactors = FALSE)

# ============================================================
# Figure 1
# Variance decomposition
# ============================================================

fig1_dat <- data.frame(
	dataset = rep(c("Cheruvu et al., 2022", "HuSkinDB"), each = 3),
	component = rep(c("Compound", "Source", "Residual"), 2),
	percent = c(
		var10$percent[var10$component == "Compound"],
		var10$percent[var10$component == "Reference"],
		var10$percent[var10$component == "Residual"],
		rep31$compound_variance_percent,
		rep31$reference_variance_percent,
		rep31$residual_variance_percent
	)
)

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
	coord_cartesian(xlim = c(0, 100), clip = "off"
	) +
	labs(
		title = "Source-associated variation in human skin permeability measurements",
		x = "Variance component (%)",
		y = NULL
	) +
	theme_manuscript(12) +
	theme(
		legend.position = "top",
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(15, 24, 12, 24)
	)

save_fig(fig1, "Figure_1_variance_decomposition", 9, 3.7)

# ============================================================
# Figure 2
# Source-blind versus known-source prediction
# ============================================================

blind50$group <- paste(blind50$dataset, blind50$specification)

blind50$group[blind50$group == "Cheruvu S1"] <- "Cheruvu, 2022 — S1"
blind50$group[blind50$group == "Cheruvu S2"] <- "Cheruvu, 2022 — S2"

blind50$stage <- factor(
	blind50$approach,
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

blind50$rmse_label <- sprintf("%.3f", blind50$RMSE)

blind50$gain_label <- ifelse(
	blind50$approach == "conventional",
	"",
	ifelse(
		blind50$improvement_vs_conventional >= 0,
		paste0("−", round(blind50$improvement_vs_conventional, 1), "% RMSE"),
		paste0("+", abs(round(blind50$improvement_vs_conventional, 1)), "% RMSE")
	)
)

ylim2 <- range(blind50$RMSE)
ylim2 <- c(ylim2[1] - 0.08, ylim2[2] + 0.08)

fig2 <- ggplot(blind50, aes(x = stage, y = RMSE, group = group)) +
	geom_line(linewidth = 0.9, color = c_grey) +
	geom_point(aes(fill = stage), shape = 21, size = 4.4, stroke = 0.6, color = "white") +
	geom_text(
		aes(label = rmse_label),
		vjust = -1.25,
		size = 3.4
	) +
	geom_text(
		aes(label = gain_label),
		vjust = 1.8,
		size = 3.05,
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
		title = "Source context accounts for most of the predictive improvement",
		subtitle = "All panels use the same RMSE scale; percentages are relative to conventional compound aggregation",
		x = NULL,
		y = "Leave-one-compound-out RMSE"
	) +
	theme_manuscript(11) +
	theme(
		legend.position = "none",
		panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.5),
		plot.margin = margin(15, 22, 22, 22),
		axis.text.x = element_text(angle = 45, hjust = 1, size = 2)
	)

save_fig(fig2, "Figure_2_source_information_prediction", 10, 4.8)

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
		size = 4.3,
		fill = "white",
		color = c_grey_dark,
		stroke = 1.1
	) +
	geom_point(
		aes(x = RMSE_calibrated),
		shape = 21,
		size = 4.5,
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
		size = 3.4,
		fontface = "bold",
		color = c_blue
	) +
	scale_x_continuous(expand = expansion(mult = c(0.03, 0.18))) +
	labs(
		title = "One shared compound calibrates a completely unseen source",
		subtitle = "Open circle = HuSkin prediction; filled circle = after one-anchor source calibration",
		x = "RMSE on SkinPiX-only compounds",
		y = NULL
	) +
	theme_manuscript(10.5) +
	theme(
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(10, 24, 8, 18)
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
		size = 3,
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
			"TRUE" = "Lower error after calibration",
			"FALSE" = "Higher error after calibration"
		)
	) +
	coord_equal(xlim = c(0, lim3b), ylim = c(0, lim3b)) +
	labs(
		title = "Calibration reduces error for most external compounds",
		subtitle = paste0(n_improved, "/", n_total, " SkinPiX-only compounds move below the identity line"),
		x = "Absolute error before calibration",
		y = "Absolute error after calibration"
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "bottom",
		legend.direction = "horizontal"
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
		linewidth = 0.9
	) +
	geom_point(
		data = p3c_ind,
		aes(x = anchors, y = RMSE_improvement_percent, color = reference_id),
		size = 2.8
	) +
	geom_line(
		data = match41,
		aes(x = anchors, y = RMSE_improvement_percent),
		color = c_blue,
		linewidth = 1.7
	) +
	geom_point(
		data = match41,
		aes(x = anchors, y = RMSE_improvement_percent),
		color = c_blue,
		size = 4
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
		size = 3.2
	) +
	scale_color_manual(values = source_cols) +
	scale_x_continuous(breaks = 1:3) +
	labs(
		title = "More shared compounds provide additional calibration information",
		subtitle = "Colored lines are the same three external sources; dark blue is their aggregate RMSE reduction",
		x = "Number of shared anchor compounds",
		y = "RMSE improvement (%)"
	) +
	theme_manuscript(10.5) +
	theme(legend.position = "bottom")

fig3 <- p3a / (p3b | p3c) +
	plot_layout(heights = c(0.9, 1.15)) +
	plot_annotation(
		title = "Calibration transfers a HuSkin-trained model to independent SkinPiX measurements",
		subtitle = "HuSkinDB was used for model fitting only; evaluation used 54 SkinPiX-only compounds from four unseen sources",
		tag_levels = "A",
		theme = theme(
			plot.title = element_text(face = "bold", size = 14),
			plot.subtitle = element_text(size = 10.5),
			plot.margin = margin(12, 16, 8, 16)
		)
	)

save_fig(fig3, "Figure_3_skinpix_transfer", 10.5, 8.4)

# ============================================================
# Figure 4
# Source-specific falsification
# ============================================================

# ------------------------------------------------------------
# A. Aggregate prediction effect
# ------------------------------------------------------------

p4a <- data.frame(
	anchors = rep(wrong50$anchors, 2),
	calibration = rep(
		c(
			"Anchors from the same held-out source",
			"Offset borrowed from a different source"
		),
		each = nrow(wrong50)
	),
	improvement = c(
		wrong50$correct_RMSE_improvement,
		wrong50$wrong_RMSE_improvement
	)
)

p4a$calibration <- factor(
	p4a$calibration,
	levels = c(
		"Anchors from the same held-out source",
		"Offset borrowed from a different source"
	)
)

plot4a <- ggplot(p4a, aes(x = anchors, y = improvement, color = calibration)) +
	geom_hline(yintercept = 0, linetype = 3, color = c_grey) +
	geom_line(linewidth = 1.5) +
	geom_point(size = 4) +
	geom_text(
		aes(label = paste0(round(improvement, 1), "%")),
		vjust = ifelse(p4a$improvement >= 0, -1, 1.5),
		fontface = "bold",
		size = 3.2,
		show.legend = FALSE
	) +
	scale_color_manual(
		values = c(
			"Anchors from the same held-out source" = c_blue,
			"Offset borrowed from a different source" = c_red
		)
	) +
	scale_x_continuous(breaks = 1:5) +
	labs(
		title = "Only source-matched calibration improves aggregate prediction",
		x = "Number of anchor compounds",
		y = "RMSE improvement vs no calibration (%)"
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "bottom",
		legend.text = element_text(size = 8.5)
	)

# ------------------------------------------------------------
# B. How consistently does each procedure help?
# ------------------------------------------------------------

prop_same <- aggregate(
	correct_RMSE_improvement ~ anchors,
	data = wrongref50,
	FUN = function(x) 100 * mean(x > 0)
)

prop_other <- aggregate(
	wrong_RMSE_improvement ~ anchors,
	data = wrongref50,
	FUN = function(x) 100 * mean(x > 0)
)

names(prop_same)[2] <- "percent"
names(prop_other)[2] <- "percent"

prop_same$calibration <- "Anchors from the same held-out source"
prop_other$calibration <- "Offset borrowed from a different source"

p4b <- rbind(prop_same, prop_other)

p4b$calibration <- factor(
	p4b$calibration,
	levels = c(
		"Anchors from the same held-out source",
		"Offset borrowed from a different source"
	)
)

plot4b <- ggplot(p4b, aes(x = factor(anchors), y = percent, fill = calibration)) +
	geom_col(
		position = position_dodge(width = 0.72),
		width = 0.64
	) +
	geom_text(
		aes(label = paste0(round(percent), "%")),
		position = position_dodge(width = 0.72),
		vjust = -0.4,
		size = 3.1,
		fontface = "bold"
	) +
	scale_fill_manual(
		values = c(
			"Anchors from the same held-out source" = c_blue,
			"Offset borrowed from a different source" = c_red
		)
	) +
	scale_y_continuous(
		limits = c(0, 108),
		breaks = seq(0, 100, 20),
		expand = c(0, 0)
	) +
	labs(
		title = "Source-matched offsets help across held-out references",
		subtitle = "Percentage of references with lower RMSE than the uncalibrated prediction",
		x = "Number of anchor compounds",
		y = "References improved (%)"
	) +
	theme_manuscript(10.5) +
	theme(legend.position = "none")

fig4 <- plot4a | plot4b +
	plot_annotation(
		title = "The calibration signal is specific to the experimental source",
		subtitle = "Replacing the held-out source offset with an offset estimated from another source reverses the benefit",
		tag_levels = "A",
		theme = theme(
			plot.title = element_text(face = "bold", size = 14),
			plot.subtitle = element_text(size = 10.5),
			plot.margin = margin(12, 16, 8, 16)
		)
	)

save_fig(fig4, "Figure_4_source_specific_falsification", 10.5, 5.4)

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

counts <- tapply(ad_plot$in_domain, ad_plot$group, function(x) paste0(sum(x), "/", length(x), " in-domain"))

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
		size = 2.25,
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
		subtitle = "Distance is normalized to the fold-specific 95% Mahalanobis threshold; values >1 exceed that distance criterion",
		x = "Mahalanobis distance / 95% threshold",
		y = NULL
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "top",
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(16, 30, 14, 30)
	)

save_fig(figs1, "Figure_S1_applicability_domain", 9.2, 5.2)

# ============================================================
# Figure S2
# HuSkin matched-source calibration
# ============================================================

s2_ind <- matchref31
s2_ind$line_type <- ifelse(s2_ind$reference_id == "HUR_48", "HUR_48", "Individual held-out source")

s2_end <- s2_ind[s2_ind$anchors == max(s2_ind$anchors), ]

figs2 <- ggplot() +
	geom_hline(yintercept = 0, linetype = 3, color = c_grey) +
	geom_line(
		data = s2_ind[s2_ind$reference_id != "HUR_48", ],
		aes(x = anchors, y = RMSE_improvement_percent, group = reference_id, color = "Individual held-out source"),
		linewidth = 0.7,
		alpha = 0.65
	) +
	geom_point(
		data = s2_ind[s2_ind$reference_id != "HUR_48", ],
		aes(x = anchors, y = RMSE_improvement_percent, color = "Individual held-out source"),
		size = 1.9,
		alpha = 0.7
	) +
	geom_line(
		data = s2_ind[s2_ind$reference_id == "HUR_48", ],
		aes(x = anchors, y = RMSE_improvement_percent, color = "HUR_48"),
		linewidth = 1,
		show.legend = TRUE
	) +
	geom_point(
		data = s2_ind[s2_ind$reference_id == "HUR_48", ],
		aes(x = anchors, y = RMSE_improvement_percent, color = "HUR_48"),
		size = 2.4,
		show.legend = FALSE
	) +
	geom_line(
		data = match31,
		aes(x = anchors, y = RMSE_improvement_percent, color = "Aggregate across all 9 sources"),
		linewidth = 1.8
	) +
	geom_point(
		data = match31,
		aes(x = anchors, y = RMSE_improvement_percent, color = "Aggregate across all 9 sources"),
		size = 4
	) +
	geom_text(
		data = match31,
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			label = paste0(round(RMSE_improvement_percent, 1), "%")
		),
		vjust = -1,
		color = c_blue,
		fontface = "bold",
		size = 3.1
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
		size = 2.8,
		show.legend = FALSE
	) +
	scale_color_manual(
		values = c(
			"Individual held-out source" = c_grey_dark,
			"Aggregate across all 9 sources" = c_blue,
			"HUR_48" = c_red
		),
		breaks = c(
			"Aggregate across all 9 sources",
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
		subtitle = "Dark blue is the cohort-level RMSE reduction across all nine sources; it is not the arithmetic mean of the individual percentages",
		x = "Number of anchor compounds",
		y = "Source-level RMSE improvement (%)"
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "top",
		plot.margin = margin(15, 42, 12, 18)
	)

save_fig(figs2, "Figure_S2_huskin_source_calibration", 9.2, 5.8)

cat(
	"\nManuscript figures created\n",
	"Figure 1: variance decomposition\n",
	"Figure 2: source-information prediction\n",
	"Figure 3: SkinPiX external transfer\n",
	"Figure 4: source-specific falsification\n",
	"Figure S1: applicability-domain diagnostics\n",
	"Figure S2: matched HuSkin calibration heterogeneity\n",
	sep = ""
)