# ============================================================
# 99b_make_figures.R
# ============================================================

library(ggplot2)
library(patchwork)

source("R/00_config.R")

cheruvu <- read.csv(
	"data/processed/01_cheruvu.csv",
	stringsAsFactors = FALSE
)

huskin <- read.csv(
	"data/processed/30_huskindb.csv",
	stringsAsFactors = FALSE
)

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

loco60 <- read.csv("R/results/60_loco_bootstrap_ci.csv", stringsAsFactors = FALSE)
var60 <- read.csv("R/results/60_variance_component_ci.csv", stringsAsFactors = FALSE)
cal60 <- read.csv("R/results/60_huskin_calibration_ci.csv", stringsAsFactors = FALSE)
wrong60 <- read.csv("R/results/60_wrong_source_ci.csv", stringsAsFactors = FALSE)

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
# Within-compound heterogeneity and variance decomposition
# ============================================================

# ------------------------------------------------------------
# Panel A. Centered repeated measurements
# ------------------------------------------------------------

prep_centered <- function(x, dataset_label) {
	x <- x[, c("compound_id", "reference_id", "log_kp")]
	x <- x[
		!is.na(x$compound_id) &
		!is.na(x$reference_id) &
		!is.na(x$log_kp),
	]

	n_measurements <- table(x$compound_id)
	repeated <- names(n_measurements[n_measurements >= 2])

	x <- x[as.character(x$compound_id) %in% repeated, ]

	x$compound_mean <- ave(
		x$log_kp,
		x$compound_id,
		FUN = mean
	)

	x$centered_log_kp <- x$log_kp - x$compound_mean

	spread <- aggregate(
		centered_log_kp ~ compound_id,
		data = x,
		FUN = function(z) diff(range(z))
	)

	spread <- spread[
		order(spread$centered_log_kp),
	]

	spread$compound_order <- seq_len(nrow(spread))

	x$compound_order <- spread$compound_order[
		match(x$compound_id, spread$compound_id)
	]

	source_key <- unique(
	x[, c("compound_id", "reference_id")]
)

source_key <- source_key[
	order(source_key$compound_id, source_key$reference_id),
]

source_key$source_color_index <- NA_integer_

compound_groups <- split(
	seq_len(nrow(source_key)),
	source_key$compound_id
)

for (idx in compound_groups) {
	n_ref <- length(idx)

	if (n_ref == 1) {
		color_index <- 1
	} else {
		color_index <- round(
			seq(
				1,
				100,
				length.out = n_ref
			)
		)
	}

	source_key$source_color_index[idx] <- color_index
}

x$source_color_index <- source_key$source_color_index[
	match(
		paste(x$compound_id, x$reference_id),
		paste(source_key$compound_id, source_key$reference_id)
	)
]

	ranges <- aggregate(
		centered_log_kp ~ compound_order,
		data = x,
		FUN = function(z) c(min = min(z), max = max(z))
	)

	ranges <- data.frame(
		compound_order = ranges$compound_order,
		ymin = ranges$centered_log_kp[, "min"],
		ymax = ranges$centered_log_kp[, "max"]
	)

	x$dataset <- dataset_label
	ranges$dataset <- dataset_label

	list(
		points = x,
		ranges = ranges
	)
}

centered_ch <- prep_centered(
	cheruvu,
	"Cheruvu"
)

centered_hu <- prep_centered(
	huskin,
	"HuSkinDB"
)

fig1a_points <- rbind(
	centered_ch$points,
	centered_hu$points
)

fig1a_ranges <- rbind(
	centered_ch$ranges,
	centered_hu$ranges
)

fig1a_points$dataset <- factor(
	fig1a_points$dataset,
	levels = c("Cheruvu", "HuSkinDB"),
	labels = c(
		"Cheruvu (n = 81 repeated compounds)",
		"HuSkinDB (n = 89 repeated compounds)"
	)
)

fig1a_ranges$dataset <- factor(
	fig1a_ranges$dataset,
	levels = c("Cheruvu", "HuSkinDB"),
	labels = c(
		"Cheruvu (n = 81 repeated compounds)",
		"HuSkinDB (n = 89 repeated compounds)"
	)
)

fig1a <- ggplot() +
	geom_hline(
		yintercept = 0,
		linetype = "dashed",
		linewidth = 0.4,
		color = "grey55"
	) +
	geom_linerange(
		data = fig1a_ranges,
		aes(
			x = compound_order,
			ymin = ymin,
			ymax = ymax
		),
		color = "grey65",
		linewidth = 0.4
	) +
	geom_point(
		data = fig1a_points,
		aes(
			x = compound_order,
			y = centered_log_kp,
			color = factor(source_color_index)
		),
		size = 1.8,
		alpha = 0.85,
		position = position_jitter(
			width = 0.12,
			height = 0
		)
	) +
	ylim(-4,4) +
	facet_wrap(
		~ dataset,
		ncol = 1,
		scales = "free_x"
	) +
	scale_color_manual(
		values = setNames(
			palette(100),
			as.character(seq_len(100))
		)
	) +
	annotate(
		"text",
		x = Inf,
		y = Inf,
		label = "Colors distinguish sources within each compound",
		hjust = 1.4,
		vjust = 1.4,
		size = 3.2,
		color = "grey35"
	) +
	labs(
		title = "Within-compound permeability measurements",
		x = "Repeated compounds, ordered by within-compound spread",
		y = expression(Observed~logK[p] - Mean~logK[p])
	) +
	theme_manuscript(12) +
	theme(
		legend.position = "none",
		axis.title.y = element_text(vjust = -10),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank(),
		plot.title = element_text(hjust = 0),
		panel.grid.major.x = element_blank(),
		panel.grid.minor.x = element_blank(),
		strip.background = element_blank(),
		strip.text = element_text(face = "bold", hjust = 0),
		plot.margin = margin(15, 24, 12, 24)
	)

# ------------------------------------------------------------
# Panel B. Variance decomposition
# ------------------------------------------------------------

fig1b_dat <- var60[, c(
	"dataset",
	"component",
	"estimate_percent"
)]

fig1b_dat$dataset <- factor(
	fig1b_dat$dataset,
	levels = c("HuSkinDB", "Cheruvu")
)

fig1b_dat$component <- factor(
	fig1b_dat$component,
	levels = c("Residual", "Source", "Compound")
)

fig1b_dat$label <- paste0(
	round(fig1b_dat$estimate_percent, 1),
	"%"
)

fig1b_dat$text_color <- ifelse(
	fig1b_dat$component == "Compound",
	"white",
	"black"
)

fig1b <- ggplot(
	fig1b_dat,
	aes(
		x = estimate_percent,
		y = dataset,
		fill = component
	)
) +
	geom_col(
		width = 0.6
	) +
	geom_text(
		aes(
			label = label,
			color = text_color
		),
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
		breaks = c(
			"Compound",
			"Source",
			"Residual"
		)
	) +
	scale_x_continuous(
		breaks = seq(0, 100, 20),
		expand = c(0, 0)
	) +
	coord_cartesian(
		xlim = c(0, 100),
		clip = "off"
	) +
	labs(
	title = "Variance contribution by component",
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

# ------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------

fig1 <- fig1a / fig1b +
	plot_layout(
		heights = c(2.0, 1)
	) +
	plot_annotation(
		tag_levels = "A"
	) &
	theme(
		plot.tag = element_text(
			face = "bold",
			size = 20,
			hjust = 3
		)
	)

save_fig(
	fig1,
	"Figure_1_variance_decomposition",
	9,
	8.5
)

# ------------------------------------------------------------
# Figure 2. Source information and prediction
# ------------------------------------------------------------

fig2 <- loco60[loco60$approach %in% c("hierarchical_source_blind", "hierarchical_known_source") & (loco60$dataset == "Cheruvu" | (loco60$dataset == "HuSkinDB" & loco60$specification == "S1")), ]
fig2$group <- paste(fig2$dataset, fig2$specification)
fig2_order <- c("Cheruvu S1", "Cheruvu S2", "HuSkinDB S1")
stopifnot(nrow(fig2) == 6L, all(table(fig2$group) == 2L))

fig2$blind <- fig2$approach == "hierarchical_source_blind"
fig2$y <- 4 - match(fig2$group, fig2_order) + ifelse(fig2$blind, 0.14, -0.14)
fig2$plot_col <- ifelse(fig2$blind, pal[7], pal[1])
fig2$plot_pch <- ifelse(fig2$blind, 17, 16)

fig2_ticks <- seq(floor(min(c(0, fig2$improvement_ci_low)) / 10) * 10, ceiling(max(c(0, fig2$improvement_ci_high)) / 10) * 10, by = 10)

pdf("R/figures/Figure_2_source_information_prediction.pdf", width = 7.6, height = 3.8)
par(mar = c(4.3, 6.7, 1.4, 0.8), mgp = c(2.2, 0.6, 0), tcl = -0.25)

plot(NA_real_, NA_real_, type = "n", xlim = range(fig2_ticks) + c(-1, 1), ylim = c(0.6, 3.8), axes = FALSE, xlab = "RMSE reduction relative to conventional QSPR (%)", ylab = "")
abline(v = fig2_ticks, col = "grey92")
abline(h = c(1.5, 2.5), col = "grey94")
abline(v = 0, col = "grey45", lwd = 1.2)

segments(fig2$improvement_ci_low, fig2$y, fig2$improvement_ci_high, fig2$y, col = fig2$plot_col, lwd = 2)
segments(fig2$improvement_ci_low, fig2$y - 0.055, fig2$improvement_ci_low, fig2$y + 0.055, col = fig2$plot_col, lwd = 1.5)
segments(fig2$improvement_ci_high, fig2$y - 0.055, fig2$improvement_ci_high, fig2$y + 0.055, col = fig2$plot_col, lwd = 1.5)
points(fig2$improvement_percent, fig2$y, pch = fig2$plot_pch, col = fig2$plot_col, cex = 1.35)

axis(1, at = fig2_ticks)
axis(2, at = c(3, 2, 1), labels = fig2_order, las = 1, tick = FALSE, cex.axis = 1.05)
legend("topright", legend = c("Source blind", "Source aware"), col = c(pal[7], pal[1]), pch = c(17, 16), horiz = TRUE, bty = "n", cex = 0.9, inset = c(0.01, 0.01))

dev.off()

# ============================================================
# Figure 3
# External calibration: HuSkinDB to SkinPiX
# ============================================================

# ------------------------------------------------------------
# A. One-anchor calibration by SkinPiX source
# ------------------------------------------------------------

p3a_dat <- ref41[ref41$anchors == 1, ]
p3a_dat <- p3a_dat[order(p3a_dat$RMSE_uncalibrated, decreasing = TRUE), ]
p3a_dat$source <- factor(p3a_dat$reference_id, levels = rev(p3a_dat$reference_id))
p3a_dat$change <- 100 * (p3a_dat$RMSE_calibrated / p3a_dat$RMSE_uncalibrated - 1)
p3a_dat$change_label <- paste0(sprintf("%.1f", abs(p3a_dat$change)), "% ", ifelse(p3a_dat$change < 0, "lower", "higher"))
p3a_dat$change_group <- ifelse(p3a_dat$change < 0, "Lower RMSE", "Higher RMSE")
p3a_dat$arrow_end <- p3a_dat$RMSE_calibrated + sign(p3a_dat$RMSE_uncalibrated - p3a_dat$RMSE_calibrated) * 0.015

p3a_points <- rbind(
	data.frame(source = p3a_dat$source, RMSE = p3a_dat$RMSE_uncalibrated, stage = "Before calibration"),
	data.frame(source = p3a_dat$source, RMSE = p3a_dat$RMSE_calibrated, stage = "After calibration")
)
p3a_points$stage <- factor(p3a_points$stage, levels = c("Before calibration", "After calibration"))

p3a <- ggplot(p3a_dat, aes(y = source)) +
	geom_segment(
		aes(x = RMSE_uncalibrated, xend = arrow_end, yend = source),
		color = c_grey,
		linewidth = 1.2,
		arrow = grid::arrow(type = "closed", length = grid::unit(0.09, "inches"))
	) +
	geom_point(
		data = p3a_points,
		aes(x = RMSE, fill = stage),
		shape = 21,
		size = 4.4,
		stroke = 1,
		color = c_grey_dark
	) +
	geom_text(
		aes(x = pmax(RMSE_uncalibrated, RMSE_calibrated), label = change_label, color = change_group),
		hjust = -0.18,
		size = 3.6,
		fontface = "bold"
	) +
	scale_fill_manual(values = c("Before calibration" = "white", "After calibration" = c_blue)) +
	scale_color_manual(values = c("Lower RMSE" = c_teal, "Higher RMSE" = c_red), guide = "none") +
	scale_x_continuous(breaks = seq(0.4, 1.2, 0.2)) +
	coord_cartesian(xlim = c(0.34, 1.32), clip = "off") +
	labs(
		title = "RMSE before and after one-anchor calibration",
		x = "RMSE on SkinPiX-only compounds",
		y = NULL
	) +
	theme_manuscript(12) +
	theme(
	    legend.position = c(0.8, 1.2),
    	legend.justification = c(0, 1),
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(8, 34, 6, 8)
	) +
	guides(fill = guide_legend(override.aes = list(shape = 21, size = 4, color = c_grey_dark)))

# ------------------------------------------------------------
# B. Compounds with reduced absolute error
# ------------------------------------------------------------

p3b_dat <- aggregate(
	cbind(abs_error_uncalibrated, abs_error_calibrated) ~ compound_id,
	data = pred41,
	mean
)

n_improved <- sum(p3b_dat$abs_error_calibrated < p3b_dat$abs_error_uncalibrated)
n_total <- nrow(p3b_dat)
n_other <- n_total - n_improved

p3b_summary <- data.frame(
	outcome = factor(
		c("Error decreased", "Error increased or unchanged"),
		levels = c("Error decreased", "Error increased or unchanged")
	),
	count = c(n_improved, n_other)
)
p3b_summary$percent <- 100 * p3b_summary$count / n_total
p3b_summary$label <- paste0(p3b_summary$count, "/", n_total, " (", round(p3b_summary$percent), "%)")

p3b <- ggplot(p3b_summary, aes(x = outcome, y = percent, fill = outcome)) +
	geom_col(width = 0.62) +
	geom_text(aes(label = label), vjust = -0.6, fontface = "bold", size = 3.8) +
	scale_fill_manual(values = c("Error decreased" = c_teal, "Error increased or unchanged" = c_red)) +
	scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), expand = c(0, 0)) +
	scale_x_discrete(labels = c("Error\ndecreased", "Error increased\nor unchanged")) +
	labs(
		title = "Change in compound-level error",
		x = NULL,
		y = "SkinPiX-only compounds (%)"
	) +
	theme_manuscript(11.5) +
	theme(
		legend.position = "none",
		axis.title.y = element_text(margin = margin(r = 5)),
		axis.line.x = element_blank(),
		axis.ticks.x = element_blank(),
		plot.margin = margin(8, 10, 8, 8)
	)

# ------------------------------------------------------------
# C. Anchor panels in the same three sources
# ------------------------------------------------------------

p3c_dat <- match41[match41$anchors %in% 1:3, ]

p3c <- ggplot(p3c_dat, aes(x = anchors, y = RMSE_improvement_percent)) +
	geom_hline(yintercept = 0, linetype = 3, color = c_grey) +
	geom_line(color = c_blue, linewidth = 1.4) +
	geom_point(color = c_blue, size = 3.8) +
	geom_text(
		aes(label = paste0(sprintf("%.1f", RMSE_improvement_percent), "%")),
		vjust = -1.1,
		fontface = "bold",
		color = c_blue,
		size = 3.8
	) +
	scale_x_continuous(breaks = 1:3) +
	scale_y_continuous(limits = c(0, 35), breaks = seq(0, 30, 10), expand = c(0, 0)) +
	coord_cartesian(xlim = c(0.9, 3.1), clip = "off") +
	labs(
		title = "Same three SkinPiX sources",
		x = "Number of anchor compounds",
		y = "Mean source-level RMSE reduction (%)"
	) +
	theme_manuscript(11.5) +
	theme(plot.margin = margin(8, 12, 8, 8),
		axis.title.y = element_text(margin = margin(r = 5)),
		plot.tag = element_text(margin = margin(r = 20)) 
		)

# ------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------

fig3 <- p3a / (free(p3b, type = "label") | p3c) +
	plot_layout(heights = c(1.05, 1), widths = c(1, 1)) +
	plot_annotation(tag_levels = "A") & # USE THE & OPERATOR HERE
	theme(
		plot.tag = element_text(size = 18, face = "bold"), # Changes A, B, C size globally
		plot.margin = margin(8, 10, 8, 10)
	)

save_fig(fig3, "Figure_3_skinpix_transfer", 10, 8)

# ============================================================
# Figure 4
# Source-specific calibration control
# ============================================================

wrong60 <- read.csv("R/results/60_wrong_source_ci.csv", stringsAsFactors = FALSE)
wrong60 <- wrong60[order(wrong60$anchors), ]

rows <- seq_len(nrow(wrong60))
x_limits <- range(c(wrong60$correct_ci_low, wrong60$correct_ci_high, wrong60$wrong_ci_low, wrong60$wrong_ci_high, 0))
x_padding <- diff(x_limits) * 0.07
x_limits <- x_limits + c(-x_padding, x_padding)

pdf("R/figures/Figure_4_source_specific_falsification.pdf", width = 8.2, height = 4.7)
par(mar = c(4.5, 8, 1.2, 1), cex.axis = 1, cex.lab = 1.1)

plot(0, 0, type = "n", xlim = x_limits, ylim = c(0.5, nrow(wrong60) + 0.5),
	axes = FALSE, xlab = "RMSE reduction relative to no calibration (%)", ylab = "")
abline(v = pretty(x_limits), col = "grey92")
abline(v = 0, col = "grey55", lwd = 1.3)

segments(wrong60$correct_ci_low, rows, wrong60$correct_ci_high, rows, col = pal[1], lwd = 1.7)
segments(wrong60$correct_ci_low, rows - 0.09, wrong60$correct_ci_low, rows + 0.09, col = pal[1], lwd = 1.7)
segments(wrong60$correct_ci_high, rows - 0.09, wrong60$correct_ci_high, rows + 0.09, col = pal[1], lwd = 1.7)
points(wrong60$correct_improvement_percent, rows, pch = 16, col = pal[1], cex = 1.35)

segments(wrong60$wrong_ci_low, rows, wrong60$wrong_ci_high, rows, col = pal[7], lwd = 1.7)
segments(wrong60$wrong_ci_low, rows - 0.09, wrong60$wrong_ci_low, rows + 0.09, col = pal[7], lwd = 1.7)
segments(wrong60$wrong_ci_high, rows - 0.09, wrong60$wrong_ci_high, rows + 0.09, col = pal[7], lwd = 1.7)
points(wrong60$wrong_improvement_percent, rows, pch = 17, col = pal[7], cex = 1.35)

axis(1, at = pretty(x_limits))
axis(2, at = rows,
	labels = paste0(wrong60$anchors, ifelse(wrong60$anchors == 1, " anchor", " anchors"),
		"  (n=", wrong60$n_references, ")"),
	las = 1, tick = FALSE)

legend("bottomright",
	legend = c("Target-source anchors", "Other-source anchors"),
	col = c(pal[1], pal[7]), pch = c(16, 17),
	bty = "n", cex = 0.9)

dev.off()

# ============================================================
# Figure S1
# Applicability-domain diagnostics
# ============================================================

ad_ch_s1 <- ad21[
	ad21$specification == "S1",
]

ad_ch_s1$group <- "Cheruvu et al., 2022 - S1"
ad_ch_s1$ratio <- ad_ch_s1$mahalanobis / qchisq(0.95, 2)

ad_ch_s2 <- ad21[
	ad21$specification == "S2",
]

ad_ch_s2$group <- "Cheruvu et al., 2022 - S2"
ad_ch_s2$ratio <- ad_ch_s2$mahalanobis / qchisq(0.95, 6)

ad_hu <- ad31
ad_hu$group <- "HuSkinDB"
ad_hu$ratio <- ad_hu$mahalanobis / qchisq(0.95, 2)

ad_sp <- domain41
ad_sp$group <- "SkinPiX transfer"
ad_sp$ratio <- ad_sp$mahalanobis / qchisq(0.95, 2)

ad_plot <- rbind(
	ad_ch_s1[c(
		"group",
		"compound_id",
		"ratio",
		"in_domain"
	)],
	ad_ch_s2[c(
		"group",
		"compound_id",
		"ratio",
		"in_domain"
	)],
	ad_hu[c(
		"group",
		"compound_id",
		"ratio",
		"in_domain"
	)],
	ad_sp[c(
		"group",
		"compound_id",
		"ratio",
		"in_domain"
	)]
)

ad_plot$ratio <- pmax(
	ad_plot$ratio,
	0.001
)

ad_plot$group <- factor(
	ad_plot$group,
	levels = c(
		"SkinPiX transfer",
		"HuSkinDB",
		"Cheruvu et al., 2022 - S2",
		"Cheruvu et al., 2022 - S1"
	)
)

counts <- tapply(
	ad_plot$in_domain,
	ad_plot$group,
	function(x) paste0(
		sum(x),
		"/",
		length(x),
		" in-domain"
	)
)

group_labels <- setNames(
	paste0(
		names(counts),
		"    ",
		counts
	),
	names(counts)
)

figs1 <- ggplot(
	ad_plot,
	aes(x = ratio, y = group)
) +
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
		breaks = c(
			0.01,
			0.1,
			1,
			10,
			100
		),
		labels = c(
			"0.01",
			"0.1",
			"1",
			"10",
			"100"
		)
	) +
	scale_y_discrete(
		labels = group_labels
	) +
	labs(
		title = "Chemical applicability-domain diagnostics",
		subtitle = "Distance is normalized to the relevant 95% Mahalanobis threshold; values >1 exceed the distance criterion",
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

save_fig(
	figs1,
	"Figure_S1_applicability_domain",
	9.2,
	5.2
)

# ============================================================
# Figure S2
# HuSkin matched-source calibration heterogeneity
# ============================================================

s2_ind <- matchref31

s2_end <- s2_ind[
	s2_ind$anchors == max(s2_ind$anchors),
]

figs2 <- ggplot() +
	geom_hline(
		yintercept = 0,
		linetype = 3,
		color = c_grey
	) +
	geom_ribbon(
		data = cal60,
		aes(
			x = anchors,
			ymin = improvement_ci_low,
			ymax = improvement_ci_high
		),
		fill = c_blue,
		alpha = 0.15,
		inherit.aes = FALSE
	) +
	geom_line(
		data = s2_ind[
			s2_ind$reference_id != "HUR_48",
		],
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			group = reference_id,
			color = "Individual held-out source"
		),
		linewidth = 0.7,
		alpha = 0.6
	) +
	geom_point(
		data = s2_ind[
			s2_ind$reference_id != "HUR_48",
		],
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			color = "Individual held-out source"
		),
		size = 1.8,
		alpha = 0.65
	) +
	geom_line(
		data = s2_ind[
			s2_ind$reference_id == "HUR_48",
		],
		aes(
			x = anchors,
			y = RMSE_improvement_percent,
			color = "HUR_48"
		),
		linewidth = 1
	) +
	geom_point(
		data = s2_ind[
			s2_ind$reference_id == "HUR_48",
		],
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
		linewidth = 1.8
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
			label = paste0(
				round(improvement_percent, 1),
				"%"
			)
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
			color = ifelse(
				reference_id == "HUR_48",
				"HUR_48",
				"Individual held-out source"
			)
		),
		hjust = -0.12,
		size = 2.7,
		show.legend = FALSE
	) +
	scale_color_manual(
		values = c(
			"Mean source-level RMSE" = c_blue,
			"Individual held-out source" = c_grey_dark,
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
		title = "Calibration benefit across the matched HuSkin sources",
		subtitle = "Blue shading shows the 95% source-bootstrap interval; individual source percentages are shown separately",
		x = "Number of anchor compounds",
		y = "RMSE improvement (%)"
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "top",
		plot.margin = margin(15, 42, 12, 18)
	)

save_fig(
	figs2,
	"Figure_S2_huskin_source_calibration",
	9.2,
	5.8
)

# ------------------------------------------------------------
# Calibration failure analysis
# ------------------------------------------------------------

summary70 <- read.csv("R/results/70_failure_source_summary.csv", stringsAsFactors = FALSE)
selection70 <- read.csv("R/results/70_failure_anchor_selections.csv", stringsAsFactors = FALSE)
calibration31 <- read.csv("R/results/31_calibration_performance.csv", stringsAsFactors = FALSE)

primary70 <- summary70[summary70$anchors == 1, ]
primary70 <- primary70[order(primary70$dataset, primary70$reference_id), ]

# Use the original paired one-anchor evaluations for the advance check.
check_rows <- calibration31[calibration31$anchors == 1, ]
check_rows$RMSE_reduction <- check_rows$RMSE_uncalibrated - check_rows$RMSE_calibrated
known70 <- aggregate(RMSE_reduction ~ reference_id, data = check_rows[check_rows$task == "Known compound", ], mean)
novel70 <- aggregate(RMSE_reduction ~ reference_id, data = check_rows[check_rows$task == "New compound", ], mean)
check70 <- merge(known70, novel70, by = "reference_id", suffixes = c("_known", "_novel"))
check70$agrees <- (check70$RMSE_reduction_known > 0) == (check70$RMSE_reduction_novel > 0)

pdf("R/figures/Figure_calibration_failure.pdf", width = 8, height = 7)
layout(
	matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE),
	heights = c(1, 1)
)
par(mar = c(4.5, 4.5, 4.5, 1))

# A. Anchor offset compared with test bias
scale_a <- function(x) asinh(x)
values_a <- c(primary70$test_mean_residual, primary70$offset_min, primary70$offset_max)
limit_a <- max(abs(values_a), na.rm = TRUE)
limit_x <- 4.2
limit_y <- 3
ticks_a <- c(-4, -2, -1, -0.5, 0, 0.5, 1, 2, 4)
ticks_a <- ticks_a[abs(ticks_a) <= limit_a]

cols_a <- ifelse(
	primary70$RMSE_calibrated > primary70$RMSE_uncalibrated,
	pal[7],
	pal[1]
)
pch_a <- ifelse(primary70$dataset == "SkinPiX", 17, 16)

plot(
	scale_a(primary70$test_mean_residual),
	scale_a(primary70$offset_mean),
	type = "n",
	axes = FALSE,
	xlim = scale_a(c(-limit_x, limit_x)),
	ylim = scale_a(c(-limit_y, limit_y)),
	xlab = "Mean uncalibrated test residual",
	ylab = "Mean anchor offset",
	main = "Does the anchor offset match test bias?"
)

wrong_sign_fill <- adjustcolor("grey50", alpha.f = 0.10)

rect(
	scale_a(-6), 0, 0, scale_a(6),
	col = wrong_sign_fill, border = NA
)
rect(
	0, scale_a(-6), scale_a(6), 0,
	col = wrong_sign_fill, border = NA
)

abline(h = 0, v = 0, col = "grey85")
abline(0, 1, lty = 2, col = "grey45")

segments(
	scale_a(primary70$test_mean_residual),
	scale_a(primary70$offset_min),
	scale_a(primary70$test_mean_residual),
	scale_a(primary70$offset_max),
	col = adjustcolor(cols_a, alpha.f = 0.55)
)
points(
	scale_a(primary70$test_mean_residual),
	scale_a(primary70$offset_mean),
	pch = pch_a,
	col = cols_a,
	cex = 0.85
)
axis(1, at = scale_a(ticks_a), labels = ticks_a)
axis(2, at = scale_a(ticks_a), labels = ticks_a, las = 1)
box()

legend(
	"topleft",
	legend = c("RMSE decreased", "RMSE increased", "SkinPiX source"),
	pch = c(16, 16, 17),
	col = c(pal[1], pal[7], "grey35"),
	bty = "n",
	cex = 0.8
)

mtext("A", side = 3, line = 2, adj = 0, cex = 1.4, font = 2)

# B. HUR_48
hur48 <- selection70[selection70$reference_id == "HUR_48", ]
hur48_limits <- range(c(hur48$offset, hur48$test_mean_residual[1], 0))
hur48_padding <- max(diff(hur48_limits) * 0.1, 0.02)

boxplot(
	offset ~ anchors,
	data = hur48,
	outline = FALSE,
	col = adjustcolor(pal[2], alpha.f = 0.2),
	xlab = "Number of anchors",
	ylab = "Estimated source offset",
	main = "HUR_48",
	ylim = hur48_limits + c(-hur48_padding, hur48_padding)
)
set.seed(2026)
stripchart(
	offset ~ anchors,
	data = hur48,
	vertical = TRUE,
	method = "jitter",
	jitter = 0.12,
	add = TRUE,
	pch = 16,
	cex = 0.55,
	col = adjustcolor(pal[1], alpha.f = 0.45)
)
abline(h = 0, col = "grey80")
abline(h = hur48$test_mean_residual[1], col = pal[7], lwd = 2, lty = 2)

for (k in sort(unique(hur48$anchors))) {
	mtext(paste0("n=", sum(hur48$anchors == k)), side = 3, at = k, line = 0.4, cex = 0.8)
}
mtext("B", side = 3, line = 2, adj = 0, cex = 1.4, font = 2)

# C. SPX_R18
spx18 <- selection70[selection70$reference_id == "SPX_R18", ]
spx18_limits <- range(c(spx18$offset, spx18$test_mean_residual[1], 0))
spx18_padding <- max(diff(spx18_limits) * 0.1, 0.02)

plot(
	1, 0,
	type = "n",
	xlim = c(0.5, 3.5),
	ylim = spx18_limits + c(-spx18_padding, spx18_padding),
	xaxt = "n",
	xlab = "Number of anchors",
	ylab = "Estimated source offset",
	main = "SPX_R18"
)
axis(1, at = 1:3)
abline(h = 0, col = "grey80")
abline(h = spx18$test_mean_residual[1], col = pal[7], lwd = 2, lty = 2)

for (k in sort(unique(spx18$anchors))) {
	offsets <- spx18$offset[spx18$anchors == k]
	x_jitter <- if (length(offsets) == 1) 0 else seq(-0.08, 0.08, length.out = length(offsets))

	segments(k - 0.13, mean(offsets), k + 0.13, mean(offsets), lwd = 2, col = "black")
	points(k + x_jitter, offsets, pch = 16, cex = 1.25, col = pal[1])
	mtext(paste0("n=", length(offsets)), side = 3, at = k, line = 0.4, cex = 0.8)
}
mtext("C", side = 3, line = 2, adj = 0, cex = 1.4, font = 2)

dev.off()

cat("Diagnostic sources: ", nrow(primary70), "\n", sep = "")


cat(
	"\nManuscript figures created\n",
	"Figure 1: variance decomposition\n",
	"Figure 2: LOCO prediction with bootstrap uncertainty\n",
	"Figure 3: SkinPiX cross-dataset transfer\n",
	"Figure 4: pairwise double-holdout falsification\n",
	"Figure S1: applicability-domain diagnostics\n",
	"Figure S2: matched HuSkin calibration heterogeneity\n",
	sep = ""
)