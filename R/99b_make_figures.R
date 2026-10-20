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

loco_ci <- read.csv(
	"R/results/60_loco_bootstrap_ci.csv",
	stringsAsFactors = FALSE
)

fig2 <- loco_ci[loco_ci$dataset %in% c("Cheruvu", "HuSkinDB"), ]
fig2 <- fig2[!(fig2$dataset == "HuSkinDB" & fig2$specification != "S1"), ]

fig2$panel <- ifelse(
	fig2$dataset == "Cheruvu" & fig2$specification == "S1",
	"Cheruvu et al., 2022 - S1\n148 compounds, 476 observations",
	ifelse(
		fig2$dataset == "Cheruvu" & fig2$specification == "S2",
		"Cheruvu et al., 2022 - S2\n148 compounds, 476 observations",
		"HuSkinDB - S1\n251 compounds, 546 observations"
	)
)

fig2$panel <- factor(
	fig2$panel,
	levels = c(
		"Cheruvu et al., 2022 - S1\n148 compounds, 476 observations",
		"Cheruvu et al., 2022 - S2\n148 compounds, 476 observations",
		"HuSkinDB - S1\n251 compounds, 546 observations"
	)
)

fig2$approach_label <- ifelse(
	fig2$approach == "conventional",
	"Aggregate\nQSPR",
	ifelse(
		fig2$approach == "hierarchical_source_blind",
		"Hierarchy\nonly",
		"Hierarchy\n+ source"
	)
)

fig2$approach_label <- factor(
	fig2$approach_label,
	levels = c("Aggregate\nQSPR", "Hierarchy\nonly", "Hierarchy\n+ source")
)

fmt_num <- function(x, digits = 3) {
	formatC(x, format = "f", digits = digits)
}

fig2$rmse_label <- fmt_num(fig2$RMSE, 3)

fig2$change_label <- ifelse(
	fig2$approach == "conventional",
	"reference",
	ifelse(
		fig2$improvement_percent >= 0,
		paste0(formatC(fig2$improvement_percent, format = "f", digits = 1), "% lower"),
		paste0(formatC(abs(fig2$improvement_percent), format = "f", digits = 1), "% higher")
	)
)

cols <- palette(7)
fig2_cols <- c(
	"Aggregate\nQSPR" = cols[2],
	"Hierarchy\nonly" = cols[4],
	"Hierarchy\n+ source" = cols[7]
)

y_min <- min(fig2$RMSE_ci_low, na.rm = TRUE)
y_max <- max(fig2$RMSE_ci_high, na.rm = TRUE)
y_pad <- 0.18 * (y_max - y_min)

panel_tags <- data.frame(
	panel = factor(
		c(
			"Cheruvu et al., 2022 - S1\n148 compounds, 476 observations",
			"Cheruvu et al., 2022 - S2\n148 compounds, 476 observations",
			"HuSkinDB - S1\n251 compounds, 546 observations"
		),
		levels = levels(fig2$panel)
	),
	approach_label = factor(
		c("Aggregate\nQSPR", "Aggregate\nQSPR", "Aggregate\nQSPR"),
		levels = levels(fig2$approach_label)
	),
	y = c(y_max + y_pad * 0.80, y_max + y_pad * 0.80, y_max + y_pad * 0.80),
	label = c("A", "B", "C"),
	stringsAsFactors = FALSE
)

p2 <- ggplot(
	fig2,
	aes(x = approach_label, y = RMSE, group = 1)
) +
	geom_line(
		linewidth = 0.7,
		colour = "grey55"
	) +
	geom_errorbar(
		aes(ymin = RMSE_ci_low, ymax = RMSE_ci_high),
		width = 0.10,
		linewidth = 0.6,
		colour = "grey35"
	) +
	geom_point(
		aes(fill = approach_label),
		shape = 21,
		size = 4.0,
		stroke = 0.5,
		colour = "black"
	) +
	geom_text(
		aes(label = rmse_label),
		vjust = -0.95,
		size = 3.4
	) +
	geom_text(
		data = fig2[fig2$approach != "conventional", ],
		aes(label = change_label),
		vjust = 1.85,
		size = 3.0
	) +
	geom_text(
		data = panel_tags,
		aes(x = approach_label, y = y, label = label),
		inherit.aes = FALSE,
		hjust = -0.45,
		vjust = 1,
		fontface = "bold",
		size = 5.2
	) +
	facet_wrap(
		~ panel,
		nrow = 1
	) +
	scale_fill_manual(
		values = fig2_cols,
		guide = "none"
	) +
	coord_cartesian(
		ylim = c(y_min - 0.02, y_max + y_pad),
		clip = "off"
	) +
	labs(
		title = "",
		subtitle = "Points show LOCO RMSE; error bars show 95% compound-bootstrap intervals",
		x = NULL,
		y = expression("LOCO RMSE in " * log[10] * K[p])
	) +
	theme_bw(base_size = 11) +
	theme(
		plot.title = element_text(face = "bold", size = 15),
		plot.subtitle = element_text(size = 10),
		strip.background = element_rect(fill = "grey94", colour = "grey75"),
		strip.text = element_text(face = "bold", size = 10.5),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		axis.text.x = element_text(size = 9),
		axis.title.y = element_text(size = 11),
		plot.margin = margin(12, 18, 12, 12)
	)

ggsave(
	"R/figures/Figure_2_source_information_prediction.pdf",
	p2,
	width = 11.2,
	height = 5.0,
	units = "in"
)


# ============================================================
# Figure 3
# External calibration: HuSkinDB to SkinPiX
# ============================================================

# ------------------------------------------------------------
# A. One-anchor calibration by SkinPiX source
# ------------------------------------------------------------

p3a_dat <- ref41[ref41$anchors == 1, ]

p3a_dat <- p3a_dat[
	order(p3a_dat$RMSE_uncalibrated, decreasing = TRUE),
]

p3a_dat$source <- factor(
	p3a_dat$reference_id,
	levels = rev(p3a_dat$reference_id)
)

p3a_dat$change <- 100 * (
	p3a_dat$RMSE_calibrated /
		p3a_dat$RMSE_uncalibrated - 1
)

p3a_dat$change_label <- paste0(
	ifelse(
		p3a_dat$change < 0,
		"RMSE reduced ",
		"RMSE increased "
	),
	sprintf("%.1f", abs(p3a_dat$change)),
	"%"
)

p3a_points <- rbind(
	data.frame(
		source = p3a_dat$source,
		RMSE = p3a_dat$RMSE_uncalibrated,
		stage = "Before calibration"
	),
	data.frame(
		source = p3a_dat$source,
		RMSE = p3a_dat$RMSE_calibrated,
		stage = "After calibration"
	)
)

p3a_points$stage <- factor(
	p3a_points$stage,
	levels = c(
		"Before calibration",
		"After calibration"
	)
)

p3a <- ggplot(p3a_dat, aes(y = source)) +
	geom_segment(
		aes(
			x = RMSE_uncalibrated,
			xend = RMSE_calibrated,
			yend = source
		),
		color = c_grey_light,
		linewidth = 2.2
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
		aes(
			x = pmax(RMSE_uncalibrated, RMSE_calibrated),
			label = change_label
		),
		hjust = -0.15,
		size = 3.4,
		fontface = "bold",
		color = c_blue
	) +
	scale_fill_manual(
		values = c(
			"Before calibration" = "white",
			"After calibration" = c_teal
		)
	) +
	scale_x_continuous(
		breaks = seq(0.4, 1.6, 0.2)
	) +
	coord_cartesian(
		xlim = c(0.38, 1.32),
		clip = "off"
	) +
	labs(
		title = "RMSE before and after one-anchor calibration",
		x = "RMSE on SkinPiX-only compounds",
		y = NULL
	) +
	theme_manuscript(11) +
	theme(
		legend.position = "bottom",
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(10, 20, 8, 12)
	) +
	guides(
		fill = guide_legend(
			override.aes = list(
				shape = 21,
				size = 4,
				color = c_grey_dark
			)
		)
	)

# ------------------------------------------------------------
# B. Proportion of compounds with reduced error
# ------------------------------------------------------------

p3b_dat <- aggregate(
	cbind(
		abs_error_uncalibrated,
		abs_error_calibrated
	) ~ compound_id,
	data = pred41,
	mean
)

n_improved <- sum(
	p3b_dat$abs_error_calibrated <
		p3b_dat$abs_error_uncalibrated
)

n_total <- nrow(p3b_dat)
n_other <- n_total - n_improved

p3b_summary <- data.frame(
	outcome = factor(
		c("Error decreased", "Error increased or unchanged"),
		levels = c(
			"Error decreased",
			"Error increased or unchanged"
		)
	),
	count = c(n_improved, n_other)
)

p3b_summary$percent <- 100 *
	p3b_summary$count / n_total

p3b_summary$label <- paste0(
	p3b_summary$count,
	"/",
	n_total,
	" (",
	round(p3b_summary$percent),
	"%)"
)

p3b <- ggplot(
	p3b_summary,
	aes(x = outcome, y = percent, fill = outcome)
) +
	geom_col(
		width = 0.62
	) +
	geom_text(
		aes(label = label),
		vjust = -0.6,
		fontface = "bold",
		size = 3.6
	) +
	scale_fill_manual(
		values = c(
			"Error decreased" = c_teal,
			"Error increased or unchanged" = c_red
		)
	) +
	scale_y_continuous(
		limits = c(0, 105),
		breaks = seq(0, 100, 25),
		expand = c(0, 0)
	) +
	scale_x_discrete(
		labels = c(
			"Error\ndecreased",
			"Error increased\nor unchanged"
		)
	) +
	labs(
		title = "Compounds with reduced prediction error",
		x = NULL,
		y = "SkinPiX-only compounds (%)"
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "none",
		axis.line.x = element_blank(),
		axis.ticks.x = element_blank(),
		plot.margin = margin(10, 12, 10, 12)
	)

# ------------------------------------------------------------
# C. Combined RMSE reduction with additional anchors
# ------------------------------------------------------------

p3c_dat <- match41[
	match41$anchors %in% 1:3,
]

p3c <- ggplot(
	p3c_dat,
	aes(
		x = anchors,
		y = RMSE_improvement_percent
	)
) +
	geom_hline(
		yintercept = 0,
		linetype = 3,
		color = c_grey
	) +
	geom_line(
		color = c_blue,
		linewidth = 1.5
	) +
	geom_point(
		color = c_blue,
		size = 4
	) +
	geom_text(
		aes(
			label = paste0(
				sprintf("%.1f", RMSE_improvement_percent),
				"%"
			)
		),
		vjust = -1.2,
		fontface = "bold",
		color = c_blue,
		size = 3.5
	) +
	scale_x_continuous(
		breaks = 1:3
	) +
	coord_cartesian(
		xlim = c(0.9,3.1),
		clip = "off"
	) +
	scale_y_continuous(
		limits = c(0, 35),
		breaks = seq(0, 30, 10),
		expand = c(0, 0)
	) +
	labs(
		title = "Calibration with additional anchors",
		subtitle = "Same three eligible sources",
		x = "Number of anchor compounds",
		y = "Combined RMSE reduction (%)"
	) +
	theme_manuscript(10.5) +
	theme(
		plot.margin = margin(10, 12, 10, 12)
	)

# ------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------

fig3 <- p3a / (p3b | p3c) +
	plot_layout(
		heights = c(1, 1.1),
		widths = c(1, 1)
	) +
	plot_annotation(
		title = "",
		subtitle = paste0(
			n_total,
			" SkinPiX-only compounds from ",
			nrow(p3a_dat),
			" unseen experimental sources"
		),
		tag_levels = "A",
		theme = theme(
			plot.title = element_text(
				face = "bold",
				size = 13
			),
			plot.subtitle = element_text(size = 10),
			plot.margin = margin(12, 14, 8, 14)
		)
	)

save_fig(
	fig3,
	"Figure_3_skinpix_transfer",
	10.5,
	7.8
)


# ============================================================
# Figure 4
# Source-specific calibration control
# ============================================================

wrong60 <- read.csv("R/results/60_wrong_source_ci.csv", stringsAsFactors = FALSE)

y_limits <- range(c(wrong60$correct_ci_low, wrong60$correct_ci_high, wrong60$wrong_ci_low, wrong60$wrong_ci_high, 0))
y_padding <- diff(y_limits) * 0.1
x_target <- wrong60$anchors - 0.12
x_other <- wrong60$anchors + 0.12

pdf("R/figures/Figure_4_source_specific_falsification.pdf", width = 8.2, height = 4.8)
par(mar = c(4.6, 5.4, 1, 1), cex.lab = 1.1, cex.axis = 1)

plot(wrong60$anchors, wrong60$correct_improvement_percent,
	type = "n", xlim = range(wrong60$anchors) + c(-0.5, 0.5),
	ylim = y_limits + c(-y_padding, y_padding),
	xaxt = "n", yaxt = "n",
	xlab = "Number of anchor compounds",
	ylab = "RMSE reduction relative to no calibration (%)")
axis(1, at = wrong60$anchors)
axis(2, las = 1)
abline(h = 0, lty = 3, col = "grey65")

segments(x_target, wrong60$correct_ci_low, x_target, wrong60$correct_ci_high, col = pal[1], lwd = 1.5)
segments(x_target - 0.04, wrong60$correct_ci_low, x_target + 0.04, wrong60$correct_ci_low, col = pal[1], lwd = 1.5)
segments(x_target - 0.04, wrong60$correct_ci_high, x_target + 0.04, wrong60$correct_ci_high, col = pal[1], lwd = 1.5)
points(x_target, wrong60$correct_improvement_percent, pch = 16, col = pal[1], cex = 1.3)

segments(x_other, wrong60$wrong_ci_low, x_other, wrong60$wrong_ci_high, col = pal[7], lwd = 1.5)
segments(x_other - 0.04, wrong60$wrong_ci_low, x_other + 0.04, wrong60$wrong_ci_low, col = pal[7], lwd = 1.5)
segments(x_other - 0.04, wrong60$wrong_ci_high, x_other + 0.04, wrong60$wrong_ci_high, col = pal[7], lwd = 1.5)
points(x_other, wrong60$wrong_improvement_percent, pch = 17, col = pal[7], cex = 1.3)

legend("topleft",
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

pdf("R/figures/Figure_calibration_failure.pdf", width = 9, height = 7)
par(mfrow = c(2, 2), mar = c(4.7, 4.8, 3.2, 1), mgp = c(2.8, 0.7, 0), cex.axis = 1.05, cex.lab = 1.05, cex.main = 1.05)

# A. Retrospective comparison with the test residual.
improved70 <- primary70$RMSE_calibrated < primary70$RMSE_uncalibrated
cols70 <- ifelse(improved70, pal[1], pal[7])
symbols70 <- ifelse(primary70$dataset == "HuSkinDB", 19, 17)
lims70 <- range(c(primary70$test_mean_residual, primary70$offset_min, primary70$offset_max))

plot(primary70$test_mean_residual, primary70$offset_mean,
	xlim = lims70, ylim = lims70, type = "n", asp = 1,
	xlab = "Mean uncalibrated test residual",
	ylab = "Mean anchor offset",
	main = "A. Offset and test bias")
abline(h = 0, v = 0, col = "grey85")
abline(0, 1, lty = 2, col = "grey45")
segments(primary70$test_mean_residual, primary70$offset_min, primary70$test_mean_residual, primary70$offset_max, col = cols70)
points(primary70$test_mean_residual, primary70$offset_mean, pch = symbols70, col = cols70)

legend("topleft",
	legend = c("HuSkinDB", "SkinPiX", "RMSE decreased", "RMSE increased"),
	col = c("grey30", "grey30", pal[1], pal[7]),
	pch = c(19, 17, 15, 15), bty = "n", cex = 0.8)

# B. Check available before evaluating unseen compounds.
check_cols70 <- ifelse(check70$agrees, pal[1], pal[7])
plot(check70$RMSE_reduction_known, check70$RMSE_reduction_novel,
	type = "n",
	xlab = "RMSE reduction on remaining shared compounds",
	ylab = "RMSE reduction on unseen compounds",
	main = "B. Do the two outcomes agree?")
abline(h = 0, v = 0, col = "grey65")
points(check70$RMSE_reduction_known, check70$RMSE_reduction_novel, pch = 19, col = check_cols70)
legend("topleft", c("Same outcome", "Different outcome"),
	col = c(pal[1], pal[7]), pch = 19, bty = "n", cex = 0.8)

# C and D. Variation across anchor panels for two failure examples.
for (ref in c("HUR_48", "SPX_R18")) {
	x <- selection70[selection70$reference_id == ref, ]
	panel_counts70 <- table(x$anchors)
	offset_limits70 <- range(c(x$offset, x$test_mean_residual[1], 0))
	offset_padding70 <- max(diff(offset_limits70) * 0.16, 0.04)

	boxplot(offset ~ anchors, data = x, outline = FALSE,
		col = adjustcolor(pal[2], alpha.f = 0.2),
		xlab = "Number of anchors", ylab = "Estimated source offset",
		main = paste(if (ref == "HUR_48") "C." else "D.", ref),
		ylim = offset_limits70 + c(-offset_padding70, offset_padding70))
	stripchart(offset ~ anchors, data = x, vertical = TRUE,
		method = "jitter", jitter = 0.12, add = TRUE,
		pch = 16, cex = 0.5, col = adjustcolor(pal[1], alpha.f = 0.45))
	abline(h = 0, col = "grey80")
	abline(h = x$test_mean_residual[1], col = pal[7], lwd = 2)
	axis(3, at = seq_along(panel_counts70), labels = paste0("n=", panel_counts70),
		tick = FALSE, line = -0.7, cex.axis = 0.7)
}

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