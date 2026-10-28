# 99b_make_figures.R

library(ggplot2)
library(patchwork)

palette <- function(n = 100) {
	grDevices::colorRampPalette(c(
		"#00496f", "#0f85a0", "#20b2aa", "#f0ead6",
		"#f1d94a", "#ed8b00", "#dd4124"
	))(n)
}

cols <- c(
	setNames(palette(7), c("blue", "cyan", "teal", "cream", "yellow", "orange", "red")),
	white = "white", black = "black", grey = "grey65", dark = "grey40",
	text = "grey35", reference = "grey45", shade = "grey50", zero = "grey55",
	grid = "grey92", separator = "grey94", baseline = "grey80", crosshair = "grey85"
)

source_cols <- setNames(palette(100), as.character(seq_len(100)))

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
	ggsave(paste0("R/figures/", filename, ".pdf"), plot = plot, width = width, height = height, units = "in")
}

draw_ci <- function(low, high, estimate, y, color, pch, cap, lwd) {
	segments(low, y, high, y, col = color, lwd = lwd)
	segments(low, y - cap, low, y + cap, col = color, lwd = 1.5)
	segments(high, y - cap, high, y + cap, col = color, lwd = 1.5)
	points(estimate, y, pch = pch, col = color, cex = 1.35)
}

# Inputs

cheruvu <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)
huskin <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)
loco60 <- read.csv("R/results/60_loco_bootstrap_ci.csv", stringsAsFactors = FALSE)
var60 <- read.csv("R/results/60_variance_component_ci.csv", stringsAsFactors = FALSE)
wrong60 <- read.csv("R/results/60_wrong_source_ci.csv", stringsAsFactors = FALSE)
cal60 <- read.csv("R/results/60_huskin_calibration_ci.csv", stringsAsFactors = FALSE)
matchref31 <- read.csv("R/results/31_calibration_matched_reference.csv", stringsAsFactors = FALSE)
match41 <- read.csv("R/results/41_skinpix_calibration_matched_3anchor.csv", stringsAsFactors = FALSE)
ref41 <- read.csv("R/results/41_skinpix_calibration_reference.csv", stringsAsFactors = FALSE)
pred41 <- read.csv("R/results/41_skinpix_one_anchor_predictions.csv", stringsAsFactors = FALSE)
domain41 <- read.csv("R/results/41_skinpix_transfer_domain.csv", stringsAsFactors = FALSE)
summary70 <- read.csv("R/results/70_failure_source_summary.csv", stringsAsFactors = FALSE)
selection70 <- read.csv("R/results/70_failure_anchor_selections.csv", stringsAsFactors = FALSE)
ad21 <- read.csv("R/results/21_applicability_domain.csv", stringsAsFactors = FALSE)
ad31 <- read.csv("R/results/31_applicability_domain.csv", stringsAsFactors = FALSE)

# 01. Within-compound variation and variance decomposition

prep_centered <- function(x, dataset_label) {
	x <- x[, c("compound_id", "reference_id", "log_kp")]
	x <- x[complete.cases(x), ]
	counts <- table(x$compound_id)
	x <- x[x$compound_id %in% names(counts[counts >= 2]), ]
	x$centered_log_kp <- x$log_kp - ave(x$log_kp, x$compound_id, FUN = mean)
	spread <- aggregate(centered_log_kp ~ compound_id, data = x, function(z) diff(range(z)))
	spread <- spread[order(spread$centered_log_kp), ]
	x$compound_order <- match(x$compound_id, spread$compound_id)
	source_key <- unique(x[c("compound_id", "reference_id")])
	source_key <- source_key[order(source_key$compound_id, source_key$reference_id), ]
	source_key$source_color_index <- NA_integer_
	for (idx in split(seq_len(nrow(source_key)), source_key$compound_id)) {
		source_key$source_color_index[idx] <- if (length(idx) == 1) 1 else round(seq(1, 100, length.out = length(idx)))
	}
	x$source_color_index <- source_key$source_color_index[
		match(paste(x$compound_id, x$reference_id), paste(source_key$compound_id, source_key$reference_id))
	]
	ranges <- aggregate(centered_log_kp ~ compound_order, data = x, function(z) c(min = min(z), max = max(z)))
	ranges <- data.frame(
		compound_order = ranges$compound_order,
		ymin = ranges$centered_log_kp[, "min"],
		ymax = ranges$centered_log_kp[, "max"]
	)
	x$dataset <- dataset_label
	ranges$dataset <- dataset_label
	list(points = x, ranges = ranges)
}

centered_ch <- prep_centered(cheruvu, "Cheruvu")
centered_hu <- prep_centered(huskin, "HuSkinDB")

fig1_points <- rbind(centered_ch$points, centered_hu$points)
fig1_ranges <- rbind(centered_ch$ranges, centered_hu$ranges)

dataset_labels <- c(
	"Cheruvu (n = 81 repeated compounds)",
	"HuSkinDB (n = 89 repeated compounds)"
)

fig1_points$dataset <- factor(fig1_points$dataset, levels = c("Cheruvu", "HuSkinDB"), labels = dataset_labels)
fig1_ranges$dataset <- factor(fig1_ranges$dataset, levels = c("Cheruvu", "HuSkinDB"), labels = dataset_labels)

fig1a <- ggplot() +
	geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = cols["zero"]) +
	geom_linerange(data = fig1_ranges, aes(x = compound_order, ymin = ymin, ymax = ymax), color = cols["grey"], linewidth = 0.4) +
	geom_point(
		data = fig1_points,
		aes(x = compound_order, y = centered_log_kp, color = factor(source_color_index)),
		size = 1.8, alpha = 0.85, position = position_jitter(width = 0.12, height = 0)
	) +
	ylim(-4, 4) +
	facet_wrap(~ dataset, ncol = 1, scales = "free_x") +
	scale_color_manual(values = source_cols) +
	annotate("text", x = Inf, y = Inf, label = "Colors distinguish sources within each compound", hjust = 1.4, vjust = 1.4, size = 3.2, color = cols["text"]) +
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
		strip.text = element_text(face = "bold", hjust = 0),
		plot.margin = margin(15, 24, 12, 24)
	)

fig1_variance <- var60[c("dataset", "component", "estimate_percent")]

fig1_variance$dataset <- factor(fig1_variance$dataset, levels = c("HuSkinDB", "Cheruvu"))
fig1_variance$component <- factor(fig1_variance$component, levels = c("Residual", "Source", "Compound"))
fig1_variance$label <- paste0(round(fig1_variance$estimate_percent, 1), "%")
fig1_variance$text_color <- ifelse(fig1_variance$component == "Compound", cols["white"], cols["black"])

fig1b <- ggplot(fig1_variance, aes(x = estimate_percent, y = dataset, fill = component)) +
	geom_col(width = 0.6) +
	geom_text(aes(label = label, color = text_color), position = position_stack(vjust = 0.5), fontface = "bold", size = 4) +
	scale_color_identity() +
	scale_fill_manual(
		values = c(Compound = unname(cols["blue"]), Source = unname(cols["teal"]), Residual = unname(cols["yellow"])),
		breaks = c("Compound", "Source", "Residual")
	) +
	scale_x_continuous(breaks = seq(0, 100, 20), expand = c(0, 0)) +
	coord_cartesian(xlim = c(0, 100), clip = "off") +
	labs(title = "Variance contribution by component", x = "Variance component (%)", y = NULL) +
	theme_manuscript(12) +
	theme(
		legend.position = "top",
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(15, 24, 12, 24)
	)

fig1 <- fig1a / fig1b +
	plot_layout(heights = c(2, 1)) +
	plot_annotation(tag_levels = "A") &
	theme(plot.tag = element_text(face = "bold", size = 20, hjust = 3))

save_fig(fig1, "01_variance_decomposition", 9, 8.5)

# 02. Source information and prediction

fig2 <- loco60[
	loco60$approach %in% c("hierarchical_source_blind", "hierarchical_known_source") &
	(loco60$dataset == "Cheruvu" | (loco60$dataset == "HuSkinDB" & loco60$specification == "S1")),
]

fig2$group <- paste(fig2$dataset, fig2$specification)
fig2_order <- c("Cheruvu S1", "Cheruvu S2", "HuSkinDB S1")

fig2$blind <- fig2$approach == "hierarchical_source_blind"
fig2$y <- 4 - match(fig2$group, fig2_order) + ifelse(fig2$blind, 0.14, -0.14)
fig2$plot_col <- ifelse(fig2$blind, cols["red"], cols["blue"])
fig2$plot_pch <- ifelse(fig2$blind, 17, 16)

fig2_ticks <- seq(
	floor(min(c(0, fig2$improvement_ci_low)) / 10) * 10,
	ceiling(max(c(0, fig2$improvement_ci_high)) / 10) * 10,
	by = 10
)

pdf("R/figures/02_source_information_prediction.pdf", width = 7.6, height = 3.8)
par(mar = c(4.3, 6.7, 1.4, 0.8), mgp = c(2.2, 0.6, 0), tcl = -0.25)

plot(
	NA_real_, NA_real_, type = "n", xlim = range(fig2_ticks) + c(-1, 1), ylim = c(0.6, 3.8),
	axes = FALSE, xlab = "RMSE reduction relative to conventional QSPR (%)", ylab = ""
)
abline(v = fig2_ticks, col = cols["grid"])
abline(h = c(1.5, 2.5), col = cols["separator"])
abline(v = 0, col = cols["reference"], lwd = 1.2)
draw_ci(fig2$improvement_ci_low, fig2$improvement_ci_high, fig2$improvement_percent, fig2$y, fig2$plot_col, fig2$plot_pch, 0.055, 2)
axis(1, at = fig2_ticks)
axis(2, at = c(3, 2, 1), labels = fig2_order, las = 1, tick = FALSE, cex.axis = 1.05)
legend("topright", legend = c("Source blind", "Source aware"), col = cols[c("red", "blue")], pch = c(17, 16), horiz = TRUE, bty = "n", cex = 0.9, inset = c(0.01, 0.01))

dev.off()

# 03. Source-specific calibration control

wrong60 <- wrong60[order(wrong60$anchors), ]

rows <- seq_len(nrow(wrong60))
x_limits <- range(c(wrong60$correct_ci_low, wrong60$correct_ci_high, wrong60$wrong_ci_low, wrong60$wrong_ci_high, 0))
x_limits <- x_limits + c(-1, 1) * diff(x_limits) * 0.07

pdf("R/figures/03_source_specific_falsification.pdf", width = 8.2, height = 4.7)
par(mar = c(4.5, 8, 1.2, 1), cex.axis = 1, cex.lab = 1.1)

plot(0, 0, type = "n", xlim = x_limits, ylim = c(0.5, nrow(wrong60) + 0.5), axes = FALSE, xlab = "RMSE reduction relative to no calibration (%)", ylab = "")

abline(v = pretty(x_limits), col = cols["grid"])
abline(v = 0, col = cols["zero"], lwd = 1.3)

draw_ci(wrong60$correct_ci_low, wrong60$correct_ci_high, wrong60$correct_improvement_percent, rows, cols["blue"], 16, 0.09, 1.7)
draw_ci(wrong60$wrong_ci_low, wrong60$wrong_ci_high, wrong60$wrong_improvement_percent, rows, cols["red"], 17, 0.09, 1.7)

axis(1, at = pretty(x_limits))
axis(
	2, at = rows,
	labels = paste0(wrong60$anchors, ifelse(wrong60$anchors == 1, " anchor", " anchors"), "  (n=", wrong60$n_references, ")"),
	las = 1, tick = FALSE
)

legend("bottomright", legend = c("Target-source anchors", "Other-source anchors"), col = cols[c("blue", "red")], pch = c(16, 17), bty = "n", cex = 0.9)

dev.off()

# 04. External transfer from HuSkinDB to SkinPiX

fig4_sources <- ref41[ref41$anchors == 1, ]

fig4_sources <- fig4_sources[order(fig4_sources$RMSE_uncalibrated, decreasing = TRUE), ]

fig4_sources$source <- factor(fig4_sources$reference_id, levels = rev(fig4_sources$reference_id))
fig4_sources$change <- 100 * (fig4_sources$RMSE_calibrated / fig4_sources$RMSE_uncalibrated - 1)
fig4_sources$change_label <- paste0(sprintf("%.1f", abs(fig4_sources$change)), "% ", ifelse(fig4_sources$change < 0, "lower", "higher"))
fig4_sources$change_group <- ifelse(fig4_sources$change < 0, "Lower RMSE", "Higher RMSE")
fig4_sources$arrow_end <- fig4_sources$RMSE_calibrated + sign(fig4_sources$RMSE_uncalibrated - fig4_sources$RMSE_calibrated) * 0.015

fig4_points <- rbind(
	data.frame(source = fig4_sources$source, RMSE = fig4_sources$RMSE_uncalibrated, stage = "Before calibration"),
	data.frame(source = fig4_sources$source, RMSE = fig4_sources$RMSE_calibrated, stage = "After calibration")
)

fig4_points$stage <- factor(fig4_points$stage, levels = c("Before calibration", "After calibration"))

fig4a <- ggplot(fig4_sources, aes(y = source)) +
	geom_segment(
		aes(x = RMSE_uncalibrated, xend = arrow_end, yend = source),
		color = cols["grey"], linewidth = 1.2,
		arrow = grid::arrow(type = "closed", length = grid::unit(0.09, "inches"))
	) +
	geom_point(data = fig4_points, aes(x = RMSE, fill = stage), shape = 21, size = 4.4, stroke = 1, color = cols["dark"]) +
	geom_text(
		aes(x = pmax(RMSE_uncalibrated, RMSE_calibrated), label = change_label, color = change_group),
		hjust = -0.18, size = 3.6, fontface = "bold"
	) +
	scale_fill_manual(values = setNames(unname(cols[c("white", "blue")]), c("Before calibration", "After calibration"))) +
	scale_color_manual(values = setNames(unname(cols[c("teal", "red")]), c("Lower RMSE", "Higher RMSE")), guide = "none") +
	scale_x_continuous(breaks = seq(0.4, 1.2, 0.2)) +
	coord_cartesian(xlim = c(0.34, 1.32), clip = "off") +
	labs(title = "RMSE before and after one-anchor calibration", x = "RMSE on SkinPiX-only compounds", y = NULL) +
	theme_manuscript(12) +
	theme(
		legend.position = c(0.8, 1.2),
		legend.justification = c(0, 1),
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(8, 34, 6, 8)
	) +
	guides(fill = guide_legend(override.aes = list(shape = 21, size = 4, color = cols["dark"])))

fig4_compounds <- aggregate(cbind(abs_error_uncalibrated, abs_error_calibrated) ~ compound_id, data = pred41, mean)

n_improved <- sum(fig4_compounds$abs_error_calibrated < fig4_compounds$abs_error_uncalibrated)
n_total <- nrow(fig4_compounds)
outcomes <- c("Error decreased", "Error increased or unchanged")

fig4_counts <- data.frame(outcome = factor(outcomes, levels = outcomes), count = c(n_improved, n_total - n_improved))
fig4_counts$percent <- 100 * fig4_counts$count / n_total
fig4_counts$label <- paste0(fig4_counts$count, "/", n_total, " (", round(fig4_counts$percent), "%)")

fig4b <- ggplot(fig4_counts, aes(x = outcome, y = percent, fill = outcome)) +
	geom_col(width = 0.62) +
	geom_text(aes(label = label), vjust = -0.6, fontface = "bold", size = 3.8) +
	scale_fill_manual(values = setNames(unname(cols[c("teal", "red")]), outcomes)) +
	scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), expand = c(0, 0)) +
	scale_x_discrete(labels = c("Error\ndecreased", "Error increased\nor unchanged")) +
	labs(title = "Change in compound-level error", x = NULL, y = "SkinPiX-only compounds (%)") +
	theme_manuscript(11.5) +
	theme(
		legend.position = "none",
		axis.title.y = element_text(margin = margin(r = 5)),
		axis.line.x = element_blank(),
		axis.ticks.x = element_blank(),
		plot.margin = margin(8, 10, 8, 8)
	)
fig4c <- ggplot(match41[match41$anchors %in% 1:3, ], aes(x = anchors, y = RMSE_improvement_percent)) +
	geom_hline(yintercept = 0, linetype = 3, color = cols["grey"]) +
	geom_line(color = cols["blue"], linewidth = 1.4) +
	geom_point(color = cols["blue"], size = 3.8) +
	geom_text(aes(label = paste0(sprintf("%.1f", RMSE_improvement_percent), "%")), vjust = -1.1, fontface = "bold", color = cols["blue"], size = 3.8) +
	scale_x_continuous(breaks = 1:3) +
	scale_y_continuous(limits = c(0, 35), breaks = seq(0, 30, 10), expand = c(0, 0)) +
	coord_cartesian(xlim = c(0.9, 3.1), clip = "off") +
	labs(title = "Same three SkinPiX sources", x = "Number of anchor compounds", y = "Mean source-level RMSE reduction (%)") +
	theme_manuscript(11.5) +
	theme(
		plot.margin = margin(8, 12, 8, 8),
		axis.title.y = element_text(margin = margin(r = 5)),
		plot.tag = element_text(margin = margin(r = 20))
	)
fig4 <- fig4a / (free(fig4b, type = "label") | fig4c) +
	plot_layout(heights = c(1.05, 1), widths = c(1, 1)) +
	plot_annotation(tag_levels = "A") &
	theme(plot.tag = element_text(size = 18, face = "bold"), plot.margin = margin(8, 10, 8, 10))

save_fig(fig4, "04_skinpix_transfer", 10, 8)

# 05. Calibration failure diagnostics

primary70 <- summary70[summary70$anchors == 1, ]
primary70 <- primary70[order(primary70$dataset, primary70$reference_id), ]

pdf("R/figures/05_calibration_failure.pdf", width = 8, height = 7)
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(1, 1))
par(mar = c(4.5, 4.5, 4.5, 1))

limit_a <- max(abs(c(primary70$test_mean_residual, primary70$offset_min, primary70$offset_max)), na.rm = TRUE)
ticks_a <- c(-4, -2, -1, -0.5, 0, 0.5, 1, 2, 4)
ticks_a <- ticks_a[abs(ticks_a) <= limit_a]
failure_cols <- ifelse(primary70$RMSE_calibrated > primary70$RMSE_uncalibrated, cols["red"], cols["blue"])
failure_pch <- ifelse(primary70$dataset == "SkinPiX", 17, 16)

plot(
	asinh(primary70$test_mean_residual), asinh(primary70$offset_mean),
	type = "n", axes = FALSE, xlim = asinh(c(-4.2, 4.2)), ylim = asinh(c(-3, 3)),
	xlab = "Mean uncalibrated test residual", ylab = "Mean anchor offset",
	main = "Does the anchor offset match test bias?"
)

wrong_sign_fill <- adjustcolor(cols["shade"], alpha.f = 0.10)

rect(asinh(-6), 0, 0, asinh(6), col = wrong_sign_fill, border = NA)
rect(0, asinh(-6), asinh(6), 0, col = wrong_sign_fill, border = NA)

abline(h = 0, v = 0, col = cols["crosshair"])
abline(0, 1, lty = 2, col = cols["reference"])

segments(
	asinh(primary70$test_mean_residual), asinh(primary70$offset_min),
	asinh(primary70$test_mean_residual), asinh(primary70$offset_max),
	col = adjustcolor(failure_cols, alpha.f = 0.55)
)

points(asinh(primary70$test_mean_residual), asinh(primary70$offset_mean), pch = failure_pch, col = failure_cols, cex = 0.85)

axis(1, at = asinh(ticks_a), labels = ticks_a)
axis(2, at = asinh(ticks_a), labels = ticks_a, las = 1)

box()

legend("topleft", legend = c("RMSE decreased", "RMSE increased", "SkinPiX source"), pch = c(16, 16, 17), col = cols[c("blue", "red", "text")], bty = "n", cex = 0.8)
mtext("A", side = 3, line = 2, adj = 0, cex = 1.4, font = 2)

hur48 <- selection70[selection70$reference_id == "HUR_48", ]
hur48_limits <- range(c(hur48$offset, hur48$test_mean_residual[1], 0))
hur48_padding <- max(diff(hur48_limits) * 0.1, 0.02)

boxplot(
	offset ~ anchors, data = hur48, outline = FALSE, col = adjustcolor(cols["cyan"], alpha.f = 0.2),
	xlab = "Number of anchors", ylab = "Estimated source offset", main = "HUR_48",
	ylim = hur48_limits + c(-hur48_padding, hur48_padding)
)

set.seed(2026)
stripchart(
	offset ~ anchors, data = hur48, vertical = TRUE, method = "jitter", jitter = 0.12,
	add = TRUE, pch = 16, cex = 0.55, col = adjustcolor(cols["blue"], alpha.f = 0.45)
)

abline(h = 0, col = cols["baseline"])
abline(h = hur48$test_mean_residual[1], col = cols["red"], lwd = 2, lty = 2)

for (k in sort(unique(hur48$anchors))) {
	mtext(paste0("n=", sum(hur48$anchors == k)), side = 3, at = k, line = 0.4, cex = 0.8)
}

mtext("B", side = 3, line = 2, adj = 0, cex = 1.4, font = 2)

spx18 <- selection70[selection70$reference_id == "SPX_R18", ]
spx18_limits <- range(c(spx18$offset, spx18$test_mean_residual[1], 0))
spx18_padding <- max(diff(spx18_limits) * 0.1, 0.02)

plot(
	1, 0, type = "n", xlim = c(0.5, 3.5), ylim = spx18_limits + c(-spx18_padding, spx18_padding),
	xaxt = "n", xlab = "Number of anchors", ylab = "Estimated source offset", main = "SPX_R18"
)

axis(1, at = 1:3)

abline(h = 0, col = cols["baseline"])
abline(h = spx18$test_mean_residual[1], col = cols["red"], lwd = 2, lty = 2)

for (k in sort(unique(spx18$anchors))) {
	offsets <- spx18$offset[spx18$anchors == k]
	x_jitter <- if (length(offsets) == 1) 0 else seq(-0.08, 0.08, length.out = length(offsets))
	segments(k - 0.13, mean(offsets), k + 0.13, mean(offsets), lwd = 2, col = cols["black"])
	points(k + x_jitter, offsets, pch = 16, cex = 1.25, col = cols["blue"])
	mtext(paste0("n=", length(offsets)), side = 3, at = k, line = 0.4, cex = 0.8)
}

mtext("C", side = 3, line = 2, adj = 0, cex = 1.4, font = 2)

dev.off()

# S1. Matched HuSkinDB source calibration

individual <- matchref31[matchref31$reference_id != "HUR_48", ]
highlight <- matchref31[matchref31$reference_id == "HUR_48", ]
endpoints <- matchref31[matchref31$anchors == max(matchref31$anchors), ]

figs1 <- ggplot() +
	geom_hline(yintercept = 0, linetype = 3, color = cols["grey"]) +
	geom_ribbon(data = cal60, aes(x = anchors, ymin = improvement_ci_low, ymax = improvement_ci_high), fill = cols["blue"], alpha = 0.15) +
	geom_line(data = individual, aes(x = anchors, y = RMSE_improvement_percent, group = reference_id, color = "Individual held-out source"), linewidth = 0.7, alpha = 0.6) +
	geom_point(data = individual, aes(x = anchors, y = RMSE_improvement_percent, color = "Individual held-out source"), size = 1.8, alpha = 0.65) +
	geom_line(data = highlight, aes(x = anchors, y = RMSE_improvement_percent, color = "HUR_48"), linewidth = 1) +
	geom_point(data = highlight, aes(x = anchors, y = RMSE_improvement_percent, color = "HUR_48"), size = 2.4, show.legend = FALSE) +
	geom_line(data = cal60, aes(x = anchors, y = improvement_percent, color = "Mean source-level RMSE"), linewidth = 1.8) +
	geom_point(data = cal60, aes(x = anchors, y = improvement_percent, color = "Mean source-level RMSE"), size = 3.8) +
	geom_text(data = cal60, aes(x = anchors, y = improvement_percent, label = paste0(round(improvement_percent, 1), "%")), vjust = -1, color = cols["blue"], fontface = "bold", size = 3) +
	geom_text(
		data = endpoints,
		aes(x = anchors, y = RMSE_improvement_percent, label = reference_id, color = ifelse(reference_id == "HUR_48", "HUR_48", "Individual held-out source")),
		hjust = -0.12, size = 2.7, show.legend = FALSE
	) +
	scale_color_manual(
		values = setNames(unname(cols[c("blue", "dark", "red")]), c("Mean source-level RMSE", "Individual held-out source", "HUR_48")),
		breaks = c("Mean source-level RMSE", "Individual held-out source", "HUR_48")
	) +
	scale_x_continuous(breaks = 1:5, limits = c(0.9, 5.8)) +
	labs(
		title = "Calibration benefit across the matched HuSkin sources",
		x = "Number of anchor compounds", y = "RMSE improvement (%)"
	) +
	theme_manuscript(10.5) +
	theme(legend.position = "top", plot.margin = margin(15, 42, 12, 18))

save_fig(figs1, "S1_huskin_source_calibration", 9.2, 5.8)

# S2. Applicability-domain diagnostics

prep_domain <- function(x, label, dimensions) {
	data.frame(
		group = label,
		compound_id = x$compound_id,
		ratio = x$mahalanobis / qchisq(0.95, dimensions),
		in_domain = x$in_domain
	)
}

ad_plot <- rbind(
	prep_domain(ad21[ad21$specification == "S1", ], "Cheruvu - S1", 2),
	prep_domain(ad21[ad21$specification == "S2", ], "Cheruvu - S2", 6),
	prep_domain(ad31, "HuSkinDB", 2),
	prep_domain(domain41, "SkinPiX transfer", 2)
)

ad_plot$ratio <- pmax(ad_plot$ratio, 0.001)
ad_plot$group <- factor(ad_plot$group, levels = c("SkinPiX transfer", "HuSkinDB", "Cheruvu - S2", "Cheruvu - S1"))

counts <- tapply(ad_plot$in_domain, ad_plot$group, function(x) paste0(sum(x), "/", length(x), " in-domain"))
group_labels <- setNames(paste0(names(counts), "    ", counts), names(counts))

figs2 <- ggplot(ad_plot, aes(x = ratio, y = group)) +
	geom_vline(xintercept = 1, linetype = 2, color = cols["red"], linewidth = 0.9) +
	geom_jitter(aes(fill = in_domain), shape = 21, height = 0.13, width = 0, size = 2.25, alpha = 0.75, color = cols["white"], stroke = 0.2) +
	scale_fill_manual(
		values = setNames(unname(cols[c("teal", "red")]), c("TRUE", "FALSE")),
		labels = c("TRUE" = "In-domain", "FALSE" = "Out-of-domain")
	) +
	scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
	scale_y_discrete(labels = group_labels) +
	labs(
		title = "Chemical applicability-domain diagnostics",
		x = "Mahalanobis distance / 95% threshold", y = NULL
	) +
	theme_manuscript(10.5) +
	theme(
		legend.position = "top",
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		plot.margin = margin(16, 30, 14, 30)
	)

save_fig(figs2, "S2_applicability_domain", 9.2, 5.2)

cat(
	"\nManuscript figures created\n",
	"01_variance_decomposition.pdf\n",
	"02_source_information_prediction.pdf\n",
	"03_source_specific_falsification.pdf\n",
	"04_skinpix_transfer.pdf\n",
	"05_calibration_failure.pdf\n",
	"S1_huskin_source_calibration.pdf\n",
	"S2_applicability_domain.pdf\n",
	"Diagnostic sources: ", nrow(primary70), "\n",
	sep = ""
)