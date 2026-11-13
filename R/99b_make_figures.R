# Draw the four analysis figures used by the pooled manuscript.
# Run from the project root after 99a.

library(ggplot2)
library(patchwork)

result_dir <- "R/results"
figure_dir <- "R/figures"
example_sources <- c("DOI:10.1016/s0928-0987(03)00176-3", "DOI:10.1038/jid.1969.9")

read_checked <- function(path, required) {
	x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
	if (!nrow(x) || !all(required %in% names(x))) stop("Empty input or missing columns: ", path)
	x
}
read_result <- function(name, required) {
	read_checked(file.path(result_dir, name), required)
}
fmt <- function(x, digits = 3) {
	formatC(x, format = "f", digits = digits)
}
close_enough <- function(x, y) {
	length(x) == length(y) && all(is.finite(x) & is.finite(y) &
		abs(x - y) <= 1e-7 * pmax(1, abs(y)))
}
cols <- c(
	blue = "#00496f", teal = "#20b2aa", yellow = "#f1d94a",
	red = "#dd4124", grey = "grey65", dark = "grey40"
)
theme_manuscript <- function(base_size = 11) {
	theme_classic(base_size = base_size) +
		theme(
			plot.title = element_text(face = "bold", margin = margin(b = 6)),
			plot.subtitle = element_text(size = base_size - 1, margin = margin(b = 9)),
			plot.tag = element_text(face = "bold", size = base_size + 3),
			strip.background = element_blank(),
			strip.text = element_text(face = "bold"),
			plot.margin = margin(12, 15, 12, 12)
		)
}
save_fig <- function(p, filename, width, height) {
	ggsave(file.path(figure_dir, paste0(filename, ".pdf")), plot = p,
		width = width, height = height, units = "in", device = grDevices::cairo_pdf)
}

# Load the primary analysis and the shared publication labels.
settings81 <- read_result("81_analysis_settings.csv", c("panel_method", "n_boot", "seed", "input_md5"))
settings85 <- read_result("85_analysis_settings.csv", "md5")
if (!identical(settings81$panel_method, "exhaustive")) stop("Rerun 81 with exhaustive panels.")
pooled_path <- "data/processed/80_pooled.csv"
pooled <- read_checked(pooled_path, c("compound_id", "publication_id"))
pooled_md5 <- unname(tools::md5sum(pooled_path))
if (!identical(settings81$input_md5, pooled_md5) || !identical(settings85$md5, pooled_md5)) {
	stop("Pooled input changed. Rerun 81, 85 and subsequent reporting scripts.")
}

key <- read_checked("R/additional_files/Additional_file_1_publication_index.csv", c("publication_label", "source_id"))
if (anyDuplicated(key$source_id) || anyDuplicated(key$publication_label)) stop("Duplicate publication labels.")
label <- function(ids) {
	j <- match(ids, key$source_id)
	if (anyNA(j)) stop("Publication labels are stale. Rerun 99a.")
	key$publication_label[j]
}

points <- read_result("85_centered_measurements.csv", c("compound_id", "publication_id", "centered_log_kp"))
variance <- read_result("85_variance_components.csv", c("component", "percent", "percent_ci_low", "percent_ci_high"))
source_cal <- read_result("81_pooled_calibration_by_source.csv", c(
	"source_id", "anchors", "RMSE_uncalibrated", "RMSE_calibrated"
))
all_cal <- read_result("81_pooled_calibration_summary_all.csv", c("anchors", "n_sources", "mean_source_RMSE_calibrated"))
matched_cal <- read_result("81_pooled_calibration_summary_matched.csv", "n_sources")
paired <- read_result("81_pooled_one_vs_five_anchors.csv", c("n_sources", "mean_RMSE_decrease"))
failure <- read_result("84_failure_source_summary.csv", c(
	"source_id", "anchors", "test_mean_residual", "offset_mean", "RMSE_calibrated"
))
panels <- read_result("84_failure_anchor_selections.csv", c("source_id", "anchors", "offset", "RMSE_calibrated"))

one <- source_cal[source_cal$anchors == 1L, , drop = FALSE]
overall_one <- all_cal[all_cal$anchors == 1L, , drop = FALSE]
if (nrow(overall_one) != 1L || nrow(paired) != 1L || nrow(one) != overall_one$n_sources ||
	!close_enough(mean(one$RMSE_calibrated), overall_one$mean_source_RMSE_calibrated)) {
	stop("Primary calibration summaries disagree. Rerun 81.")
}
if (!close_enough(sum(variance$percent), 100)) stop("Variance fractions do not sum to 100.")
one$Publication <- label(one$source_id)

# Figure 1. Repeated measurements and variance fractions.
compound_rows <- split(seq_len(nrow(points)), points$compound_id)
compound_summary <- do.call(
	rbind,
	lapply(compound_rows, function(ii) {
		data.frame(
			compound_id = points$compound_id[ii[1]],
			n_sources = length(unique(points$publication_id[ii])),
			n_measurements = length(ii),
			low = min(points$centered_log_kp[ii]),
			high = max(points$centered_log_kp[ii]),
			stringsAsFactors = FALSE
		)
	})
)
rownames(compound_summary) <- NULL
compound_summary$range <- with(compound_summary, high - low)
compound_summary$source_group <- ifelse(
	compound_summary$n_sources == 1L,
	"One source",
	"Multiple sources"
)

group_levels <- c("One source", "Multiple sources")
compound_summary$source_group <- factor(compound_summary$source_group, levels = group_levels)

# Rank compounds separately within each facet.
compound_summary$rank <- NA_integer_

for (group in group_levels) {
	ii <- which(compound_summary$source_group == group)
	ii <- ii[order(compound_summary$range[ii], compound_summary$compound_id[ii])]
	compound_summary$rank[ii] <- seq_along(ii)
}

point_index <- match(points$compound_id, compound_summary$compound_id)

plot_points <- points
plot_points$rank <- compound_summary$rank[point_index]
plot_points$source_group <- compound_summary$source_group[point_index]

facet_labels <- setNames(
	vapply(group_levels, function(group) {
		x <- compound_summary[compound_summary$source_group == group, , drop = FALSE]
		paste0(group, "\n", nrow(x), " compounds, ", sum(x$n_measurements), " measurements")
	}, character(1)),
	group_levels
)

# A. Observed variation among repeated measurements.
p1a <- ggplot() +
	geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
	geom_linerange(
		data = compound_summary,
		aes(x = rank, ymin = low, ymax = high),
		colour = "grey70",
		linewidth = 0.35
	) +
	geom_point(
		data = plot_points,
		aes(x = rank, y = centered_log_kp),
		colour = cols[["blue"]],
		size = 1.1,
		alpha = 0.65,
		position = position_jitter(width = 0.12, height = 0, seed = 2026)
	) +
	facet_grid(
		cols = vars(source_group),
		scales = "free_x",
		space = "free_x",
		labeller = labeller(source_group = facet_labels)
	) +
	scale_x_continuous(breaks = c(1, 20, 40, 60, 80), expand = expansion(mult = c(0.025, 0.025))) +
	labs(
		title = "Variation within compounds",
		subtitle = paste0(nrow(compound_summary), " compounds and ", nrow(plot_points), " repeated measurements"),
		x = "Compounds ranked by observed range within each group",
		y = expression(log[10] * K[p] ~ "minus compound mean")
	) +
	theme_manuscript() +
	theme(
		legend.position = "none",
		strip.background = element_blank(),
		strip.text = element_text(face = "bold", size = 10),
		panel.spacing.x = grid::unit(1.0, "lines"),
		axis.text.x = element_blank(),
		axis.ticks.x = element_blank()
	)

# B. Stacked variance components.
component_order <- c("Compound", "Publication", "Residual")

component_index <- match(component_order, as.character(variance$component))

if (anyNA(component_index)) stop("Missing variance component.")
variance_bar <- variance[component_index, ]
variance_bar$name_label <- c("Compound", "Source", "Residual")
variance_bar$percent_label <- paste0(fmt(variance_bar$percent, 1), " %")
variance_bar$ci_label <- paste0(
	"[",
	fmt(variance_bar$percent_ci_low, 1),
	", ",
	fmt(variance_bar$percent_ci_high, 1),
	"]"
)
variance_bar$xmax <- cumsum(variance_bar$percent)
variance_bar$xmin <- c(0, head(variance_bar$xmax, -1))
variance_bar$width <- with(variance_bar, xmax - xmin)
variance_bar$fill_colour <- c(cols[["blue"]], cols[["teal"]], cols[["yellow"]])
variance_bar$text_colour <- c("white", "black", "black")

# Place text within each rectangle.
variance_bar$x_percent <- with(variance_bar, xmin + 0.28 * width + c(0, 3.5, 1.5))
variance_bar$x_name <- with(variance_bar, xmin + 0.7 * width - c(0, 3, 0))
variance_bar$size <- c(4, 4.5, 4)
variance_bar$y_percent <- 0.28
variance_bar$y_ci <- -0.08
variance_bar$y_name <- with(variance_bar, (y_percent + y_ci) / 2)

p1b <- ggplot(variance_bar) +
	geom_rect(
		aes(
			xmin = xmin,
			xmax = xmax,
			ymin = -0.5,
			ymax = 0.75,
			fill = fill_colour
		),
		colour = "white",
		linewidth = 0.5
	) +
	geom_text(
		aes(x = x_percent, y = y_percent, label = percent_label, colour = text_colour),
		size = 4
	) +
	geom_text(
		aes(x = x_percent, y = y_ci, label = ci_label, colour = text_colour),
		size = 3
	) +
	geom_text(
		aes(
			x = x_name,
			y = y_name,
			label = name_label,
			colour = text_colour,
			size = size
		)
	) +
	scale_size_identity() +
	scale_fill_identity() +
	scale_colour_identity() +
	scale_x_continuous(
		limits = c(0, 100),
		breaks = seq(0, 100, 20),
		expand = expansion(mult = c(0, 0))
	) +
	scale_y_continuous(limits = c(-0.75, 0.8), breaks = NULL, expand = expansion(mult = c(0, 0))) +
	labs(
		title = "Pooled variance decomposition",
		subtitle = paste0(
			nrow(pooled), " measurements, ",
			length(unique(pooled$compound_id)),
			" compounds, ",
			length(unique(pooled$publication_id)),
			" sources"
		),
		x = "Variance component (% of total)",
		y = NULL
	) +
	theme_manuscript() +
	theme(
		legend.position = "none",
		axis.line.y = element_blank(),
		axis.ticks.y = element_blank(),
		axis.text.y = element_blank(),
		panel.grid = element_blank()
	)

p1 <- p1a / p1b + plot_layout(heights = c(2, 1)) + plot_annotation(tag_levels = "A")
save_fig(p1, "01_variance_decomposition", 7, 5.5)

# Figure 2. Calibration benefit from zero to five anchors.
matched_ids <- unique(source_cal$source_id[source_cal$anchors == 5L])
fig2_n_sources <- length(matched_ids)
matched <- source_cal[source_cal$source_id %in% matched_ids, , drop = FALSE]
if (fig2_n_sources != paired$n_sources || any(matched_cal$n_sources != fig2_n_sources)) {
	stop("Fixed subset does not agree with the saved summaries.")
}

panel_counts <- table(factor(matched$source_id, levels = matched_ids), factor(matched$anchors, levels = 1:5))
if (any(panel_counts != 1L) || any(!matched$anchors %in% 1:5)) {
	stop("Expected one result per source and anchor count.")
}

# Use the same uncalibrated baseline for each source at every anchor count.
baseline <- matched[matched$anchors == 1L, c("source_id", "RMSE_uncalibrated"), drop = FALSE]
baseline <- baseline[match(matched_ids, baseline$source_id), , drop = FALSE]
baseline_index <- match(matched$source_id, baseline$source_id)
if (!close_enough(matched$RMSE_uncalibrated, baseline$RMSE_uncalibrated[baseline_index])) {
	stop("Uncalibrated RMSE changes across anchor counts.")
}

# Rows are sources and columns are anchor counts. The zero-anchor gain is zero.
anchor_counts <- 0:5
gain_matrix <- vapply(anchor_counts, function(k) {
	if (k == 0L) return(rep(0, fig2_n_sources))
	x <- matched[matched$anchors == k, , drop = FALSE]
	j <- match(matched_ids, x$source_id)
	baseline$RMSE_uncalibrated - x$RMSE_calibrated[j]
}, numeric(fig2_n_sources))
colnames(gain_matrix) <- as.character(anchor_counts)
additional_gain <- gain_matrix[, "5"] - gain_matrix[, "1"]
if (!close_enough(mean(additional_gain), paired$mean_RMSE_decrease)) {
	stop("One-to-five decrease disagrees with the paired summary.")
}

# Resample sources with their complete trajectories.
set.seed(as.integer(settings81$seed) + 300L)
bootstrap_gain <- replicate(as.integer(settings81$n_boot), {
	sampled_sources <- sample.int(fig2_n_sources, size = fig2_n_sources, replace = TRUE)
	colMeans(gain_matrix[sampled_sources, , drop = FALSE])
})
mean_gain <- data.frame(
	anchors = anchor_counts,
	RMSE_decrease = colMeans(gain_matrix),
	ci_low = apply(bootstrap_gain, 1, quantile, probs = 0.025, names = FALSE),
	ci_high = apply(bootstrap_gain, 1, quantile, probs = 0.975, names = FALSE)
)

# Individual source trajectories, including zero anchors.
fig2_trajectories <- data.frame(
	source_id = rep(matched_ids, times = length(anchor_counts)),
	anchors = rep(anchor_counts, each = fig2_n_sources),
	RMSE_decrease = as.vector(gain_matrix),
	stringsAsFactors = FALSE
)

# Red identifies worsening from one to five anchors.
worsened_ids <- matched_ids[additional_gain < -1e-10]
fig2_trajectories$worsened_with_five <- fig2_trajectories$source_id %in% worsened_ids

p2 <- ggplot() +
	geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
	geom_ribbon(
		data = mean_gain,
		aes(x = anchors, ymin = ci_low, ymax = ci_high),
		fill = cols[["yellow"]],
		alpha = 0.7
	) +
	geom_line(
		data = fig2_trajectories[!fig2_trajectories$worsened_with_five, ],
		aes(x = anchors, y = RMSE_decrease, group = source_id),
		colour = cols[["blue"]],
		linewidth = 0.55,
		alpha = 0.75
	) +
	geom_line(
		data = fig2_trajectories[fig2_trajectories$worsened_with_five, ],
		aes(x = anchors, y = RMSE_decrease, group = source_id),
		colour = cols[["red"]],
		linewidth = 0.75
	) +
	geom_line(
		data = mean_gain,
		aes(x = anchors, y = RMSE_decrease),
		colour = cols[["teal"]],
		linewidth = 1.1
	) +
	geom_point(
		data = mean_gain,
		aes(x = anchors, y = RMSE_decrease),
		colour = cols[["teal"]],
		size = 2.5
	) +
	annotate(
		"text",
		x = 5.08,
		y = mean_gain$RMSE_decrease[mean_gain$anchors == 5L],
		label = "Mean",
		hjust = 0,
		colour = cols[["teal"]],
		size = 3.5
	) +
	scale_x_continuous(breaks = 0:5, expand = expansion(mult = c(0.02, 0.02))) +
	scale_y_continuous(expand = expansion(mult = c(0.08, 0))) +
	coord_cartesian(xlim = c(0, 5), ylim = c(NA, 1.0), clip = "off") +
	labs(
		title = "Calibration benefit with increasing number of anchors",
		x = "Number of anchor compounds",
		y = expression("RMSE decrease (" * log[10] * K[p] ~ "units)")
	) +
	theme_manuscript() +
	theme(legend.position = "none", plot.margin = margin(12, 48, 12, 12))

save_fig(p2, "02_anchor_number", 7, 4)

# Figure 3. Source-level RMSE changes with one anchor.
fig3_sources <- one
fig3_sources$RMSE_decrease <- with(fig3_sources, RMSE_uncalibrated - RMSE_calibrated)

fig3_mean_decrease <- mean(fig3_sources$RMSE_decrease)

fig3_mean_label <- paste0("Mean RMSE decrease = ", fmt(fig3_mean_decrease, 3))
fig3_sources$Outcome <- ifelse(
	fig3_sources$RMSE_decrease > 1e-10,
	"Improved",
	ifelse(fig3_sources$RMSE_decrease < -1e-10, "Worsened", "Unchanged")
)

source_order <- order(fig3_sources$RMSE_decrease, fig3_sources$source_id)
fig3_sources$Publication <- factor(fig3_sources$Publication, levels = rev(fig3_sources$Publication[source_order]))

p3 <- ggplot(
	fig3_sources, aes(y = Publication, colour = Outcome)
) +
	geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
	geom_vline(
		xintercept = fig3_mean_decrease,
		colour = cols[["teal"]],
		linetype = "dashed",
		linewidth = 0.8
	) +
	geom_segment(aes(x = 0, xend = RMSE_decrease, yend = Publication), linewidth = 0.6) +
	geom_point(aes(x = RMSE_decrease), size = 2.2) +
	annotate(
		"text",
		x = fig3_mean_decrease + 0.02,
		y = Inf,
		label = fig3_mean_label,
		hjust = 0,
		vjust = 1.2,
		size = 3.5,
		colour = cols[["teal"]]
	) +
	scale_colour_manual(
		values = c(Improved = cols[["blue"]], Worsened = cols[["red"]], Unchanged = cols[["grey"]])
	) +
	scale_x_continuous(
		breaks = function(limits) pretty(limits, n = 5),
		expand = expansion(mult = c(0.08, 0.08))
	) +
	scale_y_discrete(expand = expansion(add = 0.7)) +
	labs(
		title = "One-anchor calibration",
		y = "Source",
		x = expression("RMSE decrease (" * log[10] * K[p] ~ "units)")
	) +
	theme_manuscript() +
	theme(legend.position = "none", axis.text.y = element_text(size = 8))

save_fig(p3, "03_one_anchor_calibration", 6, max(5, 0.07 * nrow(fig3_sources) + 2))

# Figure 4. Two contrasting calibration failures.
f <- failure[failure$anchors == 1L, , drop = FALSE]
j <- match(f$source_id, one$source_id)
if (nrow(f) != nrow(one) || anyNA(j) ||
	!close_enough(f$RMSE_calibrated, one$RMSE_calibrated[j])) {
	stop("Failure results disagree with 81. Rerun 84.")
}

fig4_rows <- do.call(rbind, lapply(example_sources, function(id) {
	do.call(rbind, lapply(c(1L, 5L), function(k) {
		x <- panels[panels$source_id == id & panels$anchors == k, , drop = FALSE]
		z <- failure[failure$source_id == id & failure$anchors == k, , drop = FALSE]
		if (!nrow(x) || nrow(z) != 1L) stop("Missing example results for ", id)
		if (!close_enough(mean(x$offset), z$offset_mean) ||
			!close_enough(mean(x$RMSE_calibrated), z$RMSE_calibrated)) {
			stop("Example panel results disagree with the diagnostic summary.")
		}
		data.frame(
			source_id = id,
			anchors = k,
			row_y = if (k == 1L) 2 else 1,
			offset_mean = mean(x$offset),
			offset_low = min(x$offset),
			offset_high = max(x$offset),
			test_adjustment = z$test_mean_residual,
			RMSE = z$RMSE_calibrated,
			stringsAsFactors = FALSE
		)
	}))
}))

fig4_xlim <- with(fig4_rows, range(offset_low, offset_high, test_adjustment, 0))
fig4_xlim <- c(floor(fig4_xlim[1] * 10) / 10, ceiling(fig4_xlim[2] * 10) / 10)

fig4_example <- function(id, panel_tag, heading) {
	x <- fig4_rows[fig4_rows$source_id == id, , drop = FALSE]
	target <- unique(x$test_adjustment)
	if (length(target) != 1L) stop("Test benchmark changes across anchor counts.")

	p_offset <- ggplot(x, aes(y = row_y)) +
		annotate("segment", x = 0, xend = 0, y = 0.65, yend = 2.25, colour = "grey80", linewidth = 0.5) +
		annotate(
			"segment", x = target, xend = target, y = 0.65, yend = 2.25,
			colour = cols[["red"]], linetype = "dashed", linewidth = 0.65
		) +
		geom_segment(
			aes(x = offset_low, xend = offset_high, yend = row_y), colour = "grey60", linewidth = 1.0
		) +
		geom_point(aes(x = offset_mean), colour = cols[["teal"]], size = 3.8) +
		scale_x_continuous(
			limits = fig4_xlim, breaks = seq(-2, 2, by = 0.5), expand = expansion(mult = c(0.035, 0.035))
		) +
		scale_y_continuous(
			limits = c(0.5, 2.5), breaks = c(1, 2), labels = c("Five anchors", "One anchor"),
			expand = expansion(mult = c(0, 0))
		) +
		labs(
			title = paste(panel_tag, heading, paste0("(", label(id), ")")),
			x = expression("Estimated adjustment (" * log[10] * K[p] ~ "units)"),
			y = NULL
		) +
		theme_manuscript() +
		theme(
			plot.title = element_text(size = 12, face = "bold"),
			axis.text.y = element_text(size = 10),
			axis.ticks.y = element_blank(),
			axis.line.y = element_blank(),
			plot.margin = margin(8, 8, 8, 8)
		)

	x$rmse_label <- fmt(x$RMSE, 3)
	x$text_colour <- ifelse(
		x$anchors == 5L & x$RMSE > x$RMSE[x$anchors == 1L], cols[["red"]], cols[["dark"]]
	)
	p_rmse <- ggplot(x, aes(y = row_y)) +
		geom_text(aes(x = 0.5, label = rmse_label, colour = text_colour), size = 4.8, fontface = "bold") +
		scale_colour_identity() +
		scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
		scale_y_continuous(limits = c(0.5, 2.5), expand = c(0, 0)) +
		labs(title = "Mean RMSE", x = " ", y = NULL) +
		theme_void(base_size = 11) +
		theme(
			plot.title = element_text(size = 11, face = "bold", hjust = 0.5, margin = margin(b = 6)),
			axis.title.x = element_text(size = 11, margin = margin(t = 5.5)),
			plot.margin = margin(8, 8, 8, 8)
		)
	p_offset | p_rmse
}

fig4_top <- fig4_example(example_sources[1], "A", "Wrong direction") +
	plot_layout(widths = c(5, 1))

fig4_bottom <- fig4_example(example_sources[2], "B", "Correction too large") +
	plot_layout(widths = c(5, 1))

p4 <- (fig4_top / fig4_bottom) +
	plot_annotation(
		title = "Calibration offsets and prediction error",
		theme = theme(
			plot.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
			plot.margin = margin(12, 12, 10, 12)
		)
	)

save_fig(p4, "04_calibration_failures", 8, 4.9)

message("Manuscript figures written to ", figure_dir, ".")

# Remove the old figure names and unused raster exports.
obsolete_figures <- c(
	"02_one_anchor_calibration.pdf", "03_anchor_number.pdf",
	"01_variance_decomposition.png", "02_one_anchor_calibration.png", "03_anchor_number.png",
	"02_anchor_number.png", "03_one_anchor_calibration.png", "04_calibration_failures.png"
)
unlink(file.path(figure_dir, obsolete_figures))
