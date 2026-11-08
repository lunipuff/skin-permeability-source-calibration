# 99b_make_figures.R
#
# Run from the project root after 80-85 and 99a:
# Rscript R/99b_make_figures.R
#
# Uses saved results only.
# No models or bootstraps are rerun.

library(ggplot2)
library(patchwork)

result_dir <- "R/results"
figure_dir <- "R/figures"

# Set TRUE to generate PNG copies alongside the PDFs.
save_png <- FALSE

example_sources <- c(
  "DOI:10.1016/s0928-0987(03)00176-3",
  "DOI:10.1038/jid.1969.9"
)

# ------------------------------------------------------------
# Helpers and appearance
# ------------------------------------------------------------

read_checked <- function(path, required) {
  if (!file.exists(path)) {
    stop("Missing input: ", path)
  }

  x <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!nrow(x) || !all(required %in% names(x))) {
    stop("Empty input or missing columns: ", path)
  }

  x
}

read_result <- function(name, required) {
  read_checked(file.path(result_dir, name), required)
}

fmt <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

close_enough <- function(x, y) {
  all(
    is.finite(x) &
      is.finite(y) &
      abs(x - y) <= 1e-7 * pmax(1, abs(y))
  )
}

palette <- function(n) {
  grDevices::colorRampPalette(c(
    "#00496f", "#0f85a0", "#20b2aa", "#f0ead6",
    "#f1d94a", "#ed8b00", "#dd4124"
  ))(n)
}

cols <- c(
  blue = "#00496f",
  teal = "#20b2aa",
  yellow = "#f1d94a",
  red = "#dd4124",
  orange = "#ed8b00",
  grey = "grey65",
  dark = "grey40",
  grid = "grey92"
)

theme_manuscript <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        margin = margin(b = 6)
      ),
      plot.subtitle = element_text(
        size = base_size - 1,
        margin = margin(b = 9)
      ),
      plot.caption = element_text(
        size = base_size - 2,
        hjust = 0
      ),
      plot.tag = element_text(
        face = "bold",
        size = base_size + 3
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(12, 15, 12, 12)
    )
}

figure_index <- list()

save_fig <- function(p, filename, width, height, caption) {
  ggsave(
    file.path(figure_dir, paste0(filename, ".pdf")),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )

  if (save_png) {
    ggsave(
      file.path(figure_dir, paste0(filename, ".png")),
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      bg = "white"
    )
  }

  figure_index[[length(figure_index) + 1L]] <<-
    data.frame(
      file = paste0(filename, ".pdf"),
      width = width,
      height = height,
      caption = caption
    )
}

# ------------------------------------------------------------
# 1. Read pooled outputs
# ------------------------------------------------------------

# Reporting must not describe old sampled outputs as exhaustive evaluation.
panel_settings_files <- c(
  file.path(result_dir, "81_analysis_settings.csv"),
  list.files(file.path(result_dir, "sensitivity"),
    pattern = "^81_analysis_settings\\.csv$", recursive = TRUE, full.names = TRUE)
)
for (settings_file in panel_settings_files) {
  panel_settings <- read_checked(settings_file, "panel_method")
  if (nrow(panel_settings) != 1L || panel_settings$panel_method != "exhaustive") {
    stop("Rerun exhaustive calibration and sensitivity analyses before reporting: ", settings_file)
  }
}

key <- read_checked(
  "R/tables/S1_publication_index.csv",
  c("publication_label", "source_id")
)

if (
  anyDuplicated(key$source_id) ||
  anyDuplicated(key$publication_label)
) {
  stop("Duplicate publication labels.")
}

label <- function(ids) {
  j <- match(ids, key$source_id)

  if (anyNA(j)) {
    stop("Publication labels are stale. Rerun 99a.")
  }

  key$publication_label[j]
}

points <- read_result(
  "85_centered_measurements.csv",
  c(
    "compound_id", "publication_id",
    "centered_log_kp", "compound_order"
  )
)

variance <- read_result(
  "85_variance_components.csv",
  c(
    "component", "percent",
    "percent_ci_low", "percent_ci_high"
  )
)

source_cal <- read_result(
  "81_pooled_calibration_by_source.csv",
  c(
    "source_id", "anchors",
    "RMSE_uncalibrated", "RMSE_calibrated",
    "RMSE_reduction_percent"
  )
)

all_cal <- read_result(
  "81_pooled_calibration_summary_all.csv",
  c(
    "anchors", "n_sources", "n_sources_improved",
    "mean_source_RMSE_uncalibrated",
    "mean_source_RMSE_calibrated",
    "RMSE_reduction_percent",
    "reduction_ci_low", "reduction_ci_high"
  )
)

matched_cal <- read_result(
  "81_pooled_calibration_summary_matched.csv",
  c(
    "anchors", "n_sources",
    "RMSE_reduction_percent",
    "reduction_ci_low", "reduction_ci_high"
  )
)

paired <- read_result(
  "81_pooled_one_vs_five_anchors.csv",
  c(
    "n_sources", "mean_RMSE_decrease",
    "decrease_ci_low", "decrease_ci_high"
  )
)

failure <- read_result(
  "84_failure_source_summary.csv",
  c(
    "source_id", "anchors",
    "test_mean_residual",
    "offset_mean", "offset_min", "offset_max",
    "RMSE_uncalibrated", "RMSE_calibrated"
  )
)

panels <- read_result(
  "84_failure_anchor_selections.csv",
  c("source_id", "anchors", "offset")
)

sensitivity <- read_checked(
  file.path(
    result_dir,
    "sensitivity/83_sensitivity_comparison.csv"
  ),
  c(
    "scenario", "curation_policy", "log_tolerance",
    "one_anchor_reduction_percent",
    "one_anchor_ci_low", "one_anchor_ci_high",
    "one_to_five_RMSE_decrease",
    "one_to_five_ci_low", "one_to_five_ci_high"
  )
)

# ------------------------------------------------------------
# 2. Check consistency
# ------------------------------------------------------------

one <- source_cal[
  source_cal$anchors == 1,
  ,
  drop = FALSE
]

overall_one <- all_cal[
  all_cal$anchors == 1,
  ,
  drop = FALSE
]

if (
  nrow(overall_one) != 1 ||
  nrow(paired) != 1 ||
  nrow(one) != overall_one$n_sources
) {
  stop("Unexpected primary summaries.")
}

if (!close_enough(
  mean(one$RMSE_calibrated),
  overall_one$mean_source_RMSE_calibrated
)) {
  stop("Primary results disagree.")
}

if (!close_enough(sum(variance$percent), 100)) {
  stop("Variance fractions do not sum to 100.")
}

for (stage in c("81", "85")) {
  settings <- read.csv(
    file.path(
      result_dir,
      paste0(stage, "_analysis_settings.csv")
    ),
    stringsAsFactors = FALSE
  )

  hash <- if (stage == "81") {
    settings$input_md5
  } else {
    settings$md5
  }

  if (
    length(hash) != 1 ||
    is.na(hash) ||
    hash != unname(tools::md5sum(
      "data/processed/80_pooled.csv"
    ))
  ) {
    stop(
      "Pooled input changed. Rerun ",
      stage,
      " and subsequent reporting scripts."
    )
  }
}

if (!all(example_sources %in% failure$source_id)) {
  stop(
    "An example publication is absent. ",
    "Update example_sources to two publications represented in 84."
  )
}

one$Publication <- label(one$source_id)
failure$Publication <- label(failure$source_id)
source_cal$Publication <- label(source_cal$source_id)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 3. Figure 1: repeated measurements and variance fractions
# ------------------------------------------------------------

# Summarize repeated measurements by compound.
compound_rows <- split(
  seq_len(nrow(points)),
  points$compound_id
)

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

compound_summary$range <- with(
  compound_summary,
  high - low
)

compound_summary$source_group <- ifelse(
  compound_summary$n_sources == 1L,
  "One source",
  "Multiple sources"
)

group_levels <- c("One source", "Multiple sources")

compound_summary$source_group <- factor(
  compound_summary$source_group,
  levels = group_levels
)

# Rank compounds separately within each facet.
compound_summary$rank <- NA_integer_

for (group in group_levels) {
  ii <- which(compound_summary$source_group == group)

  ii <- ii[order(
    compound_summary$range[ii],
    compound_summary$compound_id[ii]
  )]

  compound_summary$rank[ii] <- seq_along(ii)
}

point_index <- match(
  points$compound_id,
  compound_summary$compound_id
)

plot_points <- points
plot_points$rank <- compound_summary$rank[point_index]
plot_points$source_group <-
  compound_summary$source_group[point_index]

facet_labels <- setNames(
  vapply(group_levels, function(group) {
    x <- compound_summary[
      compound_summary$source_group == group,
      , drop = FALSE
    ]

    paste0(
      group, "\n",
      nrow(x), " compounds, ",
      sum(x$n_measurements), " measurements"
    )
  }, character(1)),
  group_levels
)

# A. Observed variation among repeated measurements.
p1a <- ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey65",
    linewidth = 0.35
  ) +
  geom_linerange(
    data = compound_summary,
    aes(
      x = rank,
      ymin = low,
      ymax = high
    ),
    colour = "grey70",
    linewidth = 0.35
  ) +
  geom_point(
    data = plot_points,
    aes(
      x = rank,
      y = centered_log_kp
    ),
    colour = cols[["blue"]],
    size = 1.1,
    alpha = 0.65,
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 2026
    )
  ) +
  facet_grid(
    cols = vars(source_group),
    scales = "free_x",
    space = "free_x",
    labeller = labeller(source_group = facet_labels)
  ) +
  scale_x_continuous(
		  breaks = c(1, 20, 40, 60, 80),
		  expand = expansion(mult = c(0.025, 0.025))
	) +
  labs(
    title = "Variation within compounds",
    subtitle = paste0(
      nrow(compound_summary), " compounds and ",
      nrow(plot_points), " repeated measurements"
    ),
    x = "Compounds ranked by observed range within each group",
    y = expression(
  		log[10] * K[p] ~ "minus compound mean"
	)
  ) +
  theme_manuscript() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    panel.spacing.x = grid::unit(1.0, "lines"),
    axis.text.x = element_blank(),
	axis.ticks.x = element_blank(),
  )

# B. Stacked variance components.

component_order <- c("Compound", "Publication", "Residual")

component_index <- match(
  component_order,
  as.character(variance$component)
)

if (anyNA(component_index)) {
  stop("Missing variance component.")
}

variance_bar <- variance[component_index, ]

# Separate text fields.
variance_bar$name_label <- c(
  "Compound", "Source", "Residual"
)

variance_bar$percent_label <- paste0(
  fmt(variance_bar$percent, 1),
  " %"
)

variance_bar$ci_label <- paste0(
  "[",
  fmt(variance_bar$percent_ci_low, 1),
  ", ",
  fmt(variance_bar$percent_ci_high, 1),
  "]"
)

# Rectangle positions.
variance_bar$xmax <- cumsum(variance_bar$percent)
variance_bar$xmin <- c(0, head(variance_bar$xmax, -1))
variance_bar$width <- with(
  variance_bar,
  xmax - xmin
)

variance_bar$fill_colour <- c(
  cols[["blue"]],
  cols[["teal"]],
  cols[["yellow"]]
)

variance_bar$text_colour <- c(
  "white", "black", "black"
)

# Horizontal text positions within each rectangle.
# Increase these fractions to move the text right.
variance_bar$x_percent <- with(
  variance_bar,
  xmin + 0.28 * width
)

variance_bar$x_percent[2] <- variance_bar$x_percent[2] + 3.5
variance_bar$x_percent[3] <- variance_bar$x_percent[3] + 1.5

variance_bar$x_name <- with(
  variance_bar,
  xmin + 0.7 * width
)

variance_bar$size <- c(4,4.5,4)

variance_bar$x_name[2] <- variance_bar$x_name[2] - 3

# Vertical text positions.
variance_bar$y_percent <- 0.28
variance_bar$y_ci <- -0.08

# Centre the component name between the two numerical lines.
variance_bar$y_name <- with(
  variance_bar,
  (y_percent + y_ci) / 2
)

pooled_fig1 <- read_checked(
  "data/processed/80_pooled.csv",
  c("compound_id", "publication_id")
)

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
    aes(
      x = x_percent,
      y = y_percent,
      label = percent_label,
      colour = text_colour
    ),
    size = 4
  ) +
  geom_text(
    aes(
      x = x_percent,
      y = y_ci,
      label = ci_label,
      colour = text_colour
    ),
    size = 3
  ) +
  geom_text(
    aes(
      x = x_name,
      y = y_name,
      label = name_label,
      colour = text_colour,
      size = size
    ),
  ) +
  scale_size_identity() +
  scale_fill_identity() +
  scale_colour_identity() +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(-0.75, 0.8),
    breaks = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Pooled variance decomposition",
    subtitle = paste0(
      nrow(pooled_fig1), " measurements, ",
      length(unique(pooled_fig1$compound_id)),
      " compounds, ",
      length(unique(pooled_fig1$publication_id)),
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

save_fig(
  p1a / p1b +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(tag_levels = "A"),

  "01_variance_decomposition",
  7,
  5.5,

  paste(
    "Within-compound variation and estimated variance fractions.",
    "A: Measurements are centred on each compound's mean.",
    "Facets distinguish compounds measured repeatedly within one",
    "source from compounds measured across multiple sources.",
    "Points represent measurements and vertical lines show ranges.",
    "Compounds are ranked by observed range within each facet.",
    "Both facets use the same vertical scale.",
    "Observed ranges depend on measurement counts and do not",
    "by themselves establish a source effect.",
    "B: Variance fractions from the full pooled dataset using an",
    "intercept-only REML model with compound and source random",
    "intercepts, without molecular descriptors or temperature.",
    "Labels show estimated percentages with marginal 95% parametric",
	"bootstrap confidence intervals in brackets.",
    "Source denotes the original publication.",
    "Source-associated variance does not identify experimental causes."
  )
)

# ------------------------------------------------------------
# Shared paired-dot plot helper
# ------------------------------------------------------------

paired_plot <- function(
  x, before, after,
  title, subtitle, x_title
) {
  x$before <- x[[before]]
  x$after <- x[[after]]

  change <- x$after - x$before

  x$Outcome <- ifelse(
    change > 1e-10,
    "Worsened",
    ifelse(
      change < -1e-10,
      "Improved",
      "Unchanged"
    )
  )

  order_ids <- x$Publication[
    order(change, decreasing = TRUE)
  ]

  x$Publication <- factor(
    x$Publication,
    levels = rev(order_ids)
  )

  ggplot(x, aes(y = Publication)) +
    geom_segment(
      aes(
        x = before,
        xend = after,
        yend = Publication,
        colour = Outcome
      ),
      linewidth = 0.55
    ) +
    geom_point(
      aes(x = before),
      shape = 1,
      colour = cols[["dark"]],
      size = 2
    ) +
    geom_point(
      aes(
        x = after,
        colour = Outcome
      ),
      size = 2
    ) +
    scale_colour_manual(
      values = c(
        Improved = cols[["blue"]],
        Worsened = cols[["red"]],
        Unchanged = cols[["grey"]]
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_title,
      y = "Publication",
      caption = "Open points indicate before calibration. Filled points indicate after calibration."
    ) +
    theme_manuscript() +
    theme(
      axis.text.y = element_text(size = 8)
    )
}

# ------------------------------------------------------------
# 4. Figure 2: source-level RMSE changes
# ------------------------------------------------------------

fig2_sources <- one

fig2_sources$RMSE_decrease <- with(
  fig2_sources,
  RMSE_uncalibrated - RMSE_calibrated
)

fig2_mean_decrease <- mean(fig2_sources$RMSE_decrease)

fig2_mean_label <- paste0(
  "Mean RMSE decrease = ",
  fmt(fig2_mean_decrease, 3)
)

fig2_sources$Outcome <- ifelse(
  fig2_sources$RMSE_decrease > 1e-10,
  "Improved",
  ifelse(
    fig2_sources$RMSE_decrease < -1e-10,
    "Worsened",
    "Unchanged"
  )
)

source_order <- fig2_sources$Publication[
  order(
    fig2_sources$RMSE_decrease,
    fig2_sources$source_id
  )
]

fig2_sources$Publication <- factor(
  fig2_sources$Publication,
  levels = rev(source_order)
)

p2 <- ggplot(
  fig2_sources,
  aes(
    y = Publication,
    colour = Outcome
  )
) +
  geom_vline(
    xintercept = 0,
    colour = "grey60",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = fig2_mean_decrease,
    colour = cols[["teal"]],
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = RMSE_decrease,
      yend = Publication
    ),
    linewidth = 0.6
  ) +
  geom_point(
    aes(x = RMSE_decrease),
    size = 2.2
  ) +
	annotate(
	  "text",
	  x = fig2_mean_decrease + 0.02,
	  y = Inf,
	  label = fig2_mean_label,
	  hjust = 0,
	  vjust = 1.2,
	  size = 3.5,
	  colour = cols[["teal"]]
	) +
  scale_colour_manual(
    values = c(
      Improved = cols[["blue"]],
      Worsened = cols[["red"]],
      Unchanged = cols[["grey"]]
    )
  ) +
  scale_x_continuous(
    breaks = function(limits) pretty(limits, n = 5),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_y_discrete(
    expand = expansion(add = 0.7)
  ) +
  labs(
    title = "One-anchor calibration",
    y = "Source",
    x = expression(
  		"RMSE decrease (" * log[10] * K[p] ~ "units)"
		),
  ) +
  theme_manuscript() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8)
  )

save_fig(
  p2,
  "02_one_anchor_calibration",
  6,
  max(5, 0.07 * nrow(fig2_sources) + 2),

  paste(
    "One-anchor calibration for held-out sources, evaluated on",
    "compounds absent from training.",
    "RMSE decrease is uncalibrated minus calibrated RMSE.",
    "Positive values indicate improvement and negative values",
    "indicate worsening.",
    "Performance was averaged across every available single anchor",
    "within each source.",
    "Horizontal lines connect estimates to zero and are not",
    "confidence intervals.",
    "The teal dashed line shows the mean RMSE decrease across",
    "all sources, giving equal weight to each source.",
    "Sources are ordered by RMSE decrease.",
    "RMSE is expressed in log10 Kp units, with Kp in cm/h.",
    "Source identities are provided in Additional file 1."
  )
)

# ------------------------------------------------------------
# 5. Figure 3: calibration benefit from zero to five anchors
# ------------------------------------------------------------

matched_ids <- unique(
  source_cal$source_id[source_cal$anchors == 5L]
)

fig3_n_sources <- length(matched_ids)

m <- source_cal[
  source_cal$source_id %in% matched_ids,
  , drop = FALSE
]

if (
  fig3_n_sources != paired$n_sources ||
  any(matched_cal$n_sources != fig3_n_sources)
) {
  stop("Fixed subset does not agree with the saved summaries.")
}

panel_counts <- table(
  factor(m$source_id, levels = matched_ids),
  factor(m$anchors, levels = 1:5)
)

if (
  any(panel_counts != 1L) ||
  any(!m$anchors %in% 1:5)
) {
  stop("Expected one result per source and anchor count.")
}

# One row per source, with the same uncalibrated baseline
# used at every anchor count.
baseline <- m[
  m$anchors == 1L,
  c("source_id", "RMSE_uncalibrated"),
  drop = FALSE
]

baseline <- baseline[
  match(matched_ids, baseline$source_id),
  , drop = FALSE
]

baseline_index <- match(m$source_id, baseline$source_id)

if (!close_enough(
  m$RMSE_uncalibrated,
  baseline$RMSE_uncalibrated[baseline_index]
)) {
  stop("Uncalibrated RMSE changes across anchor counts.")
}

# Rows are sources and columns are anchor counts.
# Zero anchors has zero decrease by definition.
anchor_counts <- 0:5

gain_matrix <- vapply(
  anchor_counts,
  function(k) {
    if (k == 0L) {
      return(rep(0, fig3_n_sources))
    }

    x <- m[m$anchors == k, , drop = FALSE]
    j <- match(matched_ids, x$source_id)

    baseline$RMSE_uncalibrated -
      x$RMSE_calibrated[j]
  },
  numeric(fig3_n_sources)
)

rownames(gain_matrix) <- matched_ids
colnames(gain_matrix) <- as.character(anchor_counts)

# Check the additional benefit from one to five anchors.
additional_gain <- gain_matrix[, "5"] - gain_matrix[, "1"]

if (!close_enough(
  mean(additional_gain),
  paired$mean_RMSE_decrease
)) {
  stop("One-to-five decrease disagrees with the paired summary.")
}

# Source-level bootstrap.
# Each resampled source retains its complete trajectory.
fig3_settings <- read_result(
  "81_analysis_settings.csv",
  c("n_boot", "seed")
)

if (nrow(fig3_settings) != 1L) {
  stop("Expected one row of calibration settings.")
}

fig3_n_boot <- as.integer(fig3_settings$n_boot)
fig3_seed <- as.integer(fig3_settings$seed) + 300L

set.seed(fig3_seed)

bootstrap_gain <- replicate(
  fig3_n_boot,
  {
    sampled_sources <- sample.int(
      fig3_n_sources,
      size = fig3_n_sources,
      replace = TRUE
    )

    colMeans(
      gain_matrix[sampled_sources, , drop = FALSE]
    )
  }
)

mean_gain <- data.frame(
  anchors = anchor_counts,
  RMSE_decrease = colMeans(gain_matrix),
  ci_low = apply(
    bootstrap_gain, 1, quantile,
    probs = 0.025,
    names = FALSE
  ),
  ci_high = apply(
    bootstrap_gain, 1, quantile,
    probs = 0.975,
    names = FALSE
  )
)

# Individual source trajectories, including zero anchors.
fig3_trajectories <- data.frame(
  source_id = rep(matched_ids, times = length(anchor_counts)),
  anchors = rep(anchor_counts, each = fig3_n_sources),
  RMSE_decrease = as.vector(gain_matrix),
  stringsAsFactors = FALSE
)

# Red identifies worsening from ONE to FIVE anchors.
worsened_ids <- matched_ids[additional_gain < -1e-10]

fig3_trajectories$worsened_with_five <-
  fig3_trajectories$source_id %in% worsened_ids

fig3_n_improved <- sum(additional_gain > 1e-10)

failure_labels <- fig3_trajectories[
  fig3_trajectories$anchors == 5L &
    fig3_trajectories$worsened_with_five,
  , drop = FALSE
]

failure_labels$source_label <- label(
  failure_labels$source_id
)

p3 <- ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.4
  ) +
  geom_ribbon(
    data = mean_gain,
    aes(
      x = anchors,
      ymin = ci_low,
      ymax = ci_high
    ),
    fill = cols[["yellow"]],
    alpha = 0.7
  ) +
  geom_line(
    data = fig3_trajectories[
      !fig3_trajectories$worsened_with_five,
    ],
    aes(
      x = anchors,
      y = RMSE_decrease,
      group = source_id
    ),
    colour = cols[["blue"]],
    linewidth = 0.55,
    alpha = 0.75
  ) +
  geom_line(
    data = fig3_trajectories[
      fig3_trajectories$worsened_with_five,
    ],
    aes(
      x = anchors,
      y = RMSE_decrease,
      group = source_id
    ),
    colour = cols[["red"]],
    linewidth = 0.75
  ) +
  geom_line(
    data = mean_gain,
    aes(
      x = anchors,
      y = RMSE_decrease
    ),
    colour = cols[["teal"]],
    linewidth = 1.1
  ) +
  geom_point(
    data = mean_gain,
    aes(
      x = anchors,
      y = RMSE_decrease
    ),
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
  scale_x_continuous(
    breaks = 0:5,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0))
  ) +
  coord_cartesian(
    xlim = c(0, 5),
    ylim = c(NA, 1.0),
    clip = "off"
  ) +
  labs(
    title = "Calibration benefit with increasing number of anchors",
    x = "Number of anchor compounds",
    y = expression(
		  "RMSE decrease (" *
		    log[10] * K[p] ~ "units)"
		),
  ) +
  theme_manuscript() +
  theme(
    legend.position = "none",
    plot.margin = margin(12, 48, 12, 12)
  )

save_fig(
  p3,
  "03_anchor_number",
  7,
  4,

  paste(
    "Calibration benefit from zero to five anchors within the fixed subset.",
    "Zero anchors represents uncalibrated predictions.",
    "RMSE decrease is uncalibrated minus calibrated RMSE.",
    "Thin lines show individual sources.",
    "Red lines identify sources that worsened from one to five anchors.",
    "The blue line shows the mean decrease with equal source weighting.",
    "Shading shows pointwise 95% source-bootstrap confidence intervals",
    "relative to no calibration.",
    "The annotated additional decrease and its confidence interval",
    "compare one with five anchors.",
    "Performance averages all possible panels within each source",
    "and anchor count.",
    "The sources and test measurements remain fixed.",
    "RMSE is expressed in log10 Kp units, with Kp in cm/h.",
    "Source identities are provided in Additional file 1."
  )
)

# ------------------------------------------------------------
# 6. Figure 4: two contrasting calibration failures
# ------------------------------------------------------------

f <- failure[failure$anchors == 1L, , drop = FALSE]
j <- match(f$source_id, one$source_id)
if (nrow(f) != nrow(one) || anyNA(j) ||
    !close_enough(f$RMSE_calibrated, one$RMSE_calibrated[j])) {
  stop("Failure results disagree with 81. Rerun 84.")
}

fig4_ids <- example_sources
if (length(fig4_ids) != 2L) stop("Two example sources are required.")

fig4_rows <- do.call(rbind, lapply(fig4_ids, function(id) {
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

fig4_xlim <- range(fig4_rows$offset_low, fig4_rows$offset_high,
                   fig4_rows$test_adjustment, 0)
fig4_xlim <- c(floor(fig4_xlim[1] * 10) / 10,
               ceiling(fig4_xlim[2] * 10) / 10)

fig4_example <- function(id, panel_tag, heading) {
  x <- fig4_rows[fig4_rows$source_id == id, , drop = FALSE]
  target <- unique(x$test_adjustment)
  if (length(target) != 1L) stop("Test benchmark changes across anchor counts.")

  # Two rows make the one-versus-five comparison explicit.
  p_offset <- ggplot(x, aes(y = row_y)) +
    annotate("segment", x = 0, xend = 0,
             y = 0.65, yend = 2.25,
             colour = "grey80", linewidth = 0.5) +
    annotate("segment", x = target, xend = target,
             y = 0.65, yend = 2.25,
             colour = cols[["red"]], linetype = "dashed", linewidth = 0.65) +
    geom_segment(aes(x = offset_low, xend = offset_high, yend = row_y),
                 colour = "grey60", linewidth = 1.0) +
    geom_point(aes(x = offset_mean),
               colour = cols[["teal"]], size = 3.8) +
    scale_x_continuous(limits = fig4_xlim,
                       breaks = seq(-2, 2, by = 0.5),
                       expand = expansion(mult = c(0.035, 0.035))) +
    scale_y_continuous(limits = c(0.5, 2.5),
                       breaks = c(1, 2),
                       labels = c("Five anchors", "One anchor"),
                       expand = expansion(mult = c(0, 0))) +
    labs(title = paste(
			  panel_tag, heading, paste0("(", label(id), ")")),
         x = expression("Estimated adjustment (" * log[10] * K[p] ~ "units)"),
         y = NULL) +
    theme_manuscript() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.margin = margin(8, 8, 8, 8))

  # Prediction error sits beside the corresponding anchor row.
  x$rmse_label <- fmt(x$RMSE, 3)
  x$text_colour <- ifelse(x$anchors == 5L &
                           x$RMSE > x$RMSE[x$anchors == 1L],
                         cols[["red"]], cols[["dark"]])
  p_rmse <- ggplot(x, aes(y = row_y)) +
    geom_text(aes(x = 0.5, label = rmse_label, colour = text_colour),
              size = 4.8, fontface = "bold") +
    scale_colour_identity() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 2.5), expand = c(0, 0)) +
    labs(title = "Mean RMSE", x = " ", y = NULL) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5,
                                    margin = margin(b = 6)),
          axis.title.x = element_text(size = 11, margin = margin(t = 5.5)),
          plot.margin = margin(8, 8, 8, 8))

  p_offset | p_rmse
}

fig4_top <- fig4_example(
  fig4_ids[1],
  "A",
  "Wrong direction"
) +
  plot_layout(widths = c(5, 1))

fig4_bottom <- fig4_example(
  fig4_ids[2],
  "B",
  "Correction too large"
) +
  plot_layout(widths = c(5, 1))

p4 <- (fig4_top / fig4_bottom) +
  plot_annotation(
    title = "Calibration offsets and prediction error",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", margin = margin(b = 10)),
      plot.margin = margin(12, 12, 10, 12)
    )
  )

save_fig(
  p4, "04_calibration_failures", 8, 4.9,
  paste(
    "One- and five-anchor calibration in two sources.",
    "Points show mean estimated offsets and grey lines show full ranges",
    "across all possible panels, not confidence intervals.",
    "Red dashed lines show the mean uncalibrated test residuals, which minimize",
    "test RMSE among constant offsets. These benchmarks were calculated",
    "retrospectively and were not used for calibration.",
    "The numerical column shows mean calibrated RMSE across panels.",
    "Zero on the adjustment axis represents no correction.",
    "Both examples use the same adjustment scale.",
    "Panels can share compounds and are not independent experiments.",
    "Offsets and RMSE are expressed in log10 Kp units, with Kp in cm/h.",
    "Red RMSE values indicate higher error with five anchors than with one.",
    "The estimated correction points in the wrong direction in A and is too large in B.",
    "In both examples, five anchors narrow the offset range but increase mean RMSE."
  )
)



# ------------------------------------------------------------
# Supplementary figure: diagnostics across all eligible sources
# ------------------------------------------------------------

f$Outcome <- ifelse(
  f$RMSE_calibrated > f$RMSE_uncalibrated + 1e-10,
  "Worsened",
  "Improved or unchanged"
)

p4_supp <- ggplot(
  f,
  aes(
    x = test_mean_residual,
    y = offset_mean,
    colour = Outcome
  )
) +
  geom_hline(
    yintercept = 0,
    colour = "grey70",
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = 0,
    colour = "grey70",
    linewidth = 0.35
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = cols[["dark"]]
  ) +
  geom_linerange(
    aes(
      ymin = offset_min,
      ymax = offset_max
    ),
    alpha = 0.5,
    linewidth = 0.45
  ) +
  geom_point(size = 2) +
  scale_colour_manual(
    values = c(
      "Improved or unchanged" = cols[["blue"]],
      Worsened = cols[["red"]]
    )
  ) +
  labs(
    title = "Anchor offsets and retrospective test residuals",
    x = "Mean uncalibrated test residual",
    y = "Mean estimated offset"
  ) +
  theme_manuscript()

save_fig(
  p4_supp,
  "S3_offset_test_bias",
  7,
  4.5,

  paste(
    "One-anchor calibration diagnostics across eligible sources.",
    "Points show mean offsets and vertical lines show full ranges",
    "across all available single anchors.",
    "Ranges are not confidence intervals.",
    "The dashed diagonal indicates equality between the estimated",
    "offset and mean uncalibrated test residual.",
    "Colours indicate the change in mean RMSE across single anchors.",
    "Both axes are expressed in log10 Kp units.",
    "Test residuals were calculated retrospectively.",
    "Mean offsets alone do not determine mean RMSE because",
    "offsets vary among anchor choices."
  )
)

# ------------------------------------------------------------
# 7. Supplementary Figure S1: curation sensitivity
# ------------------------------------------------------------

sensitivity$Policy <- ifelse(
  sensitivity$curation_policy == "priority",
  "Database priority",
  paste0(
    "Agreement ",
    sensitivity$log_tolerance
  )
)

sensitivity$Policy <- factor(
  sensitivity$Policy,
  levels = rev(sensitivity$Policy)
)

forest <- function(
  data, estimate, low, high,
  title, x_title
) {
  data$estimate <- data[[estimate]]
  data$low <- data[[low]]
  data$high <- data[[high]]

  ggplot(data, aes(y = Policy)) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = cols["grey"]
    ) +
    geom_segment(
      aes(
        x = low,
        xend = high,
        yend = Policy
      ),
      colour = cols["blue"],
      linewidth = 0.7
    ) +
    geom_point(
      aes(x = estimate),
      colour = cols["blue"],
      size = 2.5
    ) +
    labs(
      title = title,
      x = x_title,
      y = NULL
    ) +
    theme_manuscript()
}

s1a <- forest(
  sensitivity,
  "one_anchor_reduction_percent",
  "one_anchor_ci_low",
  "one_anchor_ci_high",
  "One-anchor calibration",
  "RMSE reduction (%)"
)

s1b <- forest(
  sensitivity,
  "one_to_five_RMSE_decrease",
  "one_to_five_ci_low",
  "one_to_five_ci_high",
  "Five versus one anchor",
  "RMSE decrease (log10 Kp units)"
)

save_fig(
  (s1a | s1b) +
    plot_annotation(tag_levels = "A"),

  "S1_curation_sensitivity",
  10,
  4,

  paste(
    "Calibration effects and 95% publication-bootstrap intervals",
    "under each curation policy.",
    "Evaluation populations and test measurements",
    "can differ between policies."
  )
)

# ------------------------------------------------------------
# 8. Supplementary Figure S2: variance confidence intervals
# ------------------------------------------------------------

variance$Policy <- factor(
  as.character(variance$component),
  levels = c(
    "Residual",
    "Publication",
    "Compound"
  )
)

p_s2 <- forest(
  variance,
  "percent",
  "percent_ci_low",
  "percent_ci_high",
  "Pooled variance components",
  "Variance component (%)"
)

save_fig(
  p_s2,
  "S2_variance_intervals",
  6.5,
  3.5,

  paste(
    "Pooled variance fractions and 95% parametric-bootstrap",
    "intervals from the intercept-only crossed random-intercept model.",
    "Intervals are marginal and their bounds do not sum to 100%."
  )
)

# ------------------------------------------------------------
# 9. Save caption index
# ------------------------------------------------------------

write.csv(
  do.call(rbind, figure_index),
  file.path(figure_dir, "figure_index.csv"),
  row.names = FALSE
)

cat(
  "\nFour main and two supplementary figures ",
  "written to R/figures/.\n",
  "Captions: R/figures/figure_index.csv\n",
  sep = ""
)
