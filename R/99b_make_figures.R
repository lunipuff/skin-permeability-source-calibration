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
  "DOI:10.1016/S0928-0987(03)00176-3",
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
# 3. Figure 1: within-compound variation and variance components
# ------------------------------------------------------------

# Publication colours are assigned within each compound.
# They do not identify the same publication across compounds.

source_key <- unique(
  points[c("compound_id", "publication_id")]
)

source_key$colour_index <- NA_integer_

for (ii in split(
  seq_len(nrow(source_key)),
  source_key$compound_id
)) {
  ii <- ii[order(source_key$publication_id[ii])]

  source_key$colour_index[ii] <- if (length(ii) == 1) {
    1L
  } else {
    round(seq(1, 100, length.out = length(ii)))
  }
}

j <- match(
  paste(
    points$compound_id,
    points$publication_id,
    sep = "\r"
  ),
  paste(
    source_key$compound_id,
    source_key$publication_id,
    sep = "\r"
  )
)

points$colour_index <- factor(
  source_key$colour_index[j],
  levels = 1:100
)

ranges <- aggregate(
  centered_log_kp ~ compound_order,
  points,
  function(z) c(min = min(z), max = max(z))
)

ranges$low <- ranges$centered_log_kp[, "min"]
ranges$high <- ranges$centered_log_kp[, "max"]

p1a <- ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = cols["grey"]
  ) +
  geom_linerange(
    data = ranges,
    aes(
      x = compound_order,
      ymin = low,
      ymax = high
    ),
    colour = cols["grey"],
    linewidth = 0.35
  ) +
  geom_point(
    data = points,
    aes(
      compound_order,
      centered_log_kp,
      colour = colour_index
    ),
    size = 1.3,
    alpha = 0.8,
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 2026
    )
  ) +
  scale_colour_manual(
    values = setNames(
      palette(100),
      as.character(1:100)
    ),
    drop = FALSE
  ) +
  labs(
    title = "Within-compound variation",
    subtitle = paste(
      length(unique(points$compound_id)),
      "compounds with repeated measurements"
    ),
    x = "Compounds ordered by measurement range",
    y = expression(
      log[10] * K[p] ~ "minus compound mean"
    )
  ) +
  theme_manuscript() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

variance$component <- factor(
  variance$component,
  levels = c("Residual", "Publication", "Compound")
)

variance$label <- paste0(
  fmt(variance$percent, 1),
  "%"
)

variance$text_colour <- ifelse(
  variance$component == "Compound",
  "white",
  "black"
)

p1b <- ggplot(
  variance,
  aes(
    x = percent,
    y = "Pooled",
    fill = component
  )
) +
  geom_col(width = 0.55) +
  geom_text(
    aes(
      label = label,
      colour = text_colour
    ),
    position = position_stack(vjust = 0.5),
    size = 4,
    fontface = "bold"
  ) +
  scale_colour_identity() +
  scale_fill_manual(
    values = c(
      Compound = cols[["blue"]],
      Publication = cols[["teal"]],
      Residual = cols[["yellow"]]
    ),
    breaks = c(
      "Compound",
      "Publication",
      "Residual"
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, 20),
    expand = c(0, 0)
  ) +
  labs(
    title = "Pooled variance decomposition",
    x = "Variance component (%)",
    y = NULL
  ) +
  theme_manuscript() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

save_fig(
  p1a / p1b +
    plot_layout(heights = c(2.3, 1)) +
    plot_annotation(tag_levels = "A"),

  "01_variance_decomposition",
  9,
  7.5,

  paste(
    "Within-compound variation and pooled variance decomposition.",
    "A: Measurements are centred on each compound's mean;",
    "colours distinguish publications within a compound,",
    "not consistently across compounds. Lines show measurement ranges.",
    "B: REML intercept-only model with compound and publication",
    "random intercepts, without molecular descriptors or temperature.",
    "Bootstrap intervals are shown in Supplementary Figure S2",
    "and Table S3. Publication variance does not identify",
    "experimental causes."
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
      colour = cols["dark"],
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
      caption = "Open points: before. Filled points: after."
    ) +
    theme_manuscript() +
    theme(
      axis.text.y = element_text(size = 8)
    )
}

# ------------------------------------------------------------
# 4. Figure 2: one-anchor calibration
# ------------------------------------------------------------

subtitle2 <- paste0(
  "Mean RMSE: ",
  fmt(overall_one$mean_source_RMSE_uncalibrated),
  " to ",
  fmt(overall_one$mean_source_RMSE_calibrated),
  "; ",
  overall_one$n_sources_improved,
  "/",
  overall_one$n_sources,
  " publications improved\n",
  "Reduction: ",
  fmt(overall_one$RMSE_reduction_percent, 1),
  "% (95% CI ",
  fmt(overall_one$reduction_ci_low, 1),
  " to ",
  fmt(overall_one$reduction_ci_high, 1),
  "%)"
)

p2 <- paired_plot(
  one,
  "RMSE_uncalibrated",
  "RMSE_calibrated",
  "One-anchor calibration",
  subtitle2,
  "RMSE (log10 Kp units)"
) +
  labs(
    caption = "Open: uncalibrated. Filled: calibrated."
  )

save_fig(
  p2,
  "02_one_anchor_calibration",
  8.5,
  max(7, 0.18 * nrow(one) + 2),

  paste(
    "One-anchor calibration by held-out publication.",
    "Open points are uncalibrated RMSE and filled points",
    "calibrated RMSE, averaged over the saved random anchor selections.",
    "Red indicates increased mean RMSE.",
    "Publications are ordered by RMSE change.",
    "The overall percentage reduction is calculated from",
    "mean publication-level RMSE; its interval uses",
    "publication bootstrap resampling.",
    "Publication labels are defined in Supplementary Table S1."
  )
)

# ------------------------------------------------------------
# 5. Figure 3: fixed-subset anchor-number comparison
# ------------------------------------------------------------

matched_ids <- source_cal$source_id[
  source_cal$anchors == 5
]

m <- source_cal[
  source_cal$source_id %in% matched_ids,
  ,
  drop = FALSE
]

if (
  length(matched_ids) != paired$n_sources ||
  any(matched_cal$n_sources != length(matched_ids))
) {
  stop("Matched cohort mismatch.")
}

for (k in matched_cal$anchors) {
  x <- m[m$anchors == k, ]

  expected <- 100 * (
    1 -
      mean(x$RMSE_calibrated) /
      mean(x$RMSE_uncalibrated)
  )

  if (!close_enough(
    expected,
    matched_cal$RMSE_reduction_percent[
      matched_cal$anchors == k
    ]
  )) {
    stop("Matched summaries disagree.")
  }
}

p3a <- ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = cols["grey"]
  ) +
  geom_line(
    data = m,
    aes(
      anchors,
      RMSE_reduction_percent,
      group = source_id
    ),
    colour = "grey65",
    linewidth = 0.5,
    alpha = 0.8
  ) +
  geom_ribbon(
    data = matched_cal,
    aes(
      anchors,
      ymin = reduction_ci_low,
      ymax = reduction_ci_high
    ),
    fill = cols["blue"],
    alpha = 0.15
  ) +
  geom_line(
    data = matched_cal,
    aes(
      anchors,
      RMSE_reduction_percent
    ),
    colour = cols["blue"],
    linewidth = 1
  ) +
  geom_point(
    data = matched_cal,
    aes(
      anchors,
      RMSE_reduction_percent
    ),
    colour = cols["blue"],
    size = 2
  ) +
  scale_x_continuous(breaks = 1:5) +
  labs(
    title = "Anchor-number comparison",
    subtitle = paste(
      length(matched_ids),
      "publications at every anchor count"
    ),
    x = "Anchor compounds",
    y = "RMSE reduction (%)",
    caption = paste(
      "Grey: individual publications.",
      "Blue: overall estimate and 95% CI."
    )
  ) +
  theme_manuscript()

comparison <- merge(
  m[
    m$anchors == 1,
    c("source_id", "Publication", "RMSE_calibrated")
  ],
  m[
    m$anchors == 5,
    c("source_id", "RMSE_calibrated")
  ],
  by = "source_id",
  suffixes = c("_one", "_five")
)

if (!close_enough(
  mean(
    comparison$RMSE_calibrated_one -
      comparison$RMSE_calibrated_five
  ),
  paired$mean_RMSE_decrease
)) {
  stop("Paired anchor-number results disagree.")
}

p3b <- paired_plot(
  comparison,
  "RMSE_calibrated_one",
  "RMSE_calibrated_five",
  "Five versus one anchor",

  paste0(
    "Mean decrease: ",
    fmt(paired$mean_RMSE_decrease),
    "\n95% CI ",
    fmt(paired$decrease_ci_low),
    " to ",
    fmt(paired$decrease_ci_high)
  ),

  "RMSE (log10 Kp units)"
) +
  labs(
    caption = "Open: one anchor. Filled: five anchors."
  )

save_fig(
  (p3a | p3b) +
    plot_layout(widths = c(1.2, 1)) +
    plot_annotation(tag_levels = "A"),

  "03_anchor_number",
  11,
  6.8,

  paste(
    "Anchor-number comparison in the fixed subset supporting five anchors.",
    "A: Individual-publication percentage reductions and",
    "the percentage reduction in mean publication RMSE.",
    "The overall curve is not the mean of individual percentage reductions.",
    "Shading is a 95% publication-bootstrap interval.",
    "B: Paired calibrated RMSE with one and five anchors;",
    "red indicates higher RMSE with five.",
    "The same publications and test measurements are used.",
    "Additional anchors also change the shrinkage factor."
  )
)

# ------------------------------------------------------------
# 6. Figure 4: calibration failure diagnostics
# ------------------------------------------------------------

f <- failure[
  failure$anchors == 1,
  ,
  drop = FALSE
]

j <- match(f$source_id, one$source_id)

if (
  nrow(f) != nrow(one) ||
  anyNA(j) ||
  !close_enough(
    f$RMSE_calibrated,
    one$RMSE_calibrated[j]
  )
) {
  stop("Failure results disagree with 81. Rerun 84.")
}

f$Outcome <- ifelse(
  f$RMSE_calibrated >
    f$RMSE_uncalibrated + 1e-10,
  "Worsened",
  "Improved or unchanged"
)

highlight <- f[
  f$source_id %in% example_sources,
]

p4a <- ggplot(
  f,
  aes(
    test_mean_residual,
    offset_mean,
    colour = Outcome
  )
) +
  geom_hline(
    yintercept = 0,
    colour = cols["grey"],
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = 0,
    colour = cols["grey"],
    linewidth = 0.35
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = cols["dark"]
  ) +
  geom_linerange(
    aes(
      ymin = offset_min,
      ymax = offset_max
    ),
    alpha = 0.55,
    linewidth = 0.45
  ) +
  geom_point(size = 2) +
  geom_text(
    data = highlight,
    aes(label = Publication),
    nudge_x = 0.06,
    nudge_y = 0.07,
    hjust = 0,
    size = 3,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "Improved or unchanged" = cols[["blue"]],
      Worsened = cols[["red"]]
    )
  ) +
  labs(
    title = "Anchor offsets and test bias",
    x = "Mean uncalibrated test residual (retrospective)",
    y = "Mean estimated offset",
    caption = paste(
      "Vertical lines: sampled offset ranges,",
      "not confidence intervals."
    )
  ) +
  theme_manuscript()

offset_plot <- function(id) {
  x <- panels[
    panels$source_id == id,
    ,
    drop = FALSE
  ]

  z <- f[
    f$source_id == id,
    ,
    drop = FALSE
  ]

  if (!nrow(x) || nrow(z) != 1) {
    stop("Missing example panels: ", id)
  }

  x$anchors <- factor(
    x$anchors,
    levels = sort(unique(x$anchors))
  )

  ggplot(x, aes(anchors, offset)) +
    geom_hline(
      yintercept = 0,
      colour = cols["grey"],
      linewidth = 0.3
    ) +
    geom_hline(
      yintercept = z$test_mean_residual,
      colour = cols["red"],
      linetype = "dashed"
    ) +
    geom_boxplot(
      fill = cols["teal"],
      alpha = 0.5,
      width = 0.6,
      outlier.shape = NA
    ) +
    geom_point(
      size = 0.6,
      alpha = 0.2,
      position = position_jitter(
        width = 0.12,
        height = 0,
        seed = 2026
      )
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      fill = "white",
      size = 2
    ) +
    labs(
      title = paste("Anchor offsets:", label(id)),
      subtitle = paste(
        "Mean test residual:",
        fmt(z$test_mean_residual)
      ),
      x = "Anchor compounds",
      y = "Estimated offset",
      caption = paste(
        "Dashed red: mean test residual.",
        "Diamonds: mean offsets."
      )
    ) +
    theme_manuscript()
}

p4 <- p4a /
  (
    offset_plot(example_sources[1]) |
      offset_plot(example_sources[2])
  ) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(tag_levels = "A")

save_fig(
  p4,
  "04_calibration_failures",
  10,
  9,

  paste(
    "Retrospective calibration diagnostics.",
    "A: Mean offset and sampled range across single-anchor choices",
    "against mean uncalibrated test residual.",
    "Dashed diagonal indicates equality.",
    "Mean offsets alone do not determine mean RMSE",
    "across varying anchor choices.",
    "B-C: Sampled offset distributions for",
    paste(label(example_sources), collapse = " and "),
    "across anchor counts.",
    "Boxplots and point ranges are not confidence intervals,",
    "and sampled panels are not independent experiments.",
    "Test residuals were not used to select anchors",
    "or estimate calibration offsets."
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