# ============================================================
# PILOT: Hierarchical model of human skin permeability
#
# logKp_ij = theta_i + u_j + beta_T*T_ij + beta_H*H_ij + error_ij
#
# Input:
#   Table 1.xlsx
#
# Outputs:
#   pilot_model_comparison.csv
#   pilot_variance_components.csv
#   pilot_fixed_effects.csv
#   pilot_adjusted_compound_logKp.csv
#   pilot_reference_effects.csv
#   pilot_observation_predictions.csv
#
#   pilot_variance_components.png
#   pilot_raw_vs_adjusted_compound_logKp.png
#   pilot_adjusted_compound_logKp.png
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "ggplot2",
  "lme4",
  "lmerTest",
  "broom.mixed",
  "performance",
  "forcats",
  "see"
)

to_install <- packages[!packages %in% rownames(installed.packages())]

if (length(to_install) > 0) {
  install.packages(to_install)
}

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(performance)
library(forcats)


# ------------------------------------------------------------
# 1. Read raw Cheruvu table
# ------------------------------------------------------------

dat_raw <- read_excel("Table 1.xlsx")

names(dat_raw)

# Rename only the variables needed for this pilot
dat <- dat_raw %>%
  transmute(
    row_id = No,
    compound = as.factor(Compound),
    cas = `CAS No`,
    reference = as.factor(Reference),
    temperature_K = as.numeric(Texpi),
    skin_thickness_mm = as.numeric(`Skin thicknessj`),
    logKp = as.numeric(logkpl)
  ) %>%
  filter(
    !is.na(compound),
    !is.na(reference),
    !is.na(logKp),
    !is.na(temperature_K),
    !is.na(skin_thickness_mm)
  )

cat("Observations:", nrow(dat), "\n")
cat("Compounds:", n_distinct(dat$compound), "\n")
cat("References:", n_distinct(dat$reference), "\n")


# ------------------------------------------------------------
# 2. Center experimental variables
#
# Centering makes the intercept / compound estimates correspond
# to a realistic standardized experimental condition.
#
# Here I use dataset means for the pilot.
# Later we can choose scientifically meaningful values instead.
# ------------------------------------------------------------

T0 <- mean(dat$temperature_K, na.rm = TRUE)
H0 <- mean(dat$skin_thickness_mm, na.rm = TRUE)

cat("Reference temperature:", T0, "K\n")
cat("Reference skin thickness:", H0, "mm\n")

dat <- dat %>%
  mutate(
    temperature_c = temperature_K - T0,
    thickness_c = skin_thickness_mm - H0
  )


# ------------------------------------------------------------
# 3. Fit nested pilot models
# ------------------------------------------------------------

# Model 0:
# compound variation only
m0 <- lmer(
  logKp ~ 1 +
    (1 | compound),
  data = dat,
  REML = TRUE
)


# Model 1:
# compound + reference
m1 <- lmer(
  logKp ~ 1 +
    (1 | compound) +
    (1 | reference),
  data = dat,
  REML = TRUE
)


# Model 2:
# compound + reference + experimental conditions
#
# This is the main pilot model:
#
# logKp_ij =
#   mu +
#   theta_i +
#   u_j +
#   beta_T T_ij +
#   beta_H H_ij +
#   epsilon_ij
#
m2 <- lmer(
  logKp ~
    temperature_c +
    thickness_c +
    (1 | compound) +
    (1 | reference),
  data = dat,
  REML = TRUE
)


summary(m2)


# ------------------------------------------------------------
# 4. Model comparison
#
# Refit using ML for likelihood-based model comparisons.
# ------------------------------------------------------------

m0_ml <- update(m0, REML = FALSE)
m1_ml <- update(m1, REML = FALSE)
m2_ml <- update(m2, REML = FALSE)

comparison <- tibble(
  model = c(
    "Compound only",
    "Compound + reference",
    "Compound + reference + temperature + thickness"
  ),
  AIC = c(AIC(m0_ml), AIC(m1_ml), AIC(m2_ml)),
  BIC = c(BIC(m0_ml), BIC(m1_ml), BIC(m2_ml)),
  logLik = c(
    as.numeric(logLik(m0_ml)),
    as.numeric(logLik(m1_ml)),
    as.numeric(logLik(m2_ml))
  )
)

write.csv(
  comparison,
  "pilot_model_comparison.csv",
  row.names = FALSE
)

print(comparison)


# ------------------------------------------------------------
# 5. Variance decomposition
# ------------------------------------------------------------

vc <- as.data.frame(VarCorr(m2))

variance_components <- vc %>%
  transmute(
    component = case_when(
      grp == "compound" ~ "Compound",
      grp == "reference" ~ "Reference",
      grp == "Residual" ~ "Residual",
      TRUE ~ grp
    ),
    variance = vcov,
    SD = sdcor
  ) %>%
  mutate(
    proportion_total_variance =
      variance / sum(variance)
  )

write.csv(
  variance_components,
  "pilot_variance_components.csv",
  row.names = FALSE
)

print(variance_components)


# ------------------------------------------------------------
# 6. ICC-style variance proportions
# ------------------------------------------------------------

var_compound <- variance_components %>%
  filter(component == "Compound") %>%
  pull(variance)

var_reference <- variance_components %>%
  filter(component == "Reference") %>%
  pull(variance)

var_residual <- variance_components %>%
  filter(component == "Residual") %>%
  pull(variance)

var_total <- var_compound + var_reference + var_residual

icc_summary <- tibble(
  quantity = c(
    "Compound variance proportion",
    "Reference variance proportion",
    "Residual variance proportion"
  ),
  estimate = c(
    var_compound / var_total,
    var_reference / var_total,
    var_residual / var_total
  )
)

write.csv(
  icc_summary,
  "pilot_variance_proportions.csv",
  row.names = FALSE
)

print(icc_summary)


# ------------------------------------------------------------
# 7. Fixed experimental effects
# ------------------------------------------------------------

fixed_effects <- broom.mixed::tidy(
  m2,
  effects = "fixed",
  conf.int = TRUE
)

write.csv(
  fixed_effects,
  "pilot_fixed_effects.csv",
  row.names = FALSE
)

print(fixed_effects)


# ------------------------------------------------------------
# 8. Estimate adjusted compound-level permeability
#
# For compound i:
#
# theta_i_hat =
#    overall intercept +
#    compound random intercept
#
# Because temperature and thickness are centered,
# this corresponds to permeability at:
#
# temperature = T0
# thickness   = H0
# reference effect = 0
#
# ------------------------------------------------------------

compound_re <- ranef(m2, condVar = TRUE)$compound

compound_adjusted <- tibble(
  compound = rownames(compound_re),
  compound_random_effect = compound_re[, "(Intercept)"]
)

intercept <- fixef(m2)[["(Intercept)"]]

compound_adjusted <- compound_adjusted %>%
  mutate(
    adjusted_logKp = intercept + compound_random_effect
  )


# ------------------------------------------------------------
# 9. Add raw compound summaries
# ------------------------------------------------------------

raw_compound_summary <- dat %>%
  group_by(compound) %>%
  summarise(
    n_observations = n(),
    n_references = n_distinct(reference),
    raw_mean_logKp = mean(logKp),
    raw_sd_logKp = sd(logKp),
    raw_min_logKp = min(logKp),
    raw_max_logKp = max(logKp),
    raw_range_logKp = raw_max_logKp - raw_min_logKp,
    .groups = "drop"
  )

compound_adjusted <- compound_adjusted %>%
  left_join(
    raw_compound_summary,
    by = "compound"
  ) %>%
  arrange(adjusted_logKp)

write.csv(
  compound_adjusted,
  "pilot_adjusted_compound_logKp.csv",
  row.names = FALSE
)

print(head(compound_adjusted, 20))


# ------------------------------------------------------------
# 10. Reference effects
# ------------------------------------------------------------

reference_re <- ranef(m2)$reference

reference_effects <- tibble(
  reference = rownames(reference_re),
  reference_effect = reference_re[, "(Intercept)"]
) %>%
  arrange(reference_effect)

write.csv(
  reference_effects,
  "pilot_reference_effects.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Observation-level fitted values and residuals
# ------------------------------------------------------------

dat_output <- dat %>%
  mutate(
    fitted_logKp = fitted(m2),
    residual = resid(m2),
    abs_residual = abs(residual)
  )

write.csv(
  dat_output,
  "pilot_observation_predictions.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Figure 1: Variance decomposition
# ------------------------------------------------------------

p_var <- variance_components %>%
  ggplot(
    aes(
      x = fct_reorder(component, proportion_total_variance),
      y = proportion_total_variance
    )
  ) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = NULL,
    y = "Proportion of total variance",
    title = "Variance decomposition of reported skin permeability",
    subtitle =
      "Crossed mixed-effects model with compound and literature-reference random effects"
  ) +
  theme_classic(base_size = 12)

ggsave(
  "pilot_variance_components.png",
  p_var,
  width = 7,
  height = 4.5,
  dpi = 300
)


# ------------------------------------------------------------
# 13. Figure 2: Raw mean versus context-adjusted permeability
# ------------------------------------------------------------

p_compare <- compound_adjusted %>%
  ggplot(
    aes(
      x = raw_mean_logKp,
      y = adjusted_logKp
    )
  ) +
  geom_point(
    aes(size = n_observations),
    alpha = 0.65
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  labs(
    x = "Raw mean logKp",
    y = "Reference- and condition-adjusted logKp",
    size = "Observations",
    title = "Raw versus context-adjusted compound permeability",
    subtitle =
      "Adjusted values account for reference, temperature, and skin thickness"
  ) +
  theme_classic(base_size = 12)

ggsave(
  "pilot_raw_vs_adjusted_compound_logKp.png",
  p_compare,
  width = 6.5,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------
# 14. Figure 3:
# adjusted permeability for compounds with >= 3 measurements
#
# Useful sanity-check figure.
# ------------------------------------------------------------

plot_compounds <- compound_adjusted %>%
  filter(n_observations >= 3) %>%
  arrange(adjusted_logKp) %>%
  mutate(
    compound = factor(
      compound,
      levels = compound
    )
  )

p_compounds <- ggplot(
  plot_compounds,
  aes(
    x = adjusted_logKp,
    y = compound
  )
) +
  geom_point() +
  labs(
    x = "Adjusted logKp",
    y = NULL,
    title = "Context-adjusted compound permeability",
    subtitle =
      "Compounds with at least three permeability observations"
  ) +
  theme_classic(base_size = 10)

ggsave(
  "pilot_adjusted_compound_logKp.png",
  p_compounds,
  width = 7,
  height = max(6, nrow(plot_compounds) * 0.20),
  dpi = 300,
  limitsize = FALSE
)


# ------------------------------------------------------------
# 15. Model diagnostics
# ------------------------------------------------------------

png(
  "pilot_model_diagnostics.png",
  width = 1800,
  height = 900,
  res = 180
)

par(mfrow = c(1, 2))

plot(
  fitted(m2),
  resid(m2),
  xlab = "Fitted logKp",
  ylab = "Residual",
  main = "Residuals vs fitted"
)

abline(h = 0, lty = 2)

qqnorm(
  resid(m2),
  main = "Residual Q-Q plot"
)

qqline(
  resid(m2)
)

dev.off()


# ------------------------------------------------------------
# 16. Basic model diagnostics to console
# ------------------------------------------------------------

performance::check_singularity(m2)

performance::r2(m2)

performance::icc(m2)

performance::check_model(m2)


# ------------------------------------------------------------
# 17. Save key pilot metadata
# ------------------------------------------------------------

metadata <- tibble(
  quantity = c(
    "Observations",
    "Compounds",
    "References",
    "Reference temperature K",
    "Reference skin thickness mm"
  ),
  value = c(
    nrow(dat),
    n_distinct(dat$compound),
    n_distinct(dat$reference),
    T0,
    H0
  )
)

write.csv(
  metadata,
  "pilot_metadata.csv",
  row.names = FALSE
)