# ============================================================
# 32_compare_calibration_methods.R
# Conventional QSPR + calibration vs hierarchical + calibration
# ============================================================

library(lme4)

dat <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)
dat$compound_id <- factor(dat$compound_id)
dat$reference_id <- factor(dat$reference_id)

ctrl <- lmerControl(
	optimizer = "bobyqa",
	optCtrl = list(maxfun = 100000)
)

set.seed(2026)

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

metrics <- function(obs, pred) {
	c(
		RMSE = sqrt(mean((obs - pred)^2)),
		MAE = mean(abs(obs - pred))
	)
}

profile <- function(x) {
	p <- aggregate(
		x[c("mw", "log_kow")],
		list(compound_id = x$compound_id),
		median
	)

	y <- aggregate(
		x$log_kp,
		list(compound_id = x$compound_id),
		mean
	)

	names(y)[2] <- "mean_log_kp"

	merge(
		p,
		y,
		by = "compound_id",
		sort = FALSE
	)
}

scale_data <- function(train, test = NULL) {
	p <- profile(train)

	center <- sapply(
		p[c("mw", "log_kow")],
		mean
	)

	scale <- sapply(
		p[c("mw", "log_kow")],
		sd
	)

	for (v in c("mw", "log_kow")) {
		p[[paste0(v, "_z")]] <- (p[[v]] - center[v]) / scale[v]
		train[[paste0(v, "_z")]] <- (train[[v]] - center[v]) / scale[v]

		if (!is.null(test)) {
			test[[paste0(v, "_z")]] <- (test[[v]] - center[v]) / scale[v]
		}
	}

	list(
		train = train,
		test = test,
		profile = p
	)
}

source_variance <- function(x) {
	nj <- table(x$reference_id)
	N <- nrow(x)
	J <- length(nj)

	if (J < 2 || N <= J) {
		stop("Insufficient source-compound residuals for variance estimation")
	}

	source_mean <- aggregate(
		residual ~ reference_id,
		data = x,
		mean
	)

	source_n <- aggregate(
		residual ~ reference_id,
		data = x,
		length
	)

	names(source_n)[2] <- "n"

	source_summary <- merge(
		source_mean,
		source_n,
		by = "reference_id",
		sort = FALSE
	)

	grand <- weighted.mean(
		source_summary$residual,
		source_summary$n
	)

	ss_between <- sum(
		source_summary$n *
		(source_summary$residual - grand)^2
	)

	x <- merge(
		x,
		source_mean,
		by = "reference_id",
		suffixes = c("", "_source"),
		sort = FALSE
	)

	ss_within <- sum(
		(x$residual - x$residual_source)^2
	)

	ms_between <- ss_between / (J - 1)
	ms_within <- ss_within / (N - J)

	n0 <- (
		N - sum(source_summary$n^2) / N
	) / (J - 1)

	tau2 <- max(
		0,
		(ms_between - ms_within) / n0
	)

	list(
		source_variance = tau2,
		residual_variance = ms_within
	)
}

shrink <- function(k, tau2, sigma2) {
	denom <- sigma2 + k * tau2

	if (denom <= 0) {
		return(0)
	}

	k * tau2 / denom
}

# ------------------------------------------------------------
# Held-out-source calibration comparison
# ------------------------------------------------------------

references <- unique(as.character(dat$reference_id))
anchor_sizes <- 1:5

out <- list()
z <- 1

for (ref in references) {
	train <- dat[as.character(dat$reference_id) != ref, ]
	held <- dat[as.character(dat$reference_id) == ref, ]

	train_compounds <- unique(as.character(train$compound_id))
	shared <- intersect(
		unique(as.character(held$compound_id)),
		train_compounds
	)
	novel <- setdiff(
		unique(as.character(held$compound_id)),
		train_compounds
	)

	if (!length(shared) || !length(novel)) {
		next
	}

	sc <- scale_data(train, held)
	train <- sc$train
	held <- sc$test
	p <- sc$profile

	# --------------------------------------------------------
	# Conventional compound-aggregate QSPR
	# --------------------------------------------------------

	lm_fit <- lm(
		mean_log_kp ~ mw_z + log_kow_z,
		data = p
	)

	train$expected_conventional <- as.numeric(
		predict(
			lm_fit,
			newdata = train
		)
	)

	train$residual_conventional <- (
		train$log_kp -
		train$expected_conventional
	)

	train_res_conventional <- aggregate(
		residual_conventional ~ reference_id + compound_id,
		data = train,
		mean
	)

	names(train_res_conventional)[3] <- "residual"

	var_conventional <- source_variance(
		train_res_conventional
	)

	# --------------------------------------------------------
	# Hierarchical QSPR
	# --------------------------------------------------------

	mix_fit <- lmer(
		log_kp ~ mw_z + log_kow_z +
			(1 | compound_id) +
			(1 | reference_id),
		data = train,
		REML = TRUE,
		control = ctrl
	)

	train$expected_hierarchical <- predict(
		mix_fit,
		newdata = train,
		re.form = ~(1 | compound_id),
		allow.new.levels = TRUE
	)

	train$residual_hierarchical <- (
		train$log_kp -
		train$expected_hierarchical
	)

	train_res_hierarchical <- aggregate(
		residual_hierarchical ~ reference_id + compound_id,
		data = train,
		mean
	)

	names(train_res_hierarchical)[3] <- "residual"

	var_hierarchical <- source_variance(
		train_res_hierarchical
	)

	# --------------------------------------------------------
	# Held-out anchors and novel compounds
	# --------------------------------------------------------

	shared_rows <- held[
		as.character(held$compound_id) %in% shared,
	]

	novel_rows <- held[
		as.character(held$compound_id) %in% novel,
	]

	shared_rows$expected_conventional <- as.numeric(
		predict(
			lm_fit,
			newdata = shared_rows
		)
	)

	novel_rows$expected_conventional <- as.numeric(
		predict(
			lm_fit,
			newdata = novel_rows
		)
	)

	shared_rows$expected_hierarchical <- predict(
		mix_fit,
		newdata = shared_rows,
		re.form = ~(1 | compound_id),
		allow.new.levels = TRUE
	)

	novel_rows$expected_hierarchical <- predict(
		mix_fit,
		newdata = novel_rows,
		re.form = NA,
		allow.new.levels = TRUE
	)

	anchor_conventional <- aggregate(
		shared_rows$log_kp -
			shared_rows$expected_conventional,
		list(compound_id = shared_rows$compound_id),
		mean
	)

	anchor_hierarchical <- aggregate(
		shared_rows$log_kp -
			shared_rows$expected_hierarchical,
		list(compound_id = shared_rows$compound_id),
		mean
	)

	names(anchor_conventional)[2] <- "residual"
	names(anchor_hierarchical)[2] <- "residual"

	# --------------------------------------------------------
	# Same anchor panels for both approaches
	# --------------------------------------------------------

	for (k in anchor_sizes) {
		if (length(shared) < k) {
			next
		}

		shrink_conventional <- shrink(
			k,
			var_conventional$source_variance,
			var_conventional$residual_variance
		)

		shrink_hierarchical <- shrink(
			k,
			var_hierarchical$source_variance,
			var_hierarchical$residual_variance
		)

		for (rep in 1:100) {
			anchors <- sample(
				shared,
				k
			)

			raw_conventional <- mean(
				anchor_conventional$residual[
					as.character(anchor_conventional$compound_id) %in% anchors
				]
			)

			raw_hierarchical <- mean(
				anchor_hierarchical$residual[
					as.character(anchor_hierarchical$compound_id) %in% anchors
				]
			)

			offset_conventional <- (
				shrink_conventional *
				raw_conventional
			)

			offset_hierarchical <- (
				shrink_hierarchical *
				raw_hierarchical
			)

			uncal_conventional <- metrics(
				novel_rows$log_kp,
				novel_rows$expected_conventional
			)

			cal_conventional <- metrics(
				novel_rows$log_kp,
				novel_rows$expected_conventional +
					offset_conventional
			)

			uncal_hierarchical <- metrics(
				novel_rows$log_kp,
				novel_rows$expected_hierarchical
			)

			cal_hierarchical <- metrics(
				novel_rows$log_kp,
				novel_rows$expected_hierarchical +
					offset_hierarchical
			)

			out[[z]] <- data.frame(
				reference_id = ref,
				anchors = k,
				replicate = rep,
				method = "Conventional",
				n_test = nrow(novel_rows),
				source_variance = var_conventional$source_variance,
				residual_variance = var_conventional$residual_variance,
				shrinkage = shrink_conventional,
				offset_raw = raw_conventional,
				offset = offset_conventional,
				RMSE_uncalibrated = uncal_conventional["RMSE"],
				RMSE_calibrated = cal_conventional["RMSE"],
				MAE_uncalibrated = uncal_conventional["MAE"],
				MAE_calibrated = cal_conventional["MAE"]
			)

			z <- z + 1

			out[[z]] <- data.frame(
				reference_id = ref,
				anchors = k,
				replicate = rep,
				method = "Hierarchical",
				n_test = nrow(novel_rows),
				source_variance = var_hierarchical$source_variance,
				residual_variance = var_hierarchical$residual_variance,
				shrinkage = shrink_hierarchical,
				offset_raw = raw_hierarchical,
				offset = offset_hierarchical,
				RMSE_uncalibrated = uncal_hierarchical["RMSE"],
				RMSE_calibrated = cal_hierarchical["RMSE"],
				MAE_uncalibrated = uncal_hierarchical["MAE"],
				MAE_calibrated = cal_hierarchical["MAE"]
			)

			z <- z + 1
		}
	}
}

performance <- do.call(
	rbind,
	out
)

# ------------------------------------------------------------
# Source-level results
# Average random anchor selections within each source
# ------------------------------------------------------------

reference <- aggregate(
	cbind(
		RMSE_uncalibrated,
		RMSE_calibrated,
		MAE_uncalibrated,
		MAE_calibrated
	) ~ reference_id + anchors + method,
	data = performance,
	mean
)

reference$RMSE_improvement_percent <- 100 * (
	1 -
	reference$RMSE_calibrated /
	reference$RMSE_uncalibrated
)

reference$MAE_improvement_percent <- 100 * (
	1 -
	reference$MAE_calibrated /
	reference$MAE_uncalibrated
)

# ------------------------------------------------------------
# All eligible sources at each anchor count
# ------------------------------------------------------------

summary_all <- aggregate(
	cbind(
		RMSE_uncalibrated,
		RMSE_calibrated,
		MAE_uncalibrated,
		MAE_calibrated
	) ~ anchors + method,
	data = reference,
	mean
)

n_all <- aggregate(
	reference_id ~ anchors + method,
	data = reference,
	function(x) length(unique(x))
)

names(n_all)[3] <- "n_references"

summary_all <- merge(
	summary_all,
	n_all,
	by = c("anchors", "method"),
	sort = FALSE
)

summary_all$RMSE_improvement_percent <- 100 * (
	1 -
	summary_all$RMSE_calibrated /
	summary_all$RMSE_uncalibrated
)

summary_all$MAE_improvement_percent <- 100 * (
	1 -
	summary_all$MAE_calibrated /
	summary_all$MAE_uncalibrated
)

summary_all$cohort <- "All eligible sources"

# ------------------------------------------------------------
# Matched sources supporting all five anchor counts
# ------------------------------------------------------------

matched_refs <- unique(
	reference$reference_id[
		reference$anchors == max(anchor_sizes)
	]
)

matched <- reference[
	reference$reference_id %in% matched_refs,
]

summary_matched <- aggregate(
	cbind(
		RMSE_uncalibrated,
		RMSE_calibrated,
		MAE_uncalibrated,
		MAE_calibrated
	) ~ anchors + method,
	data = matched,
	mean
)

summary_matched$n_references <- length(matched_refs)

summary_matched$RMSE_improvement_percent <- 100 * (
	1 -
	summary_matched$RMSE_calibrated /
	summary_matched$RMSE_uncalibrated
)

summary_matched$MAE_improvement_percent <- 100 * (
	1 -
	summary_matched$MAE_calibrated /
	summary_matched$MAE_uncalibrated
)

summary_matched$cohort <- "Matched 5-anchor-capable cohort"

# ------------------------------------------------------------
# Combined summary
# ------------------------------------------------------------

summary <- rbind(
	summary_all,
	summary_matched
)

summary <- summary[
	order(
		summary$cohort,
		summary$anchors,
		summary$method
	),
]

# ------------------------------------------------------------
# Direct paired method comparison
# Positive value favors hierarchical calibration
# ------------------------------------------------------------

conv <- reference[
	reference$method == "Conventional",
	c(
		"reference_id",
		"anchors",
		"RMSE_calibrated",
		"MAE_calibrated"
	)
]

hier <- reference[
	reference$method == "Hierarchical",
	c(
		"reference_id",
		"anchors",
		"RMSE_calibrated",
		"MAE_calibrated"
	)
]

names(conv)[3:4] <- c(
	"RMSE_conventional",
	"MAE_conventional"
)

names(hier)[3:4] <- c(
	"RMSE_hierarchical",
	"MAE_hierarchical"
)

comparison <- merge(
	conv,
	hier,
	by = c("reference_id", "anchors"),
	sort = FALSE
)

comparison$RMSE_conventional_minus_hierarchical <- (
	comparison$RMSE_conventional -
	comparison$RMSE_hierarchical
)

comparison$MAE_conventional_minus_hierarchical <- (
	comparison$MAE_conventional -
	comparison$MAE_hierarchical
)

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------

write.csv(
	performance,
	"R/results/32_calibration_method_performance.csv",
	row.names = FALSE
)

write.csv(
	reference,
	"R/results/32_calibration_method_reference.csv",
	row.names = FALSE
)

write.csv(
	summary,
	"R/results/32_calibration_method_summary.csv",
	row.names = FALSE
)

write.csv(
	comparison,
	"R/results/32_calibration_method_comparison.csv",
	row.names = FALSE
)

# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

primary <- summary[
	summary$cohort == "All eligible sources" &
	summary$anchors == 1,
]

matched_end <- summary[
	summary$cohort == "Matched 5-anchor-capable cohort" &
	summary$anchors %in% c(1, 5),
]

cat(
	"\nCalibration-method comparison complete\n\n",
	"Primary one-anchor deployment cohort\n",
	paste(
		primary$method,
		": ",
		round(primary$RMSE_uncalibrated, 3),
		" -> ",
		round(primary$RMSE_calibrated, 3),
		"; ",
		round(primary$RMSE_improvement_percent, 1),
		"% improvement across ",
		primary$n_references,
		" sources",
		sep = "",
		collapse = "\n"
	),
	"\n\nMatched 5-anchor-capable cohort\n",
	paste(
		matched_end$method,
		", ",
		matched_end$anchors,
		" anchor",
		ifelse(matched_end$anchors == 1, "", "s"),
		": ",
		round(matched_end$RMSE_improvement_percent, 1),
		"% improvement",
		sep = "",
		collapse = "\n"
	),
	"\n\nPositive conventional-minus-hierarchical RMSE favors hierarchical calibration.\n",
	sep = ""
)