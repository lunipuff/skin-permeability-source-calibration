############################################################
# 12_applicability_domain_summary.R
# Applicability-domain summary for selected QSPR model
############################################################

source("R/00_config.R")

############################################################
# Load data and selected model
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

selected_model_text <- readLines(
	path_loco_cv_selected_model,
	warn = FALSE
)

selected_model_text <- selected_model_text[1]

selected_formula <- as.formula(
	selected_model_text
)

############################################################
# User settings
############################################################

outcome_col <- "logkpl"
id_col <- "compound_id"

selected_predictors <- get_formula_predictors(
	selected_formula,
	outcome_col = outcome_col
)

central_probs <- c(
	0.10,
	0.90
)

broad_probs <- c(
	0.025,
	0.975
)

high_error_threshold <- 1

predictor_labels <- c(
	MWa = "Molecular weight",
	logKowb = "Octanol-water partition coefficient",
	Mptc = "Melting point",
	LogSaqd = "Aqueous solubility",
	LogSoce = "Organic-solvent solubility",
	Texpi = "Experimental temperature",
	Skin.thicknessj = "Skin thickness",
	Hdf = "Hydrogen-bond donor count",
	Hag = "Hydrogen-bond acceptor count",
	MVh = "Molecular volume"
)

predictor_units <- c(
	MWa = "Da",
	logKowb = "logKow",
	Mptc = "K",
	LogSaqd = "log mol/mL",
	LogSoce = "log mol/mL",
	Texpi = "K",
	Skin.thicknessj = "mm",
	Hdf = "count",
	Hag = "count",
	MVh = "cm3/mol"
)

############################################################
# Check required columns
############################################################

required_cols <- unique(c(
	id_col,
	"CAS.No",
	"Compound",
	outcome_col,
	selected_predictors
))

missing_cols <- required_cols[
	!(required_cols %in% names(df))
]

if (length(missing_cols) > 0) {
	stop(
		paste(
			"Missing required columns:",
			paste(missing_cols, collapse = ", ")
		)
	)
}

############################################################
# Prepare common complete-case dataset
############################################################

model_df <- df[
	complete.cases(df[, required_cols, drop = FALSE]),
	,
	drop = FALSE
]

model_df$row_id <- seq_len(
	nrow(model_df)
)

############################################################
# Helper functions
############################################################

get_label <- function(predictor) {
	if (predictor %in% names(predictor_labels)) {
		return(
			unname(
				predictor_labels[predictor]
			)
		)
	}

	predictor
}

get_unit <- function(predictor) {
	if (predictor %in% names(predictor_units)) {
		return(
			unname(
				predictor_units[predictor]
			)
		)
	}

	""
}

make_axis_label <- function(predictor) {
	unit <- get_unit(predictor)

	if (unit == "") {
		return(
			get_label(predictor)
		)
	}

	paste0(
		get_label(predictor),
		" (",
		unit,
		")"
	)
}

safe_rbind <- function(x) {
	x <- x[!sapply(x, is.null)]

	if (length(x) == 0) {
		return(data.frame())
	}

	all_cols <- unique(
		unlist(
			lapply(x, names)
		)
	)

	x2 <- lapply(
		x,
		function(d) {
			missing_cols <- setdiff(
				all_cols,
				names(d)
			)

			if (length(missing_cols) > 0) {
				for (m in missing_cols) {
					d[[m]] <- NA
				}
			}

			d[
				,
				all_cols,
				drop = FALSE
			]
		}
	)

	do.call(
		rbind,
		x2
	)
}

############################################################
# LOCO-CV for selected model on same complete-case dataset
############################################################

run_loco_cv <- function(data,
						model_formula,
						id_col = "compound_id",
						outcome_col = "logkpl") {

	ids <- unique(
		data[[id_col]]
	)

	pred_list <- vector(
		"list",
		length(ids)
	)

	for (i in seq_along(ids)) {

		id <- ids[i]

		training <- data[
			data[[id_col]] != id,
			,
			drop = FALSE
		]

		test <- data[
			data[[id_col]] == id,
			,
			drop = FALSE
		]

		fit <- lm(
			model_formula,
			data = training
		)

		mu <- as.numeric(
			predict(
				fit,
				newdata = test
			)
		)

		pred_list[[i]] <- data.frame(
			row_id = test$row_id,
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			observed = test[[outcome_col]],
			mu = mu,
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(
		pred_list
	)

	pred_df$residual <- pred_df$observed - pred_df$mu
	pred_df$abs_error <- abs(
		pred_df$residual
	)

	pred_df
}

############################################################
# Create descriptor-domain summary
############################################################

descriptor_summary_list <- lapply(
	selected_predictors,
	function(p) {

		x <- model_df[[p]]

		data.frame(
			predictor = p,
			predictor_label = get_label(p),
			unit = get_unit(p),
			min = min(x, na.rm = TRUE),
			q2.5 = as.numeric(
				quantile(x, 0.025, na.rm = TRUE)
			),
			q10 = as.numeric(
				quantile(x, 0.10, na.rm = TRUE)
			),
			median = median(x, na.rm = TRUE),
			q90 = as.numeric(
				quantile(x, 0.90, na.rm = TRUE)
			),
			q97.5 = as.numeric(
				quantile(x, 0.975, na.rm = TRUE)
			),
			max = max(x, na.rm = TRUE),
			stringsAsFactors = FALSE
		)
	}
)

descriptor_summary <- safe_rbind(
	descriptor_summary_list
)

############################################################
# Classify domain for each observation
############################################################

classify_one_predictor <- function(x, predictor, descriptor_summary) {

	row <- descriptor_summary[
		descriptor_summary$predictor == predictor,
	]

	if (nrow(row) != 1) {
		stop(
			paste(
				"Descriptor summary not found for predictor:",
				predictor
			)
		)
	}

	out <- rep(
		"central",
		length(x)
	)

	out[
		x < row$q10 |
			x > row$q90
	] <- "broad"

	out[
		x < row$q2.5 |
			x > row$q97.5
	] <- "outside"

	out
}

domain_class_matrix <- sapply(
	selected_predictors,
	function(p) {
		classify_one_predictor(
			x = model_df[[p]],
			predictor = p,
			descriptor_summary = descriptor_summary
		)
	}
)

domain_class_matrix <- as.data.frame(
	domain_class_matrix,
	stringsAsFactors = FALSE
)

names(domain_class_matrix) <- paste0(
	selected_predictors,
	"_domain_class"
)

model_df_with_domain <- cbind(
	model_df,
	domain_class_matrix
)

domain_rank <- c(
	central = 1,
	broad = 2,
	outside = 3
)

overall_domain_class <- apply(
	domain_class_matrix,
	1,
	function(x) {
		ranks <- domain_rank[x]

		names(domain_rank)[
			max(ranks, na.rm = TRUE)
		]
	}
)

model_df_with_domain$domain_class <- overall_domain_class

model_df_with_domain$n_predictors_outside_broad <- rowSums(
	domain_class_matrix == "outside"
)

model_df_with_domain$n_predictors_outside_central <- rowSums(
	domain_class_matrix != "central"
)

############################################################
# Get LOCO-CV predictions and merge domain information
############################################################

ad_predictions <- run_loco_cv(
	data = model_df,
	model_formula = selected_formula,
	id_col = id_col,
	outcome_col = outcome_col
)

domain_cols <- c(
	"row_id",
	selected_predictors,
	paste0(selected_predictors, "_domain_class"),
	"domain_class",
	"n_predictors_outside_broad",
	"n_predictors_outside_central"
)

ad_predictions <- merge(
	ad_predictions,
	model_df_with_domain[, domain_cols],
	by = "row_id",
	all.x = TRUE,
	all.y = FALSE
)

ad_predictions$high_error <- ad_predictions$abs_error >
	high_error_threshold

############################################################
# Error summary by domain class
############################################################

summarize_error <- function(data, group_name) {
	data.frame(
		group = group_name,
		n_observations = nrow(data),
		n_compounds = length(unique(data$compound_id)),
		RMSE = rmse(data$observed, data$mu),
		MAE = mae(data$observed, data$mu),
		R2_pred = r2_pred(data$observed, data$mu),
		median_abs_error = median(data$abs_error, na.rm = TRUE),
		n_abs_error_gt_1 = sum(data$abs_error > high_error_threshold, na.rm = TRUE),
		proportion_abs_error_gt_1 = mean(data$abs_error > high_error_threshold, na.rm = TRUE),
		stringsAsFactors = FALSE
	)
}

domain_levels <- c(
	"central",
	"broad",
	"outside"
)

error_by_domain_list <- lapply(
	domain_levels,
	function(level) {
		this_data <- ad_predictions[
			ad_predictions$domain_class == level,
		]

		if (nrow(this_data) == 0) {
			return(NULL)
		}

		summarize_error(
			this_data,
			group_name = level
		)
	}
)

error_by_domain <- safe_rbind(
	error_by_domain_list
)

error_overall <- summarize_error(
	ad_predictions,
	group_name = "overall"
)

error_by_domain <- safe_rbind(
	list(
		error_overall,
		error_by_domain
	)
)

############################################################
# High-error summary table
############################################################

high_error_summary <- ad_predictions[
	ad_predictions$high_error,
]

high_error_summary <- high_error_summary[
	order(
		-high_error_summary$abs_error
	),
]

high_error_summary <- high_error_summary[
	,
	c(
		"compound_id",
		"CAS.No",
		"Compound",
		"observed",
		"mu",
		"residual",
		"abs_error",
		"domain_class",
		"n_predictors_outside_broad",
		"n_predictors_outside_central",
		selected_predictors
	)
]

############################################################
# Save tables
############################################################

write.csv(
	ad_predictions,
	path_ad_prediction_dataset,
	row.names = FALSE
)

write.csv(
	descriptor_summary,
	path_ad_descriptor_summary,
	row.names = FALSE
)

write.csv(
	error_by_domain,
	path_ad_error_by_domain,
	row.names = FALSE
)

write.csv(
	high_error_summary,
	path_ad_high_error_summary,
	row.names = FALSE
)

write.csv(
	descriptor_summary,
	path_table_ad_descriptor_summary,
	row.names = FALSE
)

write.csv(
	error_by_domain,
	path_table_ad_error_by_domain,
	row.names = FALSE
)

write.csv(
	high_error_summary,
	path_table_ad_high_error_summary,
	row.names = FALSE
)

############################################################
# Figure: error by domain class
############################################################

make_error_by_domain_plot <- function() {

	plot_df <- ad_predictions[
		ad_predictions$domain_class %in% domain_levels,
	]

	plot_df$domain_class <- factor(
		plot_df$domain_class,
		levels = domain_levels
	)

	par(
		mar = c(5, 5, 3, 2)
	)

	boxplot(
		abs_error ~ domain_class,
		data = plot_df,
		xlab = "Applicability-domain class",
		ylab = "Absolute prediction error",
		main = "Prediction error by descriptor-domain class"
	)

	stripchart(
		abs_error ~ domain_class,
		data = plot_df,
		vertical = TRUE,
		method = "jitter",
		pch = 16,
		cex = 0.7,
		add = TRUE
	)

	abline(
		h = high_error_threshold,
		lty = 2,
		lwd = 2
	)

	legend(
		"topright",
		legend = paste0(
			"High-error threshold = ",
			high_error_threshold,
			" log unit"
		),
		lty = 2,
		lwd = 2,
		bty = "n"
	)
}

pdf(
	path_fig_ad_error_by_domain_pdf,
	width = 7,
	height = 5
)

make_error_by_domain_plot()

dev.off()

png(
	path_fig_ad_error_by_domain_png,
	width = 1800,
	height = 1200,
	res = 250
)

make_error_by_domain_plot()

dev.off()

############################################################
# Figure: absolute error by selected predictor
############################################################

make_abs_error_by_predictor_plot <- function() {

	n_predictors <- length(selected_predictors)

	n_col <- 3
	n_row <- ceiling(n_predictors / n_col)

	op <- par(
		mfrow = c(n_row, n_col),
		mar = c(4.5, 4.8, 3, 1),
		oma = c(0, 0, 2.5, 0)
	)

	on.exit(
		par(op)
	)

	for (p in selected_predictors) {

		plot(
			ad_predictions[[p]],
			ad_predictions$abs_error,
			pch = 16,
			cex = 0.8,
			xlab = make_axis_label(p),
			ylab = "Absolute prediction error",
			main = get_label(p)
		)

		if (sum(complete.cases(ad_predictions[, c(p, "abs_error")])) > 5) {
			lines(
				lowess(
					ad_predictions[[p]],
					ad_predictions$abs_error
				),
				lwd = 2
			)
		}

		abline(
			h = high_error_threshold,
			lty = 2
		)

		rug(
			ad_predictions[[p]],
			side = 1
		)
	}

	if (n_predictors < n_row * n_col) {
		for (i in seq_len(n_row * n_col - n_predictors)) {
			plot.new()
		}
	}

	mtext(
		"Absolute prediction error by selected descriptor",
		outer = TRUE,
		cex = 1.2,
		font = 2
	)
}

pdf(
	path_fig_ad_abs_error_by_predictor_pdf,
	width = 11,
	height = 8.5
)

make_abs_error_by_predictor_plot()

dev.off()

png(
	path_fig_ad_abs_error_by_predictor_png,
	width = 2750,
	height = 2125,
	res = 250
)

make_abs_error_by_predictor_plot()

dev.off()

############################################################
# Console summary
############################################################

cat("\nApplicability-domain summary complete.\n")
cat("Selected model:\n")
cat(selected_model_text, "\n\n")

cat("Selected predictors:\n")
cat(paste(selected_predictors, collapse = ", "), "\n\n")

cat("Domain class counts:\n")
print(
	table(ad_predictions$domain_class)
)

cat("\nError by domain class:\n")
print(error_by_domain)

cat("\nHigh-error observations:", nrow(high_error_summary), "\n")
cat("Output table:", path_table_ad_error_by_domain, "\n")