############################################################
# 06_heteroscedastic_uncertainty_model.R
# Evaluate heteroscedastic Gaussian uncertainty models
############################################################

source("R/00_config.R")

############################################################
# Load cleaned dataset and selected mean model
############################################################

df <- read.csv(
	path_cleaned_dataset,
	stringsAsFactors = FALSE
)

df$Compound <- trimws(df$Compound)

selected_mean_model <- readLines(
	path_loco_cv_selected_model,
	warn = FALSE
)

selected_mean_model <- selected_mean_model[1]

############################################################
# User settings
############################################################

outcome_col <- "logkpl"
id_col <- "compound_id"

mean_formula <- as.formula(selected_mean_model)

selected_model_predictors <- get_formula_predictors(
	selected_mean_model,
	outcome_col = outcome_col
)

# Candidate variance models.
# These model log(sigma_i^2).
# "mu" means the fitted predicted logKp from the mean model.
variance_models <- list(
	var_mu = c("mu"),
	var_MWa = c("MWa"),
	var_Mptc = c("Mptc"),
	var_LogSaqd = c("LogSaqd"),
	var_LogSoce = c("LogSoce"),
	var_Texpi = c("Texpi"),
	var_mu_MWa = c("mu", "MWa"),
	var_mu_Mptc = c("mu", "Mptc"),
	var_mu_LogSaqd = c("mu", "LogSaqd"),
	var_mu_LogSoce = c("mu", "LogSoce"),
	var_mu_Texpi = c("mu", "Texpi"),
	var_selected_predictors = selected_model_predictors,
	var_mu_selected_predictors = c("mu", selected_model_predictors)
)

############################################################
# Check required columns
############################################################

required_vars <- unique(c(
	all.vars(mean_formula),
	id_col,
	"CAS.No",
	"Compound",
	unlist(variance_models)
))

required_vars <- setdiff(
	required_vars,
	"mu"
)

missing_vars <- setdiff(
	required_vars,
	names(df)
)

if (length(missing_vars) > 0) {
	stop(
		paste(
			"These required columns are missing from df:",
			paste(missing_vars, collapse = ", ")
		)
	)
}

############################################################
# Prepare analysis data
############################################################

dat <- df[
	,
	required_vars,
	drop = FALSE
]

dat <- dat[
	complete.cases(dat),
	,
	drop = FALSE
]

cat("Rows used:", nrow(dat), "\n")
cat("Unique compounds:", length(unique(dat[[id_col]])), "\n")
cat("Mean model:\n")
cat(selected_mean_model, "\n")

############################################################
# Helper functions
############################################################

safe_sd <- function(x) {
	s <- sd(x, na.rm = TRUE)

	if (is.na(s) || s == 0) {
		s <- 1
	}

	s
}

scale_matrix_train <- function(X) {
	X <- as.matrix(X)

	center <- rep(0, ncol(X))
	scale <- rep(1, ncol(X))

	names(center) <- colnames(X)
	names(scale) <- colnames(X)

	for (j in seq_len(ncol(X))) {
		if (colnames(X)[j] != "(Intercept)") {
			center[j] <- mean(X[, j], na.rm = TRUE)
			scale[j] <- safe_sd(X[, j])
			X[, j] <- (X[, j] - center[j]) / scale[j]
		}
	}

	list(
		X = X,
		center = center,
		scale = scale
	)
}

scale_matrix_apply <- function(X, center, scale) {
	X <- as.matrix(X)

	for (j in seq_len(ncol(X))) {
		nm <- colnames(X)[j]

		if (nm != "(Intercept)") {
			X[, j] <- (X[, j] - center[nm]) / scale[nm]
		}
	}

	X
}

make_mean_matrix <- function(formula, data) {
	model.matrix(
		formula,
		data = data
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

			d <- d[
				,
				all_cols,
				drop = FALSE
			]

			d
		}
	)

	do.call(
		rbind,
		x2
	)
}

############################################################
# Variance model matrix functions
############################################################

make_var_raw <- function(data, mu, var_terms) {
	out <- data.frame(
		row.names = seq_len(nrow(data))
	)

	for (v in var_terms) {
		if (v == "mu") {
			out$mu <- as.numeric(mu)
		} else {
			out[[v]] <- data[[v]]
		}
	}

	out
}

make_var_scaler <- function(data, mu, var_terms) {
	Zraw <- make_var_raw(
		data = data,
		mu = mu,
		var_terms = var_terms
	)

	center <- sapply(
		Zraw,
		mean,
		na.rm = TRUE
	)

	scale <- sapply(
		Zraw,
		safe_sd
	)

	list(
		center = center,
		scale = scale
	)
}

make_var_matrix <- function(data, mu, var_terms, scaler) {
	Zraw <- make_var_raw(
		data = data,
		mu = mu,
		var_terms = var_terms
	)

	for (v in names(Zraw)) {
		Zraw[[v]] <- (Zraw[[v]] - scaler$center[v]) / scaler$scale[v]
	}

	Z <- cbind(
		"(Intercept)" = 1,
		as.matrix(Zraw)
	)

	Z
}

############################################################
# Full MLE heteroscedastic model
############################################################

fit_hetero_mle <- function(data_train,
						   mean_formula,
						   var_terms,
						   outcome_col = "logkpl",
						   maxit = 3000) {

	y <- data_train[[outcome_col]]

	X0 <- make_mean_matrix(
		mean_formula,
		data_train
	)

	X_scaled <- scale_matrix_train(X0)
	X <- X_scaled$X

	lm_init <- lm.fit(
		x = X,
		y = y
	)

	beta_init <- as.numeric(lm_init$coefficients)
	beta_init[is.na(beta_init)] <- 0

	mu_init <- as.numeric(
		X %*% beta_init
	)

	resid_init <- y - mu_init

	sigma2_init <- mean(
		resid_init^2,
		na.rm = TRUE
	)

	if (is.na(sigma2_init) || sigma2_init <= 0) {
		sigma2_init <- var(
			y,
			na.rm = TRUE
		)
	}

	var_scaler <- make_var_scaler(
		data = data_train,
		mu = mu_init,
		var_terms = var_terms
	)

	Z_init <- make_var_matrix(
		data = data_train,
		mu = mu_init,
		var_terms = var_terms,
		scaler = var_scaler
	)

	gamma_init <- rep(
		0,
		ncol(Z_init)
	)

	names(gamma_init) <- colnames(Z_init)
	gamma_init[1] <- log(sigma2_init)

	par_init <- c(
		beta_init,
		gamma_init
	)

	n_beta <- length(beta_init)
	n_gamma <- length(gamma_init)

	neg_loglik <- function(par) {
		beta <- par[seq_len(n_beta)]
		gamma <- par[n_beta + seq_len(n_gamma)]

		mu <- as.numeric(
			X %*% beta
		)

		Z <- make_var_matrix(
			data = data_train,
			mu = mu,
			var_terms = var_terms,
			scaler = var_scaler
		)

		log_var <- as.numeric(
			Z %*% gamma
		)

		log_var <- pmin(
			pmax(log_var, -20),
			20
		)

		var_i <- exp(log_var)
		resid <- y - mu

		nll <- 0.5 * sum(
			log(2 * pi) +
				log_var +
				(resid^2 / var_i)
		)

		if (!is.finite(nll)) {
			nll <- 1e100
		}

		nll
	}

	opt <- optim(
		par = par_init,
		fn = neg_loglik,
		method = "BFGS",
		control = list(
			maxit = maxit,
			reltol = 1e-9
		)
	)

	beta_hat <- opt$par[seq_len(n_beta)]
	gamma_hat <- opt$par[n_beta + seq_len(n_gamma)]

	names(beta_hat) <- colnames(X)
	names(gamma_hat) <- colnames(Z_init)

	mu_hat <- as.numeric(
		X %*% beta_hat
	)

	Z_hat <- make_var_matrix(
		data = data_train,
		mu = mu_hat,
		var_terms = var_terms,
		scaler = var_scaler
	)

	log_var_hat <- as.numeric(
		Z_hat %*% gamma_hat
	)

	log_var_hat <- pmin(
		pmax(log_var_hat, -20),
		20
	)

	sigma_hat <- sqrt(
		exp(log_var_hat)
	)

	nll <- opt$value
	logLik <- -nll
	npar <- length(beta_hat) + length(gamma_hat)
	aic <- 2 * npar + 2 * nll

	list(
		beta = beta_hat,
		gamma = gamma_hat,
		mean_formula = mean_formula,
		var_terms = var_terms,
		x_center = X_scaled$center,
		x_scale = X_scaled$scale,
		var_scaler = var_scaler,
		fitted = data.frame(
			observed = y,
			mu = mu_hat,
			sigma = sigma_hat,
			diff = y - mu_hat,
			abs_diff = abs(y - mu_hat)
		),
		logLik = logLik,
		AIC = aic,
		npar = npar,
		convergence = opt$convergence,
		message = opt$message
	)
}

############################################################
# Prediction function
############################################################

predict_hetero_mle <- function(fit, newdata) {

	X0 <- make_mean_matrix(
		fit$mean_formula,
		newdata
	)

	missing_cols <- setdiff(
		names(fit$beta),
		colnames(X0)
	)

	if (length(missing_cols) > 0) {
		stop(
			paste(
				"Missing columns in prediction mean matrix:",
				paste(missing_cols, collapse = ", ")
			)
		)
	}

	X0 <- X0[
		,
		names(fit$beta),
		drop = FALSE
	]

	X <- scale_matrix_apply(
		X0,
		fit$x_center,
		fit$x_scale
	)

	mu <- as.numeric(
		X %*% fit$beta
	)

	Z <- make_var_matrix(
		data = newdata,
		mu = mu,
		var_terms = fit$var_terms,
		scaler = fit$var_scaler
	)

	Z <- Z[
		,
		names(fit$gamma),
		drop = FALSE
	]

	log_var <- as.numeric(
		Z %*% fit$gamma
	)

	log_var <- pmin(
		pmax(log_var, -20),
		20
	)

	sigma <- sqrt(
		exp(log_var)
	)

	data.frame(
		mu = mu,
		sigma = sigma,
		lower95 = mu - 1.96 * sigma,
		upper95 = mu + 1.96 * sigma
	)
}

############################################################
# Homoscedastic baseline LOCO-CV using lm
############################################################

run_homoscedastic_cv <- function(data,
								 mean_formula,
								 id_col = "compound_id",
								 outcome_col = "logkpl") {

	ids <- unique(data[[id_col]])
	pred_list <- vector("list", length(ids))

	for (i in seq_along(ids)) {
		id <- ids[i]

		train <- data[
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
			mean_formula,
			data = train
		)

		pred <- predict(
			fit,
			newdata = test
		)

		sigma <- summary(fit)$sigma

		pred_list[[i]] <- data.frame(
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			observed = test[[outcome_col]],
			mu = as.numeric(pred),
			sigma = sigma,
			lower95 = as.numeric(pred) - 1.96 * sigma,
			upper95 = as.numeric(pred) + 1.96 * sigma,
			model = "homoscedastic_lm",
			stringsAsFactors = FALSE
		)
	}

	pred_df <- safe_rbind(pred_list)

	pred_df$diff <- pred_df$observed - pred_df$mu
	pred_df$abs_diff <- abs(pred_df$diff)
	pred_df$covered95 <- pred_df$observed >= pred_df$lower95 &
		pred_df$observed <= pred_df$upper95

	summary_df <- data.frame(
		model = "homoscedastic_lm",
		RMSE = rmse(pred_df$observed, pred_df$mu),
		MAE = mae(pred_df$observed, pred_df$mu),
		R2_pred = r2_pred(pred_df$observed, pred_df$mu),
		coverage95 = mean(pred_df$covered95, na.rm = TRUE),
		mean_sigma = mean(pred_df$sigma, na.rm = TRUE),
		median_sigma = median(pred_df$sigma, na.rm = TRUE),
		cor_abs_error_sigma = cor(
			pred_df$abs_diff,
			pred_df$sigma,
			use = "complete.obs"
		),
		failed_folds = 0,
		stringsAsFactors = FALSE
	)

	list(
		prediction = pred_df,
		summary = summary_df
	)
}

############################################################
# Heteroscedastic LOCO-CV
############################################################

run_hetero_cv_one_model <- function(data,
									mean_formula,
									var_terms,
									var_model_name,
									id_col = "compound_id",
									outcome_col = "logkpl",
									maxit = 3000) {

	ids <- unique(data[[id_col]])
	pred_list <- vector("list", length(ids))
	fold_info <- data.frame()

	for (i in seq_along(ids)) {
		id <- ids[i]

		cat(
			"Running",
			var_model_name,
			"- fold",
			i,
			"of",
			length(ids),
			"\n"
		)

		train <- data[
			data[[id_col]] != id,
			,
			drop = FALSE
		]

		test <- data[
			data[[id_col]] == id,
			,
			drop = FALSE
		]

		fit <- try(
			fit_hetero_mle(
				data_train = train,
				mean_formula = mean_formula,
				var_terms = var_terms,
				outcome_col = outcome_col,
				maxit = maxit
			),
			silent = TRUE
		)

		if (inherits(fit, "try-error")) {
			warning(
				"Model failed for fold ID: ",
				id
			)

			next
		}

		pred <- try(
			predict_hetero_mle(
				fit,
				test
			),
			silent = TRUE
		)

		if (inherits(pred, "try-error")) {
			warning(
				"Prediction failed for fold ID: ",
				id
			)

			next
		}

		pred_list[[i]] <- data.frame(
			compound_id = test[[id_col]],
			CAS.No = test$CAS.No,
			Compound = test$Compound,
			observed = test[[outcome_col]],
			mu = pred$mu,
			sigma = pred$sigma,
			lower95 = pred$lower95,
			upper95 = pred$upper95,
			model = var_model_name,
			stringsAsFactors = FALSE
		)

		fold_info <- rbind(
			fold_info,
			data.frame(
				model = var_model_name,
				fold = i,
				compound_id = id,
				CAS.No = test$CAS.No[1],
				convergence = fit$convergence,
				logLik = fit$logLik,
				AIC = fit$AIC,
				npar = fit$npar,
				stringsAsFactors = FALSE
			)
		)
	}

	pred_df <- safe_rbind(pred_list)

	pred_df$diff <- pred_df$observed - pred_df$mu
	pred_df$abs_diff <- abs(pred_df$diff)
	pred_df$covered95 <- pred_df$observed >= pred_df$lower95 &
		pred_df$observed <= pred_df$upper95

	summary_df <- data.frame(
		model = var_model_name,
		RMSE = rmse(pred_df$observed, pred_df$mu),
		MAE = mae(pred_df$observed, pred_df$mu),
		R2_pred = r2_pred(pred_df$observed, pred_df$mu),
		coverage95 = mean(pred_df$covered95, na.rm = TRUE),
		mean_sigma = mean(pred_df$sigma, na.rm = TRUE),
		median_sigma = median(pred_df$sigma, na.rm = TRUE),
		cor_abs_error_sigma = cor(
			pred_df$abs_diff,
			pred_df$sigma,
			use = "complete.obs"
		),
		failed_folds = length(ids) - length(unique(pred_df[[id_col]])),
		stringsAsFactors = FALSE
	)

	list(
		prediction = pred_df,
		summary = summary_df,
		fold_info = fold_info
	)
}

run_all_hetero_cv <- function(data,
							  mean_formula,
							  variance_models,
							  id_col = "compound_id",
							  outcome_col = "logkpl",
							  maxit = 3000) {

	all_predictions <- list()
	all_summaries <- list()
	all_fold_info <- list()

	for (m in names(variance_models)) {
		res <- run_hetero_cv_one_model(
			data = data,
			mean_formula = mean_formula,
			var_terms = variance_models[[m]],
			var_model_name = m,
			id_col = id_col,
			outcome_col = outcome_col,
			maxit = maxit
		)

		all_predictions[[m]] <- res$prediction
		all_summaries[[m]] <- res$summary
		all_fold_info[[m]] <- res$fold_info
	}

	list(
		prediction = safe_rbind(all_predictions),
		summary = safe_rbind(all_summaries),
		fold_info = safe_rbind(all_fold_info)
	)
}

############################################################
# Fit full-data models
############################################################

fit_lm_full <- lm(
	mean_formula,
	data = dat
)

full_fit_summary <- data.frame(
	model = "homoscedastic_lm",
	logLik = as.numeric(logLik(fit_lm_full)),
	AIC = AIC(fit_lm_full),
	npar = length(coef(fit_lm_full)) + 1,
	convergence = NA,
	mean_sigma = summary(fit_lm_full)$sigma,
	median_sigma = summary(fit_lm_full)$sigma,
	cor_abs_error_sigma = NA_real_,
	stringsAsFactors = FALSE
)

full_fit_list <- list()

for (m in names(variance_models)) {
	cat(
		"Fitting full-data heteroscedastic model:",
		m,
		"\n"
	)

	fit_m <- fit_hetero_mle(
		data_train = dat,
		mean_formula = mean_formula,
		var_terms = variance_models[[m]],
		outcome_col = outcome_col
	)

	full_fit_list[[m]] <- fit_m

	full_fit_summary <- rbind(
		full_fit_summary,
		data.frame(
			model = m,
			logLik = fit_m$logLik,
			AIC = fit_m$AIC,
			npar = fit_m$npar,
			convergence = fit_m$convergence,
			mean_sigma = mean(fit_m$fitted$sigma),
			median_sigma = median(fit_m$fitted$sigma),
			cor_abs_error_sigma = cor(
				fit_m$fitted$abs_diff,
				fit_m$fitted$sigma,
				use = "complete.obs"
			),
			stringsAsFactors = FALSE
		)
	)
}

full_fit_summary <- full_fit_summary[
	order(full_fit_summary$AIC),
]

############################################################
# Leave-one-compound-out CV
############################################################

cat("Running homoscedastic leave-one-compound-out CV...\n")

cv_lm <- run_homoscedastic_cv(
	data = dat,
	mean_formula = mean_formula,
	id_col = id_col,
	outcome_col = outcome_col
)

cat("Running heteroscedastic leave-one-compound-out CV...\n")

cv_hetero <- run_all_hetero_cv(
	data = dat,
	mean_formula = mean_formula,
	variance_models = variance_models,
	id_col = id_col,
	outcome_col = outcome_col,
	maxit = 3000
)

cv_summary <- safe_rbind(
	list(
		cv_lm$summary,
		cv_hetero$summary
	)
)

cv_summary <- cv_summary[
	order(
		cv_summary$RMSE,
		abs(cv_summary$coverage95 - 0.95),
		-cv_summary$cor_abs_error_sigma
	),
]

############################################################
# Select best heteroscedastic model
############################################################

cv_hetero_summary <- cv_hetero$summary

cv_hetero_summary <- cv_hetero_summary[
	order(
		-cv_hetero_summary$cor_abs_error_sigma,
		abs(cv_hetero_summary$coverage95 - 0.95),
		cv_hetero_summary$RMSE
	),
]

best_hetero_name <- cv_hetero_summary$model[1]

best_pred <- cv_hetero$prediction[
	cv_hetero$prediction$model == best_hetero_name,
]

best_pred$high_error <- abs(best_pred$diff) > 1

writeLines(
	best_hetero_name,
	path_hetero_selected_model
)

############################################################
# Save outputs
############################################################

write.csv(
	full_fit_summary,
	path_hetero_full_fit_summary,
	row.names = FALSE
)

write.csv(
	cv_summary,
	path_hetero_cv_summary,
	row.names = FALSE
)

write.csv(
	cv_summary,
	path_table_hetero_cv_summary,
	row.names = FALSE
)

write.csv(
	cv_hetero$prediction,
	path_hetero_cv_predictions_all,
	row.names = FALSE
)

write.csv(
	cv_hetero$fold_info,
	path_hetero_fold_info,
	row.names = FALSE
)

write.csv(
	best_pred,
	path_hetero_cv_predictions_best,
	row.names = FALSE
)

############################################################
# Console summary
############################################################

cat("Heteroscedastic uncertainty analysis complete.\n")
cat("Mean model:\n")
cat(selected_mean_model, "\n")
cat("Selected heteroscedastic variance model:", best_hetero_name, "\n")
cat("Best heteroscedastic RMSE:", rmse(best_pred$observed, best_pred$mu), "\n")
cat("Best heteroscedastic MAE:", mae(best_pred$observed, best_pred$mu), "\n")
cat("Best heteroscedastic R2_pred:", r2_pred(best_pred$observed, best_pred$mu), "\n")
cat("Best heteroscedastic coverage95:", mean(best_pred$covered95, na.rm = TRUE), "\n")
cat("Best heteroscedastic cor_abs_error_sigma:",
	cor(
		best_pred$abs_diff,
		best_pred$sigma,
		use = "complete.obs"
	),
	"\n"
)