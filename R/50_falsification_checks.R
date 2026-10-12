# ============================================================
# 50_falsification_checks.R
# ============================================================

library(lme4)

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
set.seed(2026)

metrics <- function(obs, pred) {
	c(
		RMSE = sqrt(mean((obs - pred)^2)),
		MAE = mean(abs(obs - pred)),
		R2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
	)
}

make_profile <- function(x, descriptors) {
	p <- aggregate(x[descriptors], list(compound_id = x$compound_id), median)
	y <- aggregate(x$log_kp, list(compound_id = x$compound_id), mean)
	names(y)[2] <- "mean_log_kp"
	merge(p, y, by = "compound_id", sort = FALSE)
}

run_source_loco <- function(dat, descriptors, dataset, specification, use_temperature = FALSE) {
	compounds <- unique(as.character(dat$compound_id))
	out <- vector("list", length(compounds))

	for (i in seq_along(compounds)) {
		id <- compounds[i]
		train <- dat[as.character(dat$compound_id) != id, ]
		test <- dat[as.character(dat$compound_id) == id, ]

		p <- make_profile(train, descriptors)
		center <- sapply(p[descriptors], mean)
		scale <- sapply(p[descriptors], sd)

		for (v in descriptors) {
			p[[paste0(v, "_z")]] <- (p[[v]] - center[v]) / scale[v]
			train[[paste0(v, "_z")]] <- (train[[v]] - center[v]) / scale[v]
			test[[paste0(v, "_z")]] <- (test[[v]] - center[v]) / scale[v]
		}

		fixed <- paste(paste0(descriptors, "_z"), collapse = " + ")
		lm_fit <- lm(as.formula(paste("mean_log_kp ~", fixed)), data = p)

		if (use_temperature) {
			temp_mean <- mean(train$temperature_c)
			train$temperature_centered <- train$temperature_c - temp_mean
			test$temperature_centered <- test$temperature_c - temp_mean
			form <- as.formula(paste("log_kp ~ temperature_centered +", fixed, "+ (1 | compound_id) + (1 | reference_id)"))
		} else {
			form <- as.formula(paste("log_kp ~", fixed, "+ (1 | compound_id) + (1 | reference_id)"))
		}

		mix_fit <- lmer(form, data = train, REML = TRUE, control = ctrl)

		test_profile <- test[1, descriptors, drop = FALSE]
		for (v in descriptors) test_profile[[paste0(v, "_z")]] <- (test_profile[[v]] - center[v]) / scale[v]

		pred_conventional <- as.numeric(predict(lm_fit, newdata = test_profile))
		pred_blind <- as.numeric(predict(mix_fit, newdata = test, re.form = NA, allow.new.levels = TRUE))
		pred_known <- as.numeric(predict(mix_fit, newdata = test, re.form = ~(1 | reference_id), allow.new.levels = TRUE))

		out[[i]] <- data.frame(
			dataset = dataset,
			specification = specification,
			compound_id = id,
			observed = test$log_kp,
			conventional = pred_conventional,
			hierarchical_source_blind = pred_blind,
			hierarchical_known_source = pred_known
		)
	}

	do.call(rbind, out)
}

summarize_loco <- function(x) {
	approaches <- c("conventional", "hierarchical_source_blind", "hierarchical_known_source")

	do.call(rbind, lapply(approaches, function(a) {
		m <- metrics(x$observed, x[[a]])

		data.frame(
			dataset = x$dataset[1],
			specification = x$specification[1],
			approach = a,
			RMSE = m["RMSE"],
			MAE = m["MAE"],
			R2 = m["R2"]
		)
	}))
}

# ------------------------------------------------------------
# 1. Source-blind LOCO
# ------------------------------------------------------------

cheruvu <- read.csv("data/processed/01_cheruvu.csv", stringsAsFactors = FALSE)
cheruvu$compound_id <- factor(cheruvu$compound_id)
cheruvu$reference_id <- factor(cheruvu$reference_id)

hu <- read.csv("data/processed/30_huskindb.csv", stringsAsFactors = FALSE)
hu$compound_id <- factor(hu$compound_id)
hu$reference_id <- factor(hu$reference_id)

ch_s1 <- run_source_loco(
	cheruvu,
	c("mw", "log_kow"),
	"Cheruvu",
	"S1",
	TRUE
)

ch_s2 <- run_source_loco(
	cheruvu,
	c("mw", "log_kow", "melting_point_k", "log_saq", "h_donor", "h_acceptor"),
	"Cheruvu",
	"S2",
	TRUE
)

hu_s1 <- run_source_loco(
	hu,
	c("mw", "log_kow"),
	"HuSkinDB",
	"S1",
	FALSE
)

source_loco_predictions <- rbind(
	ch_s1,
	ch_s2,
	hu_s1
)

source_loco <- rbind(
	summarize_loco(ch_s1),
	summarize_loco(ch_s2),
	summarize_loco(hu_s1)
)

source_loco$improvement_vs_conventional <- NA_real_

for (d in unique(source_loco$dataset)) {
	for (s in unique(source_loco$specification[source_loco$dataset == d])) {
		i <- source_loco$dataset == d & source_loco$specification == s
		base <- source_loco$RMSE[i & source_loco$approach == "conventional"]
		source_loco$improvement_vs_conventional[i] <- 100 * (1 - source_loco$RMSE[i] / base)
	}
}

write.csv(source_loco, "R/results/50_source_blind_loco.csv", row.names = FALSE)
write.csv(source_loco_predictions, "R/results/50_source_blind_loco_predictions.csv", row.names = FALSE)

# ------------------------------------------------------------
# 2. Pairwise double-holdout wrong-source falsification
# ------------------------------------------------------------

references <- unique(as.character(hu$reference_id))
anchor_sizes <- 1:5
n_anchor_reps <- 100

target_refs <- references[sapply(references, function(ref) {
	held_compounds <- unique(as.character(hu$compound_id[as.character(hu$reference_id) == ref]))
	other_compounds <- unique(as.character(hu$compound_id[as.character(hu$reference_id) != ref]))
	shared <- intersect(held_compounds, other_compounds)
	novel <- setdiff(held_compounds, other_compounds)
	length(shared) > 0 && length(novel) > 0
})]

wrong_pair <- list()
z <- 1
n_failed <- 0

for (target_ref in target_refs) {
	target_all <- hu[as.character(hu$reference_id) == target_ref, ]
	other_target <- hu[as.character(hu$reference_id) != target_ref, ]

	target_compounds <- unique(as.character(target_all$compound_id))
	target_novel <- setdiff(target_compounds, unique(as.character(other_target$compound_id)))
	target_test <- target_all[as.character(target_all$compound_id) %in% target_novel, ]

	if (!nrow(target_test)) next

	for (donor_ref in setdiff(references, target_ref)) {
		donor_all <- hu[as.character(hu$reference_id) == donor_ref, ]

		train <- hu[
			!as.character(hu$reference_id) %in% c(target_ref, donor_ref),
		]

		train_compounds <- unique(as.character(train$compound_id))
		target_shared <- intersect(target_compounds, train_compounds)
		donor_shared <- intersect(unique(as.character(donor_all$compound_id)), train_compounds)

		if (!length(target_shared) || !length(donor_shared)) next

		max_k <- min(max(anchor_sizes), length(target_shared), length(donor_shared))
		if (max_k < min(anchor_sizes)) next

		p <- make_profile(train, c("mw", "log_kow"))
		center <- sapply(p[c("mw", "log_kow")], mean)
		scale <- sapply(p[c("mw", "log_kow")], sd)

		for (v in c("mw", "log_kow")) {
			train[[paste0(v, "_z")]] <- (train[[v]] - center[v]) / scale[v]
			target_all[[paste0(v, "_z")]] <- (target_all[[v]] - center[v]) / scale[v]
			donor_all[[paste0(v, "_z")]] <- (donor_all[[v]] - center[v]) / scale[v]
			target_test[[paste0(v, "_z")]] <- (target_test[[v]] - center[v]) / scale[v]
		}

		fit <- try(
			lmer(
				log_kp ~ mw_z + log_kow_z + (1 | compound_id) + (1 | reference_id),
				data = train,
				REML = TRUE,
				control = ctrl
			),
			silent = TRUE
		)

		if (inherits(fit, "try-error")) {
			n_failed <- n_failed + 1
			next
		}

		vc <- as.data.frame(VarCorr(fit))
		v_ref <- vc$vcov[vc$grp == "reference_id"]
		v_res <- vc$vcov[vc$grp == "Residual"]

		target_anchor <- target_all[as.character(target_all$compound_id) %in% target_shared, ]
		donor_anchor <- donor_all[as.character(donor_all$compound_id) %in% donor_shared, ]

		target_anchor$expected <- predict(
			fit,
			newdata = target_anchor,
			re.form = ~(1 | compound_id),
			allow.new.levels = TRUE
		)

		donor_anchor$expected <- predict(
			fit,
			newdata = donor_anchor,
			re.form = ~(1 | compound_id),
			allow.new.levels = TRUE
		)

		target_test$expected <- predict(
			fit,
			newdata = target_test,
			re.form = NA,
			allow.new.levels = TRUE
		)

		target_residual <- aggregate(
			target_anchor$log_kp - target_anchor$expected,
			list(compound_id = target_anchor$compound_id),
			mean
		)

		donor_residual <- aggregate(
			donor_anchor$log_kp - donor_anchor$expected,
			list(compound_id = donor_anchor$compound_id),
			mean
		)

		names(target_residual)[2] <- "residual"
		names(donor_residual)[2] <- "residual"

		obs <- target_test$log_kp
		base_pred <- target_test$expected

		rmse_uncalibrated <- sqrt(mean((obs - base_pred)^2))
		mae_uncalibrated <- mean(abs(obs - base_pred))

		for (k in anchor_sizes[anchor_sizes <= max_k]) {
			shrinkage <- k * v_ref / (v_res + k * v_ref)

			rmse_correct <- numeric(n_anchor_reps)
			rmse_wrong <- numeric(n_anchor_reps)
			mae_correct <- numeric(n_anchor_reps)
			mae_wrong <- numeric(n_anchor_reps)

			for (rep in seq_len(n_anchor_reps)) {
				target_selected <- sample(target_shared, k)
				donor_selected <- sample(donor_shared, k)

				target_raw <- mean(target_residual$residual[as.character(target_residual$compound_id) %in% target_selected])
				donor_raw <- mean(donor_residual$residual[as.character(donor_residual$compound_id) %in% donor_selected])

				target_offset <- shrinkage * target_raw
				donor_offset <- shrinkage * donor_raw

				rmse_correct[rep] <- sqrt(mean((obs - (base_pred + target_offset))^2))
				rmse_wrong[rep] <- sqrt(mean((obs - (base_pred + donor_offset))^2))
				mae_correct[rep] <- mean(abs(obs - (base_pred + target_offset)))
				mae_wrong[rep] <- mean(abs(obs - (base_pred + donor_offset)))
			}

			wrong_pair[[z]] <- data.frame(
				reference_id = target_ref,
				wrong_reference_id = donor_ref,
				anchors = k,
				n_test = length(obs),
				n_target_anchors = length(target_shared),
				n_donor_anchors = length(donor_shared),
				RMSE_uncalibrated = rmse_uncalibrated,
				RMSE_correct = mean(rmse_correct),
				RMSE_wrong = mean(rmse_wrong),
				MAE_uncalibrated = mae_uncalibrated,
				MAE_correct = mean(mae_correct),
				MAE_wrong = mean(mae_wrong)
			)

			z <- z + 1
		}
	}
}

wrong_pair <- do.call(rbind, wrong_pair)

wrong_pair$correct_RMSE_improvement <- 100 * (1 - wrong_pair$RMSE_correct / wrong_pair$RMSE_uncalibrated)
wrong_pair$wrong_RMSE_improvement <- 100 * (1 - wrong_pair$RMSE_wrong / wrong_pair$RMSE_uncalibrated)
wrong_pair$RMSE_wrong_minus_correct <- wrong_pair$RMSE_wrong - wrong_pair$RMSE_correct

# ------------------------------------------------------------
# Target-source summaries
# ------------------------------------------------------------

wrong_reference <- aggregate(
	cbind(
		RMSE_uncalibrated,
		RMSE_correct,
		RMSE_wrong,
		MAE_uncalibrated,
		MAE_correct,
		MAE_wrong
	) ~ reference_id + anchors,
	data = wrong_pair,
	mean
)

n_donors <- aggregate(
	wrong_reference_id ~ reference_id + anchors,
	data = wrong_pair,
	function(x) length(unique(x))
)

names(n_donors)[3] <- "n_donors"

n_test <- aggregate(
	n_test ~ reference_id + anchors,
	data = wrong_pair,
	max
)

wrong_reference <- merge(
	wrong_reference,
	n_donors,
	by = c("reference_id", "anchors"),
	sort = FALSE
)

wrong_reference <- merge(
	wrong_reference,
	n_test,
	by = c("reference_id", "anchors"),
	sort = FALSE
)

wrong_reference$correct_RMSE_improvement <- 100 * (1 - wrong_reference$RMSE_correct / wrong_reference$RMSE_uncalibrated)
wrong_reference$wrong_RMSE_improvement <- 100 * (1 - wrong_reference$RMSE_wrong / wrong_reference$RMSE_uncalibrated)
wrong_reference$RMSE_wrong_minus_correct <- wrong_reference$RMSE_wrong - wrong_reference$RMSE_correct

wrong_reference <- wrong_reference[
	order(wrong_reference$anchors, wrong_reference$reference_id),
]

# ------------------------------------------------------------
# Cohort summaries
# ------------------------------------------------------------

wrong_summary <- aggregate(
	cbind(
		RMSE_uncalibrated,
		RMSE_correct,
		RMSE_wrong,
		MAE_uncalibrated,
		MAE_correct,
		MAE_wrong
	) ~ anchors,
	data = wrong_reference,
	mean
)

wrong_summary$correct_RMSE_improvement <- 100 * (1 - wrong_summary$RMSE_correct / wrong_summary$RMSE_uncalibrated)
wrong_summary$wrong_RMSE_improvement <- 100 * (1 - wrong_summary$RMSE_wrong / wrong_summary$RMSE_uncalibrated)
wrong_summary$correct_MAE_improvement <- 100 * (1 - wrong_summary$MAE_correct / wrong_summary$MAE_uncalibrated)
wrong_summary$wrong_MAE_improvement <- 100 * (1 - wrong_summary$MAE_wrong / wrong_summary$MAE_uncalibrated)
wrong_summary$RMSE_wrong_minus_correct <- wrong_summary$RMSE_wrong - wrong_summary$RMSE_correct

n_refs <- aggregate(
	reference_id ~ anchors,
	data = wrong_reference,
	function(x) length(unique(x))
)

names(n_refs)[2] <- "n_references"

n_pairs <- aggregate(
	wrong_reference_id ~ anchors,
	data = wrong_pair,
	length
)

names(n_pairs)[2] <- "n_pairs"

wrong_summary <- merge(wrong_summary, n_refs, by = "anchors", sort = FALSE)
wrong_summary <- merge(wrong_summary, n_pairs, by = "anchors", sort = FALSE)
wrong_summary <- wrong_summary[order(wrong_summary$anchors), ]

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------

write.csv(
	wrong_pair,
	"R/results/50_wrong_source_pair.csv",
	row.names = FALSE
)

write.csv(
	wrong_reference,
	"R/results/50_wrong_source_reference.csv",
	row.names = FALSE
)

write.csv(
	wrong_summary,
	"R/results/50_wrong_source_calibration.csv",
	row.names = FALSE
)

cat(
	"\nFalsification checks complete\n\n",
	"Source-blind LOCO\n",
	paste(
		paste0(
			source_loco$dataset, " ", source_loco$specification, " ", source_loco$approach,
			": RMSE ", round(source_loco$RMSE, 3),
			" (", round(source_loco$improvement_vs_conventional, 1), "% vs conventional)"
		),
		collapse = "\n"
	),
	"\n\nPairwise double-holdout wrong-source control\n",
	paste(
		paste0(
			wrong_summary$anchors, " anchors: matched ",
			round(wrong_summary$correct_RMSE_improvement, 1),
			"% vs different-source ",
			round(wrong_summary$wrong_RMSE_improvement, 1),
			"%; RMSE difference ",
			round(wrong_summary$RMSE_wrong_minus_correct, 3),
			" across ", wrong_summary$n_references,
			" target sources and ", wrong_summary$n_pairs, " target-donor pairs"
		),
		collapse = "\n"
	),
	"\n\nDouble-holdout model fits failed: ", n_failed,
	"\n",
	sep = ""
)