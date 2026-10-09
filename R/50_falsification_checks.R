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

# ------------------------------------------------------------
# 2. Wrong-source calibration negative control
# ------------------------------------------------------------

references <- unique(as.character(hu$reference_id))
anchor_sizes <- 1:5
source_objects <- list()

for (ref in references) {
	train <- hu[as.character(hu$reference_id) != ref, ]
	held <- hu[as.character(hu$reference_id) == ref, ]

	train_compounds <- unique(as.character(train$compound_id))
	shared <- intersect(unique(as.character(held$compound_id)), train_compounds)
	novel <- setdiff(unique(as.character(held$compound_id)), train_compounds)

	if (!length(shared) || !length(novel)) next

	p <- make_profile(train, c("mw", "log_kow"))
	center <- sapply(p[c("mw", "log_kow")], mean)
	scale <- sapply(p[c("mw", "log_kow")], sd)

	for (v in c("mw", "log_kow")) {
		train[[paste0(v, "_z")]] <- (train[[v]] - center[v]) / scale[v]
		held[[paste0(v, "_z")]] <- (held[[v]] - center[v]) / scale[v]
	}

	fit <- lmer(
		log_kp ~ mw_z + log_kow_z + (1 | compound_id) + (1 | reference_id),
		data = train,
		REML = TRUE,
		control = ctrl
	)

	vc <- as.data.frame(VarCorr(fit))
	v_ref <- vc$vcov[vc$grp == "reference_id"]
	v_res <- vc$vcov[vc$grp == "Residual"]

	anchor_rows <- held[as.character(held$compound_id) %in% shared, ]
	novel_rows <- held[as.character(held$compound_id) %in% novel, ]

	anchor_rows$expected <- predict(fit, newdata = anchor_rows, re.form = ~(1 | compound_id), allow.new.levels = TRUE)
	novel_rows$expected <- predict(fit, newdata = novel_rows, re.form = NA, allow.new.levels = TRUE)

	anchor_residual <- aggregate(anchor_rows$log_kp - anchor_rows$expected, list(compound_id = anchor_rows$compound_id), mean)
	names(anchor_residual)[2] <- "residual"

	offsets <- list()
	q <- 1

	for (k in anchor_sizes) {
		if (length(shared) < k) next

		for (rep in 1:100) {
			selected <- sample(shared, k)
			offset_raw <- mean(anchor_residual$residual[as.character(anchor_residual$compound_id) %in% selected])
			shrinkage <- k * v_ref / (v_res + k * v_ref)

			offsets[[q]] <- data.frame(
				anchors = k,
				replicate = rep,
				offset = shrinkage * offset_raw
			)
			q <- q + 1
		}
	}

	source_objects[[ref]] <- list(
		reference_id = ref,
		novel = novel_rows,
		offsets = do.call(rbind, offsets)
	)
}

eligible_refs <- names(source_objects)
negative <- list()
z <- 1

for (ref in eligible_refs) {
	obj <- source_objects[[ref]]

	for (k in anchor_sizes) {
		correct_offsets <- obj$offsets[obj$offsets$anchors == k, ]
		if (!nrow(correct_offsets)) next

		wrong_refs <- eligible_refs[sapply(eligible_refs, function(x) any(source_objects[[x]]$offsets$anchors == k))]
		wrong_refs <- setdiff(wrong_refs, ref)
		if (!length(wrong_refs)) next

		for (rep in correct_offsets$replicate) {
			correct_offset <- correct_offsets$offset[correct_offsets$replicate == rep]
			wrong_ref <- sample(wrong_refs, 1)
			wrong_offsets <- source_objects[[wrong_ref]]$offsets
			wrong_offset <- wrong_offsets$offset[wrong_offsets$anchors == k & wrong_offsets$replicate == rep]

			obs <- obj$novel$log_kp
			base_pred <- obj$novel$expected

			negative[[z]] <- data.frame(
				reference_id = ref,
				wrong_reference_id = wrong_ref,
				anchors = k,
				replicate = rep,
				n_test = length(obs),
				RMSE_uncalibrated = sqrt(mean((obs - base_pred)^2)),
				RMSE_correct = sqrt(mean((obs - (base_pred + correct_offset))^2)),
				RMSE_wrong = sqrt(mean((obs - (base_pred + wrong_offset))^2)),
				MAE_uncalibrated = mean(abs(obs - base_pred)),
				MAE_correct = mean(abs(obs - (base_pred + correct_offset))),
				MAE_wrong = mean(abs(obs - (base_pred + wrong_offset)))
			)
			z <- z + 1
		}
	}
}

negative <- do.call(rbind, negative)

wrong_summary <- aggregate(
	cbind(RMSE_uncalibrated, RMSE_correct, RMSE_wrong, MAE_uncalibrated, MAE_correct, MAE_wrong) ~ anchors,
	data = negative,
	mean
)

wrong_summary$correct_RMSE_improvement <- 100 * (1 - wrong_summary$RMSE_correct / wrong_summary$RMSE_uncalibrated)
wrong_summary$wrong_RMSE_improvement <- 100 * (1 - wrong_summary$RMSE_wrong / wrong_summary$RMSE_uncalibrated)
wrong_summary$correct_MAE_improvement <- 100 * (1 - wrong_summary$MAE_correct / wrong_summary$MAE_uncalibrated)
wrong_summary$wrong_MAE_improvement <- 100 * (1 - wrong_summary$MAE_wrong / wrong_summary$MAE_uncalibrated)

n_refs <- aggregate(reference_id ~ anchors, data = negative, function(x) length(unique(x)))
names(n_refs)[2] <- "n_references"
wrong_summary <- merge(wrong_summary, n_refs, by = "anchors", sort = FALSE)
wrong_summary <- wrong_summary[order(wrong_summary$anchors), ]

wrong_reference <- aggregate(
	cbind(RMSE_uncalibrated, RMSE_correct, RMSE_wrong, MAE_uncalibrated, MAE_correct, MAE_wrong) ~ reference_id + anchors,
	data = negative,
	mean
)

wrong_reference$correct_RMSE_improvement <- 100 * (1 - wrong_reference$RMSE_correct / wrong_reference$RMSE_uncalibrated)
wrong_reference$wrong_RMSE_improvement <- 100 * (1 - wrong_reference$RMSE_wrong / wrong_reference$RMSE_uncalibrated)

write.csv(wrong_summary, "R/results/50_wrong_source_calibration.csv", row.names = FALSE)
write.csv(wrong_reference, "R/results/50_wrong_source_reference.csv", row.names = FALSE)

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
	"\n\nWrong-source calibration control\n",
	paste(
		paste0(
			wrong_summary$anchors, " anchors: correct ",
			round(wrong_summary$correct_RMSE_improvement, 1),
			"% vs wrong ", round(wrong_summary$wrong_RMSE_improvement, 1),
			"% across ", wrong_summary$n_references, " references"
		),
		collapse = "\n"
	),
	"\n",
	sep = ""
)