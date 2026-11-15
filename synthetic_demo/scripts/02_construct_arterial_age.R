# Artery-specific proteomic age and LOWESS-adjusted age gap.
# Input: one prepared row per participant; no UK Biobank records are included.
# Protein columns must already be numeric, finite, and preprocessed/imputed.
# Synthetic protein names below are placeholders, not biological annotations.

construct_arterial_age <- function(
    data,
    artery_proteins = paste0("artery_protein_", 1:5),
    seed = 123,
    train_fraction = 0.7,
    nfolds = 5,
    lowess_fraction = 2 / 3) {
  for (package in c("caret", "glmnet")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install the required R package: ", package)
    }
  }
  if (!is.character(artery_proteins) || length(artery_proteins) < 3L ||
      anyNA(artery_proteins) || anyDuplicated(artery_proteins) ||
      any(artery_proteins %in% c("participant_id", "age", "sex", "disease_free"))) {
    stop("Supply at least three unique protein column names.")
  }
  if (length(train_fraction) != 1L || !is.finite(train_fraction) ||
      train_fraction <= 0 || train_fraction >= 1) {
    stop("train_fraction must be between 0 and 1.")
  }
  if (length(nfolds) != 1L || !is.finite(nfolds) || nfolds < 3 ||
      nfolds != as.integer(nfolds)) stop("nfolds must be an integer >= 3.")
  if (length(lowess_fraction) != 1L || !is.finite(lowess_fraction) ||
      lowess_fraction <= 0 || lowess_fraction > 1) {
    stop("lowess_fraction must be in (0, 1].")
  }
  required <- c("participant_id", "age", "sex", "disease_free", artery_proteins)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }
  data <- as.data.frame(data)
  if (!nrow(data) || anyNA(data$participant_id) || anyDuplicated(data$participant_id)) {
    stop("Supply one row per participant with unique, nonmissing IDs.")
  }
  for (variable in c("age", "sex", "disease_free", artery_proteins)) {
    if (!is.numeric(data[[variable]]) || any(!is.finite(data[[variable]]))) {
      stop(variable, " must be numeric, finite, and complete. Prepare inputs upstream.")
    }
  }
  if (any(!data$sex %in% c(0, 1)) || any(!data$disease_free %in% c(0, 1))) {
    stop("sex and disease_free must be coded 0/1.")
  }

  # disease_free replaces the original check that all selected disease dates
  # are missing. It does NOT introduce a new baseline-only disease definition.
  healthy_rows <- which(data$disease_free == 1)
  if (length(healthy_rows) < 10L || length(unique(data$age[healthy_rows])) < 2L) {
    stop("Too few disease-free participants or insufficient age variation.")
  }
  predictors <- c(artery_proteins, "sex")
  # Use exactly the same numeric matrix and column order for fitting/prediction.
  # This avoids the original as.matrix/model.matrix encoding mismatch.
  x_full <- as.matrix(data[, predictors, drop = FALSE])
  x_healthy <- x_full[healthy_rows, , drop = FALSE]
  y_healthy <- data$age[healthy_rows]

  # Retain caret's age-stratified partition, the original seed, and 70/30 split.
  set.seed(seed)
  train_idx <- as.integer(caret::createDataPartition(
    y_healthy, p = train_fraction, list = FALSE
  ))
  test_idx <- setdiff(seq_along(healthy_rows), train_idx)
  if (length(train_idx) < nfolds || length(test_idx) < 2L) {
    stop("Insufficient training/test observations for the requested partition.")
  }
  x_train <- x_healthy[train_idx, , drop = FALSE]
  y_train <- y_healthy[train_idx]
  x_test <- x_healthy[test_idx, , drop = FALSE]
  y_test <- y_healthy[test_idx]
  if (length(unique(y_train)) < 2L || length(unique(y_test)) < 2L) {
    stop("Training and test sets must both have age variation.")
  }

  # The original implementation uses FIVE-fold CV, despite its ten-fold comment.
  # glmnet standardizes predictors internally; no full-cohort scaling is added.
  # Sex is penalized alongside proteins, preserving the original specification.
  cv_fit <- glmnet::cv.glmnet(
    x_train, y_train, family = "gaussian", alpha = 1,
    nfolds = nfolds, standardize = TRUE, intercept = TRUE
  )
  best_lambda <- cv_fit$lambda.min
  fit <- glmnet::glmnet(
    x_train, y_train, family = "gaussian", alpha = 1,
    lambda = best_lambda, standardize = TRUE, intercept = TRUE
  )

  # Evaluate raw predicted age in the held-out disease-free test set.
  pred_test <- as.numeric(stats::predict(fit, newx = x_test))
  test_cor <- if (stats::sd(pred_test) > 0) stats::cor(pred_test, y_test) else NA_real_
  test_rmse <- sqrt(mean((pred_test - y_test)^2))

  # Apply the training model to every participant, without refitting on all data.
  pred_full <- as.numeric(stats::predict(fit, newx = x_full))
  raw_delta <- pred_full - data$age

  # Preserve full-cohort LOWESS age-bias correction from the original script.
  # This is a descriptive full-cohort transformation, not test-set calibration.
  lowess_fit <- stats::lowess(data$age, pred_full, f = lowess_fraction, iter = 3)
  expected_age <- stats::approx(lowess_fit$x, lowess_fit$y,
                                xout = data$age, ties = mean)$y
  age_gap <- pred_full - expected_age
  gap_sd <- stats::sd(age_gap)
  if (!is.finite(gap_sd) || gap_sd <= 0) {
    stop("The age gap has zero/non-finite SD and cannot be standardized.")
  }
  gap_mean <- mean(age_gap)
  scores <- data.frame(
    participant_id = data$participant_id,
    age = data$age,
    arterial_age = pred_full,
    arterial_age_delta = raw_delta,
    arterial_age_expected = expected_age,
    arterial_age_gap = age_gap,
    arterial_age_gap_z = (age_gap - gap_mean) / gap_sd
  )
  coefficient_matrix <- as.matrix(stats::coef(fit))
  coefficients <- data.frame(
    term = rownames(coefficient_matrix),
    coefficient = as.numeric(coefficient_matrix[, 1]), row.names = NULL
  )
  partition <- data.frame(participant_id = data$participant_id,
                          subset = "outside_disease_free_subset")
  partition$subset[healthy_rows[train_idx]] <- "training"
  partition$subset[healthy_rows[test_idx]] <- "test"

  list(
    scores = scores,
    test_predictions = data.frame(
      participant_id = data$participant_id[healthy_rows[test_idx]],
      age = y_test, arterial_age = pred_test
    ),
    metrics = data.frame(n_disease_free = length(healthy_rows),
                         n_train = length(train_idx), n_test = length(test_idx),
                         test_cor = test_cor, test_rmse = test_rmse,
                         lambda_min = best_lambda),
    coefficients = coefficients, partition = partition,
    model = fit, cv_model = cv_fit, lowess_fit = lowess_fit,
    settings = list(proteins = artery_proteins, seed = seed,
                    train_fraction = train_fraction, nfolds = nfolds,
                    lowess_fraction = lowess_fraction,
                    gap_mean = gap_mean, gap_sd = gap_sd)
  )
}

# Example (run from the directory containing this script):
# source("02_construct_arterial_age.R")
# participants <- read.csv("data/demo/synthetic_participants.csv")
# age_results <- construct_arterial_age(participants)
# age_results$metrics
# head(age_results$scores)
# Later modules should join scores to other data using participant_id.
