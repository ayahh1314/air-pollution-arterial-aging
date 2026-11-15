# Mortality-weighted air pollution index (API)
# Input: a prepared participant-level data.frame; see README.md.
# Concentrations remain in their original units (micrograms/m^3).
# This module does not read or redistribute UK Biobank data.

construct_api <- function(data) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Install the survival package before running this module.")
  }

  pollutants <- c("pm25", "pm10", "no2", "nox")
  covariates <- c("age", "sex", "race", "center", "education", "income",
                  "smoking", "drinking")
  required <- c("participant_id", "followup_days", "death_event",
                pollutants, covariates)
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }
  data <- as.data.frame(data)
  if (!nrow(data) || anyNA(data$participant_id) ||
      anyDuplicated(data$participant_id)) {
    stop("Input must contain one row per participant with unique, nonmissing IDs.")
  }
  numeric_columns <- c("followup_days", "death_event", pollutants, "age", "sex")
  for (variable in numeric_columns) {
    if (!is.numeric(data[[variable]])) stop(variable, " must be numeric.")
    if (any(!is.finite(data[[variable]]) & !is.na(data[[variable]]))) {
      stop(variable, " contains non-finite values.")
    }
  }
  if (any(!is.na(data$death_event) & !data$death_event %in% c(0, 1))) {
    stop("death_event must be coded 0/1 or NA.")
  }
  if (any(!is.na(data$sex) & !data$sex %in% c(0, 1))) {
    stop("sex must be coded 0/1 or NA.")
  }
  if (any(data$followup_days <= 0, na.rm = TRUE)) {
    stop("followup_days must be positive; resolve invalid follow-up upstream.")
  }
  for (variable in pollutants) {
    if (any(data[[variable]] < 0, na.rm = TRUE)) {
      stop(variable, " contains negative concentrations.")
    }
  }

  # Treat smoking and drinking as categorical covariates, alongside
  # the original categorical adjustment variables.
  # Sex is a numeric binary covariate; age is continuous and linear.
  for (variable in c("race", "center", "education", "income",
                     "smoking", "drinking")) {
    data[[variable]] <- factor(data[[variable]])
  }

  # Fit a separate adjusted Cox model for each pollutant.
  # Each fit uses its own complete cases, as in the original code.
  # Do not impose one common complete-case sample across all four models.
  fit_model <- function(exposure, model_data) {
    formula <- stats::as.formula(paste(
      "survival::Surv(followup_days, death_event) ~",
      paste(c(exposure, covariates), collapse = " + ")
    ))
    survival::coxph(formula, data = model_data, ties = "efron",
                    na.action = stats::na.omit)
  }

  extract_result <- function(fit, exposure) {
    coefficients <- summary(fit)$coefficients
    beta <- unname(stats::coef(fit)[exposure])
    if (length(beta) != 1L || !is.finite(beta)) {
      stop("Non-estimable exposure coefficient: ", exposure)
    }
    data.frame(
      exposure = exposure,
      beta = beta,
      standard_error = unname(coefficients[exposure, "se(coef)"]),
      hazard_ratio = exp(beta),
      p_value = unname(coefficients[exposure, "Pr(>|z|)"]),
      n = fit$n,
      events = fit$nevent,
      row.names = NULL
    )
  }

  pollutant_models <- setNames(lapply(pollutants, fit_model,
                                     model_data = data), pollutants)
  weights <- do.call(rbind, lapply(pollutants, function(variable) {
    extract_result(pollutant_models[[variable]], variable)
  }))
  beta <- setNames(weights$beta, weights$exposure)
  beta_sum <- sum(beta)
  if (abs(beta_sum) <= 1e-10 * max(1, sum(abs(beta)))) {
    stop("The sum of coefficients is too close to zero to normalize the API.")
  }
  if (any(beta < 0)) {
    warning("Negative mortality weights were estimated; inspect model results. ",
            "Weights have not been truncated or replaced.")
  }

  # Preserve the original index exactly:
  # API = 4 * sum(beta_j * concentration_j) / sum(beta_j).
  # Extract by name rather than model position; retain full precision.
  # Missing pollutant concentrations yield a missing API (no partial sum).
  weights$api_weight <- beta * (length(pollutants) / beta_sum)
  data$api <- as.vector(as.matrix(data[, pollutants, drop = FALSE]) %*%
                          (beta * length(pollutants) / beta_sum))

  # Retain the original additional model assessing API and mortality.
  # This is a model in the weight-estimation sample, not external validation.
  api_model <- fit_model("api", data)
  api_association <- extract_result(api_model, "api")

  list(data = data, weights = weights, api_association = api_association,
       pollutant_models = pollutant_models, api_model = api_model)
}

# Example usage after the synthetic data have been prepared:
# source("01_construct_api.R")
# participants <- read.csv("data/demo/synthetic_participants.csv")
# api_results <- construct_api(participants)
# participants_with_api <- api_results$data
# api_results$weights
# api_results$api_association
# Save only synthetic outputs in the public demonstration repository.
