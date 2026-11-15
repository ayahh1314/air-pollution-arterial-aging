# API-related metabolomic susceptibility score (MSS), artery model only.
# Requires prepared/imputed metabolites; no participant data are included.
# Keep screening, mediation, score construction, and interaction results separate.

construct_mss <- function(data, metabolites = paste0("metabolite_", 1:10),
                          sims = 200, seed = 123,
                          exclude_centers = c("11022", "11023", "10003")) {
  if (!requireNamespace("mediation", quietly = TRUE)) {
    stop("Install the mediation package before running this module.")
  }
  covariates <- c("age", "sex", "race", "center", "income", "education",
                  "smoking", "drinking")
  outcome <- "arterial_age_gap_z"
  required <- c("participant_id", "api", outcome, covariates, metabolites)
  if (!is.character(metabolites) || !length(metabolites) ||
      anyNA(metabolites) || anyDuplicated(metabolites) ||
      any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", metabolites)) ||
      any(metabolites %in% c("participant_id", "api", outcome, covariates))) {
    stop("Use unique metabolite column names containing letters, digits, or underscores.")
  }
  if (length(sims) != 1L || !is.finite(sims) || sims < 2 || sims != floor(sims)) {
    stop("sims must be an integer >= 2.")
  }
  missing_columns <- setdiff(required, names(data))
  if (length(missing_columns)) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "))
  }
  data <- as.data.frame(data)
  if (!nrow(data) || anyNA(data$participant_id) || anyDuplicated(data$participant_id)) {
    stop("Supply one row per participant with unique nonmissing IDs.")
  }
  for (variable in c("api", outcome, "age", "sex", metabolites)) {
    if (!is.numeric(data[[variable]]) ||
        any(!is.finite(data[[variable]]) & !is.na(data[[variable]]))) {
      stop(variable, " must be numeric and cannot contain infinite values.")
    }
  }
  if (any(!is.na(data$sex) & !data$sex %in% c(0, 1))) stop("sex must be 0/1.")
  for (variable in c("race", "center", "income", "education", "smoking", "drinking")) {
    data[[variable]] <- factor(data[[variable]])
  }
  make_formula <- function(response, predictors) {
    stats::reformulate(predictors, response = response)
  }
  fit_lm <- function(response, predictors, model_data) {
    fit <- stats::lm(make_formula(response, predictors), data = model_data,
                     na.action = stats::na.omit)
    # Store the actual formula for bootstrap model refits outside this helper.
    fit$call$formula <- stats::formula(fit)
    fit
  }
  extract_term <- function(fit, term) {
    table <- summary(fit)$coefficients
    if (!term %in% rownames(table) || !is.finite(table[term, "Estimate"])) {
      stop("Non-estimable term: ", term)
    }
    data.frame(term = term, estimate = unname(table[term, "Estimate"]),
               standard_error = unname(table[term, "Std. Error"]),
               p_value = unname(table[term, "Pr(>|t|)"]),
               n = stats::nobs(fit), row.names = NULL)
  }

  # 1. API -> standardized arterial age gap.
  main_model <- fit_lm(outcome, c("api", covariates), data)
  main_association <- extract_term(main_model, "api")
  api_iqr <- unname(stats::IQR(data$api, na.rm = TRUE))
  if (!is.finite(api_iqr) || api_iqr <= 0) stop("API IQR must be positive.")
  # Do not round IQR before effect scaling; round only the presentation table.
  main_association$api_iqr <- api_iqr
  main_association$estimate_per_iqr <- main_association$estimate * api_iqr
  main_association$se_per_iqr <- main_association$standard_error * api_iqr

  # 2. Two separate marginal screening models for each metabolite.
  # The metabolite -> age-gap screening model does not adjust for API,
  # matching the supplied original screening specification.
  api_metabolite <- do.call(rbind, lapply(metabolites, function(metabolite) {
    result <- extract_term(fit_lm(metabolite, c("api", covariates), data), "api")
    result$metabolite <- metabolite
    result
  }))
  metabolite_age <- do.call(rbind, lapply(metabolites, function(metabolite) {
    result <- extract_term(fit_lm(outcome, c(metabolite, covariates), data), metabolite)
    result$metabolite <- metabolite
    result
  }))
  # BH correction separately across metabolites in each screening family.
  api_metabolite$fdr <- stats::p.adjust(api_metabolite$p_value, method = "BH")
  metabolite_age$fdr <- stats::p.adjust(metabolite_age$p_value, method = "BH")
  candidates <- metabolites[
    api_metabolite$fdr < 0.05 & metabolite_age$fdr < 0.05
  ]

  # Preserve the original stage-specific center exclusion: screening above
  # precedes this exclusion; mediation, score, and interaction follow it.
  keep <- !is.na(data$center) & !as.character(data$center) %in% exclude_centers
  analysis_data <- droplevels(data[keep, , drop = FALSE])
  result <- list(main_association = main_association,
                 api_metabolite = api_metabolite,
                 metabolite_age = metabolite_age,
                 candidates = candidates,
                 mediation = data.frame(), selection = data.frame(),
                 weights = data.frame(), scores = data.frame(),
                 interaction_results = NULL,
                 settings = list(seed = seed, sims = sims,
                                 exclude_centers = exclude_centers,
                                 n_screening = nrow(data),
                                 n_after_center_exclusion = nrow(analysis_data)))
  if (!length(candidates)) {
    result$status <- "No metabolites passed both screening FDR thresholds."
    return(result)
  }

  # 3. Statistical mediation: API -> metabolite -> arterial age gap.
  # All models for a given mediator use the same complete-case sample.
  # Retain the original unit contrast (API 0 -> 1), 200 bootstrap draws,
  # percentile confidence intervals, and no exposure-mediator interaction.
  set.seed(seed)
  mediation_rows <- lapply(candidates, function(metabolite) {
    columns <- c(outcome, "api", metabolite, covariates)
    final_data <- droplevels(analysis_data[
      stats::complete.cases(analysis_data[, columns, drop = FALSE]), , drop = FALSE
    ])
    mediator_model <- fit_lm(metabolite, c("api", covariates), final_data)
    outcome_model <- fit_lm(outcome, c(metabolite, "api", covariates), final_data)
    med_fit <- mediation::mediate(
      mediator_model, outcome_model, treat = "api", mediator = metabolite,
      control.value = 0, treat.value = 1,
      boot = TRUE, boot.ci.type = "perc", sims = sims
    )
    data.frame(
      metabolite = metabolite, n = nrow(final_data),
      total_effect = med_fit$tau.coef, total_p = med_fit$tau.p,
      total_ci_low = med_fit$tau.ci[1], total_ci_high = med_fit$tau.ci[2],
      acme = med_fit$d.avg, acme_p = med_fit$d.avg.p,
      acme_ci_low = med_fit$d.avg.ci[1], acme_ci_high = med_fit$d.avg.ci[2],
      ade = med_fit$z.avg, ade_p = med_fit$z.avg.p,
      ade_ci_low = med_fit$z.avg.ci[1], ade_ci_high = med_fit$z.avg.ci[2],
      proportion_mediated = med_fit$n.avg,
      proportion_p = med_fit$n.avg.p,
      proportion_ci_low = med_fit$n.avg.ci[1],
      proportion_ci_high = med_fit$n.avg.ci[2], row.names = NULL
    )
  })
  mediation_table <- do.call(rbind, mediation_rows)
  mediation_table$acme_fdr <- stats::p.adjust(mediation_table$acme_p, method = "BH")
  mediation_table$proportion_fdr <- stats::p.adjust(mediation_table$proportion_p,
                                                  method = "BH")

  # Original selection uses UNADJUSTED total-effect p, ACME FDR, and PM > 5%.
  # The computed proportion FDR is reported but not a selection condition.
  mediation_table$selected <- with(mediation_table,
    is.finite(total_p) & is.finite(acme_fdr) & is.finite(proportion_mediated) &
      total_p < 0.05 & acme_fdr < 0.05 & proportion_mediated > 0.05)
  selected <- mediation_table[mediation_table$selected, , drop = FALSE]
  result$mediation <- mediation_table
  result$selection <- selected
  if (!nrow(selected)) {
    result$status <- "No metabolites passed the mediation selection thresholds."
    return(result)
  }

  # 4. Retain the supplied direction rule, not abs(ACME):
  # weight = -ACME if beta_MY < 0; otherwise weight = ACME.
  selected_metabolites <- selected$metabolite
  beta_my <- metabolite_age$estimate[match(selected_metabolites,
                                          metabolite_age$metabolite)]
  acme <- selected$acme
  weight <- ifelse(beta_my < 0, -acme, acme)
  denominator <- sum(acme)
  if (!is.finite(denominator) ||
      abs(denominator) <= 1e-10 * max(1, sum(abs(acme)))) {
    stop("The sum of selected ACMEs is too close to zero to construct MSS.")
  }
  score <- as.vector(as.matrix(analysis_data[, selected_metabolites, drop = FALSE]) %*%
                       weight) / denominator
  # The supplied code does not standardize metabolites or the resulting score.
  # Keep this raw score as mss_raw; expose a separately named z score for review.
  score_sd <- stats::sd(score, na.rm = TRUE)
  if (!is.finite(score_sd) || score_sd <= 0) stop("MSS has zero/non-finite SD.")
  analysis_data$mss_raw <- score
  analysis_data$mss_z <- (score - mean(score, na.rm = TRUE)) / score_sd
  api_cut <- stats::median(analysis_data$api, na.rm = TRUE)
  mss_cut <- stats::median(score, na.rm = TRUE)
  analysis_data$api_high <- as.integer(analysis_data$api >= api_cut)
  analysis_data$mss_high <- as.integer(score >= mss_cut)
  analysis_data$joint_group <- factor(
    1L + 2L * analysis_data$api_high + analysis_data$mss_high,
    levels = 1:4, labels = c("low_api_low_mss", "low_api_high_mss",
                            "high_api_low_mss", "high_api_high_mss")
  )
  result$weights <- data.frame(metabolite = selected_metabolites,
                               beta_my = beta_my, acme = acme,
                               direction_adjusted_weight = weight,
                               normalized_weight = weight / denominator)
  result$scores <- analysis_data[, c("participant_id", "api", outcome,
    "mss_raw", "mss_z", "api_high", "mss_high", "joint_group"), drop = FALSE]

  # 5. Original joint and interaction models use raw MSS and raw API.
  joint_model <- fit_lm(outcome, c("joint_group", covariates), analysis_data)
  binary_model <- fit_lm(outcome, c("mss_high * api_high", covariates), analysis_data)
  continuous_model <- fit_lm(outcome, c("mss_raw * api", covariates), analysis_data)
  result$interaction_results <- rbind(
    extract_term(binary_model, "mss_high:api_high"),
    extract_term(continuous_model, "mss_raw:api")
  )
  result$models <- list(main = main_model, joint = joint_model,
                        binary_interaction = binary_model,
                        continuous_interaction = continuous_model)
  result$settings$api_median <- api_cut
  result$settings$mss_median <- mss_cut
  result$settings$acme_denominator <- denominator
  result$status <- "MSS constructed."
  result
}

# Example after joining Parts 1 and 2 by participant_id:
# source("03_construct_mss.R")
# mss_results <- construct_mss(participants,
#   metabolites = paste0("metabolite_", 1:10))
# mss_results$status
# mss_results$selection
# mss_results$interaction_results
