# Run from the extracted project root: source("run_demo.R")
# CSV inputs are already supplied; Python is not required to run the analysis.
required_packages <- c("survival", "caret", "glmnet", "mediation")
missing_packages <- required_packages[!vapply(required_packages,
  requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}
source("scripts/01_construct_api.R")
source("scripts/02_construct_arterial_age.R")
source("scripts/03_construct_mss.R")
participants <- read.csv("data/demo/synthetic_participants.csv",
                         stringsAsFactors = FALSE)
dir.create("results", showWarnings = FALSE)

message("Part 1: Estimate API weights and calculate API")
api_results <- construct_api(participants)
participants <- api_results$data
write.csv(api_results$weights, "results/demo_api_weights.csv", row.names = FALSE)

message("Part 2: Train the artery age model and calculate age gap")
age_results <- construct_arterial_age(participants,
  artery_proteins = paste0("artery_protein_", 1:5))
matched_rows <- match(participants$participant_id, age_results$scores$participant_id)
stopifnot(!anyNA(matched_rows))
for (variable in setdiff(names(age_results$scores), c("participant_id", "age"))) {
  participants[[variable]] <- age_results$scores[[variable]][matched_rows]
}
write.csv(age_results$metrics, "results/demo_arterial_age_metrics.csv", row.names = FALSE)
write.csv(age_results$scores, "results/demo_arterial_age_scores.csv", row.names = FALSE)

message("Part 3: Screen metabolites, fit mediation models, construct MSS")
mss_results <- construct_mss(participants,
  metabolites = paste0("metabolite_", 1:12), sims = 200, seed = 123)
write.csv(mss_results$api_metabolite, "results/demo_api_metabolite.csv", row.names = FALSE)
write.csv(mss_results$metabolite_age, "results/demo_metabolite_age.csv", row.names = FALSE)
write.csv(mss_results$mediation, "results/demo_mediation.csv", row.names = FALSE)
write.csv(mss_results$weights, "results/demo_mss_weights.csv", row.names = FALSE)
write.csv(mss_results$scores, "results/demo_mss_scores.csv", row.names = FALSE)
if (!is.null(mss_results$interaction_results)) {
  write.csv(mss_results$interaction_results,
            "results/demo_interaction_results.csv", row.names = FALSE)
}
saveRDS(list(api = api_results, arterial_age = age_results, mss = mss_results),
        "results/demo_analysis_objects.rds")
capture.output(sessionInfo(), file = "results/demo_session_info.txt")
message(mss_results$status)
message("Finished. All results are synthetic demonstrations, not manuscript findings.")
