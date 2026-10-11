# main.R - Analysis pipeline of the bachelor thesis "Light Exposure Across the Menstrual Cycle"
# Runs all scripts in order. Set the paths in 00_setup.R and LEAMC_GGIR_sleep.R first.

script_dir <- "."   # folder containing this file

scripts <- c(
  "00_setup.R",
  "00b_utils_io.R",
  "01_import_actlumus.R",
  "LEAMC_GGIR_sleep.R",
  "02_import_actigraph.R",
  "02b_import_activity_epochs.R",
  "03_import_redcap.R",
  "04_join_trim.R",
  "05_light_metrics.R",
  "05b_activity_metrics.R",
  "06_pacts.R",
  "06b_phase_assignment.R",
  "07_descriptives.R",
  "08_missing_data.R",
  "09_inferential.R",
  "09b_robustness.R",
  "10_visualisations.R",
  "10b_thesis_numbers.R",
  "revision.R",
  "12_gamm_ar1.R",
  "13_thesis_figures.R",
  "11_export.R"
)

for (f in c("model_helpers.R", "final_models.R", "phase_assignment.R", "figure_helpers.R"))
  source(file.path(script_dir, "R", f), echo = FALSE)

for (s in scripts) {
  if (s == "LEAMC_GGIR_sleep.R" && exists("RUN_GGIR") && !isTRUE(RUN_GGIR)) {
    cat("\n=====", s, "===== (skipped, RUN_GGIR = FALSE)\n")
    next
  }
  if (s == "02b_import_activity_epochs.R" && !exists("GGIR_OUT")) GGIR_OUT <- "C:/GGIR_local_output/LEATMC"   # as in LEAMC_GGIR_sleep.R
  cat("\n=====", s, "=====\n")
  tryCatch(
    source(file.path(script_dir, s), echo = FALSE),
    error = function(e) stop(sprintf("\n[stopped in %s]\n%s", s, conditionMessage(e)), call. = FALSE)
  )
}

cat("\nPipeline finished.\n")
