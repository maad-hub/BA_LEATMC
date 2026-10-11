# 12_gamm_ar1.R - Final models for both phase assignments ("A" = main analysis, "original" = sensitivity analysis):
# LMM phase contrasts, GAMMs without and with AR(1) residuals, main tests P1-P4 (P3: GAMM with AR(1)),
# Table 3, Cycle Symptom Index rows of Table C1, covariate-adjusted GAMMs.
# Uses the data after 06_pacts.R (daily_data_06); results in OUTPUT_DIR/final_models/.

FINAL_DIR <- file.path(OUTPUT_DIR, "final_models")
dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)
write_final <- function(x, name) {
  if (!is.null(x)) utils::write.csv(x, file.path(FINAL_DIR, name), row.names = FALSE, fileEncoding = "UTF-8")
}

final_dd <- lapply(phase_tabs[c("A", "original")], function(t) ac_variant_data(daily_data_06, t, date_ranges, tz = TIMEZONE))
final_models <- lapply(names(final_dd), function(v) e2_main_models(final_dd[[v]], v, activity_trimmed))
names(final_models) <- names(final_dd)
table3 <- e2_table3(final_models)
adjusted_gamms <- e2_adjusted_gamms(final_dd)

# participants per PACTS bin (as 10b_thesis_numbers.R, block "bins") for both assignments
pacts_bins <- do.call(rbind, lapply(names(final_dd), function(v) {
  d <- final_dd[[v]][!is.na(final_dd[[v]]$cycle_bin_ov), ]
  n <- tapply(d$Id, d$cycle_bin_ov, function(x) length(unique(x)))
  data.frame(variant = v, n_bins = max(d$cycle_bin_ov), bin1_n_participants = unname(n["1"]),
             last_bin_n_participants = unname(n[as.character(max(d$cycle_bin_ov))]),
             n_participants_with_bins = length(unique(d$Id)))
}))

for (v in names(final_models)) for (n in names(final_models[[v]])) write_final(final_models[[v]][[n]], paste0(v, "_", n, ".csv"))
write_final(table3, "table3_both_assignments.csv")
write_final(adjusted_gamms, "adjusted_gamm.csv")
write_final(pacts_bins, "pacts_bins_by_assignment.csv")
cat("[12_gamm_ar1] results in", FINAL_DIR, "\n")
