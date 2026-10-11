# 11_export.R - Export of data sets and workspace; copies all outputs to the sync folder

write_csv(daily_data,  file.path(OUTPUT_DIR,"daily_data.csv"))
write_csv(bin_summary, file.path(OUTPUT_DIR,"bin_summary.csv"))
write_csv(char_table,  file.path(OUTPUT_DIR,"participant_characteristics.csv"))
if (exists("EXCLUSIONS") && nrow(EXCLUSIONS) > 0) {
  write_csv(EXCLUSIONS, file.path(OUTPUT_DIR,"QC_excluded_participants.csv"))
}
if (exists("pacts_phase_diagnosis") && nrow(pacts_phase_diagnosis) > 0) {
  write_csv(pacts_phase_diagnosis, file.path(OUTPUT_DIR,"QC_pacts_phase_diagnosis.csv"))
}
if (exists("ids_dropped_pacts_info") && nrow(ids_dropped_pacts_info) > 0) {
  write_csv(ids_dropped_pacts_info, file.path(OUTPUT_DIR,"QC_pacts_dropped_participants.csv"))
}
if (exists("menses_gaps") && nrow(menses_gaps) > 0) {
  write_csv(menses_gaps, file.path(OUTPUT_DIR,"QC_menses_gaps_implausible.csv"))
}
if (exists("actigraph_coverage_gap") && nrow(actigraph_coverage_gap) > 0) {
  write_csv(actigraph_coverage_gap, file.path(OUTPUT_DIR,"QC_actigraph_coverage_gap.csv"))
}
if (exists("npcra_phase") && nrow(npcra_phase) > 0) {
  write_csv(npcra_phase, file.path(OUTPUT_DIR,"QC_npcra_phase_daycounts.csv"))
}

save(daily_data, light_trimmed, actigraph_sleep, surveys, pacts_output,
     sens_compliance_TAT, sens_compliance_MEDI, sens_compliance_lsens, sens_compliance_csymp,
     sens_strict_TAT, sens_strict_MEDI, sens_strict_lsens, sens_strict_csymp, hc_ids_strict,
     bin_summary, char_table, metrics_phase2,
     lmm_TAT, lmm_MEDI, lmm_lsens, lmm_csymp, lmm_TAT10, lmm_TAT1000,
     lmm_IS, lmm_IV, lmm_RA, lmm_TST, lmm_sleep_eff, lmm_ls_items, lmm_drsp_items,
     gamm_TAT, gamm_lsens, gamm_csymp, gamm_ls_items, gamm_drsp_items,
     main_gamm_summary, tat_threshold_gamm_summary, secondary_light_gamm_summary, secondary_light_lmm_summary, drsp_gamm_summary, ls_gamm_summary, corr_mat,
     file=file.path(OUTPUT_DIR,"LEAMC_workspace.RData"))

cat("\n[11_export] LEAMC-Analyse abgeschlossen.\n")
cat("  Lokaler Output ->", OUTPUT_DIR, "\n")

# ---- Copy outputs to SYNC_OUTPUT_DIR ----
copy_with_retry <- function(from_dir, to_dir, max_tries = 5, wait_seconds = 3) {
  dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
  for (attempt in seq_len(max_tries)) {
    ok <- tryCatch({
      file.copy(list.files(from_dir, full.names = TRUE), to_dir,
                recursive = TRUE, overwrite = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (isTRUE(all(ok))) return(invisible(TRUE))
    if (attempt < max_tries) {
      message("  [Retry ", attempt, "/", max_tries, "] Kopieren nach '", to_dir,
              "' teilweise fehlgeschlagen - vermutlich Sync-Client. Warte ",
              wait_seconds, "s und versuche erneut...")
      Sys.sleep(wait_seconds)
    } else {
      warning("Kopieren nach '", to_dir, "' nach ", max_tries, " Versuchen nicht ",
              "vollstaendig erfolgreich. Lokale Kopie unter '", from_dir,
              "' bleibt in jedem Fall vollstaendig erhalten.")
    }
  }
}

copy_with_retry(OUTPUT_DIR, SYNC_OUTPUT_DIR)
cat("  Sync-Output ->", SYNC_OUTPUT_DIR, "\n")
