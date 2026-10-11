# 09_inferential.R - Linear mixed models (phase contrasts), GAMMs (continuous cycle time) and correlations

# ---- LMMs on participant x phase means ----

fit_lmm <- function(outcome, data=metrics_phase2) {
  safe(lmer(as.formula(paste(outcome, "~ cycle_phase_2 + (1|Id)")),
            data=filter(data, !is.na(cycle_phase_2)), REML=TRUE),
       label=paste("LMM", outcome))
}

lmm_TAT   <- fit_lmm("TAT_250_hr__mean")
lmm_MEDI  <- fit_lmm("MEDI_mean__mean")
lmm_lsens <- fit_lmm("light_sensitivity_idx__mean")
lmm_csymp <- fit_lmm("cycle_symptom_idx__mean")

# additional thresholds (10 and 1000 lx)
lmm_TAT10   <- fit_lmm("TAT_10_hr__mean")
lmm_TAT1000 <- fit_lmm("TAT_1000_hr__mean")

# high-completeness subsample
sens_compliance_TAT   <- fit_lmm("TAT_250_hr__mean",          filter(metrics_phase2, Id %in% hc_ids))
sens_compliance_MEDI  <- fit_lmm("MEDI_mean__mean",            filter(metrics_phase2, Id %in% hc_ids))
sens_compliance_lsens <- fit_lmm("light_sensitivity_idx__mean",filter(metrics_phase2, Id %in% hc_ids))
sens_compliance_csymp <- fit_lmm("cycle_symptom_idx__mean",    filter(metrics_phase2, Id %in% hc_ids))

# strict subsample (both phases)
sens_strict_TAT   <- fit_lmm("TAT_250_hr__mean",          filter(metrics_phase2, Id %in% hc_ids_strict))
sens_strict_MEDI  <- fit_lmm("MEDI_mean__mean",            filter(metrics_phase2, Id %in% hc_ids_strict))
sens_strict_lsens <- fit_lmm("light_sensitivity_idx__mean",filter(metrics_phase2, Id %in% hc_ids_strict))
sens_strict_csymp <- fit_lmm("cycle_symptom_idx__mean",    filter(metrics_phase2, Id %in% hc_ids_strict))

# item-level LMMs
lmm_ls_items   <- map(setNames(paste0("ls_",1:5,"__mean"),    paste0("ls_",1:5)),    fit_lmm)
lmm_drsp_items <- map(setNames(paste0(DRSP_ITEMS,"__mean"),   DRSP_ITEMS),           fit_lmm)

# IS/IV: participants with at least MIN_NPCRA_DAYS days in both phases
lmm_IS <- safe(lmer(IS~cycle_phase_2+(1|Id),
                     data=filter(npcra_phase, !is.na(cycle_phase_2), Id %in% npcra_ok_ids),
                     REML=TRUE),"LMM IS")
lmm_IV <- safe(lmer(IV~cycle_phase_2+(1|Id),
                     data=filter(npcra_phase, !is.na(cycle_phase_2), Id %in% npcra_ok_ids),
                     REML=TRUE),"LMM IV")
lmm_RA <- safe(lmer(RA~cycle_phase_2+(1|Id),data=filter(npcra_phase,!is.na(cycle_phase_2)),REML=TRUE),"LMM RA")

# IS/IV: all participants with a phase (sensitivity analysis)
lmm_IS_full21 <- safe(lmer(IS~cycle_phase_2+(1|Id),
                            data=filter(npcra_phase, !is.na(cycle_phase_2)), REML=TRUE),"LMM IS (voll, 21)")
lmm_IV_full21 <- safe(lmer(IV~cycle_phase_2+(1|Id),
                            data=filter(npcra_phase, !is.na(cycle_phase_2)), REML=TRUE),"LMM IV (voll, 21)")

# sleep (Prediction 4)
lmm_TST      <- fit_lmm("TST_min__mean")
lmm_sleep_eff <- fit_lmm("sleep_efficiency__mean")

# ---- GAMMs: cyclic spline over PACTS cycle time, random intercept ----

fit_gamm <- function(outcome, k = 20) {
  d <- filter(daily_data, !is.na(pacts_ov), !is.na(.data[[outcome]]))
  n_unique_ov <- n_distinct(d$pacts_ov)

  if (nrow(d) == 0) {
    message("  [SKIP] GAMM ", outcome, ": 0 Zeilen nach Filter — pacts_ov oder ",
            outcome, " sind (in Kombination) durchgehend NA.")
    return(NULL)
  }
  if (n_unique_ov < k) {
    message("  [SKIP] GAMM ", outcome, ": nur ", n_unique_ov,
            " eindeutige pacts_ov-Werte, aber k=", k, " angefordert. ",
            "k reduzieren (z.B. k=", max(3, n_unique_ov - 1), ") oder mehr Daten sammeln.")
    return(NULL)
  }

  d$Id <- factor(d$Id)

  safe(bam(as.formula(paste(outcome,
              paste0("~s(pacts_ov,bs='cc',k=", k, ")+s(Id,bs='re')"))),
           data=d, method="fREML"),
       label=paste("GAMM", outcome))
}

gamm_TAT        <- fit_gamm("TAT_250_hr")
gamm_lsens      <- fit_gamm("light_sensitivity_idx")
gamm_csymp      <- fit_gamm("cycle_symptom_idx")
gamm_ls_items   <- map(setNames(paste0("ls_",1:5),paste0("ls_",1:5)), fit_gamm)
gamm_drsp_items <- map(setNames(DRSP_ITEMS, DRSP_ITEMS), fit_gamm)

# relative amplitude (daily)
gamm_RA_act     <- fit_gamm("RA_act")

# additional thresholds
gamm_TAT10   <- fit_gamm("TAT_10_hr")
gamm_TAT1000 <- fit_gamm("TAT_1000_hr")

# EDF and p of the cycle smooth
extract_gamm_smooth <- function(model, label) {
  if (is.null(model)) return(tibble(item=label, edf=NA_real_, p=NA_real_))
  st <- summary(model)$s.table
  tibble(item=label, edf=round(st["s(pacts_ov)","edf"],2),
         p=round(st["s(pacts_ov)","p-value"],4))
}

main_gamm_summary <- bind_rows(
  extract_gamm_smooth(gamm_TAT,   "TAT_250_hr"),
  extract_gamm_smooth(gamm_lsens, "light_sensitivity_idx"),
  extract_gamm_smooth(gamm_csymp, "cycle_symptom_idx")
)
cat("[09_inferential] Haupt-GAMMs — EDF & p-Wert des s(pacts_ov)-Terms:\n")
print(main_gamm_summary)

ra_act_gamm_summary <- extract_gamm_smooth(gamm_RA_act, "RA_act")
cat("[09_inferential] RA_act-GAMM (Zusatzanalyse zu Prediction 4) — EDF & p-Wert:\n")
print(ra_act_gamm_summary)

tat_threshold_gamm_summary <- bind_rows(
  extract_gamm_smooth(gamm_TAT10,   "TAT_10_hr"),
  extract_gamm_smooth(gamm_TAT,     "TAT_250_hr"),
  extract_gamm_smooth(gamm_TAT1000, "TAT_1000_hr")
)
cat("[09_inferential] TAT-Schwellenwert-Zusatzanalyse (10/250/1000 lx):\n")
print(tat_threshold_gamm_summary)
cat("[09_inferential] LMM TAT_10_hr:\n");   print(summary(lmm_TAT10))
cat("[09_inferential] LMM TAT_1000_hr:\n"); print(summary(lmm_TAT1000))

# ---- Secondary light metrics ----
SECONDARY_LIGHT_VARS <- c(
  "MEDI_median", "MEDI_max", "MEDI_sd", "MEDI_daytime", "MEDI_night",
  "first_250", "last_250",
  "M10_mean", "M10_onset", "M10_midpoint",
  "L5_mean", "L5_onset",
  "RA_light", "centroid_LE_hr"
)

gamm_secondary_light <- map(setNames(SECONDARY_LIGHT_VARS, SECONDARY_LIGHT_VARS), fit_gamm)
lmm_secondary_light  <- map(setNames(paste0(SECONDARY_LIGHT_VARS, "__mean"), SECONDARY_LIGHT_VARS), fit_lmm)

# phase effect of an LMM
extract_lmm_fixed <- function(model, label) {
  if (is.null(model)) return(tibble(item=label, b=NA_real_, t=NA_real_, df=NA_real_, p=NA_real_))
  s <- summary(model)$coefficients
  row <- s[grepl("luteal", rownames(s)), , drop=FALSE]
  if (nrow(row) == 0) return(tibble(item=label, b=NA_real_, t=NA_real_, df=NA_real_, p=NA_real_))
  tibble(item=label, b=round(row[1,"Estimate"],4), t=round(row[1,"t value"],2),
         df=round(row[1,"df"],2), p=round(row[1,"Pr(>|t|)"],4))
}

secondary_light_gamm_summary <- map_dfr(names(gamm_secondary_light),
                                         ~ extract_gamm_smooth(gamm_secondary_light[[.x]], .x)) |>
  arrange(p)
secondary_light_lmm_summary <- map_dfr(names(lmm_secondary_light),
                                        ~ extract_lmm_fixed(lmm_secondary_light[[.x]], .x)) |>
  arrange(p)

cat("[09_inferential] Sekundaer-Lichtmetriken — GAMM (EDF/p, sortiert):\n")
print(secondary_light_gamm_summary, n = Inf)
cat("[09_inferential] Sekundaer-Lichtmetriken — LMM (b/t/df/p, sortiert):\n")
print(secondary_light_lmm_summary, n = Inf)

# item GAMMs sorted by p
item_p_adjust <- function(gamm_list) {
  raw <- vapply(gamm_list, function(m)
    if (is.null(m)) NA_real_ else summary(m)$s.table["s(pacts_ov)", "p-value"], numeric(1))
  tibble(item = names(raw),
         p_holm = round(p.adjust(raw, method = "holm"), 4),
         p_bh   = round(p.adjust(raw, method = "BH"), 4))
}

drsp_gamm_summary <- map_dfr(names(gamm_drsp_items),
                              ~ extract_gamm_smooth(gamm_drsp_items[[.x]], .x)) |>
  left_join(item_p_adjust(gamm_drsp_items), by = "item") |>
  arrange(p)

ls_gamm_summary <- map_dfr(names(gamm_ls_items),
                            ~ extract_gamm_smooth(gamm_ls_items[[.x]], .x)) |>
  left_join(item_p_adjust(gamm_ls_items), by = "item") |>
  arrange(p)
cat("[09_inferential] Light-Sensitivity-Items — EDF, p-Wert (unkorrigiert), Holm- und BH-adjustiert (Familie = 5 Items):\n")
print(ls_gamm_summary)
cat("[09_inferential] DRSP-Items — EDF, p-Wert (unkorrigiert), Holm- und BH-adjustiert (Familie = 20 Items):\n")
print(drsp_gamm_summary, n = Inf)

# ---- Spearman correlations (pooled daily data) ----

corr_vars <- c("TAT_250_hr", "MEDI_mean", "MEDI_morning", "MEDI_evening",
               "light_sensitivity_idx", "cycle_symptom_idx",
               "TST_min", "sleep_efficiency", "WASO_min",
               "drsp_fatigue", "drsp_hypersomnia", "drsp_insomnia")

corr_data <- daily_data |> select(all_of(corr_vars)) |> drop_na()

cat("[09_inferential] Korrelationsmatrix:", nrow(corr_data), "vollstaendige Zeilen von",
    nrow(daily_data), "Tagen (Listwise Deletion ueber alle", length(corr_vars), "Variablen)\n")

if (nrow(corr_data) < 3) {
  na_pct <- daily_data |> select(all_of(corr_vars)) |>
    summarise(across(everything(), ~ round(100*mean(is.na(.x)), 1))) |>
    pivot_longer(everything(), names_to="variable", values_to="pct_NA") |>
    arrange(desc(pct_NA))
  message("  [WARNUNG] corr_mat: nur ", nrow(corr_data), " vollstaendige Zeilen. ",
          "% fehlend je Variable (absteigend):")
  print(na_pct, n = Inf)
  corr_mat <- matrix(NA_real_, nrow=length(corr_vars), ncol=length(corr_vars),
                      dimnames=list(corr_vars, corr_vars))
} else {
  corr_mat <- cor(corr_data, method="spearman")
}

cat("[09_inferential] Fertig — LMMs, GAMMs und Korrelationsmatrix berechnet\n")
