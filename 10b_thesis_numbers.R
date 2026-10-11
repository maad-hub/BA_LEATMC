# 10b_thesis_numbers.R - Collects the numbers reported in the results in thesis_numbers.csv

.rows <- list()

add <- function(section, key, value, note = "") {
  if (length(value) == 0) value <- NA
  v <- if (is.numeric(value)) as.character(signif(value[1], 8)) else as.character(value[1])
  .rows[[length(.rows) + 1]] <<- data.frame(section = section, key = key, value = v,
                                           note = note, stringsAsFactors = FALSE)
}

block <- function(name, expr) {
  tryCatch(expr, error = function(e) add("export_error", name, conditionMessage(e)))
}

getobj <- function(name) if (exists(name)) get(name) else NULL

# ---- LMM row: phase means, effect, CI, t, df, p, N ----
lmm_row <- function(m, label, section = "lmm") {
  tryCatch({
    if (is.null(m)) { add(section, paste0(label, ".missing"), 1, "Modell ist NULL"); return(invisible()) }
    s  <- summary(m)$coefficients
    rn <- grep("luteal", rownames(s), value = TRUE)[1]
    b  <- s[rn, "Estimate"]
    ci_method <- "profile"
    ci <- tryCatch(suppressWarnings(suppressMessages(
            confint(m, parm = rn, level = 0.95, method = "profile", quiet = TRUE))),
          error = function(e) { ci_method <<- "Wald"
            suppressWarnings(confint(m, parm = rn, level = 0.95, method = "Wald")) })
    ci <- as.numeric(ci[1, ])
    msg <- tryCatch(paste(m@optinfo$conv$lme4$messages, collapse = "; "), error = function(e) "")
    k <- function(x) paste0(label, ".", x)
    add(section, k("follicular_mean"), s["(Intercept)", "Estimate"], "Intercept = follikulaer")
    add(section, k("luteal_mean"),     s["(Intercept)", "Estimate"] + b)
    add(section, k("b"),      b)
    add(section, k("se"),     s[rn, "Std. Error"])
    add(section, k("ci_low"), ci[1], paste("95%-CI,", ci_method))
    add(section, k("ci_high"), ci[2])
    add(section, k("t"),      s[rn, "t value"])
    add(section, k("df"),     s[rn, "df"])
    add(section, k("p"),      s[rn, "Pr(>|t|)"])
    add(section, k("n_obs"),  nobs(m))
    add(section, k("n_participants"), unname(lme4::ngrps(m)["Id"]))
    add(section, k("singular"), lme4::isSingular(m))
    add(section, k("conv_message"), if (nzchar(msg)) msg else "none")
  }, error = function(e) add("export_error", paste0(section, ".", label), conditionMessage(e)))
}

# ---- GAMM row: EDF and p ----
gamm_row <- function(m, label, section = "gamm") {
  tryCatch({
    if (is.null(m)) { add(section, paste0(label, ".missing"), 1, "Modell ist NULL"); return(invisible()) }
    st <- summary(m)$s.table
    add(section, paste0(label, ".edf"), st["s(pacts_ov)", "edf"])
    add(section, paste0(label, ".p"),   st["s(pacts_ov)", "p-value"])
    add(section, paste0(label, ".n_obs"), nrow(m$model))
    add(section, paste0(label, ".n_participants"), length(unique(m$model$Id)))
  }, error = function(e) add("export_error", paste0(section, ".", label), conditionMessage(e)))
}

# ---- 0) Software versions ----
block("meta", {
  add("meta", "export_time", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  add("meta", "r_version", R.version.string)
  add("meta", "r_citation_year", as.character(citation()$year))
  for (pkg in c("GGIR", "lme4", "lmerTest", "mgcv", "LightLogR", "menstrualcycleR", "dplyr", "psych"))
    add("meta", paste0("pkg.", pkg), tryCatch(as.character(packageVersion(pkg)), error = function(e) NA))
})

# ---- 1) Sample, Table 1, Cronbach's alpha ----
block("sample", {
  all_ids <- unique(c(surveys$Id, onboarding$Id, if (exists("EXCLUSIONS")) EXCLUSIONS$Id))
  all_ids <- all_ids[!is.na(all_ids)]
  add("sample", "n_enrolled", length(all_ids))
  add("sample", "n_analytic", nrow(char_table))
  valid_ids <- pacts_phase_diagnosis |> filter(has_valid_phase, Id %in% char_table$Id) |> pull(Id)
  add("sample", "n_valid_phase", length(valid_ids))
  add("sample", "n_no_phase", nrow(char_table) - length(valid_ids))
  add("sample", "pct_valid_phase", round(100 * length(valid_ids) / nrow(char_table), 1))
  add("sample", "n_age_outside_18_35", sum(char_table$age < 18 | char_table$age > 35, na.rm = TRUE))
  add("sample", "n_bmi_outside_18_26", sum(char_table$bmi < 18 | char_table$bmi > 26, na.rm = TRUE))
  add("sample", "n_pms_diagnosed",  sum(char_table$pms_diagnosed  == 1, na.rm = TRUE))
  add("sample", "n_pmdd_diagnosed", sum(char_table$pmdd_diagnosed == 1, na.rm = TRUE))

  for (v in c("age", "height_cm", "weight_kg", "bmi", "n_cycles", "n_cycles_completed", "mean_cycle_length"))
    for (col in c("n", "mean", "sd", "median", "min", "max"))
      add("table1", paste0(v, ".", col), participant_desc[v, col])

  a_ls   <- suppressWarnings(suppressMessages(psych::alpha(tidyr::drop_na(dplyr::select(surveys, ls_1:ls_5)))$total$raw_alpha))
  a_drsp <- suppressWarnings(suppressMessages(psych::alpha(tidyr::drop_na(dplyr::select(surveys, dplyr::all_of(DRSP_ITEMS))))$total$raw_alpha))
  add("sample", "alpha_light_sensitivity", a_ls)
  add("sample", "alpha_cycle_symptom", a_drsp)

  no_ids <- setdiff(char_table$Id, valid_ids)
  gaps <- pacts_input |> filter(menses == 1L, id %in% no_ids) |> arrange(id, date) |>
    group_by(id) |> mutate(gap = as.numeric(date - lag(date))) |> ungroup() |> filter(!is.na(gap))
  if (nrow(gaps) > 0) {
    add("sample", "nophase_gap_min", min(gaps$gap)); add("sample", "nophase_gap_max", max(gaps$gap))
    for (i in seq_len(nrow(gaps))) add("nophase_gaps", paste0(gaps$id[i], ".", gaps$date[i]), gaps$gap[i])
  }
})

# ---- 2) Days per data source ----
block("day_counts", {
  valid_ids <- pacts_phase_diagnosis |> filter(has_valid_phase, Id %in% char_table$Id) |> pull(Id)
  nophase_ids <- setdiff(char_table$Id, valid_ids)
  onsets <- pacts_input |> filter(menses == 1L) |> group_by(id) |>
    summarise(erster = min(date), letzter = max(date), .groups = "drop")
  dd <- daily_data |> filter(Id %in% char_table$Id, !is.na(light_sensitivity_idx)) |>
    left_join(onsets, by = c("Id" = "id"))
  has_ph <- !is.na(dd$cycle_phase_2)
  is_valid <- dd$Id %in% valid_ids
  add("day_counts", "n_survey_days_total", nrow(dd), "Tage mit ausgefuelltem Light-Sensitivity-Modul")
  add("day_counts", "n_in_valid_window", sum(has_ph))
  add("day_counts", "pct_in_valid_window", round(100 * mean(has_ph), 1))
  add("day_counts", "n_days_no_phase_participants", sum(dd$Id %in% nophase_ids))
  add("day_counts", "n_on_or_after_last_onset", sum(is_valid & !has_ph & !is.na(dd$letzter) & dd$Date >= dd$letzter))
  add("day_counts", "n_before_first_onset",     sum(is_valid & !has_ph & !is.na(dd$erster)  & dd$Date <  dd$erster))
  add("day_counts", "n_between_onsets_no_phase", sum(is_valid & !has_ph & !is.na(dd$erster) & !is.na(dd$letzter) &
                                                     dd$Date >= dd$erster & dd$Date < dd$letzter),
      "erwartet 0")
  by_id <- dd |> group_by(Id) |> summarise(n = n(), n_win = sum(!is.na(cycle_phase_2)), .groups = "drop")
  for (i in seq_len(nrow(by_id))) {
    add("days_by_id", paste0(by_id$Id[i], ".n_survey_days"), by_id$n[i])
    add("days_by_id", paste0(by_id$Id[i], ".n_in_window"),   by_id$n_win[i])
  }
})

# ---- 3) LMMs (Predictions 1-4) ----
block("lmm_main", {
  LMM_LIST <- list(TAT_250 = "lmm_TAT", TAT_10 = "lmm_TAT10", TAT_1000 = "lmm_TAT1000",
                   MEDI_mean = "lmm_MEDI", light_sensitivity = "lmm_lsens", cycle_symptom = "lmm_csymp",
                   IS = "lmm_IS", IV = "lmm_IV", RA = "lmm_RA", TST = "lmm_TST",
                   sleep_efficiency = "lmm_sleep_eff")
  for (nm in names(LMM_LIST)) lmm_row(getobj(LMM_LIST[[nm]]), nm, "lmm")

  lmm_row(getobj("lmm_IS_full21"), "IS_full21", "lmm_npcra_full21")
  lmm_row(getobj("lmm_IV_full21"), "IV_full21", "lmm_npcra_full21")
  if (exists("npcra_ok_ids")) {
    add("npcra", "n_ok_ids", length(npcra_ok_ids))
    add("npcra", "ok_ids", paste(sort(npcra_ok_ids), collapse = " "))
  }
  if (exists("npcra_phase")) {
    add("npcra", "n_low_n_rows", sum(npcra_phase$low_n_npcra, na.rm = TRUE))
    for (i in seq_len(nrow(npcra_phase))) {
      r <- npcra_phase[i, ]
      add("npcra_days_by_id_phase", paste0(r$Id, ".", r$cycle_phase_2, ".n_days"), r$n_days)
    }
  }
})

# ---- 4) GAMMs of the main outcomes ----
block("gamm_main", {
  G <- list(TAT_250 = "gamm_TAT", TAT_10 = "gamm_TAT10", TAT_1000 = "gamm_TAT1000",
            light_sensitivity = "gamm_lsens", cycle_symptom = "gamm_csymp",
            RA_act = "gamm_RA_act")
  for (nm in names(G)) gamm_row(getobj(G[[nm]]), nm, "gamm")
})

# ---- 5) High-completeness subsample ----
block("sens_compliance", {
  S <- list(TAT_250 = "sens_compliance_TAT", MEDI_mean = "sens_compliance_MEDI",
            light_sensitivity = "sens_compliance_lsens", cycle_symptom = "sens_compliance_csymp")
  for (nm in names(S)) lmm_row(getobj(S[[nm]]), nm, "lmm_compliance_subsample")
  S2 <- list(TAT_250 = "sens_strict_TAT", MEDI_mean = "sens_strict_MEDI",
             light_sensitivity = "sens_strict_lsens", cycle_symptom = "sens_strict_csymp")
  for (nm in names(S2)) lmm_row(getobj(S2[[nm]]), nm, "lmm_compliance_strict")
  if (exists("hc_ids_strict")) {
    add("compliance", "n_strict_participants", length(hc_ids_strict))
    add("compliance", "strict_ids", paste(hc_ids_strict, collapse = " "))
  }
  if (exists("hc_ids")) {
    add("compliance", "n_high_compliance_participants", length(hc_ids))
    add("compliance", "high_compliance_ids", paste(hc_ids, collapse = " "))
  }
})

# ---- 5b) Robustness checks (09b) ----
block("robustness", {
  if (exists("robustness_table")) {
    for (i in seq_len(nrow(robustness_table))) {
      r <- robustness_table[i, ]
      key <- paste(r$finding, r$test, sep = " | ")
      for (f in intersect(c("n_obs", "n_id", "b", "edf", "t", "df", "p", "note"), names(r)))
        add("robustness", paste0(key, ".", f), r[[f]])
    }
  }
})

# ---- 6) Missing data ----
block("mnar", {
  for (nm in c("mnar_light", "mnar_sleep")) {
    m <- getobj(nm)
    if (is.null(m)) { add("mnar", paste0(nm, ".missing"), 1); next }
    s <- summary(m)$coefficients
    for (term in c("cycle_phase_2luteal", "cycle_symptom_idx")) {
      key <- paste0(nm, ".", term)
      add("mnar", paste0(key, ".b"),  s[term, "Estimate"])
      add("mnar", paste0(key, ".se"), s[term, "Std. Error"])
      add("mnar", paste0(key, ".z"),  s[term, "z value"])
      add("mnar", paste0(key, ".p"),  s[term, "Pr(>|z|)"])
      add("mnar", paste0(key, ".odds_ratio"), exp(s[term, "Estimate"]))
    }
  }
})

# ---- 7) Completeness by participant and phase ----
block("compliance_table", {
  for (i in seq_len(nrow(missing_rep))) {
    k <- paste0(missing_rep$Id[i], ".", missing_rep$cycle_phase_2[i])
    add("compliance_by_id_phase", paste0(k, ".n_days"),    missing_rep$n[i])
    add("compliance_by_id_phase", paste0(k, ".pct_light"), missing_rep$pct_light[i])
    add("compliance_by_id_phase", paste0(k, ".pct_sleep"), missing_rep$pct_sleep[i])
    add("compliance_by_id_phase", paste0(k, ".pct_survey"), missing_rep$pct_survey[i])
  }
})

# ---- 8) Secondary light metrics (Table 2) ----
block("secondary", {
  add("secondary", "bonferroni_alpha", 0.05 / length(SECONDARY_LIGHT_VARS))
  add("secondary", "n_metrics", length(SECONDARY_LIGHT_VARS))
  for (v in SECONDARY_LIGHT_VARS) {
    gamm_row(gamm_secondary_light[[v]], v, "secondary_gamm")
    lmm_row(lmm_secondary_light[[v]],  v, "secondary_lmm")
  }
})

# ---- 9) Item level (Tables 4 and 5) ----
block("items", {
  for (v in names(gamm_ls_items))   gamm_row(gamm_ls_items[[v]],   v, "gamm_ls_item")
  for (v in names(gamm_drsp_items)) gamm_row(gamm_drsp_items[[v]], v, "gamm_drsp_item")

  # Holm and Benjamini-Hochberg within each item family
  adj_rows <- function(gamm_list, section) {
    raw <- vapply(gamm_list, function(m)
      if (is.null(m)) NA_real_ else summary(m)$s.table["s(pacts_ov)", "p-value"], numeric(1))
    ph <- p.adjust(raw, method = "holm"); pb <- p.adjust(raw, method = "BH")
    for (v in names(raw)) {
      add(section, paste0(v, ".p_holm"), ph[[v]])
      add(section, paste0(v, ".p_bh"),   pb[[v]])
    }
    add(section, "family_size", length(raw))
  }
  adj_rows(gamm_ls_items,   "gamm_ls_item")
  adj_rows(gamm_drsp_items, "gamm_drsp_item")
})

# ---- 10) Correlation matrix ----
block("corr", {
  cm <- corr_mat
  for (i in seq_len(nrow(cm))) for (j in seq_len(ncol(cm)))
    if (j > i) add("corr", paste(rownames(cm)[i], colnames(cm)[j], sep = "__"), cm[i, j])
  if (exists("corr_data")) add("corr", "n_complete_rows", nrow(corr_data))
})

# ---- 11) PACTS bins ----
block("bins", {
  bp <- daily_data |> filter(!is.na(cycle_bin_ov)) |> group_by(cycle_bin_ov) |>
    summarise(n_participants = n_distinct(Id), .groups = "drop")
  add("bins", "n_bins", max(bp$cycle_bin_ov))
  add("bins", "bin1.n_participants",  bp$n_participants[bp$cycle_bin_ov == 1])
  add("bins", "last_bin.n_participants", bp$n_participants[bp$cycle_bin_ov == max(bp$cycle_bin_ov)])
  add("bins", "n_participants_with_bins", n_distinct(daily_data$Id[!is.na(daily_data$cycle_bin_ov)]))
  if (exists("drsp_bin")) {
    add("bins", "drsp_bins_low_n", sum(drsp_bin$low_n, na.rm = TRUE))
    add("bins", "drsp_bins_zero_se_n_gt1", sum(drsp_bin$n > 1 & drsp_bin$se == 0, na.rm = TRUE))
  }
})

# ---- 12) Sleep QC: nights with (almost) 0 h sleep ----
block("sleep_qc", {
  add("sleep_qc", "n_nights", nrow(actigraph_sleep))
  add("sleep_qc", "n_TST_lt_60min",  sum(actigraph_sleep$TST_min < 60,  na.rm = TRUE))
  add("sleep_qc", "n_TST_lt_120min", sum(actigraph_sleep$TST_min < 120, na.rm = TRUE))
  add("sleep_qc", "n_sleep_eff_lt_20", sum(actigraph_sleep$sleep_efficiency < 20, na.rm = TRUE))
  short <- actigraph_sleep |> filter(TST_min < 60) |> count(Id)
  for (i in seq_len(nrow(short))) add("sleep_qc_by_id", paste0(short$Id[i], ".n_TST_lt_60min"), short$n[i])
  ggir_csv <- "C:/GGIR_local_output/LEATMC/output_LEATMC/results/QC/part4_nightsummary_sleep_full.csv"
  if (file.exists(ggir_csv)) {
    ns <- read.csv(ggir_csv)
    sh <- ns[!is.na(ns$SleepDurationInSpt) & ns$SleepDurationInSpt < 1, ]
    add("sleep_qc_ggir", "n_short_nights_lt_1h", nrow(sh))
    add("sleep_qc_ggir", "mean_fraction_night_invalid_short", mean(sh$fraction_night_invalid, na.rm = TRUE))
    add("sleep_qc_ggir", "mean_fraction_night_invalid_all",   mean(ns$fraction_night_invalid, na.rm = TRUE))
    for (cc in unique(sh$cleaningcode))
      add("sleep_qc_ggir", paste0("short_nights.cleaningcode.", cc), sum(sh$cleaningcode == cc, na.rm = TRUE))
  }
})

# ---- Write ----
thesis_numbers <- dplyr::bind_rows(.rows)
if (!exists("safe_write_csv", mode = "function")) {
  safe_write_csv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    cands <- unique(c(path, sub("(\\.[A-Za-z0-9]+)$", paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S"), "\\1"), path),
                      file.path(tempdir(), basename(path))))
    for (cp in cands) {
      r <- tryCatch({ utils::write.csv(x, cp, row.names = FALSE, fileEncoding = "UTF-8"); NULL },
                    warning = function(w) conditionMessage(w), error = function(e) conditionMessage(e))
      if (is.null(r) && file.exists(cp)) { if (cp != path) message("[safe_write_csv] Ersatzdatei: ", cp); return(invisible(cp)) }
      message("[safe_write_csv] ", cp, " nicht beschreibbar: ", r)
    }
    invisible(NA_character_)
  }
}
out_file <- file.path(OUTPUT_DIR, "thesis_numbers.csv")
out_written <- safe_write_csv(thesis_numbers, out_file)
cat("[10b_thesis_numbers]", nrow(thesis_numbers), "Werte exportiert ->", out_written, "\n")
errs <- thesis_numbers[thesis_numbers$section == "export_error", ]
if (nrow(errs) > 0) {
  cat("[10b_thesis_numbers] ACHTUNG - folgende Bloecke sind fehlgeschlagen:\n")
  print(errs[, c("key", "value")], row.names = FALSE)
}
