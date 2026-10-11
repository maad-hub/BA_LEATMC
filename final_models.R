# R/final_models.R - Models of 12_gamm_ar1.R for one phase assignment
# LMM phase contrasts (as 09_inferential.R), GAMMs with and without AR(1), main tests P1-P4
# (P3: GAMM with AR(1); Holm over the five P4 outcomes), CSI rows of Table C1, Table 3,
# covariate-adjusted GAMMs.

E2_K <- 20   # basis dimension of the cyclic spline

e2_bind <- function(lst) { lst <- Filter(Negate(is.null), lst); if (!length(lst)) return(NULL); do.call(rbind, lst) }

# GAMM without and with AR(1) for one outcome
e2_gam_pair <- function(dd, var, family, extra = "") {
  g0 <- ac_try(ac_gam_noar(dd, var, E2_K, extra), paste("GAMM noAR", var))
  g1 <- ac_try(ac_gam_ar1(dd, var, E2_K, extra), paste("GAMM AR1", var))
  rows <- rbind(ac_gam_row(g0, var, "noAR"), ac_gam_row(g1, var, "AR1"))
  rows$family <- family; rows$extra <- extra
  rows
}

# phase assignment summary and association of phase with study day
e2_confound <- function(dd) {
  d <- dd[!is.na(dd$cycle_phase_2), ]
  ph <- tapply(as.character(d$cycle_phase_2), d$Id, function(x) paste(sort(unique(x)), collapse = "+"))
  m <- ac_lmm_day(d, "study_day", label = "study_day ~ phase")
  q_all <- sum(!is.na(dd$light_sensitivity_idx)); q_ph <- sum(!is.na(dd$light_sensitivity_idx) & !is.na(dd$cycle_phase_2))
  closed <- unique(d$Id[d$count_method %in% "closed_pacts"])
  counting_only <- setdiff(unique(d$Id), closed)
  data.frame(n_participants = length(unique(dd$Id)), n_with_phase = length(ph),
             n_only_luteal = sum(ph == "luteal"), n_only_follicular = sum(ph == "follicular"),
             n_both = sum(ph == "follicular+luteal"), n_days_with_phase = nrow(d),
             n_days_follicular = sum(d$cycle_phase_2 == "follicular"), n_days_luteal = sum(d$cycle_phase_2 == "luteal"),
             n_questionnaire_days = q_all, n_questionnaire_days_with_phase = q_ph,
             mean_study_day_follicular = mean(d$study_day[d$cycle_phase_2 == "follicular"], na.rm = TRUE),
             mean_study_day_luteal = mean(d$study_day[d$cycle_phase_2 == "luteal"], na.rm = TRUE),
             studyday_b = m$b, studyday_lo95 = m$lo95, studyday_hi95 = m$hi95, studyday_t = m$t,
             studyday_df = m$df, studyday_p = m$p,
             n_days_forward_imputed = sum(d$count_method %in% "forward_imputed"),
             n_days_backward_imputed = sum(d$count_method %in% "backward_imputed"),
             n_ids_closed_cycles = length(closed),
             n_ids_counting_only = length(counting_only),
             n_ids_backward = length(unique(d$Id[d$count_method %in% "backward_imputed"])),
             n_ids_forward = length(unique(d$Id[d$count_method %in% "forward_imputed"])),
             n_ids_without_phase = length(setdiff(unique(dd$Id), unique(d$Id))),
             n_counting_only_backward = length(intersect(counting_only, unique(d$Id[d$count_method %in% "backward_imputed"]))),
             n_counting_only_forward = length(intersect(counting_only, unique(d$Id[d$count_method %in% "forward_imputed"]))))
}

# ---- All models of one assignment ----
E2_GAMM_SET <- list(
  main = c("TAT_250_hr", "light_sensitivity_idx", "cycle_symptom_idx", "TAT_10_hr", "TAT_1000_hr",
           "RA_act", "TST_min", "sleep_efficiency"),
  secondary = AC_SECONDARY, ls_items = AC_LS_ITEMS, drsp_items = AC_DRSP_ITEMS, impairment = AC_IMPAIR)

e2_main_models <- function(dd, variant, activity = NULL, verbose = TRUE) {
  say <- function(...) if (verbose) cat("  [", variant, "] ", ..., "\n", sep = "")
  out <- list()
  say("LMMs on participant x phase means")
  lv <- intersect(c("TAT_250_hr", "TAT_10_hr", "TAT_1000_hr", "MEDI_mean", "light_sensitivity_idx", "cycle_symptom_idx",
                    "RA_act", "TST_min", "sleep_efficiency", "centroid_LE_hr"), names(dd))
  lmm <- e2_bind(lapply(lv, function(v) ac_lmm_pp(dd, v, label = v)))
  if (!is.null(activity) && nrow(activity)) {
    np <- ac_try(ac_npcra(activity, dd), "IS/IV")
    if (!is.null(np)) {
      ok <- attr(np, "ok_ids")
      lmm <- rbind(lmm, ac_lmm_npcra(np, "IS", ok, "IS"), ac_lmm_npcra(np, "IV", ok, "IV"))
      out$npcra_n <- data.frame(n_ok_ids = length(ok))
    }
  }
  names(lmm)[names(lmm) == "label"] <- "outcome"
  out$lmm <- cbind(variant = variant, lmm)

  say("GAMMs without and with AR(1)")
  gam <- e2_bind(unlist(lapply(names(E2_GAMM_SET), function(f)
    lapply(intersect(E2_GAMM_SET[[f]], names(dd)), function(v) e2_gam_pair(dd, v, f))), recursive = FALSE))
  gam$p_holm_family <- NA_real_
  for (f in c("secondary", "ls_items", "drsp_items", "impairment")) for (m in c("noAR", "AR1")) {
    i <- gam$family == f & gam$model == m & !is.na(gam$p)
    gam$p_holm_family[i] <- stats::p.adjust(gam$p[i], "holm")
  }
  out$gamm <- cbind(variant = variant, gam)

  say("Main tests P1-P4")
  pl <- function(o) { r <- lmm[lmm$outcome == o, ][1, ]
    data.frame(estimate = r$b, lo95 = r$lo95, hi95 = r$hi95, stat = r$t, df = r$df, p = r$p, n_obs = r$n_obs, n_id = r$n_id,
               singular = r$singular, p_other_model = NA_real_) }
  g1 <- gam[gam$outcome == "cycle_symptom_idx" & gam$model == "AR1", ]; g0 <- gam[gam$outcome == "cycle_symptom_idx" & gam$model == "noAR", ]
  main <- rbind(
    cbind(test = "P1", outcome = "TAT_250_hr", model = "LMM", pl("TAT_250_hr")),
    cbind(test = "P2", outcome = "light_sensitivity_idx", model = "LMM", pl("light_sensitivity_idx")),
    data.frame(test = "P3", outcome = "cycle_symptom_idx", model = "GAMM with AR(1)", estimate = g1$edf[1], lo95 = NA, hi95 = NA,
               stat = g1$F[1], df = g1$ref_df[1], p = g1$p[1], n_obs = g1$n[1], n_id = g1$n_id[1], singular = NA,
               p_other_model = g0$p[1]),
    cbind(test = "P4", outcome = "IS", model = "LMM (>= 5 days in both phases)", pl("IS")),
    cbind(test = "P4", outcome = "IV", model = "LMM (>= 5 days in both phases)", pl("IV")),
    cbind(test = "P4", outcome = "RA_act", model = "LMM", pl("RA_act")),
    cbind(test = "P4", outcome = "TST_min", model = "LMM", pl("TST_min")),
    cbind(test = "P4", outcome = "sleep_efficiency", model = "LMM", pl("sleep_efficiency")))
  main$p_holm5_P4 <- NA_real_; i4 <- main$test == "P4"; main$p_holm5_P4[i4] <- stats::p.adjust(main$p[i4], "holm")
  main$phi_AR1 <- ifelse(main$test == "P3", g1$phi[1], NA)
  out$confirmatory <- cbind(variant = variant, main)

  say("Cycle Symptom Index: sensitivity analyses (Table C1)")
  out$robustness_csi <- e2_robustness_csi(dd, variant)

  say("Curves for Figure 4")
  curves <- list()
  for (v in intersect(c("TAT_250_hr", "light_sensitivity_idx", "cycle_symptom_idx"), names(dd))) {
    g <- ac_try(ac_gam_ar1(dd, v, E2_K), v)
    if (!is.null(g)) { cv <- ac_try(ac_curve(g), v); if (!is.null(cv)) curves[[length(curves) + 1]] <- cbind(outcome = v, model = "AR1", cv) }
  }
  out$curves <- cbind(variant = variant, e2_bind(curves))
  out$confound <- cbind(variant = variant, e2_confound(dd))
  out
}

# ---- Sensitivity analyses of the Cycle Symptom Index with the main model (GAMM with AR(1)) ----
# participants with corrected menses data (false entries removed, onsets added, onsets merged)
e2_corrected_ids <- function(env = globalenv()) {
  g <- function(n) if (exists(n, envir = env)) get(n, envir = env) else NULL
  unique(as.character(c(g("MENSES_FALSE_POSITIVES")$Id, g("manual_menses")$Id, g("menses_merged")$id)))
}
# participants with at least thr % of days with data in each phase (strict: and data in both phases)
e2_avail_ids <- function(dd, var, thr = 80, strict = FALSE) {
  d <- dd[!is.na(dd$cycle_phase_2), ]
  a <- stats::aggregate(list(pct = !is.na(d[[var]])), by = list(Id = d$Id, ph = as.character(d$cycle_phase_2)),
                        FUN = function(z) mean(z) * 100)
  m <- stats::aggregate(pct ~ Id, data = a, FUN = min); k <- stats::aggregate(ph ~ Id, data = a, FUN = length)
  ok <- m$Id[m$pct >= thr]
  if (strict) ok <- intersect(ok, k$Id[k$ph == 2])
  ok
}
# participants outside the protocol's age and BMI ranges (as revision.R R35)
e2_age_bmi_ids <- function(dd, age_max = 35, bmi_max = 26) {
  if (!all(c("age", "bmi") %in% names(dd))) {
    if (exists("char_table") && all(c("Id", "age", "bmi") %in% names(get("char_table")))) {
      ct <- get("char_table")
      return(intersect(as.character(ct$Id[which(ct$age > age_max | ct$bmi > bmi_max)]), unique(as.character(dd$Id))))
    }
    return(character())
  }
  a <- tapply(dd$age, dd$Id, function(z) suppressWarnings(max(z, na.rm = TRUE)))
  b <- tapply(dd$bmi, dd$Id, function(z) suppressWarnings(max(z, na.rm = TRUE)))
  names(which((is.finite(a) & a > age_max) | (is.finite(b[names(a)]) & b[names(a)] > bmi_max)))
}
e2_robustness_csi <- function(dd, variant, k = E2_K, thr = if (exists("COMPL_THR")) COMPL_THR else 80) {
  v <- "cycle_symptom_idx"; cid <- e2_corrected_ids(); abi <- e2_age_bmi_ids(dd)
  one <- function(test, d, extra = "") {
    g1 <- ac_try(ac_gam_ar1(d, v, k, extra), test); g0 <- ac_try(ac_gam_noar(d, v, k, extra), test)
    a <- ac_gam_row(g1, v, "AR1"); b <- ac_gam_row(g0, v, "noAR")
    data.frame(variant = variant, test = test, n_days = a$n, n_id = a$n_id, edf_AR1 = a$edf, F_AR1 = a$F, p_AR1 = a$p,
               phi = a$phi, edf_noAR = b$edf, p_noAR = b$p, note = "")
  }
  rows <- list(
    one("Main analysis", dd),
    one("Day length and weekend added", dd, " + s(day_length, k = 5) + weekend"),
    one(paste0("At least ", thr, "% questionnaire data in each phase"), dd[dd$Id %in% e2_avail_ids(dd, v, thr), ]),
    { r <- one(paste0("Without participants outside the age/BMI range (", if (length(abi)) paste(abi, collapse = ", ") else "none found", ")"),
               dd[!dd$Id %in% abi, ])
      if (!length(abi)) r$note <- "age/BMI not found - row equals the main analysis"
      r },
    one("Without participants with corrected menses data", dd[!dd$Id %in% cid, ]),
    one("High-completeness subsample, light-logger data", dd[dd$Id %in% e2_avail_ids(dd, "TAT_250_hr", thr), ]),
    one("... restricted to participants with data in both phases", dd[dd$Id %in% e2_avail_ids(dd, "TAT_250_hr", thr, TRUE), ]))
  ids <- sort(unique(dd$Id[!is.na(dd$pacts_ov) & !is.na(dd[[v]])]))
  loo <- e2_bind(lapply(ids, function(i) { g <- ac_try(ac_gam_ar1(dd[dd$Id != i, ], v, k), paste("LOO", i))
    a <- ac_gam_row(g, v, "AR1"); data.frame(left_out = i, edf = a$edf, p = a$p) }))
  rows[[length(rows) + 1]] <- data.frame(variant = variant, test = paste0("Leave-one-participant-out refits (", nrow(loo), " refits)"),
    n_days = NA, n_id = nrow(loo), edf_AR1 = stats::median(loo$edf, na.rm = TRUE), F_AR1 = NA, p_AR1 = max(loo$p, na.rm = TRUE),
    phi = NA, edf_noAR = NA, p_noAR = NA,
    note = paste0("EDF = median, p = maximum; p < .05 in ", sum(loo$p < .05, na.rm = TRUE), " of ", nrow(loo), " refits"))
  e2_bind(rows)
}

# ---- Table 3 with both assignments ----
E2_TABLE3_ROWS <- data.frame(
  outcome = c("TAT_250_hr", "TAT_10_hr", "TAT_1000_hr", "light_sensitivity_idx", "cycle_symptom_idx",
              "IS", "IV", "RA_act", "TST_min", "sleep_efficiency"),
  role = c("P1", "Expl.", "Expl.", "P2", "P3", "P4", "P4", "P4", "P4", "P4"), stringsAsFactors = FALSE)
e2_table3 <- function(full, variants = names(full)) {
  e2_bind(lapply(variants, function(v) {
    f <- full[[v]]; if (is.null(f)) return(NULL)
    e2_bind(lapply(seq_len(nrow(E2_TABLE3_ROWS)), function(i) {
      o <- E2_TABLE3_ROWS$outcome[i]
      l <- f$lmm[f$lmm$outcome == o, ][1, ]
      g1 <- f$gamm[f$gamm$outcome == o & f$gamm$model == "AR1" & f$gamm$extra == "", ]
      g0 <- f$gamm[f$gamm$outcome == o & f$gamm$model == "noAR" & f$gamm$extra == "", ]
      cf <- f$confirmatory[f$confirmatory$outcome == o, ]
      data.frame(variant = v, role = E2_TABLE3_ROWS$role[i], outcome = o, n_obs = l$n_obs, n_id = l$n_id,
                 mean_follicular = l$mean_fol, mean_luteal = l$mean_lut,
                 gamm_AR1_edf = if (nrow(g1)) g1$edf[1] else NA, gamm_AR1_p = if (nrow(g1)) g1$p[1] else NA,
                 gamm_noAR_edf = if (nrow(g0)) g0$edf[1] else NA, gamm_noAR_p = if (nrow(g0)) g0$p[1] else NA,
                 lmm_b = l$b, lmm_lo95 = l$lo95, lmm_hi95 = l$hi95, lmm_t = l$t, lmm_df = l$df, lmm_p = l$p,
                 singular = l$singular,
                 p_holm5_P4 = if (nrow(cf) && "p_holm5_P4" %in% names(cf)) cf$p_holm5_P4[1] else NA)
    }))
  }))
}

# ---- GAMMs of the nominally significant findings adjusted for study day or day length (each with weekend),
#      and the centroid without participants outside the age/BMI ranges (Table C1) ----
E2_ADJUSTED <- list(
  outcomes = c("cycle_symptom_idx", "TAT_10_hr", "centroid_LE_hr", "RA_act"),
  extras = c(unadjusted = "", study_day_weekend = " + s(study_day, k = 5) + weekend",
             day_length_weekend = " + s(day_length, k = 5) + weekend"),
  without_age_bmi = c("centroid_LE_hr"))
e2_adjusted_gamms <- function(dd_list) {
  e2_bind(lapply(names(dd_list), function(v) {
    dd <- dd_list[[v]]
    adj <- e2_bind(lapply(E2_ADJUSTED$outcomes, function(o)
      e2_bind(lapply(names(E2_ADJUSTED$extras), function(a) {
        cat("  [adjusted] ", v, " ", o, " ", a, "\n", sep = "")
        r <- e2_gam_pair(dd, o, "adjusted", E2_ADJUSTED$extras[[a]])
        cbind(variant = v, adjustment = a, r) }))))
    abi <- e2_age_bmi_ids(dd)
    sub <- e2_bind(lapply(E2_ADJUSTED$without_age_bmi, function(o) {
      cat("  [adjusted] ", v, " ", o, " without_age_bmi\n", sep = "")
      r <- e2_gam_pair(dd[!dd$Id %in% abi, ], o, "subset")
      r$note <- if (length(abi)) paste0(length(abi), " excluded") else "age/BMI not found - equals unadjusted"
      cbind(variant = v, adjustment = "without_age_bmi", r) }))
    e2_bind(list(adj, sub)) }))
}
