# revision.R - Additional analyses: corrections (fixed boundary knots, t-based CIs), Table 2,
# confounding of cycle phase with study time, robustness checks, rest-activity and sleep,
# repeated-measures correlations. Results in OUTPUT_DIR/revision_output/; IDs (M*, R*) label the results.

# ==== Settings, packages, helpers ====

REV_DIR  <- file.path(OUTPUT_DIR, "revision_output")
dir.create(REV_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Robust writing (as in 00b_utils_io.R) ----
if (!(exists("safe_write", mode = "function") && exists("safe_write_csv", mode = "function") &&
      exists("safe_write_lines", mode = "function"))) {
  safe_write <- function(path, writer) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    cands <- unique(c(path, sub("(\\.[A-Za-z0-9]+)$", paste0("_", stamp, "\\1"), path),
                      file.path(tempdir(), basename(path))))
    why <- character()
    for (i in seq_along(cands)) {
      r <- tryCatch({ writer(cands[i]); NULL },
                    warning = function(w) conditionMessage(w), error = function(e) conditionMessage(e))
      if (is.null(r) && file.exists(cands[i])) {
        if (i > 1) message("[safe_write] '", path, "' nicht beschreibbar (", paste(unique(why), collapse = " | "),
                           "). Stattdessen geschrieben: ", cands[i])
        return(invisible(cands[i]))
      }
      why <- c(why, if (is.null(r)) "Datei nach dem Schreiben nicht vorhanden" else r)
    }
    message("[safe_write] FEHLGESCHLAGEN fuer '", path, "': ", paste(unique(why), collapse = " | "))
    invisible(NA_character_)
  }
  safe_write_csv <- function(x, path)
    safe_write(path, function(p) utils::write.csv(x, p, row.names = FALSE, fileEncoding = "UTF-8"))
  safe_write_lines <- function(lines, path)
    safe_write(path, function(p) {
      con <- file(p, open = "w", encoding = "UTF-8"); on.exit(close(con), add = TRUE); writeLines(lines, con)
    })
}

RV_SKIP        <- character(0)  # names of blocks to skip, e.g. c("R15", "R20")
RV_SEED        <- 20260930
RV_LAT         <- 48.137  # Munich
RV_LON         <- 11.575
RV_AGE_MAX     <- 35  # R35
RV_BMI_MAX     <- 26  # R35
RV_N_REP       <- 1000  # R20: repetitions
RV_N_DAYS      <- 5  # R20: days per participant and phase
RV_MIN_MIN_PER_HOUR  <- 30  # R19/R20: minimum valid minutes per hour
RV_MIN_HOURS_PER_DAY <- 16  # R20: minimum valid hours per valid day

rv_on <- function(tag) !(tag %in% RV_SKIP)

# ---- Package check ----
rv_needed  <- c("suncalc", "rmcorr")
rv_missing <- rv_needed[!vapply(rv_needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(rv_missing) > 0) {
  message("[revision] Fehlende Pakete (betroffene Abschnitte werden 'NICHT VERFUEGBAR'): ",
          paste(rv_missing, collapse = ", "),
          "\n  install.packages(c(", paste0('"', rv_missing, '"', collapse = ", "), "))")
}

# ---- Result collection ----
results <- data.frame(id = character(), text = character(), note = character(),
                      stringsAsFactors = FALSE)
add_result <- function(id, text, note = "") {
  results <<- rbind(results, data.frame(
    id = as.character(id),
    text = paste(as.character(text), collapse = " "),
    note = paste(as.character(note), collapse = " "),
    stringsAsFactors = FALSE))
  invisible(NULL)
}

# ---- Formatting ----
rv_minus <- function(s) gsub("-", "\u2212", s, fixed = TRUE)
fmt_num <- function(x, d = 2) {
  x <- ifelse(is.na(x), NA_real_, ifelse(abs(x) < 0.5 * 10^(-d), 0, x))
  rv_minus(formatC(x, format = "f", digits = d))
}
fmt_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < .001, "< .001", sub("^0", "", formatC(p, format = "f", digits = 3))))
}
fmt_peq <- function(p) {
  ifelse(is.na(p), "p = NA", ifelse(p < .001, "p < .001", paste0("p = ", fmt_p(p))))
}

# LMM coefficient with t-based 95% CI: b +/- qt(.975, df) * SE
lmm_coef <- function(m, term = "cycle_phase_2luteal") {
  s <- summary(m)$coefficients
  if (!term %in% rownames(s))
    stop("Term '", term, "' nicht im Modell. Vorhanden: ", paste(rownames(s), collapse = ", "))
  r  <- s[term, ]
  b  <- unname(r["Estimate"]); se <- unname(r["Std. Error"]); df <- unname(r["df"])
  tq <- stats::qt(.975, df)
  list(b = b, se = se, df = df, t = unname(r["t value"]), p = unname(r["Pr(>|t|)"]),
       lo = b - tq * se, hi = b + tq * se)
}
fmt_lmm <- function(model, term = "cycle_phase_2luteal", d = 2) {
  k <- lmm_coef(model, term)
  paste0("b = ", fmt_num(k$b, d), ", 95% CI [", fmt_num(k$lo, d), ", ", fmt_num(k$hi, d),
         "], t(", formatC(k$df, format = "f", digits = 1), ") = ", fmt_num(k$t, 2), ", ",
         fmt_peq(k$p))
}
lmm_note <- function(m) {
  paste0("N = ", stats::nobs(m), " Beobachtungen / ", lme4::ngrps(m)[["Id"]], " Personen",
         if (lme4::isSingular(m)) "; SINGULAERER FIT" else "")
}

# GAMM: EDF, F and p of the cycle smooth
gam_row <- function(m, smooth = "s(pacts_ov)") {
  st <- summary(m)$s.table
  if (!smooth %in% rownames(st)) stop("Smooth '", smooth, "' nicht im Modell.")
  r <- st[smooth, ]
  list(edf = unname(r["edf"]),
       F   = unname(r[grep("^F$|^Chi", names(r))[1]]),
       p   = unname(r[grep("p-value", names(r))[1]]))
}
fmt_gam <- function(model, smooth = "s(pacts_ov)") {
  k <- gam_row(model, smooth)
  paste0("EDF = ", fmt_num(k$edf, 2), ", F = ", fmt_num(k$F, 2), ", ", fmt_peq(k$p))
}
gam_note <- function(m) {
  paste0("N = ", nrow(m$model), " Tage / ", length(unique(m$model$Id)), " Personen",
         if (!is.null(attr(m, "k_used"))) paste0("; k = ", attr(m, "k_used")) else "")
}

# runs a block in its own environment; an error is recorded in the results and the script continues
rv_run <- function(ids, expr) {
  ex  <- substitute(expr)
  env <- new.env(parent = globalenv())
  ok <- tryCatch({ eval(ex, env); TRUE }, error = function(e) {
    msg <- conditionMessage(e)
    message("  [revision] FEHLER in ", paste(ids, collapse = "/"), ": ", msg)
    for (i in setdiff(ids, results$id)) add_result(i, "NICHT VERFUEGBAR", msg)
    FALSE
  })
  invisible(ok)
}
rv_skip_note <- function(ids, tag) {
  for (i in ids) add_result(i, "NICHT VERFUEGBAR", paste0("uebersprungen (RV_SKIP enthaelt ", tag, ")"))
}
rv_get <- function(n) get0(n, envir = .GlobalEnv, inherits = FALSE)

cat("\n[revision] Start - Ausgabe nach", REV_DIR, "\n")

# ==== Data preparation (study day, weekend, day length) ====
rv_dd <- daily_data
rv_dd$Id   <- as.character(rv_dd$Id)
rv_dd$Date <- as.Date(rv_dd$Date)

rv_first <- if (exists("date_ranges")) {
  as.data.frame(date_ranges)[, c("Id", "d_first")]
} else {
  aggregate(Date ~ Id, data = rv_dd, FUN = min) |> stats::setNames(c("Id", "d_first"))
}
rv_first$Id <- as.character(rv_first$Id)
rv_first$d_first <- as.Date(rv_first$d_first)
rv_dd <- merge(rv_dd, rv_first, by = "Id", all.x = TRUE, sort = FALSE)
# study day: days since the first day
rv_dd$study_day <- as.integer(rv_dd$Date - rv_dd$d_first)
rv_dd$d_first   <- NULL
rv_dd$weekend   <- as.integer(lubridate::wday(rv_dd$Date, week_start = 1) >= 6)

# day length (h) in Munich
rv_daylength_h <- function(date, lat) {
  doy  <- as.integer(format(date, "%j"))
  decl <- 23.44 * pi / 180 * sin(2 * pi * (284 + doy) / 365)
  phi  <- lat * pi / 180
  24 * acos(pmin(1, pmax(-1, -tan(phi) * tan(decl)))) / pi
}
rv_ud <- sort(unique(rv_dd$Date[!is.na(rv_dd$Date)]))
rv_dl <- tryCatch({
  st <- suncalc::getSunlightTimes(date = rv_ud, lat = RV_LAT, lon = RV_LON,
                                  keep = c("sunrise", "sunset"), tz = TIMEZONE)
  data.frame(Date = rv_ud,
             day_length = as.numeric(difftime(st$sunset, st$sunrise, units = "hours")))
}, error = function(e) {
  message("  [revision] suncalc nicht nutzbar (", conditionMessage(e),
          ") -> astronomische Naeherung")
  data.frame(Date = rv_ud, day_length = rv_daylength_h(rv_ud, RV_LAT))
})
rv_dd <- merge(rv_dd, rv_dl, by = "Date", all.x = TRUE, sort = FALSE)
rv_dd <- rv_dd[order(rv_dd$Id, rv_dd$Date), ]
rv_dd$csi_no_headache <- {
  m <- rowMeans(rv_dd[, setdiff(DRSP_ITEMS, "drsp_headache"), drop = FALSE], na.rm = TRUE)
  m[is.nan(m)] <- NA_real_; m
}
rv_valid_ids <- unique(rv_dd$Id[!is.na(rv_dd$cycle_phase_2)])
cat("[revision] rv_dd:", nrow(rv_dd), "Tage |", length(unique(rv_dd$Id)), "Personen |",
    length(rv_valid_ids), "mit gueltiger Phase\n")

# ---- Model functions ----
rv_fit_gamm <- function(outcome, data = rv_dd, k = 20, extra = "", knots_fix = TRUE,
                        ar_rho = NULL) {
  need <- unique(c(outcome, "pacts_ov", "Id", "Date",
                   all.vars(stats::as.formula(paste0("~ 1", extra)))))
  d <- data[stats::complete.cases(data[, need, drop = FALSE]), , drop = FALSE]
  d <- d[order(d$Id, d$Date), , drop = FALSE]
  if (nrow(d) < 50) stop("zu wenige Zeilen fuer ", outcome, " (", nrow(d), ")")
  nu <- length(unique(d$pacts_ov))
  k_used <- min(k, nu - 1)
  if (k_used < 4) stop("zu wenige eindeutige pacts_ov-Werte fuer ", outcome)
  d$Id <- factor(d$Id)
  d$ar_start <- !duplicated(d$Id)
  f  <- stats::as.formula(paste0(outcome, " ~ s(pacts_ov, bs = 'cc', k = ", k_used, ")",
                                 extra, " + s(Id, bs = 're')"))
  kn <- if (knots_fix) list(pacts_ov = c(-1, 1)) else NULL
  m <- if (is.null(ar_rho)) {
    mgcv::bam(f, data = d, method = "fREML", knots = kn)
  } else {
    mgcv::bam(f, data = d, method = "fREML", knots = kn, rho = ar_rho, AR.start = d$ar_start)
  }
  attr(m, "k_used") <- k_used
  m
}

# LMM on participant x phase means
rv_pp_lmm <- function(data, var, ids = NULL, phase_col = "cycle_phase_2") {
  d <- data
  d$ph <- as.character(d[[phase_col]])
  d <- d[!is.na(d$ph) & !is.na(d[[var]]), , drop = FALSE]
  if (!is.null(ids)) d <- d[d$Id %in% ids, , drop = FALSE]
  if (nrow(d) == 0) stop("keine Daten fuer ", var)
  pp <- stats::aggregate(x = list(y = d[[var]]),
                         by = list(Id = d$Id, cycle_phase_2 = d$ph), FUN = mean)
  pp$cycle_phase_2 <- factor(pp$cycle_phase_2, levels = c("follicular", "luteal"))
  if (length(unique(pp$Id)) < 5 || length(unique(pp$cycle_phase_2)) < 2)
    stop("zu wenige Personen/Phasen fuer ", var)
  lmerTest::lmer(y ~ cycle_phase_2 + (1 | Id), data = pp)
}

# Holm-adjusted p-values
rv_holm <- function(p) stats::p.adjust(p, method = "holm")

# ==== Metadata ====
rv_run("M1", {
  all_dates <- c(rv_get("light_raw")$Date, rv_get("actigraph_sleep")$Date, rv_get("surveys")$Date)
  all_dates <- all_dates[!is.na(all_dates)]
  add_result("M1", paste0("First recording date: ", format(min(all_dates)), "; last recording date: ",
                          format(max(all_dates)), " (all sources, all participants)."),
             paste0("Analysierter Zeitraum (daily_data, nach Trimmen): ",
                    format(min(rv_dd$Date, na.rm = TRUE)), " bis ", format(max(rv_dd$Date, na.rm = TRUE))))
})

rv_run(c("M2"), {
  last_sunday <- function(y, m) {
    d <- seq(as.Date(sprintf("%d-%02d-01", y, m)), by = "day", length.out = 31)
    d <- d[as.integer(format(d, "%m")) == m]
    d[max(which(as.integer(format(d, "%u")) == 7))]
  }
  yrs <- sort(unique(as.integer(format(rv_dd$Date[!is.na(rv_dd$Date)], "%Y"))))
  dst <- as.Date(unlist(lapply(yrs, function(y) c(last_sunday(y, 3), last_sunday(y, 10)))),
                 origin = "1970-01-01")
  rng <- stats::aggregate(Date ~ Id, data = rv_dd, FUN = range)
  hit <- vapply(seq_len(nrow(rng)), function(i)
    any(dst >= rng$Date[i, 1] & dst <= rng$Date[i, 2]), logical(1))
  tz_light <- attr(rv_get("light_raw")$Datetime, "tzone")
  tz_act   <- attr(rv_get("activity_epochs")$Datetime, "tzone")
  add_result("M2", paste0(
    "Timestamps are in ", TIMEZONE, " (light data tzone attribute: '", paste(tz_light, collapse = ","),
    "'; activity data: '", paste(tz_act, collapse = ","), "'). ",
    if (any(hit)) paste0("A daylight-saving transition (", paste(format(dst), collapse = ", "),
                         ") fell inside the recording window of ", sum(hit), " participant(s): ",
                         paste(rng$Id[hit], collapse = ", "), ".")
    else "No daylight-saving transition fell inside any participant's recording window."),
    paste0("Geprueft: letzter Sonntag im Maerz/Oktober der Jahre ", paste(yrs, collapse = ", ")))
})

rv_run("M5", {
  lr <- rv_get("light_raw")
  med <- tapply(as.numeric(lr$Datetime), lr$Id, function(x) stats::median(diff(sort(x))))
  add_result("M5", paste0("ActLumus sampling interval (median of timestamp differences): ",
                          stats::median(med, na.rm = TRUE), " s."),
             paste0("Pro Person: min ", min(med, na.rm = TRUE), " s, max ", max(med, na.rm = TRUE), " s"))
})

rv_run("M37", {
  ex <- rv_get("EXCLUSIONS")
  ids <- ex$Id[grepl("BMI", ex$reason)]
  s <- rv_get("surveys")
  b <- vapply(ids, function(i) { v <- s$bmi[s$Id == i]; v <- v[!is.na(v)]; if (length(v)) v[1] else NA_real_ }, numeric(1))
  add_result("M37", paste0("Participant(s) excluded for BMI above the inclusion range: ",
                           paste0(ids, " (BMI = ", fmt_num(b, 1), ")", collapse = "; "), "."),
             "Quelle: surveys$bmi (demographics) der in EXCLUSIONS gelisteten ID")
})

rv_run(c("M7", "M8", "M40"), {
  nm <- unique(c(names(rv_get("surveys")), names(rv_get("onboarding"))))
  dep <- nm[grepl("depress|sad|blue|down", nm, ignore.case = TRUE)]
  men <- nm[grepl("mens", nm, ignore.case = TRUE)]
  chr <- nm[grepl("chrono|meq|mctq|morning|evening|lark|owl", nm, ignore.case = TRUE)]
  add_result("M7", paste0("DRSP item 'felt depressed, sad, down or blue' in data: ",
                          if (length(dep)) paste("candidate columns:", paste(dep, collapse = ", "))
                          else "no matching column in surveys/onboarding. DRSP_ITEMS used: ",
                          if (!length(dep)) paste(DRSP_ITEMS, collapse = ", ") else ""),
             "Namenssuche in names(surveys) und names(onboarding); bitte gegen das Data Dictionary pruefen")
  add_result("M8", paste0("Variables documenting menses onset: ", paste(men, collapse = ", ")),
             "Spalten mit 'mens' im Namen; PACTS nutzt surveys$menses (1 = erster Tag) plus menses_extra")
  add_result("M40", paste0("Chronotype variables: ",
                           if (length(chr)) paste(chr, collapse = ", ") else "none found"),
             "Namenssuche (chrono|meq|mctq|morning|evening|lark|owl)")
  safe_write_csv(data.frame(column = nm), file.path(REV_DIR, "M7_M8_M40_column_names.csv"))
})

rv_run(c("M17", "M18", "M19"), {
  cnt <- function(nm, idcol) {
    o <- rv_get(nm); if (is.null(o)) return("object not found")
    paste0(nrow(o), " rows / ", length(unique(o[[idcol]])), " participants")
  }
  add_result("M17", paste0("Menses false-positive corrections (MENSES_FALSE_POSITIVES): ",
                           cnt("MENSES_FALSE_POSITIVES", "Id")), "")
  add_result("M18", paste0("Manually added menses entries (manual_menses): ", cnt("manual_menses", "Id")), "")
  add_result("M19", paste0("Merged duplicate menses markings (menses_merged): ", cnt("menses_merged", "id")),
             "Zeilen/Personen je Korrekturtyp; Objekte stammen aus 03_import_redcap.R / 06_pacts.R")
})

# R10: data completeness vs. study day
rv_run("R10", {
  d <- rv_dd[!is.na(rv_dd$study_day) & rv_dd$study_day >= 1 & rv_dd$study_day <= 40, ]
  d$compl <- (as.integer(!is.na(d$TAT_250_hr)) + as.integer(!is.na(d$TST_min)) +
                as.integer(!is.na(d$light_sensitivity_idx))) / 3
  ct <- suppressWarnings(stats::cor.test(d$study_day, d$compl, method = "spearman", exact = FALSE))
  add_result("R10", paste0("Spearman correlation between study day (1-40) and daily data completeness ",
                           "(share of the three sources with data on that day): rho = ",
                           fmt_num(unname(ct$estimate), 2), ", ", fmt_peq(ct$p.value), "; N = ", nrow(d),
                           " participant-days."),
             "EXPLIZITE Neu-Definition (Id x Tag, 3 Quellen: ActLumus/Actigraph/Fragebogen). Der urspruengliche Wert rho = .18 stammt aus einem Ad-hoc-Snippet; weicht der Wert hier ab, den Text anpassen")
})

# ==== Corrections: fixed boundary knots, t-based CIs ====
RV_GAMM_SET <- unique(c("TAT_250_hr", "light_sensitivity_idx", "cycle_symptom_idx",
                        "TAT_10_hr", "TAT_1000_hr", "RA_act", "centroid_LE_hr",
                        SECONDARY_LIGHT_VARS, "MEDI_mean", "MEDI_morning", "MEDI_evening",
                        paste0("ls_", 1:5), DRSP_ITEMS, "TST_min", "sleep_efficiency",
                        "drsp_impair_work", "drsp_impair_social", "drsp_impair_relations"))
RV_GAMM_SET <- RV_GAMM_SET[RV_GAMM_SET %in% names(rv_dd)]

rv_gamm_k <- list()
cat("[revision] Fitte", length(RV_GAMM_SET), "GAMMs mit festen Knoten c(-1, 1) ...\n")
for (v in RV_GAMM_SET) {
  rv_gamm_k[[v]] <- tryCatch(rv_fit_gamm(v), error = function(e) {
    message("  [revision] GAMM ", v, " fehlgeschlagen: ", conditionMessage(e)); NULL })
}

rv_old_gamm <- function(v) {
  m <- switch(v,
              TAT_250_hr = rv_get("gamm_TAT"), light_sensitivity_idx = rv_get("gamm_lsens"),
              cycle_symptom_idx = rv_get("gamm_csymp"), TAT_10_hr = rv_get("gamm_TAT10"),
              TAT_1000_hr = rv_get("gamm_TAT1000"), RA_act = rv_get("gamm_RA_act"), NULL)
  if (is.null(m) && grepl("^ls_", v)) m <- rv_get("gamm_ls_items")[[v]]
  if (is.null(m) && v %in% DRSP_ITEMS) m <- rv_get("gamm_drsp_items")[[v]]
  if (is.null(m) && v %in% SECONDARY_LIGHT_VARS) m <- rv_get("gamm_secondary_light")[[v]]
  m
}
rv_gam_stats <- function(m) {
  if (is.null(m)) return(c(edf = NA_real_, F = NA_real_, p = NA_real_))
  k <- tryCatch(gam_row(m), error = function(e) list(edf = NA_real_, F = NA_real_, p = NA_real_))
  c(edf = k$edf, F = k$F, p = k$p)
}

rv_run("KNOTS", {
  tb <- do.call(rbind, lapply(RV_GAMM_SET, function(v) {
    o <- rv_gam_stats(rv_old_gamm(v)); n <- rv_gam_stats(rv_gamm_k[[v]])
    data.frame(type = "GAMM", model = v, old_edf = o["edf"], new_edf = n["edf"],
               old_p = o["p"], new_p = n["p"],
               sig_changed = ifelse(is.na(o["p"]) | is.na(n["p"]), NA, (o["p"] < .05) != (n["p"] < .05)),
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  assign("rv_changed_gamm", tb, envir = .GlobalEnv)
  chg <- tb$model[which(tb$sig_changed %in% TRUE)]
  main <- tb[tb$model %in% c("TAT_250_hr", "light_sensitivity_idx", "cycle_symptom_idx", "RA_act", "centroid_LE_hr"), ]
  add_result("KNOTS", paste0(
    "All cyclic GAMMs were refitted with fixed boundary knots (knots = list(pacts_ov = c(-1, 1))). ",
    "Main models (old -> new): ",
    paste0(main$model, ": EDF ", fmt_num(main$old_edf), " -> ", fmt_num(main$new_edf),
           ", p ", fmt_p(main$old_p), " -> ", fmt_p(main$new_p), collapse = "; "), ". ",
    if (length(chg)) paste0("Models whose p-value crossed .05: ", paste(chg, collapse = ", "), ".")
    else "No model's p-value crossed .05."),
    "Details in changed_numbers.csv (type = GAMM). Alle Zahlen im Text/den Tabellen mit 'new_*' ersetzen, wenn die Korrektur uebernommen wird")
})

# t-based CIs for all LMM phase contrasts
rv_lmm_objects <- c("lmm_TAT", "lmm_MEDI", "lmm_lsens", "lmm_csymp", "lmm_TAT10", "lmm_TAT1000",
                    "lmm_IS", "lmm_IV", "lmm_RA", "lmm_TST", "lmm_sleep_eff",
                    "lmm_IS_full21", "lmm_IV_full21",
                    paste0("sens_compliance_", c("TAT", "MEDI", "lsens", "csymp")),
                    paste0("sens_strict_", c("TAT", "MEDI", "lsens", "csymp")))
rv_run("CI", {
  term <- "cycle_phase_2luteal"
  rows <- lapply(rv_lmm_objects, function(nm) {
    m <- rv_get(nm)
    if (is.null(m)) return(data.frame(type = "LMM_CI", model = nm, note = "Objekt fehlt", stringsAsFactors = FALSE))
    k <- lmm_coef(m, term)
    old <- tryCatch(suppressWarnings(suppressMessages(
      stats::confint(m, parm = term, level = .95, method = "profile", quiet = TRUE))),
      error = function(e) NULL)
    old_m <- "profile"
    if (is.null(old)) {
      old <- tryCatch(suppressWarnings(stats::confint(m, parm = term, level = .95, method = "Wald")),
                      error = function(e) NULL)
      old_m <- "Wald"
    }
    data.frame(type = "LMM_CI", model = nm, b = k$b, se = k$se, df = k$df, p = k$p,
               old_ci_low = if (is.null(old)) NA_real_ else as.numeric(old[1, 1]),
               old_ci_high = if (is.null(old)) NA_real_ else as.numeric(old[1, 2]),
               old_ci_method = if (is.null(old)) NA_character_ else old_m,
               new_ci_low = k$lo, new_ci_high = k$hi, new_ci_method = "t (Satterthwaite-df)",
               note = lmm_note(m), stringsAsFactors = FALSE)
  })
  rv_changed_ci <- dplyr::bind_rows(rows)
  show <- rv_changed_ci[rv_changed_ci$model %in% c("lmm_TAT", "lmm_lsens", "lmm_csymp", "lmm_RA"), ]
  add_result("CI", paste0(
    "95% confidence intervals of all categorical LMM phase contrasts were recomputed as t-based intervals ",
    "(b +/- qt(.975, df) * SE, Satterthwaite df). Examples (old profile CI -> new t CI): ",
    paste0(show$model, ": [", fmt_num(show$old_ci_low), ", ", fmt_num(show$old_ci_high), "] -> [",
           fmt_num(show$new_ci_low), ", ", fmt_num(show$new_ci_high), "]", collapse = "; "), "."),
    "Details in changed_numbers.csv (type = LMM_CI). Sekundaermetriken: siehe table2_secondary_metrics.csv")
  allchg <- dplyr::bind_rows(rv_get("rv_changed_gamm"), rv_changed_ci)
  safe_write_csv(allchg, file.path(REV_DIR, "changed_numbers.csv"))
})

# ---- R13 / Table 2: 16 secondary metrics, Holm over 16 LMM p-values ----
RV_SEC16 <- c(setdiff(SECONDARY_LIGHT_VARS, "M10_midpoint"), "MEDI_mean", "MEDI_morning", "MEDI_evening")
rv_run("R13", {
  rows <- lapply(RV_SEC16, function(v) {
    ml <- tryCatch(fit_lmm(paste0(v, "__mean")), error = function(e) NULL)
    mg <- rv_gamm_k[[v]]
    k <- if (!is.null(ml)) lmm_coef(ml) else NULL
    g <- rv_gam_stats(mg)
    data.frame(variable = v,
               b = if (is.null(k)) NA else k$b, ci_low = if (is.null(k)) NA else k$lo,
               ci_high = if (is.null(k)) NA else k$hi, t = if (is.null(k)) NA else k$t,
               df = if (is.null(k)) NA else k$df, p_lmm = if (is.null(k)) NA else k$p,
               gamm_edf = g["edf"], gamm_p = g["p"],
               singular = if (is.null(ml)) NA else lme4::isSingular(ml), stringsAsFactors = FALSE,
               row.names = NULL)
  })
  t2 <- do.call(rbind, rows)
  t2$p_holm_16 <- rv_holm(t2$p_lmm)
  t2$sig_nominal <- t2$p_lmm < .05
  t2$sig_holm <- t2$p_holm_16 < .05
  t2 <- t2[order(t2$p_lmm), ]
  safe_write_csv(t2, file.path(REV_DIR, "table2_secondary_metrics.csv"))
  assign("rv_table2", t2, envir = .GlobalEnv)
  add_result("R13", paste0(
    "Table 2 now contains 16 metrics (M10 midpoint dropped; daily mean, morning and evening mEDI added), ",
    "sorted by LMM p with Holm-adjusted p over the 16 LMM contrasts. Nominal p < .05: ",
    paste(t2$variable[which(t2$sig_nominal)], collapse = ", "), ". Holm p < .05: ",
    if (any(t2$sig_holm, na.rm = TRUE)) paste(t2$variable[which(t2$sig_holm)], collapse = ", ") else "none", "."),
    "Bonferroni-Schwelle bei 16 Tests: alpha = .0031; Tabelle: table2_secondary_metrics.csv")
  for (v in c("MEDI_mean", "MEDI_morning", "MEDI_evening")) {
    ml <- tryCatch(fit_lmm(paste0(v, "__mean")), error = function(e) NULL)
    add_result(paste0("R13-", v),
               paste0(if (!is.null(ml)) paste0("LMM: ", fmt_lmm(ml)) else "LMM nicht berechenbar",
                      "; GAMM: ", if (!is.null(rv_gamm_k[[v]])) fmt_gam(rv_gamm_k[[v]]) else "nicht berechenbar",
                      "; Holm p (16): ", fmt_p(t2$p_holm_16[t2$variable == v])),
               paste0(if (!is.null(ml)) lmm_note(ml) else "", " | ",
                      if (!is.null(rv_gamm_k[[v]])) gam_note(rv_gamm_k[[v]]) else ""))
  }
})

# ---- R24: GAMMs of daily TST and sleep efficiency ----
rv_run("R24", {
  for (v in c("TST_min", "sleep_efficiency")) {
    m <- rv_gamm_k[[v]]
    if (is.null(m)) { add_result(paste0("R24-", v), "NICHT VERFUEGBAR", "GAMM nicht berechenbar"); next }
    add_result(paste0("R24-", v), paste0("GAMM (cyclic spline, ovulation-centered cycle time): ", fmt_gam(m), "."),
               gam_note(m))
  }
})

# ---- R29: impairment items, Holm within the family ----
rv_run("R29", {
  v3 <- c("drsp_impair_work", "drsp_impair_social", "drsp_impair_relations")
  st <- t(vapply(v3, function(v) rv_gam_stats(rv_gamm_k[[v]]), numeric(3)))
  ph <- rv_holm(st[, "p"])
  for (i in seq_along(v3)) {
    m <- rv_gamm_k[[v3[i]]]
    add_result(paste0("R29-", sub("drsp_impair_", "", v3[i])),
               if (is.null(m)) "NICHT VERFUEGBAR" else paste0(fmt_gam(m), "; Holm-adjusted p (3 impairment items) = ", fmt_p(ph[i])),
               if (is.null(m)) "GAMM nicht berechenbar" else gam_note(m))
  }
})

# ---- Item GAMMs: Holm within the LS items (5) and the DRSP items (20) ----
rv_run("HOLM_ITEMS", {
  mk <- function(vs, fam) {
    st <- t(vapply(vs, function(v) rv_gam_stats(rv_gamm_k[[v]]), numeric(3)))
    ost <- t(vapply(vs, function(v) rv_gam_stats(rv_old_gamm(v)), numeric(3)))
    data.frame(family = fam, item = vs, edf = st[, "edf"], p_new_knots = st[, "p"],
               p_old = ost[, "p"], p_holm = rv_holm(st[, "p"]),
               p_bh_fdr = stats::p.adjust(st[, "p"], method = "BH"), row.names = NULL)
  }
  tb <- rbind(mk(paste0("ls_", 1:5), "Light Sensitivity items (5)"), mk(DRSP_ITEMS, "DRSP items (20)"))
  tb <- tb[order(tb$family, tb$p_new_knots), ]
  safe_write_csv(tb, file.path(REV_DIR, "item_gamm_holm.csv"))
  s5 <- tb[tb$family == "Light Sensitivity items (5)" & tb$p_new_knots < .05, ]
  s20 <- tb[tb$family == "DRSP items (20)" & tb$p_new_knots < .05, ]
  add_result("HOLM_ITEMS", paste0(
    "Item-level GAMMs, Holm-adjusted within family. Light Sensitivity items nominally significant: ",
    if (nrow(s5)) paste0(s5$item, " (", fmt_peq(s5$p_new_knots), ", Holm ", fmt_peq(s5$p_holm), ")", collapse = "; ") else "none",
    ". DRSP items nominally significant: ",
    if (nrow(s20)) paste0(s20$item, " (", fmt_peq(s20$p_new_knots), ", Holm ", fmt_peq(s20$p_holm), ")", collapse = "; ") else "none", "."),
    "ZUSATZ, nicht im Prompt; p mit festen Knoten gerechnet (p_old = bisherige Werte). Tabelle: item_gamm_holm.csv")
})

# ==== Confounding of cycle phase with study time ====
rv_run(c("R6", "R7", "R8", "R9"), {
  d <- rv_dd[!is.na(rv_dd$cycle_phase_2), ]
  ph_by_id <- tapply(as.character(d$cycle_phase_2), d$Id, function(x) paste(sort(unique(x)), collapse = "+"))
  n_lut <- sum(ph_by_id == "luteal"); n_fol <- sum(ph_by_id == "follicular"); n_both <- sum(ph_by_id == "follicular+luteal")
  add_result("R6", paste0(n_lut, " of ", length(ph_by_id), " participants with an assigned phase contributed only luteal-phase data (",
                          n_fol, " only follicular, ", n_both, " both phases)."), "Basis: Tage mit cycle_phase_2 != NA")
  mf <- mean(d$study_day[d$cycle_phase_2 == "follicular"], na.rm = TRUE)
  ml <- mean(d$study_day[d$cycle_phase_2 == "luteal"], na.rm = TRUE)
  add_result("R7", paste0("Mean study day in the follicular phase: ", fmt_num(mf, 1), "."), "")
  add_result("R8", paste0("Mean study day in the luteal phase: ", fmt_num(ml, 1), "."), "")
  m <- lmerTest::lmer(study_day ~ cycle_phase_2 + (1 | Id), data = d[!is.na(d$study_day), ])
  add_result("R9", paste0("Linear mixed model of study day on phase: ", fmt_lmm(m, d = 1), " (days)."), lmm_note(m))
})

# R32: missingness of actigraphy data with study day
rv_run("R32", {
  d <- rv_dd
  d$is_missing_sleep <- as.integer(is.na(d$TST_min))
  d$study_day10 <- d$study_day / 10
  d <- d[stats::complete.cases(d[, c("is_missing_sleep", "cycle_phase_2", "cycle_symptom_idx", "study_day10", "Id")]), ]
  ctl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  m <- lme4::glmer(is_missing_sleep ~ cycle_phase_2 + cycle_symptom_idx + study_day10 + (1 | Id),
                   data = d, family = binomial, control = ctl)
  s <- summary(m)$coefficients
  f <- function(term) {
    b <- s[term, "Estimate"]; se <- s[term, "Std. Error"]
    paste0("OR = ", fmt_num(exp(b), 2), ", 95% CI [", fmt_num(exp(b - 1.96 * se), 2), ", ",
           fmt_num(exp(b + 1.96 * se), 2), "], ", fmt_peq(s[term, "Pr(>|z|)"]))
  }
  add_result("R32", paste0("Missingness of actigraphy (no TST on a day), logistic mixed model: luteal vs follicular ",
                           f("cycle_phase_2luteal"), "; cycle symptom index (per point) ", f("cycle_symptom_idx"),
                           "; study day (per 10 days) ", f("study_day10"), "."),
             paste0("N = ", nrow(d), " Tage / ", length(unique(d$Id)), " Personen; Wald-KIs; Phase: Referenz = follicular"))
})

# ==== Robustness ====
# sleep period windows from the GGIR night summary
rv_spt_windows <- function(sl) {
  sl <- as.data.frame(sl)
  sl <- sl[!is.na(sl$sleep_onset_hr) & !is.na(sl$sleep_offset_hr) & !is.na(sl$Date), ]
  base <- as.numeric(as.POSIXct(paste(format(as.Date(sl$Date), "%Y-%m-%d"), "00:00:00"), tz = TIMEZONE))
  on_h <- sl$sleep_onset_hr
  on_h <- ifelse(on_h < 12, on_h + 24, on_h)  # onset before 12:00 = after midnight
  st <- base + on_h * 3600
  en <- base + sl$sleep_offset_hr * 3600
  en <- ifelse(en <= st, en + 86400, en)  # wake-up on the next day
  data.frame(Id = as.character(sl$Id), Date = as.Date(sl$Date), start = st, end = en,
             stringsAsFactors = FALSE)
}
rv_mask_spt <- function(df, win) {
  inside <- rep(FALSE, nrow(df))
  for (id in unique(df$Id)) {
    w <- win[win$Id == id, ]
    if (nrow(w) == 0) next
    w <- w[order(w$start), ]
    ii <- which(df$Id == id)
    tt <- as.numeric(df$Datetime[ii])
    k <- findInterval(tt, w$start)
    inside[ii] <- k > 0 & tt < w$end[pmax(k, 1)]
  }
  inside
}
rv_win <- tryCatch(rv_spt_windows(actigraph_sleep), error = function(e) NULL)

# R15: remove light within the GGIR sleep periods and recompute the metrics
if (rv_on("R15")) {
  rv_run(c("R15", "R15-diag"), {
    if (is.null(rv_win) || nrow(rv_win) == 0) stop("keine Schlaffenster aus actigraph_sleep ableitbar")
    dur <- (rv_win$end - rv_win$start) / 3600
    add_result("R15-diag", paste0("SPT windows reconstructed from GGIR night summary: ", nrow(rv_win),
                                  " nights; duration median ", fmt_num(stats::median(dur), 1), " h (range ",
                                  fmt_num(min(dur), 1), "-", fmt_num(max(dur), 1), " h)."),
               paste0("PLAUSIBILITAETSCHECK: sleep_onset_hr Bereich ",
                      paste(round(range(actigraph_sleep$sleep_onset_hr, na.rm = TRUE), 1), collapse = "-"),
                      ", sleep_offset_hr Bereich ",
                      paste(round(range(actigraph_sleep$sleep_offset_hr, na.rm = TRUE), 1), collapse = "-"),
                      ". Wirken Dauern unplausibel (<3 h oder >14 h), stimmt die Zeitkonvention nicht"))
    lt <- light_trimmed
    lt$in_spt <- rv_mask_spt(lt, rv_win)
    share <- mean(lt$in_spt)
    lt2 <- lt[!lt$in_spt, ]
    cat("[revision] R15: ", round(100 * share, 1), "% der Licht-Epochen liegen in Schlaffenstern\n", sep = "")
    daily2 <- lt2 |>
      dplyr::group_by(Id, Date) |>
      dplyr::summarise(
        n_epochs_wake = dplyr::n(),
        MEDI_mean = mean(MEDI, na.rm = TRUE),
        TAT_250_hr = tryCatch(as.numeric(duration_above_threshold(MEDI, Datetime, threshold = 250,
                                                                  as.df = FALSE)) / 3600,
                              error = function(e) NA_real_),
        centroid_LE_hr = tryCatch({
          cle <- centroidLE(MEDI, Datetime)
          if (!is.null(cle) && !is.na(cle)) lubridate::hour(cle) + lubridate::minute(cle) / 60 +
            lubridate::second(cle) / 3600 else NA_real_
        }, error = function(e) NA_real_),
        .groups = "drop")
    d2 <- merge(as.data.frame(daily2), rv_dd[, c("Id", "Date", "cycle_phase_2")], by = c("Id", "Date"))
    assign("rv_light_awake_daily", d2, envir = .GlobalEnv)
    safe_write_csv(d2, file.path(REV_DIR, "R15_light_awake_daily.csv"))
    m_c <- rv_pp_lmm(d2, "centroid_LE_hr"); m_c0 <- rv_pp_lmm(rv_dd, "centroid_LE_hr")
    add_result("R15", paste0("After setting light data within GGIR sleep-period windows to missing (",
                             fmt_num(100 * share, 1), "% of epochs removed), the centroid of light exposure still differed by phase: ",
                             fmt_lmm(m_c), " (original, all-day data: ", fmt_lmm(m_c0), ")."),
               lmm_note(m_c))
    for (v in c("TAT_250_hr", "MEDI_mean")) {
      mv <- tryCatch(rv_pp_lmm(d2, v), error = function(e) NULL)
      if (!is.null(mv)) add_result(paste0("R15-", v), paste0(v, " (waking hours only): ", fmt_lmm(mv)), lmm_note(mv))
    }
  })
} else rv_skip_note(c("R15", "R15-diag"), "R15")

# R16 / R33 / R34: high-completeness subsample (hc_ids as in 08; strict = both phases)
rv_run(c("R16", "R33", "R34"), {
  hc_def <- rv_get("hc_ids")
  if (is.null(hc_def)) stop("hc_ids fehlt (08_missing_data.R nicht gelaufen?)")
  hc_strict <- rv_get("hc_ids_strict")
  strict_src <- "hc_ids_strict aus 08_missing_data.R"
  if (is.null(hc_strict)) {
    mr <- rv_get("missing_rep")
    if (is.null(mr)) stop("weder hc_ids_strict noch missing_rep vorhanden")
    mn <- tapply(mr$pct_light, as.character(mr$Id), min)
    kk <- tapply(mr$pct_light, as.character(mr$Id), length)
    hc_strict <- names(mn)[kk[names(mn)] == 2 & mn >= COMPL_THR]
    strict_src <- "hc_ids_strict fehlte in der Session und wurde hier aus missing_rep nachgerechnet (beide Phasen, jeweils >= COMPL_THR %)"
  }
  hc <- list(`as defined` = hc_def, strict = hc_strict)
  per_set <- function(fun) {
    paste(vapply(names(hc), function(nm) {
      ids <- hc[[nm]]
      if (length(ids) == 0) return(paste0(nm, " (0 participants): not computable"))
      tryCatch(paste0(nm, " (", length(ids), " participants): ", fun(ids)),
               error = function(e) paste0(nm, ": not computable (", conditionMessage(e), ")"))
    }, character(1)), collapse = " | ")
  }
  add_result("R16", paste0("Centroid, high-completeness subsample - ",
                           per_set(function(ids) fmt_lmm(rv_pp_lmm(rv_dd, "centroid_LE_hr", ids = ids)))),
             paste0("as defined = hc_ids: >= COMPL_THR % ActLumus-Tage in jeder VORHANDENEN Phase (08_missing_data.R); strict: ", strict_src))
  add_result("R33", paste0("Relative amplitude, high-completeness subsample - ",
                           per_set(function(ids) fmt_lmm(rv_pp_lmm(rv_dd, "RA_act", ids = ids)))),
             "dieselbe Teilstichprobe wie R16 (ActLumus-basiert definiert)")
  add_result("R34", paste0("Cycle Symptom Index GAMM, high-completeness subsample - ",
                           per_set(function(ids) {
                             m <- rv_fit_gamm("cycle_symptom_idx", data = rv_dd[rv_dd$Id %in% ids, ])
                             paste0(fmt_gam(m), " [", gam_note(m), "]")
                           })),
             "dieselbe Teilstichprobe wie R16")
})

# R35: without participants aged over 35 or with BMI over 26
rv_run("R35", {
  ex <- char_table$Id[which(char_table$age > RV_AGE_MAX | char_table$bmi > RV_BMI_MAX)]
  keep <- setdiff(unique(rv_dd$Id), ex)
  lmm_vars <- c("TAT_250_hr", "light_sensitivity_idx", "cycle_symptom_idx", "centroid_LE_hr",
                "RA_act", "TST_min", "sleep_efficiency")
  rows <- list()
  for (v in lmm_vars) for (sub in c("full", "restricted")) {
    ids <- if (sub == "full") NULL else keep
    m <- tryCatch(rv_pp_lmm(rv_dd, v, ids = ids), error = function(e) NULL)
    if (!is.null(m)) { k <- lmm_coef(m)
      rows[[length(rows) + 1]] <- data.frame(analysis = v, model = "LMM", sample = sub,
        n_ids = lme4::ngrps(m)[["Id"]], b_or_edf = k$b, p = k$p, text = fmt_lmm(m), stringsAsFactors = FALSE) }
  }
  for (v in c("IS", "IV")) for (sub in c("full", "restricted")) {
    ids <- if (sub == "full") npcra_ok_ids else intersect(npcra_ok_ids, keep)
    npd <- as.data.frame(npcra_phase)
    npd <- npd[!is.na(npd$cycle_phase_2) & npd$Id %in% ids & !is.na(npd[[v]]), ]
    m <- tryCatch(lmerTest::lmer(stats::as.formula(paste(v, "~ cycle_phase_2 + (1 | Id)")), data = npd), error = function(e) NULL)
    if (!is.null(m)) { k <- lmm_coef(m)
      rows[[length(rows) + 1]] <- data.frame(analysis = v, model = "LMM", sample = sub,
        n_ids = lme4::ngrps(m)[["Id"]], b_or_edf = k$b, p = k$p, text = fmt_lmm(m), stringsAsFactors = FALSE) }
  }
  for (v in c("TAT_250_hr", "light_sensitivity_idx", "cycle_symptom_idx", "centroid_LE_hr", "RA_act")) for (sub in c("full", "restricted")) {
    dat <- if (sub == "full") rv_dd else rv_dd[rv_dd$Id %in% keep, ]
    m <- tryCatch(rv_fit_gamm(v, data = dat), error = function(e) NULL)
    if (!is.null(m)) { k <- gam_row(m)
      rows[[length(rows) + 1]] <- data.frame(analysis = v, model = "GAMM", sample = sub,
        n_ids = length(unique(m$model$Id)), b_or_edf = k$edf, p = k$p, text = fmt_gam(m), stringsAsFactors = FALSE) }
  }
  tb <- do.call(rbind, rows)
  safe_write_csv(tb, file.path(REV_DIR, "R35_sensitivity_age_bmi.csv"))
  r <- tb[tb$sample == "restricted", ]
  add_result("R35", paste0("Excluding ", length(ex), " participant(s) with age > ", RV_AGE_MAX, " or BMI > ", RV_BMI_MAX,
                           " (", paste(ex, collapse = ", "), "): ",
                           paste0(r$analysis, " [", r$model, "] ", r$text, collapse = "; "), "."),
             "Vergleich mit voller Stichprobe: R35_sensitivity_age_bmi.csv (Spalte sample)")
})

# ==== Rest-activity and sleep ====
rv_hourly <- NULL
rv_run(c("R19", "R20"), {
  a <- activity_trimmed
  a$hour_start <- as.POSIXct(floor(as.numeric(a$Datetime) / 3600) * 3600, origin = "1970-01-01", tz = TIMEZONE)
  h <- a |>
    dplyr::group_by(Id, hour_start) |>
    dplyr::summarise(ACC = mean(ACC, na.rm = TRUE), n_min = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n_min >= RV_MIN_MIN_PER_HOUR)
  h <- as.data.frame(h)
  h$Date <- as.Date(h$hour_start, tz = TIMEZONE)
  h$hr <- as.POSIXlt(h$hour_start, tz = TIMEZONE)$hour
  h$th <- as.numeric(h$hour_start) / 3600
  h <- merge(h, rv_dd[, c("Id", "Date", "cycle_phase_2")], by = c("Id", "Date"), all.x = TRUE, sort = FALSE)
  assign("rv_hourly", h[order(h$Id, h$hour_start), ], envir = .GlobalEnv)
})

# IS/IV (Witting formulas) from hourly values
rv_is_iv <- function(x, hh, th, min_n = 48) {
  ok <- !is.na(x); x <- x[ok]; hh <- hh[ok]; th <- th[ok]
  n <- length(x)
  if (n < min_n) return(c(IS = NA_real_, IV = NA_real_))
  xb <- mean(x); den <- sum((x - xb)^2)
  if (den == 0) return(c(IS = NA_real_, IV = NA_real_))
  hm <- tapply(x, hh, mean); hn <- tapply(x, hh, length)
  IS <- sum(hn * (hm - xb)^2) / den
  o <- order(th); x <- x[o]; th <- th[o]
  adj <- diff(th) == 1
  m <- sum(adj)
  IV <- if (m < 24) NA_real_ else (n * sum(diff(x)[adj]^2)) / (m * den)
  c(IS = IS, IV = IV)
}

# R19: IS/IV from hourly data
rv_run("R19", {
  if (is.null(rv_hourly)) stop("stuendliche Aktivitaetsdaten nicht verfuegbar")
  h <- rv_hourly[!is.na(rv_hourly$cycle_phase_2), ]
  g <- split(h, list(h$Id, h$cycle_phase_2), drop = TRUE)
  rows <- lapply(g, function(x) {
    r <- rv_is_iv(x$ACC, x$hr, x$th)
    data.frame(Id = x$Id[1], cycle_phase_2 = as.character(x$cycle_phase_2[1]), IS_h = r["IS"], IV_h = r["IV"],
               n_days = length(unique(x$Date)), row.names = NULL)
  })
  tb <- do.call(rbind, rows)
  tb$cycle_phase_2 <- factor(tb$cycle_phase_2, levels = c("follicular", "luteal"))
  tb <- merge(tb, npcra_phase[, c("Id", "cycle_phase_2", "IS", "IV")], by = c("Id", "cycle_phase_2"), all.x = TRUE)
  safe_write_csv(tb, file.path(REV_DIR, "R19_IS_IV_hourly.csv"))
  ok <- tb[tb$Id %in% npcra_ok_ids, ]
  mk <- function(var, dd) lmerTest::lmer(stats::as.formula(paste(var, "~ cycle_phase_2 + (1 | Id)")), data = dd[!is.na(dd[[var]]), ])
  m_is <- mk("IS_h", ok); m_iv <- mk("IV_h", ok)
  cr <- stats::cor(tb$IS, tb$IS_h, use = "complete.obs")
  cr2 <- stats::cor(tb$IV, tb$IV_h, use = "complete.obs")
  add_result("R19", paste0("IS and IV recomputed from hourly aggregated activity (mean ENMO per clock hour, >= ", RV_MIN_MIN_PER_HOUR,
                           " valid minutes per hour; participants with >= ", MIN_NPCRA_DAYS, " days in both phases): IS ",
                           fmt_lmm(m_is), "; IV ", fmt_lmm(m_iv), ". Agreement with the 1-minute-based values (person x phase): r = ",
                           fmt_num(cr, 2), " (IS), r = ", fmt_num(cr2, 2), " (IV)."),
             paste0("IS/IV hier mit eigener Witting-Implementierung (rv_is_iv); ", lmm_note(m_is)))
})

# R20: equal-length windows (5 random valid days per participant and phase)
if (rv_on("R20")) {
  rv_run("R20", {
    if (is.null(rv_hourly)) stop("stuendliche Aktivitaetsdaten nicht verfuegbar")
    h <- rv_hourly[!is.na(rv_hourly$cycle_phase_2), ]
    ra <- rv_dd[, c("Id", "Date", "RA_act")]
    grp <- split(h, list(h$Id, h$cycle_phase_2), drop = TRUE)
    prep <- lapply(grp, function(x) {
      cnt <- tapply(x$hr, x$Date, length)
      days <- sort(as.Date(names(cnt)[cnt >= RV_MIN_HOURS_PER_DAY]))
      if (length(days) < RV_N_DAYS) return(NULL)
      H <- matrix(NA_real_, length(days), 24)
      for (i in seq_along(days)) { xx <- x[x$Date == days[i], ]; H[i, xx$hr + 1] <- xx$ACC }
      r <- ra$RA_act[match(paste(x$Id[1], days), paste(ra$Id, ra$Date))]
      list(Id = x$Id[1], phase = as.character(x$cycle_phase_2[1]), days = as.numeric(days), H = H, ra = r)
    })
    prep <- Filter(Negate(is.null), prep)
    set.seed(RV_SEED)
    acc <- lapply(prep, function(p) c(IS = 0, IV = 0, RA = 0, nIS = 0, nIV = 0, nRA = 0))
    for (rep in seq_len(RV_N_REP)) {
      for (g in seq_along(prep)) {
        p <- prep[[g]]
        idx <- sort(sample.int(length(p$days), RV_N_DAYS))
        Hs <- p$H[idx, , drop = FALSE]
        x <- as.vector(t(Hs)); hh <- rep(0:23, times = length(idx))
        th <- as.vector(t(outer(p$days[idx] * 24, 0:23, "+")))
        r <- rv_is_iv(x, hh, th, min_n = 24)
        ra_m <- mean(p$ra[idx], na.rm = TRUE)
        a <- acc[[g]]
        if (!is.na(r["IS"])) { a["IS"] <- a["IS"] + r["IS"]; a["nIS"] <- a["nIS"] + 1 }
        if (!is.na(r["IV"])) { a["IV"] <- a["IV"] + r["IV"]; a["nIV"] <- a["nIV"] + 1 }
        if (!is.nan(ra_m) && !is.na(ra_m)) { a["RA"] <- a["RA"] + ra_m; a["nRA"] <- a["nRA"] + 1 }
        acc[[g]] <- a
      }
    }
    tb <- do.call(rbind, lapply(seq_along(prep), function(g) {
      a <- acc[[g]]
      data.frame(Id = prep[[g]]$Id, cycle_phase_2 = prep[[g]]$phase,
                 IS = a["IS"] / a["nIS"], IV = a["IV"] / a["nIV"], RA = a["RA"] / a["nRA"], row.names = NULL)
    }))
    tb$cycle_phase_2 <- factor(tb$cycle_phase_2, levels = c("follicular", "luteal"))
    both <- names(which(table(tb$Id) == 2)); tb <- tb[tb$Id %in% both, ]
    safe_write_csv(tb, file.path(REV_DIR, "R20_equal_window_means.csv"))
    mk <- function(var) lmerTest::lmer(stats::as.formula(paste(var, "~ cycle_phase_2 + (1 | Id)")), data = tb[!is.na(tb[[var]]), ])
    m1 <- mk("IS"); m2 <- mk("IV"); m3 <- mk("RA")
    add_result("R20", paste0("Equal-length windows (", RV_N_DAYS, " randomly drawn valid days per participant and phase, ", RV_N_REP,
                             " repetitions, means per participant and phase; ", length(both), " participants with both phases): IS ",
                             fmt_lmm(m1), "; IV ", fmt_lmm(m2), "; RA ", fmt_lmm(m3), "."),
               paste0("gueltiger Tag = >= ", RV_MIN_HOURS_PER_DAY, " Stunden mit Daten; IS/IV mit rv_is_iv (stuendlich); ", lmm_note(m1)))
  })
} else rv_skip_note("R20", "R20")

# R22: M10 and L5 activity levels
rv_run("R22", {
  m10 <- rv_pp_lmm(rv_dd, "M10_act_mean"); l5 <- rv_pp_lmm(rv_dd, "L5_act_mean")
  k1 <- lmm_coef(m10); k2 <- lmm_coef(l5)
  dir <- if (k1$p < .05 && k2$p < .05 && k1$b > 0 && k2$b < 0) "both (M10 rises and L5 falls)"
         else if (k1$p < .05 && k1$b > 0) "M10 rises"
         else if (k2$p < .05 && k2$b < 0) "L5 falls"
         else if (k1$p < .05 || k2$p < .05) "a significant change in the opposite direction or in one measure only (see estimates)"
         else "neither M10 nor L5 changes significantly"
  add_result("R22", paste0("M10 activity level (most active 10 h): ", fmt_lmm(m10), "; L5 activity level (least active 5 h): ",
                           fmt_lmm(l5), "."), lmm_note(m10))
  add_result("R22-D", paste0("Result for the interpretation: ", dir, " (nominal p < .05)."),
             "Einheit: ENMO (mg) wie von GGIR geliefert")
})

# R23: sleep midpoint (hours since 12:00)
rv_run("R23", {
  if (is.null(rv_win)) stop("keine Schlaffenster")
  w <- rv_win
  w$sleep_mid_h12 <- ((w$start + w$end) / 2 - (as.numeric(as.POSIXct(paste(format(w$Date, "%Y-%m-%d"), "12:00:00"), tz = TIMEZONE)))) / 3600
  d <- merge(rv_dd, w[, c("Id", "Date", "sleep_mid_h12")], by = c("Id", "Date"), all.x = TRUE, sort = FALSE)
  m <- rv_pp_lmm(d, "sleep_mid_h12"); k <- lmm_coef(m)
  mc <- tryCatch(rv_pp_lmm(rv_dd, "centroid_LE_hr"), error = function(e) NULL)
  kc <- if (is.null(mc)) NULL else lmm_coef(mc)
  same <- if (is.null(kc)) "centroid contrast not available" else
    if (sign(k$b) == sign(kc$b)) "same direction as the centroid shift (luteal vs follicular)" else "opposite direction to the centroid shift"
  add_result("R23", paste0("Sleep midpoint (hours after 12:00 of the night's start date; positive = later): ", fmt_lmm(m), "."), lmm_note(m))
  add_result("R23-D", paste0("The sleep midpoint shifts in the ", same, " (sleep midpoint b = ", fmt_num(k$b),
                             " h", if (!is.null(kc)) paste0("; centroid b = ", fmt_num(kc$b), " h") else "", ")."),
             "Richtung = Vorzeichen der luteal-vs-follikulaer-Differenz, unabhaengig von Signifikanz")
})

# ==== Associations: repeated-measures correlations ====
rv_rmcorr <- function(d, x, y) {
  d <- d[stats::complete.cases(d[, c("Id", x, y)]), c("Id", x, y)]
  d$Id <- factor(d$Id)
  keep <- names(which(table(d$Id) >= 3)); d <- d[d$Id %in% keep, ]
  d$Id <- droplevels(d$Id); d <- as.data.frame(d)
  out <- tryCatch(rmcorr::rmcorr(participant = "Id", measure1 = x, measure2 = y, dataset = d),
                  error = function(e) eval(bquote(rmcorr::rmcorr(participant = .(as.name("Id")),
                                                                  measure1 = .(as.name(x)),
                                                                  measure2 = .(as.name(y)), dataset = d))))
  list(r = unname(out$r), lo = unname(out$CI[1]), hi = unname(out$CI[2]), p = unname(out$p),
       df = unname(out$df), n = nrow(d), n_id = nlevels(d$Id))
}
fmt_rm <- function(o) paste0("r_rm = ", fmt_num(o$r, 2), ", 95% CI [", fmt_num(o$lo, 2), ", ", fmt_num(o$hi, 2), "], ",
                             fmt_peq(o$p), " (df = ", o$df, ")")

rv_rm_vars <- c(LSI = "light_sensitivity_idx", CSI = "cycle_symptom_idx", TAT250 = "TAT_250_hr",
                Centroid = "centroid_LE_hr", RA = "RA_act", TST = "TST_min", SE = "sleep_efficiency")
rv_rm_mat <- NULL
rv_run(c("R25", "R26", "R27", "R28"), {
  if (!requireNamespace("rmcorr", quietly = TRUE))
    stop("Paket 'rmcorr' fehlt: install.packages(\"rmcorr\")")
  nm <- names(rv_rm_vars)
  prs <- utils::combn(nm, 2, simplify = FALSE)
  rm_first_err <- NULL
  rows <- lapply(prs, function(p) {
    o <- tryCatch(rv_rmcorr(rv_dd, rv_rm_vars[[p[1]]], rv_rm_vars[[p[2]]]),
                  error = function(e) { if (is.null(rm_first_err)) rm_first_err <<- conditionMessage(e); NULL })
    if (is.null(o)) return(data.frame(var1 = p[1], var2 = p[2], note = "nicht berechenbar", stringsAsFactors = FALSE))
    data.frame(var1 = p[1], var2 = p[2], r_rm = o$r, ci_low = o$lo, ci_high = o$hi, p = o$p, df = o$df,
               n_obs = o$n, n_ids = o$n_id, stringsAsFactors = FALSE)
  })
  tb <- dplyr::bind_rows(rows)
  if (!("p" %in% names(tb)))
    stop("rmcorr fuer ALLE Variablenpaare fehlgeschlagen; erste Fehlermeldung: ", rm_first_err)
  tb$p_holm <- rv_holm(tb$p)
  safe_write_csv(tb, file.path(REV_DIR, "R25_rmcorr_pairs.csv"))
  M <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm)); diag(M) <- 1
  for (i in seq_len(nrow(tb))) { M[tb$var1[i], tb$var2[i]] <- tb$r_rm[i]; M[tb$var2[i], tb$var1[i]] <- tb$r_rm[i] }
  assign("rv_rm_mat", M, envir = .GlobalEnv); assign("rv_rm_tab", tb, envir = .GlobalEnv)
  sig <- tb[which(tb$p_holm < .05), ]
  add_result("R25", paste0("Repeated-measures correlations (within-person) among LSI, CSI, TAT250, centroid, RA, TST and sleep efficiency: ",
                           nrow(sig), " of ", nrow(tb), " pairs significant after Holm correction",
                           if (nrow(sig)) paste0(" (", paste0(sig$var1, "-", sig$var2, " r_rm = ", fmt_num(sig$r_rm, 2), collapse = "; "), ")") else "", "."),
             "alle Paare: R25_rmcorr_pairs.csv; Matrix fuer Abbildung 3")
  o26 <- rv_rmcorr(rv_dd, "light_sensitivity_idx", "cycle_symptom_idx")
  add_result("R26", paste0("Light Sensitivity Index - Cycle Symptom Index (repeated-measures correlation): ", fmt_rm(o26), "."),
             paste0(o26$n, " Tage / ", o26$n_id, " Personen"))
  o27 <- rv_rmcorr(rv_dd, "light_sensitivity_idx", "csi_no_headache")
  add_result("R27", paste0("Same, with a Cycle Symptom Index excluding the headache item: ", fmt_rm(o27), "."),
             "Grund: LS-Item 4 (Kopfschmerz/Augenbelastung) ueberlappt inhaltlich mit DRSP-Headache")
  o28 <- rv_rmcorr(rv_dd, "light_sensitivity_idx", "TAT_250_hr")
  add_result("R28", paste0("Light Sensitivity Index - TAT250 (repeated-measures correlation): ", fmt_rm(o28), "."),
             paste0(o28$n, " Tage / ", o28$n_id, " Personen"))
})

# ==== Write result files ====
safe_write_csv(results, file.path(REV_DIR, "results_for_placeholders.csv"))
rv_esc <- function(x) gsub("\n", "<br>", gsub("|", "\\|", x, fixed = TRUE), fixed = TRUE)
rv_md <- c("| ID | Text | Hinweis |", "|---|---|---|",
           paste0("| ", rv_esc(results$id), " | ", rv_esc(results$text), " | ", rv_esc(results$note), " |"))
rv_md_file <- safe_write_lines(rv_md, file.path(REV_DIR, "results_for_placeholders.md"))
safe_write_lines(utils::capture.output(utils::sessionInfo()), file.path(REV_DIR, "sessionInfo.txt"))

n_bad <- sum(results$text == "NICHT VERFUEGBAR")
cat("\n[revision] Fertig:", nrow(results), "Ergebniszeilen,", n_bad, "davon 'NICHT VERFUEGBAR'.\n")
cat("  ->", rv_md_file, "\n")
if (n_bad > 0) {
  cat("  Betroffene IDs und Fehlermeldungen:\n")
  print(results[results$text == "NICHT VERFUEGBAR", c("id", "note")], row.names = FALSE)
}
