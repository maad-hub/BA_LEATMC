# 09b_robustness.R - Robustness checks for the nominally significant findings
# (Cycle Symptom Index trajectory, centroid of light exposure, relative amplitude)

options(dplyr.summarise.inform = FALSE)

if (!exists("safe_write_csv", mode = "function")) {
  safe_write_csv <- function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    cands <- unique(c(path,
                      sub("(\\.[A-Za-z0-9]+)$", paste0("_", stamp, "\\1"), path),
                      file.path(tempdir(), basename(path))))
    why <- character()
    for (i in seq_along(cands)) {
      r <- tryCatch({ utils::write.csv(x, cands[i], row.names = FALSE, fileEncoding = "UTF-8"); NULL },
                    warning = function(w) conditionMessage(w),
                    error   = function(e) conditionMessage(e))
      if (is.null(r) && file.exists(cands[i])) {
        if (i > 1) message("[safe_write_csv] '", path, "' nicht beschreibbar (", paste(unique(why), collapse = " | "),
                           "). Stattdessen geschrieben: ", cands[i])
        return(invisible(cands[i]))
      }
      why <- c(why, if (is.null(r)) "Datei nach dem Schreiben nicht vorhanden" else r)
    }
    message("[safe_write_csv] FEHLGESCHLAGEN fuer '", path, "': ", paste(unique(why), collapse = " | "))
    invisible(NA_character_)
  }
}

# ---- Helpers ----
daylength_h <- function(date, lat = 48.137) {
  doy  <- as.integer(format(date, "%j"))
  decl <- 23.44 * pi / 180 * sin(2 * pi * (284 + doy) / 365)
  phi  <- lat * pi / 180
  24 * acos(pmin(1, pmax(-1, -tan(phi) * tan(decl)))) / pi
}

dd <- daily_data |>
  mutate(dl   = daylength_h(Date),
         dl_c = dl - mean(dl, na.rm = TRUE))

# participants with corrected menses data
corrected_ids <- unique(c(
  if (exists("MENSES_FALSE_POSITIVES")) MENSES_FALSE_POSITIVES$Id,
  if (exists("manual_menses"))          manual_menses$Id,
  if (exists("menses_merged"))          menses_merged$id))
cat("[09b_robustness] Personen mit korrigierten Menses-Daten:",
    paste(sort(corrected_ids), collapse = " "), "\n")

# participant x phase means
pp_of <- function(var, data = dd) {
  data |> filter(!is.na(cycle_phase_2)) |>
    group_by(Id, cycle_phase_2) |>
    summarise(y = mean(.data[[var]], na.rm = TRUE),
              n_days = sum(!is.na(.data[[var]])), .groups = "drop") |>
    filter(is.finite(y))
}

# subsample by data availability (as hc_ids in 08)
avail_ids <- function(var, strict = FALSE) {
  a <- dd |> filter(!is.na(cycle_phase_2)) |>
    group_by(Id, cycle_phase_2) |>
    summarise(pct = mean(!is.na(.data[[var]])) * 100, .groups = "drop") |>
    group_by(Id) |> summarise(m = min(pct), k = n(), .groups = "drop") |>
    filter(m >= COMPL_THR)
  if (strict) a <- filter(a, k == 2)
  a$Id
}

# LMM on participant x phase means
fit_pp <- function(pp) {
  if (is.null(pp) || nrow(pp) == 0 || n_distinct(pp$Id) < 5 ||
      n_distinct(pp$cycle_phase_2) < 2) return(tibble())
  m <- suppressMessages(suppressWarnings(
    lmerTest::lmer(y ~ cycle_phase_2 + (1 | Id), data = pp)))
  s <- summary(m)$coefficients["cycle_phase_2luteal", ]
  tibble(n_obs = nrow(pp), n_id = n_distinct(pp$Id),
         b = unname(s["Estimate"]), t = unname(s["t value"]),
         df = unname(s["df"]), p = unname(s["Pr(>|t|)"]),
         note = if (lme4::isSingular(m)) "singular fit" else "")
}

# day-level LMM
fit_day <- function(var, adjust = FALSE, data = dd) {
  d <- data |> filter(!is.na(cycle_phase_2), !is.na(.data[[var]]))
  f <- paste(var, "~ cycle_phase_2", if (adjust) "+ dl_c", "+ (1 | Id)")
  m <- suppressMessages(suppressWarnings(lmerTest::lmer(as.formula(f), data = d)))
  s <- summary(m)$coefficients["cycle_phase_2luteal", ]
  tibble(n_obs = nrow(d), n_id = n_distinct(d$Id),
         b = unname(s["Estimate"]), t = unname(s["t value"]),
         df = unname(s["df"]), p = unname(s["Pr(>|t|)"]),
         note = if (adjust) "adjustiert fuer Tageslaenge" else "Tageswerte, nicht aggregiert")
}

# paired Wilcoxon test
wil_row <- function(pp) {
  w <- pp |> select(Id, cycle_phase_2, y) |>
    tidyr::pivot_wider(names_from = cycle_phase_2, values_from = y) |>
    tidyr::drop_na(follicular, luteal)
  if (nrow(w) < 5) return(tibble())
  tibble(n_obs = 2 * nrow(w), n_id = nrow(w),
         b = median(w$luteal - w$follicular),
         p = suppressWarnings(wilcox.test(w$luteal, w$follicular,
                                          paired = TRUE, exact = FALSE)$p.value),
         note = "b = Median der Differenz (luteal - follikulaer)")
}

# GAMM as in 09
gamm_fit <- function(outcome, data = dd, k = 20, extra = "") {
  d <- data |> filter(!is.na(pacts_ov), !is.na(.data[[outcome]])) |>
    mutate(Id = factor(Id))
  if (nrow(d) < 50 || n_distinct(d$pacts_ov) < k) return(tibble())
  m <- suppressMessages(suppressWarnings(
    bam(as.formula(paste0(outcome, " ~ s(pacts_ov, bs='cc', k=", k, ")", extra,
                          " + s(Id, bs='re')")),
        data = d, method = "fREML")))
  st <- summary(m)$s.table
  tibble(n_obs = nrow(d), n_id = n_distinct(d$Id),
         edf = unname(st["s(pacts_ov)", "edf"]),
         p   = unname(st["s(pacts_ov)", "p-value"]))
}

loo_note <- function(l, stat, left) {
  bad <- l$p >= .05
  sprintf("%d Fits; %s %.3f bis %.3f; p max %.4f; p<.05 in %d/%d%s | b/EDF = Median, p = Maximum",
          nrow(l), stat, min(l[[stat]]), max(l[[stat]]), max(l$p),
          sum(!bad), nrow(l),
          if (any(bad)) paste0("; p>=.05 ohne: ", paste(left[bad], collapse = " ")) else "")
}

res <- list()
add_row <- function(finding, test, x) {
  res[[length(res) + 1]] <<-
    if (is.null(x) || nrow(x) == 0) tibble(finding = finding, test = test, note = "nicht berechenbar")
    else bind_cols(tibble(finding = finding, test = test), x)
}

# ---- (1) Cycle Symptom Index trajectory ----
F1 <- "1 Symptomverlauf (CSI, GAMM)"
add_row(F1, "Gesamtstichprobe [erwartet: EDF 2.33, p = .030]", gamm_fit("cycle_symptom_idx"))
for (kk in c(6, 10))
  add_row(F1, paste0("Spline-Basis k = ", kk), gamm_fit("cycle_symptom_idx", k = kk))
add_row(F1, "adjustiert fuer Tageslaenge", gamm_fit("cycle_symptom_idx", extra = " + dl_c"))
add_row(F1, "ohne Personen mit korrigierten Menses-Daten",
        gamm_fit("cycle_symptom_idx", data = filter(dd, !Id %in% corrected_ids)))
add_row(F1, paste0("nur Personen mit >=", COMPL_THR, "% Fragebogen-Tage je Phase"),
        gamm_fit("cycle_symptom_idx", data = filter(dd, Id %in% avail_ids("cycle_symptom_idx"))))

somatic <- c("drsp_breast_tender", "drsp_bloated", "drsp_food_cravings",
             "drsp_appetite", "drsp_fatigue")
dd_red <- dd |> mutate(csi_ohne_somatisch =
            rowMeans(pick(all_of(setdiff(DRSP_ITEMS, somatic))), na.rm = TRUE))
add_row(F1, "Index OHNE die 5 einzeln signifikanten Items (explorativ)",
        gamm_fit("csi_ohne_somatisch", data = dd_red))

ids_g <- unique(dd$Id[!is.na(dd$cycle_symptom_idx) & !is.na(dd$pacts_ov)])
l <- map_dfr(ids_g, \(i) gamm_fit("cycle_symptom_idx", data = filter(dd, Id != i)) |>
               mutate(left_out = i))
add_row(F1, "Leave-one-participant-out (GAMM)",
        tibble(n_id = n_distinct(l$left_out), edf = median(l$edf), p = max(l$p),
               note = loo_note(l, "edf", l$left_out)))

# ---- (2) Centroid of light exposure ----
F2 <- "2 Centroid Lichtexposition (LMM)"
v  <- "centroid_LE_hr"
pp <- pp_of(v)
add_row(F2, "Gesamtstichprobe [erwartet: b = -0.70, t(13.4) = -3.79, p = .0022]", fit_pp(pp))
l <- map_dfr(unique(pp$Id), \(i) fit_pp(filter(pp, Id != i)) |> mutate(left_out = i))
add_row(F2, "Leave-one-participant-out (LMM)",
        tibble(n_id = n_distinct(l$left_out), b = median(l$b), p = max(l$p),
               note = loo_note(l, "b", l$left_out)))
add_row(F2, "nur Person x Phase mit >=5 Tagen", fit_pp(filter(pp, n_days >= 5)))
add_row(F2, paste0("ActLumus-Daten >=", COMPL_THR, "% der Tage je Phase mit Daten"),
        fit_pp(filter(pp, Id %in% avail_ids("TAT_250_hr"))))
add_row(F2, "wie oben, strikt (beide Phasen)",
        fit_pp(filter(pp, Id %in% avail_ids("TAT_250_hr", strict = TRUE))))
add_row(F2, "ohne Personen mit korrigierten Menses-Daten",
        fit_pp(filter(pp, !Id %in% corrected_ids)))
add_row(F2, "Tagesebene, ohne Adjustierung", fit_day(v))
add_row(F2, "Tagesebene + Tageslaenge", fit_day(v, adjust = TRUE))
add_row(F2, "gepaarter Wilcoxon", wil_row(pp))
add_row(F2, "GAMM Gesamtstichprobe [erwartet: EDF 1.80, p = .022]", gamm_fit(v))
add_row(F2, "GAMM + Tageslaenge", gamm_fit(v, extra = " + dl_c"))

# ---- (3) Relative amplitude ----
F3 <- "3 Relative Amplitude (RA_act, LMM)"
v  <- "RA_act"
pp <- pp_of(v)
add_row(F3, "Gesamtstichprobe [erwartet: b = 0.08, t(19.3) = 2.59, p = .018; N 36/20]", fit_pp(pp))
l <- map_dfr(unique(pp$Id), \(i) fit_pp(filter(pp, Id != i)) |> mutate(left_out = i))
add_row(F3, "Leave-one-participant-out (LMM)",
        tibble(n_id = n_distinct(l$left_out), b = median(l$b), p = max(l$p),
               note = loo_note(l, "b", l$left_out)))
add_row(F3, "nur Person x Phase mit >=5 Tagen", fit_pp(filter(pp, n_days >= 5)))
add_row(F3, paste0("Actigraph-Daten >=", COMPL_THR, "% der Tage je Phase mit Daten"),
        fit_pp(filter(pp, Id %in% avail_ids("RA_act"))))
add_row(F3, "wie oben, strikt (beide Phasen)",
        fit_pp(filter(pp, Id %in% avail_ids("RA_act", strict = TRUE))))
add_row(F3, "ohne Personen mit korrigierten Menses-Daten",
        fit_pp(filter(pp, !Id %in% corrected_ids)))
add_row(F3, "Tagesebene, ohne Adjustierung", fit_day(v))
add_row(F3, "Tagesebene + Tageslaenge", fit_day(v, adjust = TRUE))
add_row(F3, "gepaarter Wilcoxon", wil_row(pp))
add_row(F3, "GAMM Gesamtstichprobe [erwartet: EDF 1.57, p = .053]", gamm_fit(v))

# ---- Output ----
robustness_table <- bind_rows(res) |>
  mutate(across(any_of(c("b", "t", "df", "edf")), ~ round(.x, 3)),
         p = signif(p, 3))
cat("\n[09b_robustness] Robustheitstests:\n")
print(robustness_table |> select(any_of(c("finding", "test", "n_obs", "n_id", "b",
                                          "edf", "t", "df", "p", "note"))),
      n = Inf, width = Inf)
rb_file <- safe_write_csv(robustness_table, file.path(OUTPUT_DIR, "robustness_table.csv"))
cat("[09b_robustness] gespeichert:", rb_file, "\n")
