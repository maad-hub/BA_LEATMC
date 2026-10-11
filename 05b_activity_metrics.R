# 05b_activity_metrics.R - Daily rest-activity metrics from ENMO (M10, L5, relative amplitude)

EMPTY_ACTIVITY <- tibble(Id=character(), Date=as.Date(character()),
  M10_act_mean=numeric(), M10_act_onset=numeric(), M10_act_midpoint=numeric(),
  L5_act_mean=numeric(), L5_act_onset=numeric(), RA_act=numeric(),
  n_epochs_act=integer())

if (nrow(activity_trimmed) == 0) {
  warning("Keine Aktivitaets-Epochen (activity_trimmed leer) — M10_act/L5_act/RA_act werden NA. ",
          "Wurde LEAMC_GGIR_sleep.R mit Part 5 neu durchlaufen? Siehe 02b_import_activity_epochs.R.")
  activity_daily <- EMPTY_ACTIVITY
} else {

  cat("[05b_activity_metrics] Starte Berechnung fuer", n_distinct(activity_trimmed$Id),
      "Teilnehmerinnen x", n_distinct(activity_trimmed$Date), "Tage...\n")
  t0 <- Sys.time()

  activity_daily_raw <- activity_trimmed |>
    ungroup() |>
    group_by(Id, Date) |>
    reframe(
      n_epochs_act = n(),
      safe_bright_dark(ACC, Datetime, "brightest", "10 hours"),
      safe_bright_dark(ACC, Datetime, "darkest", "5 hours", loop = TRUE),
      .groups = "drop"
    )

  cat("[05b_activity_metrics] Rohberechnung fertig nach",
      round(difftime(Sys.time(), t0, units="secs"), 1), "Sekunden\n")

  activity_daily <- activity_daily_raw |>
    mutate(
      M10_act_mean_raw = find_col(activity_daily_raw, "brightest", "mean"),
      M10_act_onset     = to_hr(find_col(activity_daily_raw, "brightest", "onset")),
      M10_act_midpoint  = to_hr(find_col(activity_daily_raw, "brightest", "midpoint")),
      L5_act_mean_raw   = find_col(activity_daily_raw, "darkest", "mean"),
      L5_act_onset      = to_hr(find_col(activity_daily_raw, "darkest", "onset")),
      M10_act_mean      = as.numeric(M10_act_mean_raw),
      L5_act_mean       = as.numeric(L5_act_mean_raw),
      RA_act            = if_else(!is.na(M10_act_mean) & !is.na(L5_act_mean) & (M10_act_mean+L5_act_mean) > 0,
                                   (M10_act_mean - L5_act_mean) / (M10_act_mean + L5_act_mean), NA_real_)
    ) |>
    select(Id, Date, M10_act_mean, M10_act_onset, M10_act_midpoint,
           L5_act_mean, L5_act_onset, RA_act, n_epochs_act)
}

cat("[05b_activity_metrics] activity_daily:", nrow(activity_daily), "Zeilen |",
    n_distinct(activity_daily$Id), "Teilnehmerinnen\n")

daily_data <- daily_data |>
  full_join(activity_daily, by=c("Id","Date"))
