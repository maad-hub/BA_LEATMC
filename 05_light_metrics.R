# 05_light_metrics.R - Daily light-exposure metrics (LightLogR)

EMPTY_LIGHT <- tibble(Id=character(), Date=as.Date(character()),
  TAT_250_hr=numeric(), TAT_1000_hr=numeric(), TAT_10_hr=numeric(),
  MEDI_mean=numeric(), MEDI_median=numeric(), MEDI_max=numeric(), MEDI_sd=numeric(),
  MEDI_morning=numeric(), MEDI_evening=numeric(), MEDI_daytime=numeric(), MEDI_night=numeric(),
  first_250=numeric(), last_250=numeric(),
  M10_mean=numeric(), M10_onset=numeric(), M10_midpoint=numeric(),
  L5_mean=numeric(), L5_onset=numeric(), RA_light=numeric(),
  centroid_LE_hr=numeric(), brown_compliance_pct=numeric(), n_epochs=integer())

# ---- Helper: bright_dark_period() with an NA row on error ----
safe_bright_dark <- function(Light.vector, Time.vector, period, timespan, loop = FALSE) {
  prefix <- paste0(period, "_", gsub(" ", "", timespan))
  tryCatch(
    bright_dark_period(Light.vector = Light.vector, Time.vector = Time.vector,
                        period = period, timespan = timespan, loop = loop, as.df = TRUE),
    error = function(e) {
      tibble(
        !!paste0(prefix, "_mean")     := NA_real_,
        !!paste0(prefix, "_onset")    := as.POSIXct(NA),
        !!paste0(prefix, "_midpoint") := as.POSIXct(NA),
        !!paste0(prefix, "_offset")   := as.POSIXct(NA)
      )
    }
  )
}

# ---- Helper: timing_above_threshold() ----
safe_timing <- function(Light.vector, Time.vector, threshold) {
  tryCatch(
    timing_above_threshold(Light.vector = Light.vector, Time.vector = Time.vector,
                            threshold = threshold, as.df = TRUE),
    error = function(e) {
      tibble(
        !!paste0("first_timing_above_", threshold) := as.POSIXct(NA),
        !!paste0("last_timing_above_", threshold)  := as.POSIXct(NA)
      )
    }
  )
}

# decimal hours
to_hr <- function(x) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(NA_real_)
  hour(x) + minute(x)/60 + second(x)/3600
}

if (nrow(light_trimmed) == 0) {
  warning("Keine ActLumus-Daten — LightLogR-Metriken werden NA")
  light_daily <- EMPTY_LIGHT
} else {

  cat("[05_light_metrics] Starte Berechnung fuer", n_distinct(light_trimmed$Id),
      "Teilnehmerinnen x", n_distinct(light_trimmed$Date), "Tage",
      "— das kann je nach Datenmenge mehrere Minuten dauern...\n")

  t0 <- Sys.time()

  light_daily_raw <- light_trimmed |>
    ungroup() |>
    group_by(Id, Date) |>
    reframe(
      # mEDI summary statistics
      MEDI_mean    = mean(MEDI, na.rm=TRUE),
      MEDI_median  = median(MEDI, na.rm=TRUE),
      MEDI_max     = max(MEDI, na.rm=TRUE),
      MEDI_sd      = sd(MEDI, na.rm=TRUE),
      MEDI_morning = mean(MEDI[hour(Datetime) %in% 6:9],   na.rm=TRUE),  # phase-advancing window
      MEDI_evening = mean(MEDI[hour(Datetime) %in% 19:22], na.rm=TRUE),  # phase-delaying window
      MEDI_daytime = mean(MEDI[hour(Datetime) %in% 8:19],  na.rm=TRUE),
      MEDI_night   = mean(MEDI[hour(Datetime) %in% 0:4],   na.rm=TRUE),  # QC: should be close to 0
      n_epochs     = n(),

      # time above 10, 250 and 1000 lx mEDI
      TAT_250_hr  = tryCatch(as.numeric(duration_above_threshold(
                     MEDI, Datetime, threshold=250,  as.df=FALSE))/3600, error=function(e) NA_real_),
      TAT_1000_hr = tryCatch(as.numeric(duration_above_threshold(
                     MEDI, Datetime, threshold=1000, as.df=FALSE))/3600, error=function(e) NA_real_),
      TAT_10_hr   = tryCatch(as.numeric(duration_above_threshold(
                     MEDI, Datetime, threshold=10,   as.df=FALSE))/3600, error=function(e) NA_real_),

      # first and last time above 250 lx
      safe_timing(MEDI, Datetime, threshold = 250),

      # brightest 10 h and darkest 5 h
      safe_bright_dark(MEDI, Datetime, "brightest", "10 hours"),
      safe_bright_dark(MEDI, Datetime, "darkest", "5 hours", loop = TRUE),

      # centroid of light exposure
      centroid_LE_hr = tryCatch({
        cle <- centroidLE(MEDI, Datetime)
        if (!is.null(cle) && !is.na(cle)) hour(cle)+minute(cle)/60+second(cle)/3600 else NA_real_
      }, error=function(e) NA_real_),

      # Brown et al. (2022) recommendations
      brown_compliance_pct = tryCatch({
        db <- tibble(MEDI=MEDI, Datetime=Datetime) |> Brown_cut() |> Brown2reference()
        mean(db$Reference.check, na.rm=TRUE) * 100
      }, error = function(e) NA_real_),

      .groups = "drop"
    )

  cat("[05_light_metrics] Rohberechnung fertig nach",
      round(difftime(Sys.time(), t0, units="secs"), 1), "Sekunden\n")

  # ---- Derived values ----
  find_col <- function(df, ...) {
    pat <- paste(c(...), collapse=".*")
    hit <- grep(pat, names(df), value = TRUE, ignore.case = TRUE)
    if (length(hit) == 0) return(rep(NA, nrow(df)))
    df[[hit[1]]]
  }

  light_daily <- light_daily_raw |>
    mutate(
      first_250    = to_hr(find_col(light_daily_raw, "first", "timing")),
      last_250     = to_hr(find_col(light_daily_raw, "last",  "timing")),
      M10_mean_raw = find_col(light_daily_raw, "brightest", "mean"),
      M10_onset    = to_hr(find_col(light_daily_raw, "brightest", "onset")),
      M10_midpoint = to_hr(find_col(light_daily_raw, "brightest", "midpoint")),
      L5_mean_raw  = find_col(light_daily_raw, "darkest", "mean"),
      L5_onset     = to_hr(find_col(light_daily_raw, "darkest", "onset")),
      M10_mean     = as.numeric(M10_mean_raw),
      L5_mean      = as.numeric(L5_mean_raw),
      RA_light     = if_else(!is.na(M10_mean) & !is.na(L5_mean) & (M10_mean+L5_mean) > 0,
                              (M10_mean - L5_mean) / (M10_mean + L5_mean), NA_real_)
    ) |>
    select(Id, Date, TAT_250_hr, TAT_1000_hr, TAT_10_hr,
           MEDI_mean, MEDI_median, MEDI_max, MEDI_sd,
           MEDI_morning, MEDI_evening, MEDI_daytime, MEDI_night,
           first_250, last_250,
           M10_mean, M10_onset, M10_midpoint, L5_mean, L5_onset, RA_light,
           centroid_LE_hr, brown_compliance_pct, n_epochs)
}

cat("[05_light_metrics] light_daily:", nrow(light_daily), "Zeilen |",
    n_distinct(light_daily$Id), "Teilnehmerinnen\n")

daily_data <- daily_data_raw |>
  full_join(light_daily, by=c("Id","Date"))
