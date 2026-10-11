# LEAMC_GGIR_sleep.R - Sleep detection from ActiGraph raw data (.gt3x) with GGIR;
# part 5 additionally exports the activity time series (ENMO) used for IS/IV

suppressPackageStartupMessages({
  library(GGIR)
  library(tidyverse)
  library(lubridate)
})

# ---- Paths ----
if (!exists("BASE_DIR")) BASE_DIR <- "C:/Users/astri/LRZ Sync+Share/LEATMC"  # adjust if needed
if (!exists("TIMEZONE")) TIMEZONE <- "Europe/Berlin"
if (!exists("validate_id_format")) validate_id_format <- function(ids, source_label) {
  ids <- unique(ids[!is.na(ids)])
  bad <- ids[!grepl("^P\\d{3}$", ids)]
  if (length(bad) > 0) {
    warning("[", source_label, "] ID(s) entsprechen nicht dem Format P0xx: ",
            paste(bad, collapse=", "))
  }
  invisible(bad)
}

# .gt3x files are stored in the participant folders (P0xx) below BASE_DIR
RAW_DIR <- BASE_DIR

# GGIR output is written to a local (not synchronised) folder
GGIR_OUT <- "C:/GGIR_local_output/LEATMC"  # adjust if needed

# ---- Memory settings: no parallel processing and small batches on machines with little RAM ----
GGIR_PARALLEL <- if (!exists("GGIR_PARALLEL")) TRUE else GGIR_PARALLEL

GGIR_BATCH_SIZE <- if (!exists("GGIR_BATCH_SIZE")) 3 else GGIR_BATCH_SIZE

dir.create(GGIR_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Step 1: run GGIR in batches ----

n_files_total <- length(list.files(RAW_DIR, pattern = "\\.gt3x$", recursive = TRUE))
cat("[GGIR] ", n_files_total, " Rohdateien gefunden, werden in Batches von ",
    GGIR_BATCH_SIZE, " Datei(en) verarbeitet.\n", sep = "")

batch_starts <- seq(1, n_files_total, by = GGIR_BATCH_SIZE)
n_batches    <- length(batch_starts)
batch_times  <- numeric(0)

for (i in seq_along(batch_starts)) {
  start_idx <- batch_starts[i]
  end_idx   <- min(start_idx + GGIR_BATCH_SIZE - 1, n_files_total)
  cat("\n[GGIR] ===== Batch ", i, "/", n_batches, ": Datei ", start_idx, " bis ", end_idx,
      " von ", n_files_total, " =====\n", sep = "")
  batch_t0 <- Sys.time()

  safe(GGIR(
    datadir    = RAW_DIR,
    outputdir  = GGIR_OUT,
    studyname  = "LEATMC",
    # parts 1-5 (part 5: activity time series for IS/IV)
    mode       = c(1, 2, 3, 4, 5),
    overwrite  = FALSE,  # skip files that were already processed
    do.parallel = GGIR_PARALLEL,
    f0 = start_idx,
    f1 = end_idx,

    # general
    desiredtz       = TIMEZONE,
    sensor.location = "wrist",

    # sleep detection (see Methods)
    HASPT.algo = "HDCZA",  # sleep period time window
    HASIB.algo = "vanHees2015",  # sleep/wake within the sleep period

    # participant ID from the folder name
    storefolderstructure = TRUE,

    # part 5: export of the activity time series (aggregated to 1 min in 02b_import_activity_epochs.R)
    windowsizes         = c(60, 900, 3600),
    acc.metric           = "ENMO",
    save_ms5rawlevels    = TRUE,
    save_ms5raw_format   = "csv",

    # reports
    do.report    = c(4),
    visualreport = FALSE
  ), label = paste("GGIR Batch", start_idx, "-", end_idx))

  batch_elapsed <- as.numeric(difftime(Sys.time(), batch_t0, units = "mins"))
  batch_times   <- c(batch_times, batch_elapsed)
  avg_time      <- mean(batch_times)
  batches_left  <- n_batches - i
  eta_min       <- avg_time * batches_left

  cat("[GGIR] Batch ", i, "/", n_batches, " fertig nach ", round(batch_elapsed, 1),
      " Minuten (Durchschnitt bisher: ", round(avg_time, 1), " Min/Batch).\n", sep = "")
  if (batches_left > 0) {
    cat("[GGIR] Geschaetzte Restzeit: ~", round(eta_min, 0), " Minuten (~",
        round(eta_min / 60, 1), " Std.) fuer die verbleibenden ", batches_left,
        " Batch(es). Diese Schaetzung wird nach jedem Batch genauer.\n", sep = "")
  }
}

cat("\n[GGIR] Alle Batches abgearbeitet. Ergebnisse in:", GGIR_OUT, "\n")

# ---- Step 2: read part4_nightsummary_sleep_full.csv ----

night_files <- list.files(GGIR_OUT, pattern = "nightsummary_sleep_full\\.csv$",
                          recursive = TRUE, full.names = TRUE)

if (length(night_files) == 0) {
  stop("Keine 'part4_nightsummary_sleep_full.csv' gefunden unter ", GGIR_OUT,
       ". Pruef manuell mit list.files('", GGIR_OUT, "', recursive=TRUE), ",
       "ob GGIR Part 4 ueberhaupt durchlaufen hat (Konsolen-Output oben pruefen).")
}
cat("[GGIR] Gefundene Night-Summary-Datei(en):\n"); print(night_files)

nightsummary <- read_csv(night_files[1], show_col_types = FALSE)
cat("[GGIR] Roh-Zeilen (Naechte x Teilnehmerinnen):", nrow(nightsummary), "\n")
cat("[GGIR] Spalten:", paste(names(nightsummary), collapse=", "), "\n")

# ---- Step 3: participant ID from the file name ----

filename_col <- grep("filename", names(nightsummary), ignore.case = TRUE, value = TRUE)[1]
if (is.na(filename_col)) {
  stop("Keine 'filename'-Spalte in der Night-Summary gefunden. Verfuegbare Spalten: ",
       paste(names(nightsummary), collapse = ", "))
}
cat("[GGIR] Verwende Spalte '", filename_col, "' fuer die Teilnehmerinnen-ID.\n", sep = "")

nightsummary <- nightsummary |>
  mutate(Id = str_extract(.data[[filename_col]], "P\\d{3}"))

if (any(is.na(nightsummary$Id))) {
  bad_examples <- head(unique(nightsummary[[filename_col]][is.na(nightsummary$Id)]))
  warning(sum(is.na(nightsummary$Id)), " Zeile(n) ohne erkennbares P0xx-Muster im Dateinamen. ",
          "Beispielwerte: ", paste(bad_examples, collapse = ", "))
}

cat("[GGIR] Erkannte Teilnehmerinnen-IDs:", paste(sort(unique(nightsummary$Id)), collapse=" "), "\n")
validate_id_format(nightsummary$Id, "GGIR filename-based ID")

# ---- Step 4: check units and convert ----

detect_and_convert_to_minutes <- function(x, label) {
  med <- median(x, na.rm = TRUE)
  if (is.na(med)) return(x)
  if (med <= 16) {
    cat("  [Einheit] ", label, ": Median = ", round(med,2),
        " -> als STUNDEN erkannt, x60 zu Minuten konvertiert\n", sep = "")
    return(x * 60)
  } else {
    cat("  [Einheit] ", label, ": Median = ", round(med,2),
        " -> bereits als MINUTEN erkannt, unveraendert\n", sep = "")
    return(x)
  }
}

cat("[GGIR] Pruefe Einheiten der Schlafmetriken:\n")
sleep_data <- nightsummary |>
  mutate(
    TST_min  = detect_and_convert_to_minutes(SleepDurationInSpt, "SleepDurationInSpt (TST)"),
    WASO_min = detect_and_convert_to_minutes(WASO, "WASO"),
    sleep_onset_hr  = sleeponset,
    sleep_offset_hr = wakeup,
    # sleep efficiency = sleep duration / sleep period duration
    sleep_efficiency = if_else(SptDuration > 0,
                               (SleepDurationInSpt / SptDuration) * 100,
                               NA_real_),
    n_awakenings = number_of_awakenings,
    # calendar_date format differs between GGIR versions
    Date = as.Date(suppressWarnings(
      lubridate::parse_date_time(calendar_date, orders = c("ymd", "dmy", "mdy"))
    ))
  )

# check date parsing
date_parse_failed <- is.na(sleep_data$Date) & !is.na(sleep_data$calendar_date)
if (any(date_parse_failed)) {
  warning(sum(date_parse_failed), " Datum(s) konnten in keinem der probierten Formate ",
          "(ymd/dmy/mdy) geparst werden. Beispiel-Rohwerte: ",
          paste(head(sleep_data$calendar_date[date_parse_failed]), collapse=", "))
}

actigraph_sleep <- sleep_data |>
  filter(!is.na(Id), !is.na(Date)) |>
  select(Id, Date, sleep_onset_hr, sleep_offset_hr, TST_min, WASO_min,
         sleep_efficiency, n_awakenings) |>
  distinct(Id, Date, .keep_all = TRUE) |>
  arrange(Id, Date)

# ---- Step 5: quality control and save ----

cat("\n[GGIR] ===== Zusammenfassung =====\n")
if (nrow(actigraph_sleep) == 0) {
  warning("actigraph_sleep hat 0 Zeilen nach dem Filtern auf gueltige Id/Date. ",
          "Pruef die Warnungen weiter oben (P0xx-Muster / Datumsparsing).")
} else {
  cat("Naechte gesamt:  ", nrow(actigraph_sleep), "\n")
  cat("Teilnehmerinnen: ", n_distinct(actigraph_sleep$Id), "\n")
  cat("Datumsbereich:   ", as.character(min(actigraph_sleep$Date)), "bis",
      as.character(max(actigraph_sleep$Date)), "\n")
  cat("TST_min  - Median:", round(median(actigraph_sleep$TST_min, na.rm=TRUE), 1),
      "| Bereich:", round(min(actigraph_sleep$TST_min, na.rm=TRUE), 1), "-",
      round(max(actigraph_sleep$TST_min, na.rm=TRUE), 1), "\n")
  cat("(Realistisch waeren TST_min ca. 300-600 [5-10h]. Falls die Zahlen ganz\n",
      " anders aussehen, nochmal die Einheiten-Erkennung in Schritt 4 pruefen.)\n", sep = "")
}

saveRDS(actigraph_sleep, file.path(BASE_DIR, "actigraph_sleep_ggir.rds"))
cat("\n[GGIR] Gespeichert:", file.path(BASE_DIR, "actigraph_sleep_ggir.rds"), "\n")
cat("Jetzt kannst du 02_import_actigraph.R (bzw. main.R) erneut laufen lassen.\n")
