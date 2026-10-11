# 02b_import_activity_epochs.R - Import of the GGIR part-5 activity time series (ENMO) for IS/IV

# locate ms5.outraw below GGIR_OUT
ms5raw_candidates <- list.files(GGIR_OUT, pattern = "^ms5\\.outraw$",
                                 recursive = TRUE, full.names = TRUE, include.dirs = TRUE)
ms5raw_dir <- if (length(ms5raw_candidates) > 0) ms5raw_candidates[1] else
  file.path(GGIR_OUT, "output_LEATMC", "meta", "ms5.outraw")

if (!dir.exists(ms5raw_dir)) {
  meta_dirs <- list.files(GGIR_OUT, pattern = "^meta$", recursive = TRUE,
                           full.names = TRUE, include.dirs = TRUE)
  found_str <- if (length(meta_dirs) > 0) {
    paste(unlist(lapply(meta_dirs, list.files)), collapse = ", ")
  } else "kein 'meta'-Ordner unter GGIR_OUT gefunden"
  warning("Kein ms5.outraw-Ordner gefunden unter '", ms5raw_dir, "' - wurde ",
          "LEAMC_GGIR_sleep.R mit Part 5 (mode enthaelt 5) und ",
          "save_ms5rawlevels=TRUE neu durchlaufen? activity_epochs bleibt leer, ",
          "IS/IV werden dadurch spaeter NA statt eines falschen Licht-basierten Werts. ",
          "Tatsaechlich vorhanden in meta/: ", found_str)
  activity_epochs <- tibble(Id = character(), Datetime = as.POSIXct(character()),
                             ACC = numeric(), Date = as.Date(character()))
} else {
  files_act <- list.files(ms5raw_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  cat("[02b_import_activity_epochs]", length(files_act), "Part-5-Zeitreihen-Dateien gefunden.\n")

  if (length(files_act) == 0) {
    warning("ms5.outraw-Ordner existiert, enthaelt aber keine .csv-Dateien. ",
            "save_ms5raw_format wurde vermutlich nicht auf 'csv' gesetzt.")
    activity_epochs <- tibble(Id = character(), Datetime = as.POSIXct(character()),
                               ACC = numeric(), Date = as.Date(character()))
  } else {
    read_one_activity_file <- function(path) {
      id <- str_extract(basename(path), "P\\d{3}")
      if (is.na(id)) {
        warning("Keine P0xx-ID im Dateinamen erkannt: ", basename(path), " - uebersprungen.")
        return(NULL)
      }
      d <- safe(read_csv(path, show_col_types = FALSE), paste("Activity-Epochen", basename(path)))
      if (is.null(d) || nrow(d) == 0) return(NULL)

      # time stamp: 'timestamp' or 'timenum' (seconds since 1970-01-01 UTC)
      if ("timestamp" %in% names(d)) {
        dt <- as.POSIXct(d$timestamp, tz = TIMEZONE)
      } else if ("timenum" %in% names(d)) {
        dt <- as.POSIXct(d$timenum, origin = "1970-01-01", tz = "UTC")
        attr(dt, "tzone") <- TIMEZONE
      } else {
        warning("Weder 'timestamp' noch 'timenum' in ", basename(path), " gefunden - uebersprungen.")
        return(NULL)
      }
      if (!"ACC" %in% names(d)) {
        warning("Keine 'ACC'-Spalte in ", basename(path),
                " - verfuegbare Spalten: ", paste(names(d), collapse=", "))
        return(NULL)
      }

      # remove epochs flagged as invalid (non-wear) by GGIR
      if ("invalidepoch" %in% names(d)) {
        keep <- d$invalidepoch == 0
        d  <- d[keep, , drop = FALSE]
        dt <- dt[keep]
      }

      tibble(Id = id, Datetime = dt, ACC = d$ACC) |>
        filter(!is.na(Datetime), !is.na(ACC))
    }

    activity_list <- lapply(files_act, read_one_activity_file)
    activity_epochs_raw <- bind_rows(activity_list) |>
      mutate(Id = standardize_id(Id))

    if (nrow(activity_epochs_raw) > 0) validate_id_format(activity_epochs_raw$Id, "GGIR Part 5 filename-based ID")

    # aggregate to 1-min epochs
    activity_epochs <- activity_epochs_raw |>
      mutate(Datetime = as.POSIXct(floor(as.numeric(Datetime) / 60) * 60,
                                    origin = "1970-01-01", tz = TIMEZONE)) |>
      group_by(Id, Datetime) |>
      summarise(ACC = mean(ACC, na.rm = TRUE), .groups = "drop") |>
      mutate(Date = as.Date(Datetime, tz = TIMEZONE))

    cat("[02b_import_activity_epochs] Aktivitaets-Epochen:",
        n_distinct(activity_epochs$Id), "Teilnehmerinnen |",
        nrow(activity_epochs_raw), "Roh-Epochen ->", nrow(activity_epochs),
        "1-Minuten-Epochen (auf 60s aggregiert)\n")
  }
}
