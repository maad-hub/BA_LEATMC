# 01_import_actlumus.R - Import of the ActLumus light-logger data

files_al <- list.files(BASE_DIR, pattern = "\\.txt$|\\.zip$",
                        recursive = TRUE, full.names = TRUE)
files_al <- files_al[grepl("P\\d{3}[/\\\\]", files_al)]

if (length(files_al) == 0) {
  warning("Keine ActLumus-Dateien gefunden — Lichtmetriken werden NA")
  light_raw <- tibble(Id = character(), Datetime = as.POSIXct(character()),
                      MEDI = numeric(), Date = as.Date(character()))
} else {
  light_raw <- map(files_al, \(f) safe({
    pid <- basename(dirname(f))
    # standard ActLumus import
    result <- tryCatch(
      import$ActLumus(filename = f, tz = TIMEZONE,
                      manual.id = pid,
                      remove_duplicates = TRUE,
                      auto.plot = FALSE, silent = TRUE),
      error = function(e) NULL
    )
    # fallback: ActLogModel 2.0.0 export (_Report.txt, semicolon-separated)
    if (is.null(result)) {
      lines_raw <- readLines(f, encoding = "latin1")
      header_line <- which(grepl("^DATE/TIME", lines_raw))[1]
      if (!is.na(header_line)) {
        df <- read.csv2(f, skip = header_line - 1, sep = ";",
                        header = TRUE, check.names = FALSE,
                        fileEncoding = "latin1") |>
          as_tibble() |>
          rename(Datetime = `DATE/TIME`,
                 MEDI     = `MELANOPIC EDI`) |>
          mutate(
            Datetime = as.POSIXct(Datetime, format = "%d/%m/%Y %H:%M:%S",
                                   tz = TIMEZONE),
            MEDI     = as.numeric(MEDI),
            Id       = pid
          ) |>
          select(Id, Datetime, MEDI) |>
          filter(!is.na(Datetime), !is.na(MEDI))
        result <- df
      }
    }
    result
  }, label = basename(f))) |> bind_rows() |>
    # P012: remove epochs recorded before the study start
    filter(!(Id == "P012" &
               Datetime < as.POSIXct("2026-06-09 00:00:00", tz = TIMEZONE))) |>
    mutate(Date = as_date(Datetime),
           Id   = standardize_id(Id))

  cat("[01_import_actlumus] ActLumus:", n_distinct(light_raw$Id), "Teilnehmerinnen |",
      nrow(light_raw), "Epochen\n")
  validate_id_format(light_raw$Id, "ActLumus folder-based ID")
}
