# 02_import_actigraph.R - Load the GGIR sleep results (from LEAMC_GGIR_sleep.R)

ggir_rds <- file.path(BASE_DIR, "actigraph_sleep_ggir.rds")

# read with retries (the file may briefly be locked by the sync client)
read_rds_with_retry <- function(path, max_tries = 5, wait_seconds = 3) {
  for (attempt in seq_len(max_tries)) {
    result <- tryCatch(readRDS(path), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    if (attempt < max_tries) {
      message("  [Retry ", attempt, "/", max_tries, "] Konnte '", path,
              "' noch nicht lesen (", conditionMessage(result),
              ") - vermutlich kurzzeitig durch Sync-Client gesperrt. ",
              "Warte ", wait_seconds, "s und versuche erneut...")
      Sys.sleep(wait_seconds)
    } else {
      stop("Konnte '", path, "' nach ", max_tries, " Versuchen nicht lesen: ",
           conditionMessage(result),
           ". Falls das wiederholt passiert: LRZ Sync+Share waehrend des ",
           "gesamten main.R-Laufs pausieren, nicht nur waehrend GGIR.")
    }
  }
}

if (file.exists(ggir_rds)) {
  actigraph_sleep <- read_rds_with_retry(ggir_rds) |>
    mutate(Id = standardize_id(Id))
  cat("[02_import_actigraph] Actigraph:", n_distinct(actigraph_sleep$Id), "Teilnehmerinnen |",
      nrow(actigraph_sleep), "Naechte\n")
} else {
  warning("GGIR-RDS nicht gefunden — zuerst LEAMC_GGIR_sleep.R ausfuehren. Sleep = NA.")
  actigraph_sleep <- tibble(Id=character(), Date=as.Date(character()),
    sleep_onset_hr=numeric(), sleep_offset_hr=numeric(), TST_min=numeric(),
    WASO_min=numeric(), sleep_efficiency=numeric(), n_awakenings=integer())
}
