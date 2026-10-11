# 00b_utils_io.R - Write result files without stopping the pipeline
# (falls back to a time-stamped copy if the target file is locked)

safe_write <- function(path, writer) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  cands <- unique(c(path,
                    sub("(\\.[A-Za-z0-9]+)$", paste0("_", stamp, "\\1"), path),
                    file.path(tempdir(), basename(path))))
  reasons <- character()
  for (i in seq_along(cands)) {
    res <- tryCatch({ writer(cands[i]); NULL },
                    warning = function(w) conditionMessage(w),
                    error   = function(e) conditionMessage(e))
    if (is.null(res) && file.exists(cands[i])) {
      if (i > 1) {
        message("[safe_write] '", path, "' war nicht beschreibbar (",
                paste(unique(reasons), collapse = " | "), "). ",
                "Vermutlich ist die Datei in Excel/einem Viewer geoeffnet. ",
                "Stattdessen geschrieben: ", cands[i])
      }
      return(invisible(cands[i]))
    }
    reasons <- c(reasons, if (is.null(res)) "Datei nach dem Schreiben nicht vorhanden" else res)
  }
  message("[safe_write] FEHLGESCHLAGEN fuer '", path, "': ",
          paste(unique(reasons), collapse = " | "),
          "\n  Pruefen: dir.exists('", dirname(path), "'); Datei in Excel schliessen; ",
          "Schreibrechte/Sync-Ordner.")
  invisible(NA_character_)
}

safe_write_csv <- function(x, path) {
  safe_write(path, function(p) utils::write.csv(x, p, row.names = FALSE, fileEncoding = "UTF-8"))
}

safe_write_lines <- function(lines, path) {
  safe_write(path, function(p) {
    con <- file(p, open = "w", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines(lines, con)
  })
}
