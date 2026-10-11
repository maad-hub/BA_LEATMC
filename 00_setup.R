# 00_setup.R - Packages, paths, settings and helper functions

suppressPackageStartupMessages({
  library(LightLogR)
  library(menstrualcycleR)
  library(tidyverse)
  library(lubridate)
  library(janitor)
  library(patchwork)
  library(psych)
  library(lme4)
  library(lmerTest)  # load after lme4 (Satterthwaite df and p-values)
  library(mgcv)
  library(corrplot)
  library(gt)
  library(ggbeeswarm)
})

# ---- Paths (use forward slashes) ----
BASE_DIR    <- "C:/Users/astri/LRZ Sync+Share/LEATMC"

# most recent REDCap export in BASE_DIR
redcap_candidates <- list.files(BASE_DIR, pattern = "^LightExposureAcrossT_DATA_.*\\.csv$",
                                 full.names = TRUE)
if (length(redcap_candidates) == 0) {
  stop("Keine REDCap-Export-Datei (Muster 'LightExposureAcrossT_DATA_*.csv') in ",
       BASE_DIR, " gefunden. Pruef, ob der Export dort liegt und der Dateiname ",
       "noch mit 'LightExposureAcrossT_DATA_' beginnt.")
}
REDCAP_FILE <- redcap_candidates[which.max(file.info(redcap_candidates)$mtime)]
cat("[00_setup] REDCap-Datei:", basename(REDCAP_FILE), "\n")
if (length(redcap_candidates) > 1) {
  message("  Hinweis: ", length(redcap_candidates), " passende Dateien gefunden, ",
          "neueste (", basename(REDCAP_FILE), ") wird verwendet. Alte Exporte ggf. aufraeumen.")
}

# outputs are written locally and copied to SYNC_OUTPUT_DIR at the end (11_export.R)
OUTPUT_DIR      <- "C:/GGIR_local_output/LEAMC_analysis_output"
SYNC_OUTPUT_DIR <- file.path(BASE_DIR, "output")
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
TAB_DIR    <- file.path(OUTPUT_DIR, "tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

TIMEZONE  <- "Europe/Berlin"
N_BINS    <- 28L  # PACTS bins for descriptive figures
COMPL_THR <- 80  # completeness threshold (%) for the high-completeness subsamples
MIN_BIN_N <- 5  # bins with fewer observations are flagged (low_n)
MIN_NPCRA_DAYS <- 5  # minimum valid days per participant and phase for IS/IV

# TRUE = run GGIR (LEAMC_GGIR_sleep.R); FALSE = use the saved GGIR results
RUN_GGIR <- TRUE

# ---- Helper: error-safe call (returns NULL on error) ----
safe <- function(expr, label = "") {
  tryCatch(expr, error = function(e) {
    message("  [SKIP] ", label, ": ", conditionMessage(e))
    NULL
  })
}

# ---- Helper: standardise participant IDs ----
standardize_id <- function(x) toupper(trimws(as.character(x)))

# ---- Helper: check ID format (P + three digits) ----
validate_id_format <- function(ids, source_label) {
  ids <- unique(ids[!is.na(ids)])
  bad <- ids[!grepl("^P\\d{3}$", ids)]
  if (length(bad) > 0) {
    warning("[", source_label, "] ", length(bad), " ID(s) entsprechen nicht dem Format ",
            "P0xx (nur Ziffern): ", paste(bad, collapse=", "),
            " — vermutlich Tippfehler (z.B. Buchstabe 'O' statt Ziffer '0'). ",
            "Bitte Rohdaten pruefen; diese IDs werden NICHT automatisch korrigiert ",
            "und joinen daher vermutlich gegen nichts.")
  }
  invisible(bad)
}

cat("[00_setup] OK — BASE_DIR:", BASE_DIR, "\n")
