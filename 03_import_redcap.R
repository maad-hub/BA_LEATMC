# 03_import_redcap.R - Import and scoring of the REDCap/MyCap questionnaires

redcap_raw <- read_csv(REDCAP_FILE, show_col_types = FALSE)

# ---- Split modules ----
baseline_I    <- redcap_raw |> filter(redcap_repeat_instrument == "module_epyppx1r4kemxk69")
baseline_II   <- redcap_raw |> filter(redcap_repeat_instrument == "module_xxj4o7cv7zj8qagq")
daily_drsp    <- redcap_raw |> filter(redcap_repeat_instrument == "module_34ornc7vxo0ga7rd")
daily_light   <- redcap_raw |> filter(redcap_repeat_instrument == "module_yorhe884t5obmknp")
daily_csd     <- redcap_raw |> filter(redcap_repeat_instrument == "module_ucee05sjr6b0ksog")
menses_module <- redcap_raw |> filter(redcap_repeat_instrument == "module_9f4nubx1u2v9hc56")

# ---- Merge duplicate record IDs (app reinstalled) ----
DUPLICATE_ID_MAP <- c(
  "jk8dOc3UQCvv0CvwWPaeW" = "CFbulWqT3O70oO0qzRkZb",
  "3MXuloCQPUl6ML8qdvmMF" = "Fw8EKNwkEYPAOa6BLxEOz",  # P022
  "L9YQG_w9RVCrXOdjRvi1z" = "zf0xDM_OzNzsYSNQAQq7L",  # P029
  "7XSyY3htZ90-kbuqaR5ZG" = "cV2Q53dVm7km6NhU9jtHW",  # P031
  "7bsKbcwCkEkgTs9q8nsIL" = "4u3QQ2qARESMfzb9EiujG"  # P024
)

dupe_keys <- names(DUPLICATE_ID_MAP)[duplicated(names(DUPLICATE_ID_MAP))]
if (length(dupe_keys) > 0) {
  warning("DUPLICATE_ID_MAP enthaelt mehrfach vergebene Quell-IDs: ",
          paste(dupe_keys, collapse=", "), " — nur der erste Eintrag greift!")
}

fix_id <- function(df) mutate(df, field_record_id = if_else(
  field_record_id %in% names(DUPLICATE_ID_MAP),
  unname(DUPLICATE_ID_MAP[field_record_id]),
  field_record_id))

baseline_I    <- fix_id(baseline_I)
baseline_II   <- fix_id(baseline_II)
daily_drsp    <- fix_id(daily_drsp)
daily_light   <- fix_id(daily_light)
daily_csd     <- fix_id(daily_csd)
menses_module <- fix_id(menses_module)

# ---- Map REDCap record IDs to study IDs ----
id_map <- baseline_I |>
  filter(!is.na(field_wjxd11lihf59oetp)) |>
  group_by(field_record_id) |> slice(1) |> ungroup() |>
  transmute(redcap_id = as.character(field_record_id),
            Id        = standardize_id(field_wjxd11lihf59oetp))

# ---- Correct mistyped study IDs ----
ID_TYPO_CORRECTIONS <- c("PO10" = "P010")

n_typo_fixed <- sum(id_map$Id %in% names(ID_TYPO_CORRECTIONS))
if (n_typo_fixed > 0) {
  message("[03_import_redcap] ", n_typo_fixed, " ID-Tippfehler korrigiert: ",
          paste(names(ID_TYPO_CORRECTIONS), "->", ID_TYPO_CORRECTIONS, collapse=", "))
}
id_map <- id_map |>
  mutate(Id = if_else(Id %in% names(ID_TYPO_CORRECTIONS),
                       unname(ID_TYPO_CORRECTIONS[Id]), Id))

# ---- Remove test submissions ----
TEST_SUBMISSION_IDS <- c("TEST_ID_01", "123", "TEST_RECOVERYID")

n_test_removed <- sum(id_map$Id %in% TEST_SUBMISSION_IDS)
if (n_test_removed > 0) {
  message("[03_import_redcap] ", n_test_removed, " Test-/Dummy-Submission(en) aus id_map ",
          "entfernt: ", paste(id_map$Id[id_map$Id %in% TEST_SUBMISSION_IDS], collapse=", "))
}
id_map <- id_map |> filter(!Id %in% TEST_SUBMISSION_IDS)

get_id <- function(x) id_map$Id[match(as.character(x), id_map$redcap_id)]

validate_id_format(id_map$Id, "id_map / REDCap self-entered ID")

cat("[03_import_redcap] REDCap:", nrow(id_map), "Studien-IDs:",
    paste(sort(id_map$Id), collapse=" "), "\n")

n_baseline_II_dupes <- baseline_II |> count(field_record_id) |> filter(n > 1) |> nrow()
if (n_baseline_II_dupes > 0) {
  warning(n_baseline_II_dupes, " Teilnehmerin(nen) mit mehrfacher baseline_II-Submission ",
          "gefunden — dedupliziert auf die vollstaendigste Zeile (siehe demographics-Join).")
}

# ---- Demographics (baseline I and II) and AUDIT score ----
n_baseline_I_dupes <- baseline_I |> count(field_record_id) |> filter(n > 1) |> nrow()
if (n_baseline_I_dupes > 0) {
  warning(n_baseline_I_dupes, " Teilnehmerin(nen) mit mehrfacher baseline_I-Submission ",
          "gefunden (z.B. nach Duplicate-ID-Merge oben) — dedupliziert auf die ",
          "vollstaendigste Zeile (siehe demographics-Konstruktion).")
}

demographics <- baseline_I |>
  # keep the most complete submission per participant
  group_by(field_record_id) |>
  arrange(desc(!is.na(field_7b644otbqsidyi9g)), .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  transmute(
    Id        = get_id(field_record_id),
    age       = as.numeric(field_7b644otbqsidyi9g),
    weight_kg = as.numeric(field_po2c5f7omf6nh4ld),
    height_cm = as.numeric(field_md7qtlp9eznwmgx3),
    bmi       = weight_kg / (height_cm / 100)^2,
    audit_raw_1  = as.integer(field_ovjne2zlb93du7sx) - 1L,
    audit_raw_2  = as.integer(field_fsenqmh3zif630rg) - 1L,
    audit_raw_3  = as.integer(field_un0gb4teoelxi5hh) - 1L,
    audit_raw_4  = as.integer(field_j4mglxyebr26j3yy) - 1L,
    audit_raw_5  = as.integer(field_6hmizwzhcian9z6y) - 1L,
    audit_raw_6  = as.integer(field_whd1ri29ffal2jlb) - 1L,
    audit_raw_7  = as.integer(field_0katid9hjl78nub0) - 1L,
    audit_raw_8  = as.integer(field_ufqwn26x6d8y06ud) - 1L,
    audit_9  = case_when(as.integer(field_patvl08dti1ijm2r)==1~0L,
                          as.integer(field_patvl08dti1ijm2r)==2~2L,
                          as.integer(field_patvl08dti1ijm2r)==3~4L, TRUE~NA_integer_),
    audit_10 = case_when(as.integer(field_d5t7eog789ixjtwm)==1~0L,
                          as.integer(field_d5t7eog789ixjtwm)==2~2L,
                          as.integer(field_d5t7eog789ixjtwm)==3~4L, TRUE~NA_integer_)
  ) |>
  rowwise() |>
  mutate(audit_total = sum(c(audit_raw_1, audit_raw_2, audit_raw_3, audit_raw_4,
                                   audit_raw_5, audit_raw_6, audit_raw_7, audit_raw_8,
                                   audit_9, audit_10), na.rm=TRUE)) |>
  ungroup() |>
  # AUDIT >= 8: hazardous use (WHO)
  mutate(audit_hazardous = audit_total >= 8) |>
  left_join(
    # keep one baseline II submission per participant
    baseline_II |>
      group_by(field_record_id) |>
      arrange(desc(!is.na(field_rf0t3lbc0uu6xhqt)), .by_group = TRUE) |>
      slice(1) |>
      ungroup() |>
      transmute(
        Id             = get_id(field_record_id),
        pms_diagnosed  = as.integer(field_rf0t3lbc0uu6xhqt),
        pmdd_diagnosed = as.integer(field_za66hcbp19ye12ef)),
    by = "Id")

# ---- Daily DRSP (20 symptom items, 3 impairment items; 1-6) ----
surveys_drsp <- daily_drsp |>
  transmute(
    Id   = get_id(field_record_id),
    Date = as.Date(as.POSIXct(field_response_time_1, tz=TIMEZONE), tz=TIMEZONE),
    menses                = as.integer(field_r2ta1s949d2gnsm3),  # first day of menses
    currently_menstruating= as.integer(field_2zeg2auur05jelan),
    menstrual_flow        = as.numeric(field_r4nzelaw5eb5rjhz),
    drsp_hopeless         = as.numeric(field_r6h3sbd3hawmqj58),
    drsp_worthless        = as.numeric(field_htzi3ffxlbmci6e2),
    drsp_anxious          = as.numeric(field_kvdln1tw9ha9jiwf),
    drsp_mood_swings      = as.numeric(field_39hu6drf7qfkioke),
    drsp_rejection        = as.numeric(field_27rffesa3xwfz11g),
    drsp_angry            = as.numeric(field_crues42m75t407e5),
    drsp_conflicts        = as.numeric(field_0cwuxa63lf1hwo1h),
    drsp_less_interest    = as.numeric(field_a8xdtyf516c1l14j),
    drsp_concentration    = as.numeric(field_t959twxfccwalubq),
    drsp_fatigue          = as.numeric(field_si8ujsaoahpvk7wn),
    drsp_appetite         = as.numeric(field_ns53ei3uzpg6gwfz),
    drsp_food_cravings    = as.numeric(field_6lkm954pb1gymmaa),
    drsp_hypersomnia      = as.numeric(field_q2x7f5g3mob5wre5),
    drsp_insomnia         = as.numeric(field_umk6fyylaivvlv8v),
    drsp_overwhelmed      = as.numeric(field_60bi2i0l09p8a21a),
    drsp_out_of_control   = as.numeric(field_xb8nw1ohgcioyzvc),
    drsp_breast_tender    = as.numeric(field_s1f4qwxghabpeemi),
    drsp_bloated          = as.numeric(field_0csxp3mg8cdwga73),
    drsp_headache         = as.numeric(field_0z11t5jlzpvp26ge),
    drsp_joint_pain       = as.numeric(field_3dsgrvdciuaa4n5u),
    # impairment items are not part of the symptom index
    drsp_impair_work      = as.numeric(field_ihw96rrdp9glhyni),
    drsp_impair_social    = as.numeric(field_lhlziqlvp2zf544z),
    drsp_impair_relations = as.numeric(field_4oesr5yrs5yt8bek)
  )

# ---- Daily light sensitivity (5 items; 1-6) ----
surveys_light <- daily_light |>
  transmute(
    Id   = get_id(field_record_id),
    Date = as.Date(as.POSIXct(field_response_time_2, tz=TIMEZONE), tz=TIMEZONE),
    ls_1 = as.numeric(field_7iy80i2xxvcksmoi),  # more sensitive to light than usual
    ls_2 = as.numeric(field_s4irmsvzggo80nox),  # normal indoor lighting uncomfortable
    ls_3 = as.numeric(field_4x8hurz9bn93o8pl),  # bright light worsened mood/irritability
    ls_4 = as.numeric(field_o22l7z6bvzzdigbj),  # light triggered headache/eye strain
    ls_5 = as.numeric(field_wummjf5so9ytapm9)  # avoided or adjusted lighting
  )

# ---- Daily sleep diary (naps, substances, sleep medication) ----
surveys_csd <- daily_csd |>
  transmute(
    Id               = get_id(field_record_id),
    Date             = as.Date(as.POSIXct(field_response_time_3, tz=TIMEZONE), tz=TIMEZONE),
    nap_count        = as.numeric(field_4kuc71kwk4y6ds8z),
    nap_duration_min = as.numeric(field_7ba9hhqpipltmbg3),
    alcohol_drinks   = as.numeric(field_zyajcs7dwine8t0t),
    caffeine_drinks  = as.numeric(field_j7exzdvscdbuyfk9),
    sleep_medication = as.integer(field_gl8owd5jxmn67atf)
  )

# ---- Menses event module (fills gaps in the daily menses item) ----
menses_extra <- menses_module |>
  filter(!is.na(field_191a05wk8fiib93u)) |>
  transmute(Id = get_id(field_record_id),
            Date = as.Date(as.POSIXct(field_191a05wk8fiib93u, tz=TIMEZONE), tz=TIMEZONE),
            menses_extra = 1L) |>
  distinct(Id, Date, .keep_all = TRUE)

# ---- Onboarding: first day of the last menstrual period (used only for PACTS) ----
onboarding <- baseline_I |>
  filter(!is.na(field_kiimkjcx38pvon78)) |>
  group_by(field_record_id) |> slice(1) |> ungroup() |>
  transmute(Id = get_id(field_record_id),
            Date = as.Date(field_kiimkjcx38pvon78),
            menses = 1L)

# ---- Manual corrections: confirmed menses onsets missing in the raw data ----
manual_menses <- tribble(
  ~Id,    ~Date,                 ~note,
  "P018", as.Date("2026-07-24"), "Manuell bestaetigter Blutungsbeginn, fehlte in Rohdaten (laengerer Zyklus)",
  "P010", as.Date("2026-07-07"), "currently_menstruating=1 am 07./08.07., 'erster Tag'-Flag nie gesetzt; 07.07. ist erster verfuegbarer Tag nach Erhebungsluecke (01.07. noch 0) - gibt plausiblen ~29-Tage-Zyklus ab Onboarding-Anker 08.06.",
  "P028", as.Date("2026-08-21"), "currently_menstruating=1, 'erster Tag'-Flag nie gesetzt; Blutungsbeginn nach blutungsfreiem Vortag"
) |>
  mutate(Id = standardize_id(Id))

# ---- Manual corrections: confirmed false menses-onset entries ----
MENSES_FALSE_POSITIVES <- tribble(
  ~Id,    ~Date,                 ~note,
  "P007", as.Date("2026-07-07"), "Isolierter Einzeltag mit currently_menstruating=1 (Tag davor/danach beide 0) - Spotting, faelschlich als 'erster Tag der Menses' markiert. Echter Zyklusbeginn: 15.07. (mehrtaegige Blutung 15.-17.07.).",
  "P005", as.Date("2026-06-14"), "Zweiter von zwei aufeinanderfolgenden Tagen (13./14.06.) jeweils als 'erster Tag' markiert - erzeugte kuenstlichen 1-Tage-Zyklus in cyclenum. Echter Zyklusbeginn: 13.06.",
  "P019", as.Date("2026-07-28"), "Zweiter von zwei aufeinanderfolgenden Tagen (27./28.07.) jeweils als 'erster Tag' markiert - erzeugte kuenstlichen 1-Tage-Zyklus in cyclenum. Echter Zyklusbeginn: 27.07.",
  "P009", as.Date("2026-06-19"), "Widerspruechlich: menses=1 aber currently_menstruating=0. Der Folgetag (20.06.) hat menses=1 UND currently_menstruating=1 - echter Zyklusbeginn.",
  "P009", as.Date("2026-07-20"), "Widerspruechlich: menses=1 aber currently_menstruating=0. Der Vortag (19.07.) hat menses=1 UND currently_menstruating=1 - echter Zyklusbeginn.",
  "P017", as.Date("2026-07-16"), "Widerspruechlich: menses=1 aber currently_menstruating=0. Der Folgetag (17.07.) hat menses=1 UND currently_menstruating=1, gefolgt von durchgehender Blutung 18./19.07. - echter Zyklusbeginn. Erzeugte kuenstlichen 1-Tage-Zyklus in cyclenum (faelschlich n_cycles=3 statt 2)."
) |>
  mutate(Id = standardize_id(Id))

validate_id_format(MENSES_FALSE_POSITIVES$Id, "MENSES_FALSE_POSITIVES")
if (nrow(MENSES_FALSE_POSITIVES) > 0) {
  message("[03_import_redcap] ", nrow(MENSES_FALSE_POSITIVES), " Menses-Fehlmarkierung(en) werden entfernt:")
  print(MENSES_FALSE_POSITIVES, n = Inf)
}

validate_id_format(manual_menses$Id, "manual_menses")

if (nrow(manual_menses) > 0) {
  message("[03_import_redcap] ", nrow(manual_menses), " manuelle Menses-Korrektur(en) aktiv:")
  print(manual_menses, n = Inf)
}

# ---- De-duplicate: one row per Id and Date in each module ----
dedupe_daily <- function(df, label) {
  n_before <- nrow(df)
  dup_by_id <- df |> count(Id, Date) |> filter(n > 1) |> count(Id, wt = n - 1, name = "n_dup")

  # prefer the row with menses = 1
  if ("menses" %in% names(df)) {
    df <- df |> arrange(Id, Date, desc(!is.na(menses) & menses == 1))
  }
  df <- df |> distinct(Id, Date, .keep_all = TRUE)
  n_dup <- n_before - nrow(df)
  if (n_dup > 0) {
    message("  [Dedupe] ", label, ": ", n_dup, " doppelte Id+Date-Eintraege gefunden ",
            "und auf die jeweils bevorzugte Zeile reduziert (menses=1 falls vorhanden, ",
            "sonst die erste) (", n_before, " -> ", nrow(df), " Zeilen).")
    dup_rate <- n_dup / n_before
    if (dup_rate > 0.10) {
      n_ids_affected <- nrow(dup_by_id)
      message("    [WARNUNG] Duplikat-Rate ", round(dup_rate * 100, 1), "% ist ungewoehnlich hoch ",
              "(betrifft ", n_ids_affected, " ID(s)). Aufschluesselung nach Id (absteigend):")
      print(dup_by_id |> arrange(desc(n_dup)), n = Inf)
      message("    Sind viele/alle IDs etwa gleich stark betroffen -> eher strukturelles ",
              "Export-/Import-Problem (z.B. ueberlappende Exportzeitraeume) pruefen, statt ",
              "einzelne Doppel-Submissions anzunehmen.")
    }
  }
  df
}

surveys_drsp  <- dedupe_daily(surveys_drsp,  "surveys_drsp")
surveys_light <- dedupe_daily(surveys_light, "surveys_light")
surveys_csd   <- dedupe_daily(surveys_csd,   "surveys_csd")

# ---- Join the daily modules; the daily menses item has priority over the event module ----
DRSP_ITEMS <- c("drsp_hopeless","drsp_worthless","drsp_anxious","drsp_mood_swings",
                "drsp_rejection","drsp_angry","drsp_conflicts","drsp_less_interest",
                "drsp_concentration","drsp_fatigue","drsp_appetite","drsp_food_cravings",
                "drsp_hypersomnia","drsp_insomnia","drsp_overwhelmed","drsp_out_of_control",
                "drsp_breast_tender","drsp_bloated","drsp_headache","drsp_joint_pain")

surveys <- surveys_drsp |>
  full_join(menses_extra, by=c("Id","Date")) |>
  full_join(surveys_light, by=c("Id","Date")) |>
  full_join(surveys_csd,   by=c("Id","Date")) |>
  left_join(demographics,  by="Id") |>
  mutate(menses = if_else(!is.na(menses_extra) & menses_extra == 1L, 1L, menses)) |>
  select(-menses_extra) |>
  filter(!is.na(Id)) |>
  rowwise() |>
  mutate(
    light_sensitivity_idx = mean(c(ls_1,ls_2,ls_3,ls_4,ls_5), na.rm=TRUE),
    cycle_symptom_idx     = mean(c_across(all_of(DRSP_ITEMS)),  na.rm=TRUE)
  ) |>
  ungroup()

cat("[03_import_redcap] Surveys:", nrow(surveys), "Zeilen |", n_distinct(surveys$Id), "Teilnehmerinnen\n")
