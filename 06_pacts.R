# 06_pacts.R - Ovulation-centred cycle time (PACTS, menstrualcycleR), phase labels and exclusions

# menses onsets: onboarding (LMP), manual corrections and daily surveys; ovulation not confirmed (ovtoday = 0),
# PACTS imputes ovulation as next onset - 15 days (-1 = menses onset, 0 = ovulation, +1 = next onset)
pacts_input <- bind_rows(
  onboarding    |> transmute(id=Id, date=Date, menses=1L, ovtoday=0L),
  manual_menses |> transmute(id=Id, date=Date, menses=1L, ovtoday=0L),
  surveys       |> transmute(id=Id, date=Date, menses=as.integer(menses), ovtoday=0L)
) |>
  filter(!is.na(id), !is.na(date), !is.na(menses)) |>
  arrange(id, date) |>
  group_by(id, date) |> slice(1) |> ungroup()

# remove confirmed false onset entries
if (exists("MENSES_FALSE_POSITIVES") && nrow(MENSES_FALSE_POSITIVES) > 0) {
  pacts_input <- pacts_input |>
    left_join(MENSES_FALSE_POSITIVES |> transmute(id=Id, date=Date, is_false_positive=TRUE),
              by=c("id","date")) |>
    mutate(menses = if_else(!is.na(is_false_positive), 0L, menses)) |>
    select(-is_false_positive)
}

# ---- Onsets less than 7 days apart are one cycle start (earliest onset without currently_menstruating = 0 is kept) ----
MENSES_MERGE_GAP_DAYS <- 7

menses_onsets <- pacts_input |> filter(menses == 1L) |>
  left_join(surveys |> distinct(Id, Date, .keep_all = TRUE) |>
              select(id = Id, date = Date, currently_menstruating),
            by = c("id", "date")) |>
  arrange(id, date) |>
  group_by(id) |>
  mutate(gap_prev = as.numeric(date - lag(date)),
         cluster  = cumsum(is.na(gap_prev) | gap_prev >= MENSES_MERGE_GAP_DAYS)) |>
  group_by(id, cluster) |>
  mutate(n_cl      = n(),
         ok        = is.na(currently_menstruating) | currently_menstruating != 0L,
         keep_date = if (any(ok)) min(date[ok]) else min(date)) |>
  ungroup()

menses_merged <- menses_onsets |> filter(n_cl > 1, date != keep_date)
if (nrow(menses_merged) > 0) {
  message("[06_pacts] ", nrow(menses_merged), " Menses-Termin(e) mit < ", MENSES_MERGE_GAP_DAYS,
          " Tagen Abstand zum vorigen zusammengefuehrt (frueheren, nicht widerspruechlichen Termin behalten):")
  print(menses_merged |> select(id, verworfen = date, behalten = keep_date, currently_menstruating), n = Inf)
  pacts_input <- pacts_input |>
    left_join(menses_merged |> transmute(id, date, .merged = TRUE), by = c("id", "date")) |>
    mutate(menses = if_else(!is.na(.merged), 0L, menses)) |>
    select(-.merged)
} else {
  message("[06_pacts] Keine Menses-Termine mit < ", MENSES_MERGE_GAP_DAYS, " Tagen Abstand gefunden.")
}

# plausibility check of the intervals between onsets
menses_gaps <- pacts_input |>
  filter(menses == 1L) |>
  arrange(id, date) |>
  group_by(id) |>
  mutate(gap_days = as.numeric(date - lag(date))) |>
  ungroup() |>
  filter(!is.na(gap_days), gap_days < 21 | gap_days > 40) |>
  mutate(gap_type = if_else(gap_days < 21, "zu_kurz", "zu_lang"))

if (nrow(menses_gaps) > 0) {
  gaps_short <- filter(menses_gaps, gap_type == "zu_kurz")
  gaps_long  <- filter(menses_gaps, gap_type == "zu_lang")
  if (nrow(gaps_short) > 0) {
    warning(nrow(gaps_short), " Menses-Onset(s) liegen < 21 Tage nach dem ",
            "vorherigen Onset derselben Person — biologisch nicht als eigener ",
            "Zyklus plausibel. Betroffene IDs: ",
            paste(unique(gaps_short$id), collapse=", "),
            ". Vermutlich Spotting faelschlich als 'erster Tag der Menses' markiert. ",
            "Rohdaten fuer diese IDs pruefen (siehe menses_gaps).")
  }
  if (nrow(gaps_long) > 0) {
    warning(nrow(gaps_long), " Menses-Onset(s) liegen > 40 Tage nach dem ",
            "vorherigen Onset derselben Person — vermutlich wurde dazwischen ",
            "eine echte Menses nicht geloggt. Betroffene IDs: ",
            paste(unique(gaps_long$id), collapse=", "),
            ". PACTS berechnet daraus trotzdem einen 'vollstaendigen' Zyklus ",
            "mit unplausibler Laenge. Rohdaten fuer diese IDs pruefen ",
            "(siehe menses_gaps) - ggf. manual_menses ergaenzen (siehe P018-Beispiel oben).")
  }
}

# participants with at least two onsets
ids_ok <- pacts_input |> filter(menses==1) |> count(id) |>
  filter(n>=2) |> pull(id)

# PACTS per participant; failures are reported
pacts_output <- map(unique(filter(pacts_input, id %in% ids_ok)$id), function(pid) {
  safe(pacts_scaling(
    data=filter(pacts_input, id==pid),
    id=id, date=date, menses=menses, ovtoday=ovtoday,
    lower_cyclength_bound=21, upper_cyclength_bound=40
  ), label=pid)
}) |> bind_rows()

cat("[06_pacts] PACTS:", n_distinct(pacts_output$id[!is.na(pacts_output$cyclic_time_imp_ov)]),
    "Teilnehmerinnen mit gueltiger Zykluszeit\n")

# participants dropped by PACTS
ids_dropped <- setdiff(ids_ok, unique(pacts_output$id))
if (length(ids_dropped) > 0) {
  drop_info <- pacts_input |>
    filter(id %in% ids_dropped, menses == 1L) |>
    arrange(id, date) |>
    group_by(id) |>
    summarise(n_menses = n(),
              menses_dates = paste(date, collapse=", "),
              .groups="drop")
  warning(length(ids_dropped), " Teilnehmerin(nen) hatten genug Menses-Eintraege fuer einen ",
          "PACTS-Versuch, tauchen aber in KEINER Zeile der PACTS-Ausgabe auf (vermutlich ",
          "einziger/aller Zyklen ausserhalb 21-40 Tage): ", paste(ids_dropped, collapse=", "),
          ". Details in ids_dropped_pacts_info.")
  print(drop_info, n = Inf)
  ids_dropped_pacts_info <- drop_info
}

# ---- Rename and merge into daily_data ----
daily_data <- daily_data |>
  left_join(
    pacts_output |> mutate(Id=as.character(id)) |>
      select(Id, Date=date,
             cycle_number   = cyclenum,
             ovtoday_impute,
             pacts_ov       = cyclic_time_imp_ov,  # ovulation-centred (-1 to +1)
             pacts_men      = cyclic_time_impute,  # menses-centred
             cycle_length   = mcyclength,
             luteal_length,
             cycle_incomplete),
    by=c("Id","Date")
  ) |>
  mutate(
    # cycle bins (descriptive figures only)
    cycle_bin_ov = cut(pacts_ov, breaks=seq(-1,1,length.out=N_BINS+1),
                       labels=seq_len(N_BINS), include.lowest=TRUE) |> as.integer(),
    cycle_phase_2 = case_when(pacts_ov <  0 ~ "follicular",
                               pacts_ov >= 0 ~ "luteal",
                               TRUE ~ NA_character_) |>
                   factor(levels=c("follicular","luteal")),
    # four-phase label from PACTS thresholds
    cycle_phase_4 = case_when(
      pacts_ov >= -1    & pacts_ov < -0.64 ~ "menstrual",
      pacts_ov >= -0.64 & pacts_ov <  0    ~ "follicular",
      pacts_ov == 0                         ~ "ovulatory",
      pacts_ov >   0                        ~ "luteal",
      TRUE                                  ~ NA_character_) |>
      factor(levels=c("menstrual","follicular","ovulatory","luteal"))
  )

# ---- Exclusions ----
EXCLUSIONS <- tribble(
  ~Id,    ~reason,
  "P031", "BMI-Einschlusskriterium verletzt (31.6, Bereich 18-26)",
  "P030", "frueher Studienabbruch (0 Survey-Zeilen, 0 Actigraph-Naechte)",
  "P027", "frueher Studienabbruch (5 Survey-Zeilen, 0 Actigraph-Naechte, 0 Menses)"
)
EXCLUDED_PARTICIPANT_IDS <- EXCLUSIONS$Id

n_excluded <- sum(daily_data$Id %in% EXCLUDED_PARTICIPANT_IDS, na.rm=TRUE)
if (n_excluded > 0) {
  message("[06_pacts] ", length(unique(daily_data$Id[daily_data$Id %in% EXCLUDED_PARTICIPANT_IDS])),
          " Teilnehmerin(nen) vor der Analyse ausgeschlossen (", n_excluded, " Zeilen entfernt):")
  print(EXCLUSIONS, n = Inf)
}
daily_data <- daily_data |> filter(!Id %in% EXCLUDED_PARTICIPANT_IDS)

# ---- Diagnosis: participants without a valid PACTS phase ----
all_ids_with_menses <- pacts_input |> filter(menses == 1L) |> distinct(id) |> pull(id)

pacts_phase_diagnosis <- tibble(Id = all_ids_with_menses) |>
  filter(!Id %in% EXCLUDED_PARTICIPANT_IDS) |>
  rowwise() |>
  mutate(
    n_menses = sum(pacts_input$id == Id & pacts_input$menses == 1L),
    has_valid_phase = Id %in% unique(daily_data$Id[!is.na(daily_data$cycle_phase_2)]),
    reason = case_when(
      has_valid_phase                    ~ "OK \u2014 hat gueltige Phase",
      n_menses < 2                       ~ "zu wenig Menses-Eintraege (<2)",
      !(Id %in% unique(pacts_output$id)) ~ "von PACTS komplett verworfen (0 Zeilen, z.B. Zyklus ausserhalb 21-40 Tage)",
      TRUE                                ~ "PACTS liefert Zeilen, aber cyclic_time_imp_ov bleibt NA (z.B. laufender/unvollstaendiger letzter Zyklus)"
    )
  ) |>
  ungroup()

cat("[06_pacts] Phasenzuordnung-Diagnose (Grund je Teilnehmerin):\n")
print(pacts_phase_diagnosis |> count(reason), n = Inf)
