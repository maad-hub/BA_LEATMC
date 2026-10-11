# 07_descriptives.R - Descriptive statistics

# ---- Participant x phase means ----

ALL_METRICS <- intersect(c(
  "TAT_250_hr","TAT_1000_hr","TAT_10_hr",
  "MEDI_mean","MEDI_median","MEDI_max","MEDI_sd",
  "MEDI_morning","MEDI_evening","MEDI_daytime","MEDI_night",
  "first_250","last_250","M10_mean","M10_onset","M10_midpoint",
  "L5_mean","L5_onset","RA_light","centroid_LE_hr","brown_compliance_pct",
  "light_sensitivity_idx","ls_1","ls_2","ls_3","ls_4","ls_5",
  "cycle_symptom_idx", DRSP_ITEMS,
  "drsp_impair_work","drsp_impair_social","drsp_impair_relations",
  "menstrual_flow","currently_menstruating",
  "TST_min","sleep_efficiency","WASO_min","n_awakenings",
  "sleep_onset_hr","sleep_offset_hr",
  "nap_count","nap_duration_min","alcohol_drinks","caffeine_drinks","sleep_medication"
), names(daily_data))

metrics_phase2 <- daily_data |>
  filter(!is.na(cycle_phase_2)) |>
  group_by(Id, cycle_phase_2) |>
  summarise(n_days=n(),
            across(all_of(ALL_METRICS),
                   list(mean=~mean(.x,na.rm=TRUE), sd=~sd(.x,na.rm=TRUE)),
                   .names="{.col}__{.fn}"),
            .groups="drop")

# bin summary for group figures (mean with 95% CI)
bin_summary <- daily_data |>
  filter(!is.na(cycle_bin_ov)) |>
  group_by(cycle_bin=cycle_bin_ov) |>
  summarise(n=n(),
            TAT_250_mean       = mean(TAT_250_hr, na.rm=TRUE),
            TAT_250_sd         = sd(TAT_250_hr, na.rm=TRUE),
            cycle_symptom_mean = mean(cycle_symptom_idx, na.rm=TRUE),
            cycle_symptom_sd   = sd(cycle_symptom_idx, na.rm=TRUE),
            MEDI_mean          = mean(MEDI_mean, na.rm=TRUE),
            MEDI_sd            = sd(MEDI_mean, na.rm=TRUE),
            RA_light_mean      = mean(RA_light, na.rm=TRUE),
            .groups="drop") |>
  mutate(TAT_250_ci = 1.96 * TAT_250_sd / sqrt(n),
         low_n = n < MIN_BIN_N)

low_n_bins <- bin_summary |> filter(low_n)
if (nrow(low_n_bins) > 0) {
  message("[07_descriptives] ", nrow(low_n_bins), " PACTS-Bin(s) mit n < ", MIN_BIN_N,
          " (instabiler Mittelwert, kein Ausreisserfilter fuer MEDI vorhanden):")
  print(select(low_n_bins, cycle_bin, n, MEDI_mean, MEDI_sd, TAT_250_mean), n = Inf)
}

# Cronbach's alpha
cat("[07_descriptives] Cronbach a - Light Sensitivity:",
    round(safe(psych::alpha(drop_na(select(surveys,ls_1:ls_5)))$total$raw_alpha, "a ls"), 3), "\n")
cat("[07_descriptives] Cronbach a - Cycle Symptom    :",
    round(safe(psych::alpha(drop_na(select(surveys,all_of(DRSP_ITEMS))))$total$raw_alpha, "a drsp"), 3), "\n")

# ---- Rest-activity: IS/IV from the activity time series, RA from daily RA_act ----

if (nrow(activity_trimmed) > 0) {
  npcra_overall <- activity_trimmed |>
    group_by(Id) |>
    summarise(
      IS = safe(interdaily_stability( Light.vector=ACC, Datetime.vector=Datetime), "IS"),
      IV = safe(intradaily_variability(Light.vector=ACC, Datetime.vector=Datetime), "IV"),
      .groups="drop")
} else {
  npcra_overall <- tibble(Id=character(), IS=numeric(), IV=numeric())
}
npcra_overall <- npcra_overall |>
  left_join(daily_data |> filter(!is.na(RA_act)) |>
              group_by(Id) |> summarise(RA=mean(RA_act,na.rm=TRUE),.groups="drop"),
            by="Id")

daily_data <- daily_data |>
  left_join(npcra_overall |> rename(IS_overall=IS, IV_overall=IV, RA_overall=RA), by="Id")

npcra_phase_days <- activity_trimmed |>
  left_join(select(daily_data,Id,Date,cycle_phase_2), by=c("Id","Date")) |>
  filter(!is.na(cycle_phase_2)) |>
  distinct(Id, cycle_phase_2, Date) |>
  count(Id, cycle_phase_2, name="n_days")

npcra_phase <- daily_data |>
  filter(!is.na(cycle_phase_2)) |>
  group_by(Id, cycle_phase_2) |>
  summarise(RA=mean(RA_act,na.rm=TRUE), .groups="drop") |>
  left_join(
    activity_trimmed |>
      left_join(select(daily_data,Id,Date,cycle_phase_2), by=c("Id","Date")) |>
      filter(!is.na(cycle_phase_2)) |>
      group_by(Id, cycle_phase_2) |>
      summarise(
        IS=safe(interdaily_stability( Light.vector=ACC, Datetime.vector=Datetime),"IS"),
        IV=safe(intradaily_variability(Light.vector=ACC, Datetime.vector=Datetime),"IV"),
        .groups="drop"),
    by=c("Id","cycle_phase_2")) |>
  left_join(npcra_phase_days, by=c("Id","cycle_phase_2")) |>
  mutate(low_n_npcra = n_days < MIN_NPCRA_DAYS)

n_low_npcra <- sum(npcra_phase$low_n_npcra, na.rm=TRUE)
cat("[07_descriptives] npcra_phase: ", n_low_npcra, " von ", nrow(npcra_phase),
    " Id\u00d7Phase-Zeilen mit < ", MIN_NPCRA_DAYS, " Tagen (low_n_npcra=TRUE):\n", sep="")
print(npcra_phase |> filter(low_n_npcra) |> select(Id, cycle_phase_2, n_days, IS, IV), n = Inf)

# participants with at least MIN_NPCRA_DAYS days in both phases (IS/IV analyses)
npcra_ok_ids <- npcra_phase |>
  group_by(Id) |>
  summarise(ok = all(!low_n_npcra) && n() == 2, .groups="drop") |>
  filter(ok) |> pull(Id)
cat("[07_descriptives] npcra_ok_ids: ", length(npcra_ok_ids),
    " von 21 Teilnehmerinnen mit >= ", MIN_NPCRA_DAYS,
    " Tagen in BEIDEN Phasen (fuer lmm_IS/lmm_IV):\n", sep="")
print(sort(npcra_ok_ids))

# ---- Data completeness ----

missing_rep <- daily_data |>
  filter(!is.na(cycle_phase_2)) |>
  group_by(Id, cycle_phase_2) |>
  summarise(n=n(),
            pct_light  = mean(!is.na(TAT_250_hr))  * 100,
            pct_sleep  = mean(!is.na(TST_min))      * 100,
            pct_survey = mean(!is.na(light_sensitivity_idx)) * 100,
            .groups="drop")

cat("[07_descriptives] Unter", COMPL_THR, "% ActLumus-Compliance:",
    sum(missing_rep |> group_by(Id) |>
          summarise(m=min(pct_light),.groups="drop") |>
          pull(m) < COMPL_THR), "Teilnehmerinnen\n")

# ---- QC: actigraphy coverage ----

actigraph_raw_n    <- actigraph_sleep |> count(Id, name="n_raw_nights")
actigraph_visible_n <- daily_data |> filter(!is.na(TST_min)) |> count(Id, name="n_visible_in_daily_data")

actigraph_coverage_gap <- actigraph_raw_n |>
  full_join(actigraph_visible_n, by="Id") |>
  mutate(across(starts_with("n_"), ~ replace_na(.x, 0)),
         gap = n_raw_nights - n_visible_in_daily_data) |>
  filter(gap > 0) |>
  arrange(desc(gap))

if (nrow(actigraph_coverage_gap) > 0) {
  message("[07_descriptives] Actigraph-Abdeckungsluecke bei ", nrow(actigraph_coverage_gap),
          " ID(s) — GGIR hat mehr Naechte verarbeitet, als in daily_data sichtbar sind ",
          "(Ursache: survey-verankerter Join in 04_join_trim.R, kein Datenverlust bei GGIR selbst). ",
          "Details in actigraph_coverage_gap / QC_actigraph_coverage_gap.csv.")
  print(actigraph_coverage_gap, n = Inf)
}

actigraph_fully_missing <- actigraph_coverage_gap |> filter(n_visible_in_daily_data == 0)
if (nrow(actigraph_fully_missing) > 0) {
  warning(nrow(actigraph_fully_missing), " ID(s) haben rohe GGIR-Naechte, tauchen aber ",
          "0x in daily_data auf (", paste(actigraph_fully_missing$Id, collapse=", "),
          ") — pruef auf ID-Mismatch (siehe validate_id_format-Warnungen oben), nicht auf ",
          "reine Survey-Compliance-Luecke.")
}

# ---- Participant characteristics ----

cycle_length_per_cycle <- daily_data |>
  filter(!is.na(cycle_number), !is.na(cycle_length)) |>
  filter(is.na(cycle_incomplete) | cycle_incomplete == 0) |>
  filter(cycle_length >= 21, cycle_length <= 40) |>
  distinct(Id, cycle_number, cycle_length) |>
  group_by(Id) |>
  summarise(mean_cycle_length = mean(cycle_length, na.rm=TRUE), .groups="drop")

n_cycles_completed_tbl <- daily_data |>
  filter(!is.na(cycle_number), is.na(cycle_incomplete) | cycle_incomplete == 0) |>
  distinct(Id, cycle_number) |>
  count(Id, name = "n_cycles_completed")

char_table <- daily_data |>
  group_by(Id) |>
  summarise(age=first(na.omit(age)), height_cm=first(na.omit(height_cm)),
            weight_kg=first(na.omit(weight_kg)), bmi=first(na.omit(bmi)),
            n_cycles = { v <- cycle_number[!is.na(cycle_number)];
                         if (length(v) == 0) NA_real_ else max(v) },
            .groups="drop") |>
  left_join(n_cycles_completed_tbl, by="Id") |>
  left_join(cycle_length_per_cycle, by="Id") |>
  left_join(distinct(select(surveys, Id, pms_diagnosed, pmdd_diagnosed)), by="Id")

participant_desc <- psych::describe(char_table |> select(where(is.numeric)))
