# 08_missing_data.R - Missing-data models and high-completeness subsamples

# ---- Missingness by cycle phase and symptoms ----

daily_data <- daily_data |>
  mutate(is_missing_light = as.integer(is.na(TAT_250_hr)),
         is_missing_sleep = as.integer(is.na(TST_min)))

mnar_light <- safe(glmer(
  is_missing_light ~ cycle_phase_2 + cycle_symptom_idx + (1 | Id),
  data = daily_data, family = binomial), "MNAR light")

mnar_sleep <- safe(glmer(
  is_missing_sleep ~ cycle_phase_2 + cycle_symptom_idx + (1 | Id),
  data = daily_data, family = binomial), "MNAR sleep")

# ---- High-completeness subsample ----

hc_ids <- missing_rep |> group_by(Id) |>
  summarise(m=min(pct_light),.groups="drop") |>
  filter(m>=COMPL_THR) |> pull(Id)

cat("[08_missing_data] High-Compliance-Subsample:", length(hc_ids), "Teilnehmerinnen\n")

# strict version: participants with data in both phases
hc_ids_strict <- missing_rep |> group_by(Id) |>
  summarise(m = min(pct_light), k = n(), .groups = "drop") |>
  filter(k == 2, m >= COMPL_THR) |> pull(Id)
cat("[08_missing_data] davon mit Daten in BEIDEN Phasen (strikt):", length(hc_ids_strict),
    "| nur 1 Phase:", length(setdiff(hc_ids, hc_ids_strict)), "\n")
