# 04_join_trim.R - Join the three data sources and trim the first and last day

daily_data_raw <- surveys |>
  full_join(actigraph_sleep, by=c("Id","Date"))

# common first and last day per participant across all sources
date_ranges <- bind_rows(
  surveys         |> select(Id, Date),
  light_raw       |> select(Id, Date),
  actigraph_sleep |> select(Id, Date)
) |>
  group_by(Id) |>
  summarise(d_first=min(Date), d_last=max(Date), .groups="drop")

daily_data_raw <- daily_data_raw |>
  left_join(date_ranges, by="Id") |>
  filter(Date > d_first, Date < d_last) |>
  select(-d_first, -d_last)

light_trimmed <- light_raw |>
  left_join(date_ranges, by="Id") |>
  filter(Date > d_first, Date < d_last) |>
  select(-d_first, -d_last)

# activity time series: same date ranges
if (nrow(activity_epochs) > 0) {
  activity_trimmed <- activity_epochs |>
    left_join(date_ranges, by="Id") |>
    filter(Date > d_first, Date < d_last) |>
    select(-d_first, -d_last)
} else {
  activity_trimmed <- activity_epochs
}

cat("[04_join_trim] activity_trimmed:", nrow(activity_trimmed), "Aktivitaets-Epochen |",
    n_distinct(activity_trimmed$Id), "Teilnehmerinnen\n")

cat("[04_join_trim] Nach Trimmen:", nrow(daily_data_raw), "Teilnehmerin-Tage |",
    n_distinct(daily_data_raw$Id), "Teilnehmerinnen\n")
