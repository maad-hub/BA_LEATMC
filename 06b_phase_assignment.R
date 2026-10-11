# 06b_phase_assignment.R - Extended cycle-phase assignment used for the main analyses
# "A" = extended assignment (closed cycles plus forward/backward counting from prospective onsets),
# "original" = closed cycles of 21-40 days only (as 06_pacts.R; sensitivity analysis).
# daily_data receives the phase columns of "A"; the state after 06_pacts.R is kept as daily_data_06.

phase_ctx <- pa_build_ctx()
phase_tabs <- list(original = pa_assign(phase_ctx, "original"), A = pa_assign(phase_ctx, "A"))

# check: the own PACTS implementation reproduces the assignment of 06_pacts.R
local({
  o <- phase_tabs$original[phase_tabs$original$in_days, ]
  k <- match(paste(as.character(daily_data$Id), as.Date(daily_data$Date)), paste(o$Id, o$Date))
  same <- function(x, y) (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
  ok <- c(pacts = all(same(daily_data$pacts_ov, o$pacts_ov[k])),
          phase = all(same(as.character(daily_data$cycle_phase_2), as.character(o$cycle_phase_2[k]))),
          bin = all(same(daily_data$cycle_bin_ov, o$cycle_bin_ov[k])))
  if (!all(ok)) stop("[06b] original assignment does not reproduce daily_data: ", paste(names(ok)[!ok], collapse = ", "))
})

daily_data_06 <- daily_data
daily_data <- pa_apply_assignment(daily_data, phase_tabs$A)
cat("[06b] participants with a cycle phase: original", length(unique(daily_data_06$Id[!is.na(daily_data_06$cycle_phase_2)])),
    "| extended", length(unique(daily_data$Id[!is.na(daily_data$cycle_phase_2)])), "\n")
