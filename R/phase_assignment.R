# R/phase_assignment.R - Cycle-phase assignment
#   "original": closed cycles of 21-40 days only (PACTS, reproduces 06_pacts.R)
#   "A":        extended assignment (main analysis): additionally forward counting of follicular days
#               from the last prospective onset and, for participants whose LMP interval exceeds the
#               upper bound, backward counting of luteal days (15 days before each prospective onset)
# PACTS as in menstrualcycleR 0.1.5 (cyclic_time_imp_ov, ovtoday = 0):
#   ov = next onset - 15; follicular (onset <= t <= ov): -1 * (1 - (t - onset) / (ov - onset));
#   luteal (ov < t < next onset): 1 + (-1 * (1 - (t - ov) / 15)); ovulation day = 0

PA_LUTEAL_DAYS <- 15L

# ---- Labels (as in 06_pacts.R) ----
pa_label_bin <- function(p, n_bins = 28L) {
  as.integer(cut(p, breaks = seq(-1, 1, length.out = n_bins + 1),
                 labels = seq_len(n_bins), include.lowest = TRUE))
}
pa_label_phase2 <- function(p) {
  factor(dplyr::case_when(p < 0 ~ "follicular", p >= 0 ~ "luteal", TRUE ~ NA_character_),
         levels = c("follicular", "luteal"))
}
pa_label_phase4 <- function(p) {
  factor(dplyr::case_when(
    p >= -1    & p < -0.64 ~ "menstrual",
    p >= -0.64 & p <  0    ~ "follicular",
    p == 0                 ~ "ovulatory",
    p >  0                 ~ "luteal",
    TRUE                   ~ NA_character_),
    levels = c("menstrual", "follicular", "ovulatory", "luteal"))
}

# ---- Menses onsets with source and correction status ----
# entries dated before the start of monitoring count as retrospective, not prospective
pa_anchor_table <- function(pacts_input, onboarding, prospective_dates,
                            manual = NULL, false_pos = NULL, merged = NULL,
                            monitoring_start = NULL) {
  pin <- pacts_input |> dplyr::transmute(Id = as.character(id), Date = as.Date(date),
                                         menses_final = as.integer(menses))
  lmp <- onboarding |> dplyr::filter(!is.na(Id), !is.na(Date)) |>
    dplyr::transmute(Id = as.character(Id), Date = as.Date(Date), is_lmp = TRUE) |>
    dplyr::distinct(Id, Date, .keep_all = TRUE)
  pro <- prospective_dates |> dplyr::filter(!is.na(Id), !is.na(Date)) |>
    dplyr::transmute(Id = as.character(Id), Date = as.Date(Date), is_prospective = TRUE) |>
    dplyr::distinct(Id, Date, .keep_all = TRUE)
  flag <- function(x, nm) {
    if (is.null(x) || nrow(x) == 0) return(tibble::tibble(Id = character(), Date = as.Date(character()), !!nm := logical()))
    names(x)[tolower(names(x)) == "id"] <- "Id"; names(x)[tolower(names(x)) == "date"] <- "Date"
    x |> dplyr::transmute(Id = as.character(Id), Date = as.Date(Date), !!nm := TRUE) |>
      dplyr::distinct(Id, Date, .keep_all = TRUE)
  }
  cand <- dplyr::bind_rows(lmp |> dplyr::select(Id, Date), pro |> dplyr::select(Id, Date),
                           pin |> dplyr::filter(menses_final == 1L) |> dplyr::select(Id, Date)) |>
    dplyr::distinct()
  ms <- if (is.null(monitoring_start)) tibble::tibble(Id = character(), start = as.Date(character())) else
    monitoring_start |> dplyr::transmute(Id = as.character(Id), start = as.Date(start)) |> dplyr::distinct(Id, .keep_all = TRUE)
  cand |>
    dplyr::left_join(ms, by = "Id") |>
    dplyr::left_join(lmp, by = c("Id", "Date")) |>
    dplyr::left_join(pro, by = c("Id", "Date")) |>
    dplyr::left_join(pin, by = c("Id", "Date")) |>
    dplyr::left_join(flag(manual, "is_added"), by = c("Id", "Date")) |>
    dplyr::left_join(flag(false_pos, "is_removed"), by = c("Id", "Date")) |>
    dplyr::left_join(flag(merged, "is_merged"), by = c("Id", "Date")) |>
    dplyr::mutate(
      dplyr::across(c(is_lmp, is_prospective, is_added, is_removed, is_merged), ~ dplyr::coalesce(.x, FALSE)),
      used = !is.na(menses_final) & menses_final == 1L,
      before_monitoring = !is.na(start) & Date < start,
      source = dplyr::case_when((is_prospective | is_added) & !before_monitoring ~ "prospective",
                                is_lmp ~ "retrospective_LMP",
                                is_prospective | is_added ~ "retrospective_entry",
                                TRUE ~ "unknown"),
      correction = dplyr::case_when(is_removed ~ "removed", is_merged ~ "merged",
                                    is_added ~ "added", TRUE ~ "original")) |>
    dplyr::select(Id, Date, source, correction, used, is_lmp, before_monitoring) |>
    dplyr::arrange(Id, Date)
}

# cycles between consecutive onsets of one participant (closed / open)
pa_cycles <- function(anchor_dates, anchor_source, last_date, rule = c(21, 40)) {
  o <- order(anchor_dates); a <- anchor_dates[o]; s <- anchor_source[o]
  n <- length(a)
  if (n == 0) return(data.frame())
  nxt <- c(a[-1], as.Date(NA))
  L <- as.numeric(nxt - a)
  data.frame(onset = a, onset_source = s, next_onset = nxt,
             status = ifelse(is.na(nxt), "open", "closed"),
             length_days = L,
             in_rule = !is.na(L) & L >= rule[1] & L <= rule[2],
             days_observed_open = ifelse(is.na(nxt), as.numeric(last_date - a) + 1, NA_real_),
             stringsAsFactors = FALSE)
}

# ---- PACTS for closed cycles (same steps as menstrualcycleR) ----
pa_pacts_closed <- function(dates, anchors, rule = c(21, 40)) {
  dates <- as.Date(dates); anchors <- sort(unique(as.Date(anchors)))
  raw <- rep(NA_real_, length(dates))
  cyc <- rep(NA_integer_, length(dates))
  if (length(anchors) >= 2) {
    for (k in seq_len(length(anchors) - 1)) {
      a <- anchors[k]; b <- anchors[k + 1]; L <- as.numeric(b - a)
      if (L < rule[1] || L > rule[2]) next
      ov <- b - PA_LUTEAL_DAYS
      idx <- which(dates >= a & dates < b)
      t <- dates[idx]
      fol <- t <= ov
      v <- numeric(length(t))
      v[fol]  <- -1 * (1 - (as.numeric(t[fol] - a) / as.numeric(ov - a)))
      v[!fol] <- 1 + (-1 * (1 - (as.numeric(t[!fol] - ov) / (14 + 1))))
      raw[idx] <- v
      cyc[idx] <- k
    }
  }
  value <- raw
  # onset day with a scaled previous day = +1 (end of the previous cycle)
  i_on <- which(dates %in% anchors)
  prev_ok <- i_on > 1
  i_on <- i_on[prev_ok]
  hit <- i_on[!is.na(raw[i_on - 1]) & dates[i_on - 1] == dates[i_on] - 1]
  value[hit] <- 1
  # ovulation day = 0
  ov_days <- as.Date(vapply(seq_len(max(0, length(anchors) - 1)), function(k) {
    L <- as.numeric(anchors[k + 1] - anchors[k])
    if (L < rule[1] || L > rule[2]) NA_real_ else as.numeric(anchors[k + 1] - PA_LUTEAL_DAYS)
  }, numeric(1)), origin = "1970-01-01")
  value[dates %in% ov_days[!is.na(ov_days)]] <- 0
  list(value = value, raw = raw, cycle_k = cyc)
}

# ---- Assignments ----
PA_VARIANTS <- list(
  original = list(type = "original", cycle_rule = c(21, 40)),
  A        = list(type = "A",        cycle_rule = c(21, 40))
)
PA_DEFAULTS <- list(buffer_days = 3L, lhat_include_lmp = TRUE)

# cycle length for forward counting: mean of the participant's closed in-rule cycles, otherwise sample median
pa_lhat_table <- function(cycles_all, ids, rule, include_lmp) {
  pool <- cycles_all[cycles_all$status == "closed" & cycles_all$in_rule &
                       cycles_all$Id %in% ids &
                       (include_lmp | cycles_all$onset_source == "prospective"), , drop = FALSE]
  med <- if (nrow(pool)) stats::median(pool$length_days) else NA_real_
  ind <- if (nrow(pool)) stats::aggregate(length_days ~ Id, data = pool, FUN = mean) else
    data.frame(Id = character(), length_days = numeric())
  out <- data.frame(Id = ids, stringsAsFactors = FALSE)
  out$L_individual <- ind$length_days[match(out$Id, ind$Id)]
  out$n_cycles_individual <- vapply(out$Id, function(i) sum(pool$Id == i), integer(1))
  attr(out, "sample_median") <- med
  attr(out, "n_pool") <- nrow(pool)
  out
}

# phase assignment for all days in ctx$days (one row per participant and calendar day)
pa_assign <- function(ctx, variant, spec = PA_VARIANTS[[variant]], opts = PA_DEFAULTS) {
  type <- spec$type; rule <- spec$cycle_rule
  stopifnot(type %in% c("original", "A"))
  ids <- ctx$ids
  anc <- ctx$anchors[ctx$anchors$used & ctx$anchors$Id %in% ids, , drop = FALSE]

  pin <- ctx$pacts_input
  days <- ctx$days[ctx$days$Id %in% ids, c("Id", "Date"), drop = FALSE]

  cycles_all <- do.call(rbind, lapply(ids, function(i) {
    a <- anc[anc$Id == i, , drop = FALSE]
    if (!nrow(a)) return(NULL)
    last_d <- max(c(days$Date[days$Id == i], pin$date[pin$id == i]), na.rm = TRUE)
    cbind(Id = i, pa_cycles(a$Date, a$source, last_d, rule), stringsAsFactors = FALSE)
  }))
  if (is.null(cycles_all)) cycles_all <- data.frame(Id = character(), onset = as.Date(character()),
    onset_source = character(), next_onset = as.Date(character()), status = character(),
    length_days = numeric(), in_rule = logical(), days_observed_open = numeric())
  include_lmp <- isTRUE(opts$lhat_include_lmp)
  lhat <- pa_lhat_table(cycles_all, ids, rule, include_lmp)
  med  <- attr(lhat, "sample_median")

  # participants whose interval from the recalled LMP exceeds the upper bound
  li_ids <- unique(cycles_all$Id[cycles_all$status == "closed" &
                                   cycles_all$onset_source != "prospective" &
                                   cycles_all$length_days > rule[2]])

  res <- vector("list", length(ids))
  for (j in seq_along(ids)) {
    i <- ids[j]
    a_i <- anc[anc$Id == i, , drop = FALSE]; a_i <- a_i[order(a_i$Date), , drop = FALSE]
    d_i <- sort(unique(days$Date[days$Id == i]))
    g_lo <- min(c(d_i, pin$date[pin$id == i], a_i$Date), na.rm = TRUE)
    g_hi <- max(c(d_i, pin$date[pin$id == i], a_i$Date), na.rm = TRUE)
    grid <- seq(g_lo, g_hi, by = "day")
    n <- length(grid)
    val <- rep(NA_real_, n); meth <- rep(NA_character_, n); src <- rep(NA_character_, n)

    pc <- pa_pacts_closed(grid, a_i$Date, rule)
    val <- pc$raw
    meth[!is.na(val)] <- "closed_pacts"
    # anchor source: follicular days belong to the opening onset, ovulation and luteal days to the closing onset
    k_ok <- which(!is.na(pc$cycle_k))
    kk <- pc$cycle_k[k_ok]
    src[k_ok] <- ifelse(val[k_ok] < 0, a_i$source[kk], a_i$source[kk + 1])

    is_pro <- a_i$source == "prospective"
    cyc_i <- cycles_all[cycles_all$Id == i, , drop = FALSE]
    scaled_after <- if (nrow(cyc_i)) cyc_i$status == "closed" & cyc_i$in_rule else logical(0)

    do_backward <- type == "A" & i %in% li_ids
    do_forward <- type == "A"

    # (1) backward counting: luteal days before prospective onsets
    bw_anchor <- if (do_backward) is_pro else rep(FALSE, nrow(a_i))
    for (k in which(bw_anchor)) {
      a <- a_i$Date[k]; ov <- a - PA_LUTEAL_DAYS
      prev_a <- if (k > 1) a_i$Date[k - 1] else as.Date(-Inf)
      idx <- which(grid >= ov & grid < a & grid > prev_a & is.na(val))
      t <- grid[idx]
      val[idx] <- 1 + (-1 * (1 - (as.numeric(t - ov) / (14 + 1))))
      meth[idx] <- "backward_imputed"; src[idx] <- "prospective"
    }

    # (2) onset day = +1 if the previous day is scaled; ovulation days of closed cycles = 0
    pre <- val
    i_on <- which(grid %in% a_i$Date); i_on <- i_on[i_on > 1]
    hit <- i_on[!is.na(pre[i_on - 1])]
    val[hit] <- 1; meth[hit] <- meth[hit - 1]
    src[hit] <- a_i$source[match(grid[hit], a_i$Date)]
    if (nrow(cyc_i)) {
      ovd <- cyc_i$next_onset[scaled_after] - PA_LUTEAL_DAYS
      val[grid %in% ovd] <- 0
    }

    # (3) forward counting: follicular days after prospective onsets without a scaled following cycle,
    #     up to buffer_days before the estimated ovulation (onset + round(L_hat) - 15)
    if (do_forward && nrow(a_i)) {
      fw_anchor <- which(is_pro & !scaled_after)
      L_ind <- lhat$L_individual[lhat$Id == i]
      for (k in fw_anchor) {
        L_hat <- if (is.na(L_ind)) med else L_ind
        if (is.na(L_hat)) next
        a <- a_i$Date[k]
        ov_hat <- a + floor(L_hat + 0.5) - PA_LUTEAL_DAYS
        nxt <- if (k < nrow(a_i)) a_i$Date[k + 1] else as.Date(Inf)
        idx <- which(grid > a & grid <= ov_hat - opts$buffer_days & grid < nxt & is.na(val))
        t <- grid[idx]
        val[idx] <- -1 * (1 - (as.numeric(t - a) / as.numeric(ov_hat - a)))
        meth[idx] <- "forward_imputed"; src[idx] <- "prospective"
        ia <- which(grid == a)
        if (length(ia) && is.na(val[ia]) && length(idx)) {
          val[ia] <- -1; meth[ia] <- "forward_imputed"; src[ia] <- "prospective"
        }
      }
    }

    out <- data.frame(Id = i, Date = grid, pacts_ov = val, count_method = meth,
                      anchor_source = src, stringsAsFactors = FALSE)
    out$in_days <- out$Date %in% d_i
    res[[j]] <- out
  }
  tab <- do.call(rbind, res)
  if (is.null(tab)) tab <- data.frame(Id = character(), Date = as.Date(character()), pacts_ov = numeric(),
                                      count_method = character(), anchor_source = character(), in_days = logical())

  tab$cycle_phase_2 <- pa_label_phase2(tab$pacts_ov)
  tab$cycle_phase_4 <- pa_label_phase4(tab$pacts_ov)
  tab$cycle_bin_ov  <- pa_label_bin(tab$pacts_ov)

  tab$long_interval <- tab$Id %in% li_ids
  tab$variant <- variant
  attr(tab, "lhat") <- lhat
  attr(tab, "cycles") <- cycles_all
  attr(tab, "long_interval_ids") <- li_ids
  attr(tab, "spec") <- spec
  tab
}

# ---- Input from the session after 06_pacts.R ----
pa_build_ctx <- function(env = globalenv()) {
  g <- function(n) if (exists(n, envir = env)) get(n, envir = env) else NULL
  need <- c("pacts_input", "onboarding", "surveys", "manual_menses", "daily_data", "date_ranges")
  miss <- need[!vapply(need, exists, logical(1), envir = env)]
  if (length(miss)) stop("[pa_build_ctx] objects missing (run 00-06 first): ", paste(miss, collapse = ", "))
  dd <- g("daily_data"); sv <- g("surveys"); mm <- g("manual_menses")
  ids <- sort(unique(as.character(dd$Id[!is.na(dd$Id)])))
  prospective <- dplyr::bind_rows(
    sv |> dplyr::filter(!is.na(menses), menses == 1L) |> dplyr::select(Id, Date),
    mm |> dplyr::select(Id, Date))
  anchors <- pa_anchor_table(g("pacts_input"), g("onboarding"), prospective, manual = mm,
                             false_pos = g("MENSES_FALSE_POSITIVES"), merged = g("menses_merged"),
                             monitoring_start = g("date_ranges") |> dplyr::transmute(Id, start = d_first))
  days <- dd |> dplyr::transmute(Id = as.character(Id), Date = as.Date(Date)) |>
    dplyr::filter(!is.na(Id), !is.na(Date)) |> dplyr::distinct()
  list(pacts_input = g("pacts_input") |> dplyr::mutate(id = as.character(id), date = as.Date(date)),
       anchors = anchors, days = days, ids = ids)
}

# replaces the four phase columns of daily_data with those of an assignment
pa_apply_assignment <- function(daily, tab) {
  t <- tab[tab$in_days, , drop = FALSE]
  k <- match(paste(as.character(daily$Id), as.Date(daily$Date)), paste(t$Id, t$Date))
  daily$pacts_ov      <- t$pacts_ov[k]
  daily$cycle_phase_2 <- factor(as.character(t$cycle_phase_2[k]), levels = c("follicular", "luteal"))
  daily$cycle_phase_4 <- factor(as.character(t$cycle_phase_4[k]), levels = c("menstrual", "follicular", "ovulatory", "luteal"))
  daily$cycle_bin_ov  <- as.integer(t$cycle_bin_ov[k])
  daily
}
