# R/model_helpers.R - Data preparation and model helpers for 12_gamm_ar1.R
# LMMs on participant x phase means and GAMMs as in 09_inferential.R / revision.R (fixed boundary knots);
# GAMMs with AR(1) residuals: mgcv::gamm with nlme::corAR1 on the study day, so missing days remain gaps.

AC_DRSP_ITEMS <- c("drsp_hopeless","drsp_worthless","drsp_anxious","drsp_mood_swings",
                   "drsp_rejection","drsp_angry","drsp_conflicts","drsp_less_interest",
                   "drsp_concentration","drsp_fatigue","drsp_appetite","drsp_food_cravings",
                   "drsp_hypersomnia","drsp_insomnia","drsp_overwhelmed","drsp_out_of_control",
                   "drsp_breast_tender","drsp_bloated","drsp_headache","drsp_joint_pain")
AC_IMPAIR <- c("drsp_impair_work", "drsp_impair_social", "drsp_impair_relations")
AC_LS_ITEMS <- paste0("ls_", 1:5)
AC_SECONDARY <- c("MEDI_median", "MEDI_max", "MEDI_sd", "MEDI_daytime", "MEDI_night",
                  "first_250", "last_250", "M10_mean", "M10_onset", "L5_mean", "L5_onset",
                  "RA_light", "centroid_LE_hr", "MEDI_mean", "MEDI_morning", "MEDI_evening")
AC_LAT <- 48.137; AC_LON <- 11.575   # Munich

ac_try <- function(expr, label = "") tryCatch(expr, error = function(e) {
  message("  [skip] ", label, ": ", conditionMessage(e)); NULL })

# ---- Daily data with the phase columns of one assignment, study day, weekend and day length ----
ac_variant_data <- function(daily, tab = NULL, date_ranges, tz = "Europe/Berlin") {
  dd <- as.data.frame(daily)
  dd$Id <- as.character(dd$Id); dd$Date <- as.Date(dd$Date)
  dd <- dd[!is.na(dd$Id) & !is.na(dd$Date), , drop = FALSE]
  if (!is.null(tab)) {
    t <- tab[tab$in_days, , drop = FALSE]
    k <- match(paste(dd$Id, dd$Date), paste(t$Id, t$Date))
    dd$pacts_ov      <- t$pacts_ov[k]
    dd$cycle_phase_2 <- factor(as.character(t$cycle_phase_2[k]), levels = c("follicular", "luteal"))
    dd$cycle_phase_4 <- t$cycle_phase_4[k]
    dd$cycle_bin_ov  <- t$cycle_bin_ov[k]
    dd$count_method  <- t$count_method[k]
  } else if (!"count_method" %in% names(dd)) {
    dd$count_method <- ifelse(is.na(dd$pacts_ov), NA_character_, "closed_pacts")
  }
  dr <- as.data.frame(date_ranges); dr$Id <- as.character(dr$Id)
  dd$study_day <- as.integer(dd$Date - as.Date(dr$d_first[match(dd$Id, dr$Id)]))
  dd$weekend <- as.integer(as.POSIXlt(dd$Date)$wday %in% c(0, 6))
  ud <- sort(unique(dd$Date))
  dl <- tryCatch({
    st <- suncalc::getSunlightTimes(date = ud, lat = AC_LAT, lon = AC_LON, keep = c("sunrise", "sunset"), tz = tz)
    as.numeric(difftime(st$sunset, st$sunrise, units = "hours"))
  }, error = function(e) {
    doy <- as.integer(format(ud, "%j")); decl <- 23.44 * pi / 180 * sin(2 * pi * (284 + doy) / 365)
    24 * acos(pmin(1, pmax(-1, -tan(AC_LAT * pi / 180) * tan(decl)))) / pi
  })
  dd$day_length <- dl[match(dd$Date, ud)]
  dd <- dd[order(dd$Id, dd$Date), , drop = FALSE]
  dd
}

# ---- LMM on participant x phase means ----
ac_pp <- function(dd, var, ids = NULL) {
  d <- dd[!is.na(dd$cycle_phase_2) & !is.na(dd[[var]]) & !is.nan(dd[[var]]), , drop = FALSE]
  if (!is.null(ids)) d <- d[d$Id %in% ids, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  pp <- stats::aggregate(list(y = d[[var]]), by = list(Id = d$Id, cycle_phase_2 = as.character(d$cycle_phase_2)), FUN = mean)
  nd <- stats::aggregate(list(n_days = d[[var]]), by = list(Id = d$Id, cycle_phase_2 = as.character(d$cycle_phase_2)), FUN = length)
  pp <- merge(pp, nd, by = c("Id", "cycle_phase_2"))
  pp$cycle_phase_2 <- factor(pp$cycle_phase_2, levels = c("follicular", "luteal"))
  pp
}
ac_lmm_fit <- function(pp, f = y ~ cycle_phase_2 + (1 | Id)) {
  if (is.null(pp) || length(unique(pp$Id)) < 3 || length(unique(pp$cycle_phase_2)) < 2) return(NULL)
  tryCatch(suppressMessages(suppressWarnings(lmerTest::lmer(f, data = pp, REML = TRUE))),
           error = function(e) { message("  [skip] LMM could not be fitted: ", conditionMessage(e)); NULL })
}
# phase contrast with t-based CIs (Satterthwaite df)
ac_lmm_row <- function(m, term = "cycle_phase_2luteal", label = "") {
  if (is.null(m)) return(data.frame(label = label, b = NA_real_, se = NA_real_, df = NA_real_, t = NA_real_, p = NA_real_,
                                    lo95 = NA_real_, hi95 = NA_real_, lo90 = NA_real_, hi90 = NA_real_,
                                    mean_fol = NA_real_, mean_lut = NA_real_, n_obs = NA_integer_, n_id = NA_integer_,
                                    n_both = NA_integer_, singular = NA))
  s <- summary(m)$coefficients
  if (!term %in% rownames(s)) return(ac_lmm_row(NULL, term, label))
  b <- s[term, "Estimate"]; se <- s[term, "Std. Error"]; df <- s[term, "df"]
  dat <- m@frame
  idv <- if ("Id" %in% names(dat)) dat$Id else NA
  data.frame(label = label, b = b, se = se, df = df, t = s[term, "t value"], p = s[term, "Pr(>|t|)"],
             lo95 = b - stats::qt(.975, df) * se, hi95 = b + stats::qt(.975, df) * se,
             lo90 = b - stats::qt(.95, df) * se,  hi90 = b + stats::qt(.95, df) * se,
             mean_fol = s["(Intercept)", "Estimate"], mean_lut = s["(Intercept)", "Estimate"] + b,
             n_obs = stats::nobs(m), n_id = length(unique(idv)),
             n_both = sum(table(idv) == 2), singular = lme4::isSingular(m))
}
ac_lmm_pp <- function(dd, var, ids = NULL, label = var) {
  pp <- ac_pp(dd, var, ids)
  ac_lmm_row(ac_lmm_fit(pp), label = label)
}
# day-level LMM (used for study day ~ phase)
ac_lmm_day <- function(dd, var, extra = "", label = var) {
  vars <- unique(c(var, "Id", "cycle_phase_2", all.vars(stats::as.formula(paste("~ 1", extra)))))
  d <- dd[stats::complete.cases(dd[, vars, drop = FALSE]), , drop = FALSE]
  m <- tryCatch(suppressMessages(suppressWarnings(lmerTest::lmer(
    stats::as.formula(paste(var, "~ cycle_phase_2", extra, "+ (1 | Id)")), data = d))), error = function(e) NULL)
  r <- ac_lmm_row(m, label = label)
  if (!is.null(m)) { r$n_obs <- nrow(d); r$n_id <- length(unique(d$Id)); r$n_both <- NA }
  r
}

# ---- GAMMs ----
ac_gam_data <- function(dd, var, extra = "", need_sd = FALSE) {
  vars <- unique(c(var, "pacts_ov", "Id", "Date", if (need_sd) "study_day",
                   all.vars(stats::as.formula(paste("~ 1", extra)))))
  d <- dd[stats::complete.cases(dd[, vars, drop = FALSE]), , drop = FALSE]
  d <- d[order(d$Id, d$Date), , drop = FALSE]
  d$Id <- factor(d$Id)
  d
}
ac_smooth_stats <- function(g, smooth = "s(pacts_ov)") {
  st <- summary(g)$s.table
  if (!smooth %in% rownames(st)) return(c(edf = NA, ref_df = NA, F = NA, p = NA))
  r <- st[smooth, ]
  c(edf = unname(r["edf"]), ref_df = unname(r["Ref.df"]),
    F = unname(r[grep("^F$|^Chi", names(r))[1]]), p = unname(r[grep("p-value", names(r))[1]]))
}
# without AR(1): cyclic spline with fixed boundary knots, random intercept as s(Id, bs = "re") (bam, fREML)
ac_gam_noar <- function(dd, var, k = 20, extra = "") {
  d <- ac_gam_data(dd, var, extra)
  if (nrow(d) < 50) return(NULL)
  nu <- length(unique(d$pacts_ov)); k_used <- min(k, nu - 1)
  if (k_used < 4) return(NULL)
  f <- stats::as.formula(paste0(var, " ~ s(pacts_ov, bs = 'cc', k = ", k_used, ")", extra, " + s(Id, bs = 're')"))
  m <- mgcv::bam(f, data = d, method = "fREML", knots = list(pacts_ov = c(-1, 1)))
  attr(m, "k_used") <- k_used; attr(m, "n") <- nrow(d); attr(m, "n_id") <- length(unique(d$Id))
  attr(m, "data") <- d
  m
}
# with AR(1): mgcv::gamm, random intercept, nlme::corAR1(form = ~ study_day | Id)
ac_gam_ar1 <- function(dd, var, k = 20, extra = "") {
  d <- ac_gam_data(dd, var, extra, need_sd = TRUE)
  if (nrow(d) < 50) return(NULL)
  nu <- length(unique(d$pacts_ov)); k_used <- min(k, nu - 1)
  if (k_used < 4) return(NULL)
  f <- stats::as.formula(paste0(var, " ~ s(pacts_ov, bs = 'cc', k = ", k_used, ")", extra))
  fit <- function(ctrl) mgcv::gamm(f, data = d, random = list(Id = ~1),
                                   correlation = nlme::corAR1(form = ~ study_day | Id),
                                   knots = list(pacts_ov = c(-1, 1)), method = "REML", control = ctrl)
  m <- tryCatch(fit(nlme::lmeControl(returnObject = TRUE)), error = function(e) NULL)
  note <- ""
  if (is.null(m)) {
    m <- tryCatch(fit(nlme::lmeControl(opt = "optim", returnObject = TRUE, msMaxIter = 500)), error = function(e) NULL)
    note <- "optim"
  }
  if (is.null(m)) return(NULL)
  phi <- tryCatch(unname(stats::coef(m[["lme"]]$modelStruct$corStruct, unconstrained = FALSE)), error = function(e) NA_real_)
  attr(m, "phi") <- phi; attr(m, "k_used") <- k_used; attr(m, "note") <- note
  attr(m, "n") <- nrow(d); attr(m, "n_id") <- length(unique(d$Id)); attr(m, "data") <- d
  m
}
ac_gam_row <- function(m, label, model) {
  if (is.null(m)) return(data.frame(outcome = label, model = model, edf = NA, ref_df = NA, F = NA, p = NA,
                                    phi = NA, n = NA, n_id = NA, k = NA, note = "not estimable"))
  g <- if (!is.null(m[["gam"]])) m[["gam"]] else m
  s <- ac_smooth_stats(g)
  data.frame(outcome = label, model = model, edf = s[["edf"]], ref_df = s[["ref_df"]], F = s[["F"]], p = s[["p"]],
             phi = if (is.null(attr(m, "phi"))) NA else attr(m, "phi"),
             n = attr(m, "n"), n_id = attr(m, "n_id"), k = attr(m, "k_used"),
             note = if (is.null(attr(m, "note"))) "" else attr(m, "note"))
}
# population curve with pointwise and simultaneous 95% bands (Figure 4)
ac_curve <- function(m, n = 401) {
  g <- if (!is.null(m[["gam"]])) m[["gam"]] else m
  nd <- data.frame(pacts_ov = seq(-1, 1, length.out = n))
  if ("Id" %in% all.vars(stats::formula(g))) nd$Id <- levels(g$model$Id)[1]
  Xp <- stats::predict(g, newdata = nd, type = "lpmatrix")
  keep <- !grepl("^s\\(Id\\)", colnames(Xp))
  fit <- as.numeric(Xp[, keep, drop = FALSE] %*% stats::coef(g)[keep])
  sm <- grep("^s\\(pacts_ov\\)", colnames(Xp))
  Vp <- g$Vp[sm, sm, drop = FALSE]
  se_sm <- sqrt(rowSums((Xp[, sm, drop = FALSE] %*% Vp) * Xp[, sm, drop = FALSE]))
  Vfull <- g$Vp[keep, keep, drop = FALSE]
  se <- sqrt(rowSums((Xp[, keep, drop = FALSE] %*% Vfull) * Xp[, keep, drop = FALSE]))
  set.seed(20261008)
  sims <- MASS::mvrnorm(5000, rep(0, length(sm)), Vp)
  dev <- Xp[, sm, drop = FALSE] %*% t(sims)
  crit <- stats::quantile(apply(abs(dev / pmax(se_sm, 1e-12)), 2, max), .95)
  smooth <- as.numeric(Xp[, sm, drop = FALSE] %*% stats::coef(g)[sm])
  sims_f <- MASS::mvrnorm(5000, rep(0, sum(keep)), Vfull)
  dev_f <- Xp[, keep, drop = FALSE] %*% t(sims_f)
  crit_f <- stats::quantile(apply(abs(dev_f / pmax(se, 1e-12)), 2, max), .95)
  data.frame(pacts_ov = nd$pacts_ov, fit = fit, lwr = fit - 1.96 * se, upr = fit + 1.96 * se,
             lwr_sim = fit - crit_f * se, upr_sim = fit + crit_f * se,
             centred = smooth - mean(smooth), c_lwr_pw = smooth - mean(smooth) - 1.96 * se_sm,
             c_upr_pw = smooth - mean(smooth) + 1.96 * se_sm,
             c_lwr_sim = smooth - mean(smooth) - crit * se_sm, c_upr_sim = smooth - mean(smooth) + crit * se_sm)
}

# ---- IS/IV/RA per participant and phase (as in 07_descriptives.R) ----
ac_npcra <- function(activity, dd, min_days = 5) {
  ph <- dd[!is.na(dd$cycle_phase_2), c("Id", "Date", "cycle_phase_2")]
  a <- merge(as.data.frame(activity)[, c("Id", "Datetime", "ACC", "Date")], ph, by = c("Id", "Date"))
  a <- a[order(a$Id, a$Datetime), ]
  g <- split(a, list(a$Id, as.character(a$cycle_phase_2)), drop = TRUE)
  isiv <- do.call(rbind, lapply(g, function(x) data.frame(
    Id = x$Id[1], cycle_phase_2 = as.character(x$cycle_phase_2[1]),
    IS = tryCatch(LightLogR::interdaily_stability(Light.vector = x$ACC, Datetime.vector = x$Datetime), error = function(e) NA_real_),
    IV = tryCatch(LightLogR::intradaily_variability(Light.vector = x$ACC, Datetime.vector = x$Datetime), error = function(e) NA_real_),
    n_days = length(unique(x$Date)))))
  ra <- stats::aggregate(list(RA = dd$RA_act[!is.na(dd$cycle_phase_2)]),
                         by = list(Id = dd$Id[!is.na(dd$cycle_phase_2)],
                                   cycle_phase_2 = as.character(dd$cycle_phase_2[!is.na(dd$cycle_phase_2)])),
                         FUN = function(z) mean(z, na.rm = TRUE))
  out <- merge(ra, isiv, by = c("Id", "cycle_phase_2"), all.x = TRUE)
  out$cycle_phase_2 <- factor(out$cycle_phase_2, levels = c("follicular", "luteal"))
  ok <- tapply(!is.na(out$n_days) & out$n_days >= min_days, out$Id, function(z) length(z) == 2 && all(z))
  attr(out, "ok_ids") <- names(ok)[ok %in% TRUE]
  out
}
ac_lmm_npcra <- function(np, var, ids = NULL, label = var) {
  d <- np[!is.na(np$cycle_phase_2) & !is.na(np[[var]]) & !is.nan(np[[var]]), , drop = FALSE]
  if (!is.null(ids)) d <- d[d$Id %in% ids, , drop = FALSE]
  d$y <- d[[var]]
  ac_lmm_row(ac_lmm_fit(d), label = label)
}
