# R/figure_helpers.R - Colours, labels and helpers for 13_thesis_figures.R (Okabe-Ito palette, 300 dpi)

E2_OI <- c(blue = "#0072B2", sky = "#56B4E9", green = "#009E73", orange = "#E69F00",
           vermillion = "#D55E00", purple = "#CC79A7", yellow = "#F0E442", black = "#000000", grey = "#999999")
E2_DOMAINS <- list(   # as DRSP_DOMAINS in 10_visualisations.R
  "Mood & Affect"           = c("drsp_hopeless", "drsp_worthless", "drsp_anxious", "drsp_mood_swings",
                                "drsp_rejection", "drsp_angry", "drsp_conflicts"),
  "Cognitive & Behavioural" = c("drsp_less_interest", "drsp_concentration", "drsp_overwhelmed", "drsp_out_of_control"),
  "Physical"                = c("drsp_fatigue", "drsp_appetite", "drsp_food_cravings", "drsp_hypersomnia",
                                "drsp_insomnia", "drsp_breast_tender", "drsp_bloated", "drsp_headache", "drsp_joint_pain"),
  "Impairment"              = c("drsp_impair_work", "drsp_impair_social", "drsp_impair_relations"))
E2_DOMAIN_COLS <- c("Mood & Affect" = E2_OI[["vermillion"]], "Cognitive & Behavioural" = E2_OI[["blue"]],
                    "Physical" = E2_OI[["green"]], "Impairment" = "#555555")
E2_ITEM_LABELS <- c(
  drsp_hopeless = "Hopeless", drsp_worthless = "Worthless", drsp_anxious = "Anxious", drsp_mood_swings = "Mood swings",
  drsp_rejection = "Rejection sensitivity", drsp_angry = "Angry/irritable", drsp_conflicts = "Conflicts",
  drsp_less_interest = "Less interest", drsp_concentration = "Concentration", drsp_overwhelmed = "Overwhelmed",
  drsp_out_of_control = "Out of control", drsp_fatigue = "Fatigue", drsp_appetite = "Appetite",
  drsp_food_cravings = "Food cravings", drsp_hypersomnia = "Hypersomnia", drsp_insomnia = "Insomnia",
  drsp_breast_tender = "Breast tenderness", drsp_bloated = "Bloating", drsp_headache = "Headache",
  drsp_joint_pain = "Joint pain", drsp_impair_work = "Impairment: work", drsp_impair_social = "Impairment: social",
  drsp_impair_relations = "Impairment: relationships",
  ls_1 = "More sensitive to light", ls_2 = "Indoor lighting uncomfortable", ls_3 = "Bright light worsened mood",
  ls_4 = "Light triggered headache/eye strain", ls_5 = "Avoided or adjusted lighting")
E2_PACTS_BREAKS <- c(-1, -0.5, 0, 0.5, 1)
E2_PACTS_LABELS <- c("-1\nmenses onset", "-0.5", "0\nest. ovulation", "0.5", "1\nnext onset")

e2_theme <- function(base = 10) ggplot2::theme_minimal(base_size = base) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"),
                 legend.position = "bottom")
e2_save <- function(p, file, w, h) {
  ggplot2::ggsave(paste0(file, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
  tryCatch(ggplot2::ggsave(paste0(file, ".pdf"), p, width = w, height = h), error = function(e) NULL)
  invisible(file)
}

# item means per PACTS bin with 95% CI
e2_bin_means <- function(dd, items, min_n = 5) {
  d <- dd[!is.na(dd$cycle_bin_ov), c("cycle_bin_ov", intersect(items, names(dd)))]
  e2_bind(lapply(intersect(items, names(dd)), function(it) {
    x <- d[!is.na(d[[it]]), c("cycle_bin_ov", it)]; if (!nrow(x)) return(NULL)
    s <- stats::aggregate(x[[it]], by = list(bin = x$cycle_bin_ov), FUN = function(z) c(m = mean(z), s = stats::sd(z), n = length(z)))
    s <- data.frame(bin = s$bin, mean = s$x[, "m"], sd = s$x[, "s"], n = s$x[, "n"])
    s$se <- s$sd / sqrt(s$n); s$lwr <- s$mean - 1.96 * s$se; s$upr <- s$mean + 1.96 * s$se
    s$low_n <- s$n < min_n; s$item <- it; s
  }))
}
e2_fig_items <- function(bm, file, ncol, ylim, colour_by_domain = FALSE, title = NULL) {
  if (is.null(bm) || !nrow(bm)) return(NULL)
  ord <- c(unlist(E2_DOMAINS, use.names = FALSE), paste0("ls_", 1:5))
  ord <- ord[ord %in% bm$item]
  bm$label <- factor(E2_ITEM_LABELS[bm$item], levels = E2_ITEM_LABELS[ord])
  bm$domain <- if (colour_by_domain) factor(vapply(bm$item, function(i) names(E2_DOMAINS)[vapply(E2_DOMAINS, function(d) i %in% d, TRUE)][1], ""),
                                             levels = names(E2_DOMAINS)) else factor("all")
  cols <- if (colour_by_domain) E2_DOMAIN_COLS else c(all = E2_OI[["blue"]])
  bm$lwr <- pmax(bm$lwr, ylim[1]); bm$upr <- pmin(bm$upr, ylim[2])
  p <- ggplot2::ggplot(bm, ggplot2::aes(bin, mean, colour = domain, fill = domain)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), alpha = 0.2, colour = NA) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(shape = low_n), size = 1.2) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c(`FALSE` = "n >= 5", `TRUE` = "n < 5 (unstable)"), name = NULL) +
    ggplot2::geom_vline(xintercept = 14.5, linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(values = cols, guide = if (colour_by_domain) "legend" else "none", name = NULL) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_x_continuous(breaks = c(1, 7, 14.5, 21, 28), labels = c("1", "7", "ov.", "21", "28"), limits = c(1, 28)) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::facet_wrap(~label, ncol = ncol, labeller = ggplot2::label_wrap_gen(22)) +
    ggplot2::labs(x = "PACTS bin (28 bins; dashed line = estimated ovulation)", y = "Mean score (1-6) with 95% CI", title = title) +
    e2_theme(9)
  n_rows <- ceiling(length(unique(bm$label)) / ncol)
  e2_save(p, file, 2.1 * ncol + 0.5, 1.7 * n_rows + 1.2)
}

# exclusion reasons (EXCLUSIONS from 06_pacts.R) as English labels with counts
e2_exclusion_text <- function(ex) {
  if (is.null(ex) || !nrow(ex)) return(character())
  lab <- ifelse(grepl("abbruch", ex$reason, ignore.case = TRUE), "discontinued shortly after enrolment",
                ifelse(grepl("BMI", ex$reason), "BMI above the inclusion range", ex$reason))
  tb <- table(factor(lab, levels = unique(lab)))
  paste0(names(tb), " (n = ", as.integer(tb), ")")
}
