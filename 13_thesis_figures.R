# 13_thesis_figures.R - Figures 1-7 of the thesis (300 dpi PNG and PDF) in OUTPUT_DIR/final_models/figures/
# Figure B1 (pooled Spearman correlations) is fig_corr.png from 10_visualisations.R.

FIG_FINAL_DIR <- file.path(FINAL_DIR, "figures")
dir.create(FIG_FINAL_DIR, recursive = TRUE, showWarnings = FALSE)
FF <- "DejaVu Sans"   # font with subscript digits and the >= sign (Figures 1 and 2)

# ---- Figure 1: conceptual model ----
bx <- data.frame(x = c(0, 3.2, 3.2, 6.4, 3.2), y = c(2, 3.3, 1.9, 1.9, 0.3),
                 lab = c("Menstrual cycle\n(estradiol, progesterone)", "Cycle-related symptoms\n(P3)",
                         "Behaviour: sleep-wake timing,\nactivity, time outdoors\n(P4)", "Personal light\nexposure\n(P1)",
                         "Subjective light\nsensitivity (P2)"))
ar <- data.frame(x = c(1.25, 3.2, 4.55, 1.25), y = c(2.3, 2.95, 1.9, 1.6),
                 xend = c(2.05, 3.2, 5.55, 2.15), yend = c(3.1, 2.45, 1.9, 0.55))
pm <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = ar, ggplot2::aes(x, y, xend = xend, yend = yend),
                        arrow = grid::arrow(length = grid::unit(0.16, "cm"), type = "closed"), colour = "grey35", linewidth = 0.5) +
  ggplot2::geom_label(data = bx, ggplot2::aes(x, y, label = lab), size = 3.1, family = FF,
                      label.padding = grid::unit(0.45, "lines"), fill = c("grey96", "white", "white", "white", "white")) +
  ggplot2::scale_x_continuous(limits = c(-1.3, 7.6)) + ggplot2::scale_y_continuous(limits = c(-0.4, 3.9)) +
  ggplot2::theme_void(base_size = 10)
e2_save(pm, file.path(FIG_FINAL_DIR, "fig1_concept_model"), 7.5, 3.6)

# ---- Figure 2: participant flow ----
cf <- final_models$A$confound
nid <- function(o) { x <- table3$n_id[table3$variant == "A" & table3$outcome == o]; if (length(x)) x[1] else NA }
n_ex <- nrow(EXCLUSIONS)
boxes <- data.frame(x = 0, y = c(8, 6.4, 4.6, 2.9), lab = c(
  paste0("Enrolled (n = ", cf$n_participants + n_ex, ")"),
  paste0("Analysed (n = ", cf$n_participants, ")"),
  paste0("Cycle phase assigned (n = ", cf$n_with_phase, ")\n", cf$n_ids_closed_cycles,
         " with a closed cycle of 21-40 days (PACTS)\n", cf$n_ids_counting_only,
         " without such a cycle, by counting from a prospective onset:\nluteal days backward (n = ",
         cf$n_counting_only_backward, "), follicular days forward (n = ", cf$n_counting_only_forward, ")"),
  paste0("Days with a cycle phase: ", cf$n_days_with_phase, " (follicular ", cf$n_days_follicular, ", luteal ",
         cf$n_days_luteal, ")\nQuestionnaire days with a cycle phase: ", cf$n_questionnaire_days_with_phase, " of ",
         cf$n_questionnaire_days, " (", round(100 * cf$n_questionnaire_days_with_phase / cf$n_questionnaire_days), "%)")))
side <- data.frame(x = 2.3, y = 7.2, lab = paste(c(paste0("Excluded (n = ", n_ex, ")"), e2_exclusion_text(EXCLUSIONS)), collapse = "\n"))
samp <- data.frame(x = c(-2.2, 0, 2.2), y = 0.9, lab = c(
  paste0("P1 light logger\n(TAT₂₅₀): n = ", nid("TAT_250_hr")),
  paste0("P2/P3 questionnaires\nLSI: n = ", nid("light_sensitivity_idx"), "; CSI: n = ", nid("cycle_symptom_idx")),
  paste0("P4 actigraphy\nRA/TST/SE: n = ", nid("RA_act"), "\nIS/IV (≥ 5 valid days\nin both phases): n = ",
         final_models$A$npcra_n$n_ok_ids)))
arr <- grid::arrow(length = grid::unit(0.14, "cm"))
pf <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = data.frame(y0 = c(7.75, 6.15, 3.85), y1 = c(6.65, 5.35, 3.2)),
                        ggplot2::aes(x = 0, xend = 0, y = y0, yend = y1), arrow = arr, colour = "grey40") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1.0, y = 7.2, yend = 7.2), arrow = arr, colour = "grey40") +
  ggplot2::geom_segment(data = data.frame(x1 = c(-2.2, 0, 2.2)), ggplot2::aes(x = 0, xend = x1, y = 2.55, yend = 1.55),
                        arrow = arr, colour = "grey40") +
  ggplot2::geom_label(data = boxes, ggplot2::aes(x, y, label = lab), size = 3, family = FF, label.padding = grid::unit(0.4, "lines")) +
  ggplot2::geom_label(data = side, ggplot2::aes(x, y, label = lab), size = 2.8, family = FF, label.padding = grid::unit(0.35, "lines"), fill = "grey96") +
  ggplot2::geom_label(data = samp, ggplot2::aes(x, y, label = lab), size = 2.7, family = FF, label.padding = grid::unit(0.35, "lines")) +
  ggplot2::scale_x_continuous(limits = c(-3.4, 3.7)) + ggplot2::scale_y_continuous(limits = c(0.2, 8.4)) +
  ggplot2::theme_void(base_size = 10)
e2_save(pf, file.path(FIG_FINAL_DIR, "fig2_flow"), 8, 6.2)

# ---- Figure 3: data availability per participant across PACTS bins (extended assignment) ----
OKABE <- c(Present = "#0072B2", `Partly present` = "#56B4E9", Missing = "#D55E00",
           `No day mapped to bin` = "#C8C8C8", `Outside recording window` = "#F0F0F0")
j <- daily_data_06 |>
  dplyr::transmute(Id = as.character(Id), Date = as.Date(Date),
                   has_light = !is.na(TAT_250_hr), has_sleep = !is.na(TST_min), has_quest = !is.na(light_sensitivity_idx)) |>
  dplyr::left_join(phase_tabs$A |> dplyr::filter(in_days) |> dplyr::select(Id, Date, cycle_bin_ov, count_method),
                   by = c("Id", "Date"))
j <- j[!is.na(j$cycle_bin_ov), ]
long <- rbind(
  data.frame(Id = j$Id, bin = j$cycle_bin_ov, source = "Light logger (ActLumus)", v = as.integer(j$has_light)),
  data.frame(Id = j$Id, bin = j$cycle_bin_ov, source = "Actigraphy (ActiGraph)", v = as.integer(j$has_sleep)),
  data.frame(Id = j$Id, bin = j$cycle_bin_ov, source = "Daily questionnaire", v = as.integer(j$has_quest)))
agg <- stats::aggregate(v ~ Id + bin + source, data = long, FUN = mean)
grid_av <- expand.grid(Id = sort(unique(agg$Id)), bin = seq_len(28), source = unique(long$source), stringsAsFactors = FALSE)
g <- merge(grid_av, agg, by = c("Id", "bin", "source"), all.x = TRUE)
rng <- stats::aggregate(bin ~ Id, data = long, FUN = function(x) c(min(x), max(x)))
g$bmin <- rng$bin[match(g$Id, rng$Id), 1]; g$bmax <- rng$bin[match(g$Id, rng$Id), 2]
g$status <- ifelse(is.na(g$v), ifelse(g$bin < g$bmin | g$bin > g$bmax, "Outside recording window", "No day mapped to bin"),
                   ifelse(g$v == 1, "Present", ifelse(g$v == 0, "Missing", "Partly present")))
g$status <- factor(g$status, levels = names(OKABE))
g$x <- -1 + (g$bin - 0.5) * 2 / 28
g$source <- factor(g$source, levels = c("Light logger (ActLumus)", "Actigraphy (ActiGraph)", "Daily questionnaire"))
mk <- j |> dplyr::filter(count_method %in% c("forward_imputed", "backward_imputed")) |>
  dplyr::distinct(Id, bin = cycle_bin_ov, count_method)
mk <- merge(mk, data.frame(source = levels(g$source)), by = NULL)
mk$source <- factor(mk$source, levels = levels(g$source))
mk$x <- -1 + (mk$bin - 0.5) * 2 / 28
mk$marker <- factor(ifelse(mk$count_method == "forward_imputed", "Forward-counted (follicular)", "Backward-counted (luteal)"),
                    levels = c("Forward-counted (follicular)", "Backward-counted (luteal)"))
p <- ggplot(g, aes(x = x, y = Id)) +
  geom_tile(aes(fill = status), width = 2 / 28 * 0.94, height = 0.86) +
  geom_vline(xintercept = 0, colour = "grey35", linewidth = 0.3) +
  scale_fill_manual(values = OKABE, drop = FALSE, name = NULL) +
  scale_y_discrete(limits = rev(sort(unique(g$Id)))) +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1", "-0.5", "0", "0.5", "1"), expand = c(0.005, 0)) +
  facet_wrap(~source, ncol = 1) +
  labs(x = "Cycle time (PACTS; -1 = menses onset, 0 = estimated ovulation, +1 = next menses onset)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.box = "vertical",
        strip.text = element_text(size = 8.5, hjust = 0), axis.text.y = element_text(size = 6.5),
        legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 2, order = 1))
if (nrow(mk)) p <- p + geom_point(data = mk, aes(x = x, y = Id, shape = marker), size = 0.9, colour = "#1A1A1A", inherit.aes = FALSE) +
  scale_shape_manual(values = c(`Forward-counted (follicular)` = 16, `Backward-counted (luteal)` = 2), drop = TRUE, name = NULL) +
  guides(shape = guide_legend(order = 2, override.aes = list(size = 2)))
h <- 7 + 0.42 * length(unique(g$Id)) * 3
ggsave(file.path(FIG_FINAL_DIR, "fig3_data_availability.png"), p, width = 16, height = min(h, 26), units = "cm", dpi = 300, bg = "white")
ggsave(file.path(FIG_FINAL_DIR, "fig3_data_availability.pdf"), p, width = 16, height = min(h, 26), units = "cm")

# ---- Figure 4: model-estimated trajectories (GAMM with AR(1); pointwise and simultaneous 95% bands) ----
cv <- final_models$A$curves
cv <- cv[cv$model == "AR1" & cv$outcome %in% c("TAT_250_hr", "light_sensitivity_idx", "cycle_symptom_idx"), ]
lv <- c(TAT_250_hr = 'bold(TAT[250]~"(h/day)")', light_sensitivity_idx = 'bold("Light Sensitivity Index")',
        cycle_symptom_idx = 'bold("Cycle Symptom Index")')
cv$panel <- factor(lv[cv$outcome], levels = lv)
p <- ggplot2::ggplot(cv, ggplot2::aes(pacts_ov, fit)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr_sim, ymax = upr_sim), fill = E2_OI[["sky"]], alpha = 0.35) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr), fill = E2_OI[["blue"]], alpha = 0.35) +
  ggplot2::geom_line(colour = E2_OI[["blue"]], linewidth = 0.9) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_y", labeller = ggplot2::label_parsed) +
  ggplot2::scale_x_continuous(breaks = E2_PACTS_BREAKS, labels = E2_PACTS_LABELS, limits = c(-1, 1)) +
  ggplot2::labs(x = "Ovulation-centred cycle time (PACTS)", y = "Model-estimated value") + e2_theme()
e2_save(p, file.path(FIG_FINAL_DIR, "fig4_trajectories"), 6.5, 2.4 * 3 + 0.6)

# ---- Figure 5: repeated-measures correlations (revision.R R25; * = Holm p < .05) ----
rm_pairs <- utils::read.csv(file.path(REV_DIR, "R25_rmcorr_pairs.csv"))
lab_fun <- function(x) { y <- as.list(x); y[x == "TAT250"] <- list(quote(TAT[250])); as.expression(y) }
lvls <- c("LSI", "CSI", "TAT250", "Centroid", "RA", "TST", "SE")
d <- rbind(data.frame(a = rm_pairs$var1, b = rm_pairs$var2, r = rm_pairs$r_rm, sig = rm_pairs$p_holm < .05),
           data.frame(a = rm_pairs$var2, b = rm_pairs$var1, r = rm_pairs$r_rm, sig = rm_pairs$p_holm < .05))
d$a <- factor(d$a, levels = lvls); d$b <- factor(d$b, levels = rev(lvls))
d <- d[as.integer(d$a) < (length(lvls) + 1 - as.integer(d$b)), ]
d$lab <- paste0(sub("^(-?)0\\.", "\\1.", formatC(d$r, format = "f", digits = 2)), ifelse(d$sig %in% TRUE, "*", ""))
p <- ggplot2::ggplot(d, ggplot2::aes(a, b, fill = r)) + ggplot2::geom_tile(colour = "white") +
  ggplot2::geom_text(ggplot2::aes(label = lab), size = 3.2) +
  ggplot2::scale_fill_gradient2(low = E2_OI[["blue"]], mid = "white", high = E2_OI[["vermillion"]], midpoint = 0,
                                limits = c(-1, 1), name = expression(r[rm])) +
  ggplot2::scale_x_discrete(labels = lab_fun) + ggplot2::scale_y_discrete(labels = lab_fun) +
  ggplot2::labs(x = NULL, y = NULL) +
  e2_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "right")
e2_save(p, file.path(FIG_FINAL_DIR, "fig5_rmcorr"), 6, 4.8)

# ---- Figures 6 and 7: item means per PACTS bin (common y-axis) ----
bl <- e2_bin_means(final_dd$A, AC_LS_ITEMS)
bd <- e2_bin_means(final_dd$A, c(AC_DRSP_ITEMS, AC_IMPAIR))
yl <- c(1, max(2, ceiling(2 * max(c(bl$upr[!bl$low_n], bd$upr[!bd$low_n]), na.rm = TRUE)) / 2))
e2_fig_items(bl, file.path(FIG_FINAL_DIR, "fig6_LS_items"), 3, yl)
e2_fig_items(bd, file.path(FIG_FINAL_DIR, "fig7_DRSP_items"), 5, yl, colour_by_domain = TRUE)
cat("[13_thesis_figures] figures in", FIG_FINAL_DIR, "\n")
