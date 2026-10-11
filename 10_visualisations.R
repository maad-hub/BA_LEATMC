# 10_visualisations.R - Figures of the initial analysis; fig_corr.png is Figure B1
LOW_N_ALPHA    <- c(`FALSE`=0.35, `TRUE`=0.12)
LOW_N_LINETYPE <- c(`FALSE`="solid", `TRUE`="dashed")
LOW_N_LABELS   <- c(`FALSE`=paste0("n \u2265 ", MIN_BIN_N), `TRUE`=paste0("n < ", MIN_BIN_N))

PHASE_COLS  <- c(menstrual="#c0392b",follicular="#2980b9",
                  ovulatory="#27ae60",luteal="#e67e22")
OV_VLINE    <- geom_vline(xintercept=0, linetype="dashed",colour="grey50",linewidth=0.6)
OV_VLINE_BIN <- geom_vline(xintercept=N_BINS/2, linetype="dashed",colour="grey50",linewidth=0.6)

# ---- GAMM smooth helper ----
gamm_smooth <- function(model, n=300) {
  if (is.null(model)) return(NULL)
  rng <- range(daily_data$pacts_ov, na.rm=TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) return(NULL)
  nd  <- tibble(pacts_ov=seq(rng[1],rng[2],length.out=n),
                Id=levels(factor(daily_data$Id))[1])
  p   <- predict(model, newdata=nd, exclude="s(Id)", se.fit=TRUE)
  nd |> mutate(fit=p$fit, lwr=fit-1.96*p$se.fit, upr=fit+1.96*p$se.fit)
}

# ---- GAMM plot helper ----
plot_gamm <- function(smooth_df, colour, title, ylab) {
  if (is.null(smooth_df))
    return(ggplot()+labs(title=paste(title,"(not available)"))+theme_minimal())
  ggplot(smooth_df, aes(x=pacts_ov,y=fit)) +
    geom_ribbon(aes(ymin=lwr,ymax=upr), fill=colour, alpha=0.2) +
    geom_line(colour=colour, linewidth=1.3) + OV_VLINE +
    annotate("text",x=-0.05,y=Inf,label="\u2190 Follicular",hjust=1,vjust=1.5,size=3.5,colour="grey40") +
    annotate("text",x= 0.05,y=Inf,label="Luteal phase \u2192",   hjust=0,vjust=1.5,size=3.5,colour="grey40") +
    labs(title=title, subtitle="GAMM \u00b1 95% CI  |  0 = Ovulation",
         x="PACTS (ovulation-centered)", y=ylab) +
    theme_minimal(base_size=13)
}

# ---- GAMM smooths ----
fig1 <- plot_gamm(gamm_smooth(gamm_TAT),   "steelblue",
                   "Light Exposure Across the Menstrual Cycle", "TAT 250 lux mEDI (h/day)")
fig4 <- plot_gamm(gamm_smooth(gamm_lsens), "#e74c3c",
                   "Light Sensitivity Across the Cycle", "Light Sensitivity Index (1-6)")
figA <- plot_gamm(gamm_smooth(gamm_csymp), "#c0392b",
                   "Symptom Severity Across the Cycle", "Cycle Symptom Index (1-6)")

ggsave(file.path(FIG_DIR,"fig1_GAMM_TAT250.png"), fig1, width=9, height=5, dpi=300)
ggsave(file.path(FIG_DIR,"fig4_GAMM_lsens.png"),  fig4, width=9, height=5, dpi=300)
ggsave(file.path(FIG_DIR,"figA_GAMM_symptoms.png"),figA, width=9, height=5, dpi=300)

# ---- Individual and group trajectories of TAT250 ----
fig2 <- ggplot() +
  geom_line(data=filter(daily_data,!is.na(cycle_bin_ov),!is.na(TAT_250_hr)),
            aes(x=cycle_bin_ov,y=TAT_250_hr,group=Id),
            alpha=0.2, colour="steelblue", linewidth=0.4) +
  geom_ribbon(data=bin_summary,
              aes(x=cycle_bin,ymin=TAT_250_mean-TAT_250_ci,ymax=TAT_250_mean+TAT_250_ci,
                  alpha=low_n),
              fill="steelblue") +
  geom_line(data=bin_summary, aes(x=cycle_bin,y=TAT_250_mean,linetype=low_n),
            colour="steelblue", linewidth=1.3) +
  scale_alpha_manual(values=LOW_N_ALPHA, guide="none") +
  scale_linetype_manual(values=LOW_N_LINETYPE, labels=LOW_N_LABELS, name=NULL) +
  OV_VLINE_BIN +
  labs(title="Light Exposure Trajectories",
       subtitle=paste0("Dashed/lighter = bin with n < ", MIN_BIN_N, " (unstable mean)"),
       x="PACTS bin",y="TAT 250 lux mEDI (h/day)") +
  theme_minimal(base_size=13) + theme(legend.position="bottom")
ggsave(file.path(FIG_DIR,"fig2_spaghetti_TAT.png"), fig2, width=10, height=5, dpi=300)

# ---- TAT250 by four-phase label ----
fig3 <- daily_data |> filter(!is.na(cycle_phase_4)) |>
  ggplot(aes(x=cycle_phase_4,y=TAT_250_hr,fill=cycle_phase_4)) +
  geom_boxplot(outlier.shape=NA,alpha=0.8,width=0.5) +
  geom_beeswarm(alpha=0.25,size=0.8,cex=0.6) +
  scale_fill_manual(values=PHASE_COLS,guide="none") +
  labs(title="TAT 250 lux by Cycle Phase",x=NULL,y="h/day") +
  theme_minimal(base_size=13)
ggsave(file.path(FIG_DIR,"fig3_boxplot_TAT.png"), fig3, width=7, height=5, dpi=300)

# ---- Combined panel ----
fig5 <- (fig1+fig4)/(fig2+fig3) +
  plot_annotation(title="LEAMC \u2013 Light Across the Menstrual Cycle")
ggsave(file.path(FIG_DIR,"fig5_overview.png"), fig5, width=16, height=10, dpi=300)

# ---- Example participant (LightLogR gg_days) ----
if (nrow(light_trimmed) > 0) {
  example_id <- unique(light_trimmed$Id)[1]
  example_data <- light_trimmed |>
    filter(Id == example_id) |>
    filter(Date <= min(Date) + days(9))  # first 10 days

  ggsave(file.path(FIG_DIR,"fig6_gg_days.png"),
         example_data |> gg_days(y.axis=MEDI, y.axis.label="Melanopic EDI (lux)"),
         width=14, height=8, dpi=300)
}

# ---- DRSP items by domain ----
DRSP_DOMAINS <- list(
  "Mood & Affect"           = c("drsp_hopeless","drsp_worthless","drsp_anxious",
                                 "drsp_mood_swings","drsp_rejection","drsp_angry","drsp_conflicts"),
  "Cognitive & Behavioural" = c("drsp_less_interest","drsp_concentration",
                                 "drsp_overwhelmed","drsp_out_of_control"),
  "Physical"                = c("drsp_fatigue","drsp_appetite","drsp_food_cravings",
                                 "drsp_hypersomnia","drsp_insomnia","drsp_breast_tender",
                                 "drsp_bloated","drsp_headache","drsp_joint_pain"),
  "Impairment"              = c("drsp_impair_work","drsp_impair_social","drsp_impair_relations")
)
DRSP_LABELS <- c(
  drsp_hopeless="Hopeless", drsp_worthless="Worthless", drsp_anxious="Anxious",
  drsp_mood_swings="Mood swings", drsp_rejection="Rejection sens.",
  drsp_angry="Angry/irritable", drsp_conflicts="Conflicts",
  drsp_less_interest="Less interest", drsp_concentration="Concentration",
  drsp_overwhelmed="Overwhelmed", drsp_out_of_control="Out of control",
  drsp_fatigue="Fatigue", drsp_appetite="Appetite", drsp_food_cravings="Food cravings",
  drsp_hypersomnia="Hypersomnia", drsp_insomnia="Insomnia",
  drsp_breast_tender="Breast tender.", drsp_bloated="Bloating",
  drsp_headache="Headache", drsp_joint_pain="Joint pain",
  drsp_impair_work="Work impair.", drsp_impair_social="Social impair.",
  drsp_impair_relations="Relat. impair."
)
all_drsp <- unlist(DRSP_DOMAINS, use.names=FALSE)
DOMAIN_COLS <- c("Mood & Affect"="#c0392b","Cognitive & Behavioural"="#8e44ad",
                  "Physical"="#e67e22","Impairment"="#2c3e50")

drsp_bin <- daily_data |>
  filter(!is.na(cycle_bin_ov)) |>
  select(Id, cycle_bin_ov, all_of(all_drsp)) |>
  pivot_longer(all_of(all_drsp), names_to="item", values_to="score") |>
  mutate(label  = factor(DRSP_LABELS[item], levels=DRSP_LABELS[all_drsp]),
         domain = map_chr(item, ~names(DRSP_DOMAINS)[map_lgl(DRSP_DOMAINS, \(d) .x %in% d)]) |>
                  factor(levels=names(DRSP_DOMAINS))) |>
  filter(!is.na(score)) |>
  group_by(label, domain, cycle_bin_ov) |>
  summarise(mean=mean(score,na.rm=TRUE), se=sd(score,na.rm=TRUE)/sqrt(n()), n=n(), .groups="drop") |>
  mutate(low_n = n < MIN_BIN_N)

figB <- ggplot(drsp_bin, aes(x=cycle_bin_ov,y=mean,colour=domain,fill=domain)) +
  geom_ribbon(aes(ymin=mean-1.96*se,ymax=mean+1.96*se,alpha=low_n),colour=NA) +
  geom_line(aes(linetype=low_n),linewidth=0.8) + OV_VLINE_BIN +
  scale_colour_manual(values=DOMAIN_COLS,guide="none") +
  scale_fill_manual(values=DOMAIN_COLS,guide="none") +
  scale_alpha_manual(values=LOW_N_ALPHA, guide="none") +
  scale_linetype_manual(values=LOW_N_LINETYPE, labels=LOW_N_LABELS, name=NULL) +
  facet_wrap(~label, ncol=5, labeller=label_wrap_gen(16)) +
  scale_y_continuous(limits=c(1,6),breaks=c(1,3,6)) +
  labs(title="DRSP Items Across the Menstrual Cycle",
       subtitle=paste0("Group mean \u00b1 95% CI | Dashed = n < ", MIN_BIN_N),
       x="PACTS bin", y="Mean score (1-6)") +
  theme_minimal(base_size=9) +
  theme(strip.text=element_text(size=7.5), panel.grid.minor=element_blank(),
        legend.position="bottom")
ggsave(file.path(FIG_DIR,"figB_DRSP_domain.png"), figB, width=16, height=12, dpi=300)

# ---- Light-sensitivity items ----
LS_LABELS <- c(ls_1="More sensitive",ls_2="Indoor uncomfortable",
                ls_3="Bright→mood",ls_4="Headache/eye",ls_5="Avoided light")
ls_bin <- daily_data |> filter(!is.na(cycle_bin_ov)) |>
  select(Id,cycle_bin_ov,ls_1:ls_5) |>
  pivot_longer(ls_1:ls_5,names_to="item",values_to="score") |>
  mutate(label=factor(LS_LABELS[item],levels=LS_LABELS)) |>
  filter(!is.na(score)) |>
  group_by(label,cycle_bin_ov) |>
  summarise(mean=mean(score,na.rm=TRUE),se=sd(score,na.rm=TRUE)/sqrt(n()),n=n(),.groups="drop") |>
  mutate(low_n = n < MIN_BIN_N)

figC <- ggplot(ls_bin,aes(x=cycle_bin_ov,y=mean)) +
  geom_ribbon(aes(ymin=mean-1.96*se,ymax=mean+1.96*se,alpha=low_n),fill="#e74c3c") +
  geom_line(aes(linetype=low_n),colour="#e74c3c",linewidth=0.9) + OV_VLINE_BIN +
  scale_alpha_manual(values=LOW_N_ALPHA, guide="none") +
  scale_linetype_manual(values=LOW_N_LINETYPE, labels=LOW_N_LABELS, name=NULL) +
  facet_wrap(~label,ncol=3,labeller=label_wrap_gen(22)) +
  scale_y_continuous(limits=c(1,6),breaks=c(1,3,6)) +
  labs(title="Light Sensitivity Items",
       subtitle=paste0("Dashed = n < ", MIN_BIN_N),
       x="PACTS bin",y="Mean score (1-6)") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(), legend.position="bottom")
ggsave(file.path(FIG_DIR,"figC_LS_items.png"), figC, width=12, height=7, dpi=300)

# ---- Menses check ----
figD <- daily_data |> filter(!is.na(cycle_bin_ov)) |>
  select(cycle_bin_ov, menstrual_flow, currently_menstruating) |>
  pivot_longer(-cycle_bin_ov) |>
  mutate(name=recode(name,menstrual_flow="Flow (0-3)",
                     currently_menstruating="Currently menstruating")) |>
  filter(!is.na(value)) |>
  group_by(name,cycle_bin_ov) |> summarise(m=mean(value,na.rm=TRUE),.groups="drop") |>
  ggplot(aes(x=cycle_bin_ov,y=m,colour=name)) +
  geom_line(linewidth=1) + OV_VLINE_BIN +
  facet_wrap(~name,ncol=1,scales="free_y") +
  scale_colour_manual(values=c("Flow (0-3)"="#c0392b",
                                "Currently menstruating"="#8e44ad"),guide="none") +
  labs(title="Menses Sanity Check (peak on RIGHT = correct: cycle end)",
       x="PACTS bin",y=NULL) + theme_minimal(base_size=12)
ggsave(file.path(FIG_DIR,"figD_menses_sanity.png"), figD, width=9, height=6, dpi=300)

# ---- Correlation matrix (Figure B1) ----
if (all(is.na(corr_mat))) {
  message("[10_visualisations] fig_corr.png uebersprungen — corr_mat enthaelt ",
          "keine gueltigen Werte. Siehe Warnung/Tabelle aus 09_inferential.R, ",
          "welche Variable(n) durchgehend NA sind.")
} else {
  CORR_LABELS <- c(
    TAT_250_hr             = "TAT 250 (h)",
    MEDI_mean              = "Mean MEDI",
    MEDI_morning           = "Morning MEDI",
    MEDI_evening           = "Evening MEDI",
    light_sensitivity_idx  = "Light Sensitivity",
    cycle_symptom_idx      = "Cycle Symptoms",
    TST_min                = "Total Sleep Time",
    sleep_efficiency       = "Sleep Efficiency",
    WASO_min               = "WASO",
    drsp_fatigue           = "Fatigue",
    drsp_hypersomnia       = "Hypersomnia",
    drsp_insomnia          = "Insomnia"
  )
  corr_mat_plot <- corr_mat
  dimnames(corr_mat_plot) <- list(
    CORR_LABELS[rownames(corr_mat)],
    CORR_LABELS[colnames(corr_mat)]
  )

  png(file.path(FIG_DIR,"fig_corr.png"),width=2400,height=2400,res=300)
  corrplot(corr_mat_plot,method="color",type="upper",tl.cex=0.75,
           addCoef.col="black",number.cex=0.55,
           col=colorRampPalette(c("#2166ac","white","#d6604d"))(200),
           title="Spearman Correlations",mar=c(0,0,2,0))
  dev.off()
}

# ---- Missingness heatmap ----

missing_long <- daily_data |> filter(!is.na(cycle_phase_2)) |>
  mutate(light=as.integer(!is.na(TAT_250_hr)),
         sleep=as.integer(!is.na(TST_min)),
         survey=as.integer(!is.na(light_sensitivity_idx))) |>
  select(Id, cycle_bin_ov, light, sleep, survey) |>
  pivot_longer(c(light,sleep,survey), names_to="source", values_to="value") |>
  mutate(source=recode(source,light="ActLumus",sleep="Actigraph",survey="Daily questionnaire")) |>
  group_by(Id, cycle_bin_ov, source) |>
  summarise(value = mean(value), .groups = "drop")

missing_complete <- missing_long |>
  complete(Id, cycle_bin_ov, source, fill=list(value=NA_real_)) |>
  mutate(status = case_when(
    is.na(value) ~ "outside",
    value == 1   ~ "present",
    value == 0   ~ "missing",
    TRUE         ~ "partial"
  ))

ggsave(file.path(FIG_DIR,"fig_missing.png"),
  ggplot(missing_complete, aes(x=cycle_bin_ov,y=Id,fill=status)) +
    geom_tile(colour="white",linewidth=0.2) +
    scale_fill_manual(
      values=c(present="#2ecc71", partial="#f5b041", missing="#e74c3c", outside="grey90"),
      breaks=c("present","partial","missing","outside"),
      labels=c(present="Present", partial="Partly present (bin spans >1 day)",
               missing="Missing", outside="Outside valid cycle window"),
      name=NULL) +
    facet_wrap(~source,ncol=1) +
    labs(title="Data Availability",x="PACTS bin",y=NULL) +
    theme_minimal(base_size=10)+theme(legend.position="bottom"),
  width=12, height=10, dpi=300)

cat("[10_visualisations] Alle Abbildungen gespeichert in:", FIG_DIR, "\n")
