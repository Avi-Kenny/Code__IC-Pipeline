# Choose analysis
# "Janssen" "Moderna" "AMP" "AZD1222" "Janssen (partA)" "Profiscov"
# "HVTN 705 (primary)" "HVTN 705 (all)" "RV144" "HVTN 705 (second)"
# "HVTN 705 (compare RV144)"
analysis <- "Janssen (partA)"
if (analysis=="Moderna") {
  inds <- c(5,7,9)
} else if (analysis=="Janssen (partA)") {
  inds <- c(1:64)
  # inds <- c(1:58)
} else if (analysis=="HVTN 705 (all)") {
  inds <- c(9,15:22,25:30,33:37,39)
}

df_med <- data.frame(
  "i"=integer(),
  "nde_est"=double(),"nde_se"=double(),"nde_lo"=double(),"nde_up"=double(),
  "nie_est"=double(),"nie_se"=double(),"nie_lo"=double(),"nie_up"=double(),
  "pm_est"=double(),"pm_se"=double(),"pm_lo"=double(),"pm_up"=double()
)

for (i in inds) {
  filename = paste0("rds/", analysis, " objs/ests_med_", i, ".rds")
  if (!file.exists(filename)) {
    next
  }
  ests <- readRDS(filename)
  e_nde <- as.list(ests[ests$effect=="NDE",][2:5])
  e_nie <- as.list(ests[ests$effect=="NIE",][2:5])
  e_pm <- as.list(ests[ests$effect=="PM",][2:5])
  e_all <- c(e_nde, e_nie, e_pm)
  e_all <- lapply(e_all, function(x) { round(x,3) })
  df_med[nrow(df_med)+1,] <- c(list(i=i), e_all)

}
write.table(
  df_med,
  file = paste0("mediation_results - ",analysis,".csv"),
  sep = ",",
  row.names = F
)

# !!!!! Forest plot
if (F) {

  df_med2 <- df_med
  df_merge <- cfg2$map
  df_merge$i <- c(1:nrow(df_merge))
  df_merge %<>% subset(select=c("i", "marker", "dataset", "cr2_COR"))
  df_med2 %<>% left_join(df_merge, by="i")
  df_med2 %<>% dplyr::filter(dataset %in% c(1,4,7,10))
  df_med2 %<>% dplyr::filter(marker<=4)
  df_med2 %<>% dplyr::rename("severity"=cr2_COR)
  levs_severity <- c("Moderate to Severe", "Severe", "Moderate")
  levs_marker <- c("Anti Spike IgG", "Anti RBD IgG", "Pseudovirus-nAb ID50",
                   "Phagocytic Score")
  levs_dataset <- c("Pooled", "LA", "SA", "NA")

  df_med2 %<>% dplyr::mutate(
    marker = factor(dplyr::case_when(
      marker==1 ~ levs_marker[1],
      marker==2 ~ levs_marker[2],
      marker==3 ~ levs_marker[3],
      marker==4 ~ levs_marker[4]
    ), levels=rev(levs_marker)),
    dataset = factor(dplyr::case_when(
      dataset==1 ~ levs_dataset[1],
      dataset==4 ~ levs_dataset[4],
      dataset==7 ~ levs_dataset[2],
      dataset==10 ~ levs_dataset[3]
    ), levels=levs_dataset),
    severity = factor(dplyr::case_when(
      severity==1 ~ levs_severity[1],
      severity==2 ~ levs_severity[2],
      severity==3 ~ levs_severity[3]
    ), levels=levs_severity)
  )

  forest_plot <- function(df_med, title, var_estimate, var_ci_lo, var_ci_up,
                          var_y, facet_var_x, facet_var_y, xlim) {

    # Extract data
    df_plot <- subset(df_med,select="i")
    df_plot$var_estimate <- df_med[[var_estimate]]
    df_plot$var_ci_lo <- df_med[[var_ci_lo]]
    df_plot$var_ci_up <- df_med[[var_ci_up]]
    df_plot$var_y <- df_med[[var_y]]
    df_plot$facet_var_x <- df_med[[facet_var_x]]
    df_plot$facet_var_y <- df_med[[facet_var_y]]

    # Plot
    p <- ggplot(df_plot, aes(x=var_estimate, y=var_y)) +
      geom_vline(xintercept=0, linetype="dashed", color="grey") +
      geom_point() +
      geom_errorbarh(aes(xmin=var_ci_lo, xmax=var_ci_up, height=0.5),
                     linewidth=0.4) +
      facet_grid(rows=dplyr::vars(facet_var_y),
                 cols=dplyr::vars(facet_var_x)) +
      labs(y=NULL, x=NULL, title=title) +
      coord_cartesian(xlim=xlim) +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor = element_line(linewidth=0.2),
            panel.grid.major = element_line(linewidth=0.2))
    return(p)

  }

  titles <- c("Proportion Mediated",
              "Natural Direct Effect",
              "Natural Indirect Effect")
  # titles <- c("Proportion Mediated (t_0=101)",
  #             "Natural Direct Effect (t_0=101)",
  #             "Natural Indirect Effect (t_0=101)")
  p_pm <- forest_plot(df_med2, title=titles[1],
                      var_estimate="pm_est", var_ci_lo="pm_lo",
                      var_ci_up="pm_up", var_y="marker",
                      facet_var_x="severity", facet_var_y="dataset",
                      xlim=c(-2,2))
  p_nde <- forest_plot(df_med2, title=titles[2],
                       var_estimate="nde_est", var_ci_lo="nde_lo",
                       var_ci_up="nde_up", var_y="marker",
                       facet_var_x="severity", facet_var_y="dataset",
                       xlim=c(-1,1))
  p_nie <- forest_plot(df_med2, title=titles[3],
                       var_estimate="nie_est", var_ci_lo="nie_lo",
                       var_ci_up="nie_up", var_y="marker",
                       facet_var_x="severity", facet_var_y="dataset",
                       xlim=c(-1,1))

  # Export 10x6
  ggsave(filename = paste0("Figures + Tables/", cfg2$analysis,
                           " plots/plot_pm.pdf"),
         plot=p_pm, device="pdf", width=10, height=6)
  ggsave(filename = paste0("Figures + Tables/", cfg2$analysis,
                           " plots/plot_nde.pdf"),
         plot=p_nde, device="pdf", width=10, height=6)
  ggsave(filename = paste0("Figures + Tables/", cfg2$analysis,
                           " plots/plot_nie.pdf"),
         plot=p_nie, device="pdf", width=10, height=6)

}

