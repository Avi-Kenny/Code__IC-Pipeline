###########################.
##### Qbins estimator #####
###########################.

if (F) {
  
  if ("Qbins" %in% cfg2$estimators$cr) {
    
    # Obtain estimates
    s_orig <- dat$v$s[!is.na(dat$v$s)]
    s_grid <- seq(from=min(s_orig), to=max(s_orig), length.out=101)
    ests <- est_curve(
      dat_orig = dat$v,
      estimator = "Qbins",
      params = cfg2$params,
      points = s_grid,
      dir = "decr",
      return_extra = c()
    )
    
    if (flags$save_data_objs) {
      saveRDS(ests, paste0(cfg2$analysis," plots/ests_q_",cfg2$tid,".rds"))
    }
    
    run_cve <- as.logical("CVE" %in% cfg2$plots)
    ests2 <- process_ests(ests, s_grid, run_cve=run_cve,
                          lab_risk="Risk, Qbins", lab_cve="CVE, Qbins")
    plot_data_risk <- rbind(plot_data_risk, ests2$risk)
    if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
    
  }
  
}



#################################################################.
##### Debugging: Grenander variance scale factor components #####
#################################################################.

if (F) {
  
  print("Grenander variance scale factor components")
  print("deriv_r_Mn")
  print(ests$deriv_r_Mn(seq(0,1,0.05)))
  print("f_s_n")
  print(ests$f_s_n(seq(0,1,0.05)))
  print("gamma_n")
  print(ests$gamma_n(seq(0,1,0.05)))
  print("deriv_r_Mn*f_s_n")
  print(ests$deriv_r_Mn(seq(0,1,0.05))*ests$f_s_n(seq(0,1,0.05)))
  print("deriv_r_Mn*gamma_n")
  print(ests$deriv_r_Mn(seq(0,1,0.05))*ests$gamma_n(seq(0,1,0.05)))
  print("f_s_n*gamma_n")
  print(ests$f_s_n(seq(0,1,0.05))*ests$gamma_n(seq(0,1,0.05)))
  print("deriv_r_Mn*f_s_n*gamma_n")
  print(ests$deriv_r_Mn(seq(0,1,0.05))*ests$f_s_n(seq(0,1,0.05))*
          ests$gamma_n(seq(0,1,0.05)))
  
}



#####################################.
##### Debugging: Grenander misc #####
#####################################.

if (F) {
  
  # Debugging: Examine intermediate objects
  
  x_vals <- ests$Phi_n(ests$grid)
  inds <- !base::duplicated(x_vals)
  x_vals <- x_vals[inds]
  y_vals <- -1 * ests$Gamma_os_n(ests$grid[inds])
  GCM_f <- approxfun(x=ests$gcm$x.knots, y=ests$gcm$y.knots,
                     method="linear", rule=1)
  df_A <- data.frame(x=x_vals, y=y_vals, y2=GCM_f(x_vals))
  grid2 <- round(seq(0,1,0.02),2)
  df_B <- data.frame(x=grid2, y=ests$r_Mn_Gr(grid2))
  plot1 <- ggplot(df_A, aes(x=x,y=y)) + geom_point(alpha=0.4) + geom_line(aes(y=y2)) +
    labs(title="x=Phi(S), y=Gamma_n(S)")
  plot2 <- ggplot(df_B, aes(x=x,y=y)) + geom_line() + labs(title="r_Mn estimates")
  
  ggsave(
    filename = paste0(cfg2$analysis," plots/debug_A_",cfg2$tid,".pdf"),
    plot=plot1, device="pdf", width=6, height=4
  )
  ggsave(
    filename = paste0(cfg2$analysis," plots/debug_B_",cfg2$tid,".pdf"),
    plot=plot2, device="pdf", width=6, height=4
  )
  
  # grid <- round(seq(0,1,0.01),2)
  # gcm <- approxfun(x=ests$gcm$x.knots, y=ests$gcm$y.knots, ties="ordered")
  # omega_n <- Vectorize(function(s) {
  #   ests$omega_n(x=c(0,0),s,y=100,delta=0)
  # })
  # etastar_n <- Vectorize(function(s) {
  #   ests$etastar_n(s,x=c(0,0))
  # })
  # Q_n <- Vectorize(function(s) {
  #   ests$Q_n(t=cfg2$t_0, x=c(0,0), s)
  # })
  # 
  # int_data <- data.frame(
  #   x = rep(grid,3),
  #   y = c(ests$Psi_n(grid), gcm(grid), ests$dGCM(grid)),
  #   which = rep(c("Psi_n (-1*Theta_os_n)","gcm","dGCM (-1*r_Mn)"), each=101)
  # )
  # plot1 <- ggplot(int_data, aes(x=x, y=y, color=which)) +
  #   geom_line() +
  #   theme(legend.position="bottom")
  # 
  # ggsave(
  #   filename = paste0(cfg2$analysis," plots/debug_",cfg2$tid,".pdf"),
  #   plot=plot1, device="pdf", width=6, height=4
  # )
  
  if (F) {
    
    # Q_n: conditional survival function (as a function of S)
    # !!!!! Continue
    # as.data.frame(cbind(x1=rep(0.2,n), x2=rep(1,n)))
    int_data2 <- data.frame(
      x = grid,
      y = Q_n(grid)
      # which = rep(c("",""),each=101)
    )
    ggplot(int_data2, aes(x=x, y=y)) + # color=which
      geom_line() +
      theme(legend.position="bottom")
    
    # # omega_n
    # int_data3 <- data.frame(
    #   x = grid,
    #   y = omega_n(grid)
    #   # which = rep(c("",""),each=101)
    # )
    # ggplot(int_data3, aes(x=x, y=y)) + # color=which
    #   geom_line() +
    #   theme(legend.position="bottom")
    
    # # etastar_n
    # int_data4 <- data.frame(
    #   x = grid,
    #   y = etastar_n(grid)
    #   # which = rep(c("",""),each=101)
    # )
    # ggplot(int_data4, aes(x=x, y=y)) + # color=which
    #   geom_line() +
    #   theme(legend.position="bottom") +
    #   labs(title="etastar_n")
    
  }
  
}



######################################.
##### Summary stats / DQA checks #####
######################################.

if (F) {
  
  library(ggfortify)
  
  # Alias vectors
  ind_tx <- df_tx[[cfg2$v$event]]
  ind_ct <- df_ct[[cfg2$v$event]]
  time_tx <- df_tx[[cfg2$v$time]]
  time_ct <- df_ct[[cfg2$v$time]]
  
  # Number of cases in each group
  num_case_tx <- sum(ind_tx)
  num_case_ct <- sum(ind_ct)
  num_case_tx_t_0 <- sum(ind_tx[time_tx<=cfg2$t_0])
  num_case_ct_t_0 <- sum(ind_ct[time_ct<=cfg2$t_0])
  num_atrisk_tx <- length(ind_tx)
  num_atrisk_ct <- length(ind_ct)
  print(paste0("Number of cases in vaccine group: ", num_case_tx))
  print(paste0("Number of cases in control group: ", num_case_ct))
  print(paste0("Number of cases by day ", cfg2$t_0, " in vaccine group: ",
               num_case_tx_t_0))
  print(paste0("Number of cases by day ", cfg2$t_0, " in control group: ",
               num_case_ct_t_0))
  print(paste0("Number at-risk in vaccine group: ", num_atrisk_tx))
  print(paste0("Number at-risk in control group: ", num_atrisk_ct))
  print(paste0("Naive P(COVID by day ", cfg2$t_0, ") in vaccine group: ",
               round(num_case_tx_t_0/num_atrisk_tx,3)))
  print(paste0("Naive P(COVID by day ", cfg2$t_0, ") in control group: ",
               round(num_case_ct_t_0/num_atrisk_ct,3)))
  print(paste0("Naive vaccine efficacy: ",
               round(1 - (num_case_tx_t_0/num_atrisk_tx) /
                       (num_case_ct_t_0/num_atrisk_ct),3)))
  
  # Fraction of point mass at edge
  s <- dat$v$s
  round(sum(s==min(s,na.rm=T),na.rm=T) / sum(!is.na(s)), 3)
  
  # Distribution of event times (Ph2=0 vs. Ph2=1)
  if (!is.na(cfg2$v$ph2)) {
    ggplot(
      data.frame(
        x = time_tx[which(ind_tx==1)],
        ph2 = df_tx[[cfg2$v$ph2]][which(ind_tx==1)]
      ),
      aes(x=x, fill=ph2)
    ) +
      facet_wrap(~ph2) +
      # geom_vline(xintercept=c(138,195), linetype="dashed", color="#333333") +
      geom_histogram() +
      labs(title="Distribution of event times, by Ph2 indicator", x="Time")
  }
  
  # Distribution of event times (Tx vs. Ct)
  ggplot(
    data.frame(
      x = c(time_tx[which(ind_tx==1)], time_ct[which(ind_ct==1)]),
      which = c(rep("Tx",num_case_tx), rep("Ct",num_case_ct))
    ),
    aes(x=x, fill=which)
  ) +
    facet_wrap(~which) +
    geom_vline(xintercept=195, linetype="dashed", color="grey") +
    geom_histogram() +
    labs("Distribution of event times")
  
  # Distribution of censoring times (Tx vs. Ct)
  ggplot(
    data.frame(
      x = c(time_tx[which(ind_tx==0)], time_ct[which(ind_ct==0)]),
      which = c(rep("Tx",num_atrisk_tx-num_case_tx),
                rep("Ct",num_atrisk_ct-num_case_ct))
    ),
    aes(x=x, fill=which)
  ) +
    facet_wrap(~which) +
    geom_vline(xintercept=195, linetype="dashed", color="grey") +
    geom_histogram() +
    labs("Distribution of event times")
  
  # Treatment group survival (entire cohort)
  survfit(
    formula(paste0("Surv(",cfg2$v$time,",",cfg2$v$event,")~1")),
    data = df_tx
  ) %>% autoplot()
  
  # Treatment group survival (subcohort)
  survfit(
    formula(paste0("Surv(",cfg2$v$time,",",cfg2$v$event,")~1")),
    data = df_tx,
    weights = wt.D29start1
  ) %>% autoplot()
  
}



##################.
##### Header #####
##################.

if (F) {
  
  # ...
  
}



##################.
##### Header #####
##################.

if (F) {
  
  # ...
  
}



##################.
##### Header #####
##################.

if (F) {
  
  # ...
  
}



################.
##### Misc #####
################.

if (F) {
  
  # Data processing
  {
    # Read in raw data
    df_raw_1 <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_Bli",
                                "nded_Phase_Data/adata/janssen_pooled_real_dat",
                                "a_processed_with_riskscore.csv"))
    df_raw_adcp <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_",
                                   "Blinded_Phase_Data/adata/janssen_pooled_re",
                                   "alADCP_data_processed_with_riskscore.csv"))
    df_raw_psv <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_B",
                                  "linded_Phase_Data/adata/janssen_pooled_real",
                                  "PsV_data_processed_with_riskscore.csv"))
    
    # Save datasets
    saveRDS(get(paste0("dat_orig_",d)),
            file=paste0("Janssen data/dat_orig_",d,"_Janssen.rds"))
    saveRDS(get(paste0("df_ph1_",d)),
            file=paste0("Janssen data/df_ph1_",d,"_Janssen.rds"))
    saveRDS(get(paste0("df_ct_",d)),
            file=paste0("Janssen data/df_ct_",d,"_Janssen.rds"))
    saveRDS(get(paste0("df_tx_",d)),
            file=paste0("Janssen data/df_tx_",d,"_Janssen.rds"))
    
    # Read datasets
    assign(paste0("dat_orig_",d),
           readRDS(paste0("Janssen data/dat_orig_",d,"_Janssen.rds")))
    assign(paste0("df_ph1_",d),
           readRDS(paste0("Janssen data/df_ph1_",d,"_Janssen.rds")))
    assign(paste0("df_ct_",d),
           readRDS(paste0("Janssen data/df_ct_",d,"_Janssen.rds")))
    assign(paste0("df_tx_",d),
           readRDS(paste0("Janssen data/df_tx_",d,"_Janssen.rds")))
    
  }
  
  # Use Kaplan-Meier to calculate survival and efficacy with CI
  {
    # Calculate control group survival (KM; with SE)
    srv_ct <- survfit(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
      data = df_ct
    )
    rate_ct <- 1 - srv_ct$surv[which.min(abs(srv_ct$time-cfg2$t_0))]
    ci_lo_ct <- 1 - srv_ct$upper[which.min(abs(srv_ct$time-cfg2$t_0))]
    ci_up_ct <- 1 - srv_ct$lower[which.min(abs(srv_ct$time-cfg2$t_0))]
    var_ct <- ((ci_up_ct-ci_lo_ct)/3.92)^2
    
    # Calculate treatment group survival (KM; with SE)
    srv_tx <- survfit(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
      data = df_tx
    )
    rate_tx <- 1 - srv_tx$surv[which.min(abs(srv_tx$time-cfg2$t_0))]
    ci_lo_tx <- 1 - srv_tx$upper[which.min(abs(srv_tx$time-cfg2$t_0))]
    ci_up_tx <- 1 - srv_tx$lower[which.min(abs(srv_tx$time-cfg2$t_0))]
    var_tx <- ((ci_up_tx-ci_lo_tx)/3.92)^2
    
    # Calculate overall vaccine efficacy (KM; delta method; with SE+CI)
    ve_overall <- 1 - (rate_tx/rate_ct)
    ve_se <- sqrt(rate_ct^-2*var_tx + rate_tx^2*rate_ct^-4*var_ct)
    ve_overall_lo <- ve_overall - 1.96*ve_se
    ve_overall_hi <- ve_overall + 1.96*ve_se
    print(paste0("Overall VE: ",round(100*ve_overall,1),"% (",
                 round(100*ve_overall_lo,1),"% -- ",round(100*ve_overall_hi,1),
                 "%)"))
    
    # Calculate overall vaccine efficacy (KM; within subcohort)
    srv_tx_sub <- survfit(coxph(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
      data = dplyr::filter(df_tx, ph2.D29start1==1),
      weights = wt.D29start1
    ))
    rate_tx_sub <- 1 - srv_tx_sub$surv[which.min(abs(srv_tx_sub$time-cfg2$t_0))]
    ve_subcohort <- 1 - (rate_tx_sub/rate_ct)
    print(paste0("Overall VE (subcohort): ", round(100*ve_subcohort,1), "%"))
    
    # Plot KM curve
    library(ggfortify)
    survfit(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~Trt,
      data = df_ph1_1
    ) %>% autoplot()
    
  }
  
  # Check marker quantiles and edge mass
  {
    # Alias markers
    a <- list(
      dat_orig_1$a_list[[1]],
      dat_orig_1$a_list[[2]],
      dat_orig_adcp$a_list[[3]]
    )
    a <- lapply(a, function(x) { x[!is.na(x)] }) # Remove NA values
    a <- lapply(a, function(x) { 10^x }) # Re-express on natural scale
    
    # Check that min marker values equal one-half the positivity cutoff or LOD
    lod <- c(
      2*min(a[[1]], na.rm=T), # 10.84237 = PosCutoff/2
      2*min(a[[2]], na.rm=T), # 14.08585 = PosCutoff/2
      2*min(a[[3]], na.rm=T)  # 11.57 = LOD/2
    )
    
    # Check marker quantiles
    quantile(a[[1]], probs=seq(0,1,0.1))
    quantile(a[[2]], probs=seq(0,1,0.1))
    quantile(a[[3]], probs=seq(0,1,0.1))
    
    # Check percent mass at left edge
    sum(a[[1]]==min(a[[1]]))/length(a[[1]])
    sum(a[[2]]==min(a[[2]]))/length(a[[2]])
    sum(a[[3]]==min(a[[3]]))/length(a[[3]])
    
    # Check percent mass at right edge
    sum(a[[1]]==max(a[[1]]))/length(a[[1]])
    sum(a[[2]]==max(a[[2]]))/length(a[[2]])
    sum(a[[3]]==max(a[[3]]))/length(a[[3]])
  }
  
  # Extra pieces to overlay Cox model gcomp estimator in plot
  {
    theta_ests_gcomp <- ests$gcomp(p_grid)
    ests_gcomp <- cve(theta_ests_gcomp)
    ests_gcomp <- ifelse(which,ests_gcomp,NA)
    
    plot_data_2 <- data.frame(
      x = c(x1,rep(s_grid,2),x1,x2),
      y = c(ests_cve[1],ests_cve,ests_gcomp,rep(ve_overall,2)),
      which = c("Controlled VE",rep(c("Controlled VE","Cox model"),
                                    each=length(p_grid)),
                rep("Overall VE",2)),
      ci_lo = c(ci_lo[1],ci_lo,ests_gcomp,rep(ve_ci[1],2)),
      ci_up = c(ci_up[1],ci_up,ests_gcomp,rep(ve_ci[2],2))
    )
    plot_2 <- plot_1 %+% plot_data_2
    suppressMessages({
      plot_2 <- plot_2 +
        scale_color_manual(values=c("darkblue","purple","darkgrey")) +
        scale_fill_manual(values=c("darkblue","purple","darkgrey"))
    })
    name_2 <- paste0(cfg2$analysis," plots/plot_w_Cox_",cfg2$tid,".pdf")
    ggsave(filename=name_2, plot=plot_2, device="pdf", width=6, height=4)
    
  }
  
}
