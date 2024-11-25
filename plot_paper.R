library(dplyr)
library(data.table)
library(ggplot2)
source("run_functions.R")
source("plot_utils.R")

tot_results <- read_results(filename="raw_results_long", folder_path="/projects/ec194/economic_evaluation_results/raw_results", max_n=200)

res <- tot_results$res
betas <- tot_results$betas
results <- tot_results$results

gg_color <- scale_color_brewer(palette = "Set1") 
gg_fill <- scale_fill_brewer(palette = "Set1")




# Make table
tab <- make_table(res, results, c(rep(1.8, 9), rep(2, 9), rep(2.8,9)), rep(rep(c(1,7,15), each=3), 3), rep(c("best", "lockdown", "nothing"), 9), list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(),list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(),list(), list(), list(), list(), list(), list(), list(), list(), list(), list()), rep("Main case", 27))
f <- file("table.tex")
write(knitr::kable(tab, "latex",digits=c(1,0,1,0,0,0,0,0,0,0,2,0,0,0,0,2) ), f)
close(f)


# Plot examples strategies
plot_one_strat(betas %>% select_def() %>% filter(R==2 & sev==6), results %>% select_def() %>% filter(R==2 & sev==6))
ggsave("baseline_R2_sev6.png", bg="white", width=16, height=12)

plot_one_strat(betas %>% select_def(l=1095) %>% filter(R==1.4 & sev==15), results %>% select_def(l=1095) %>% filter(R==1.4 & sev==15))
ggsave("l1095_R14_sev15.png", bg="white", width=16, height=12)

plot_one_strat(betas %>% select_def(import=210) %>% filter(R==2.5 & sev==3), results %>% select_def(import=210) %>% filter(R==2.5 & sev==3))
ggsave("import210_R25_sev3.png", bg="white", width=16, height=12)

plot_one_strat(betas %>% select_def() %>% filter(R==2 & sev==12), results %>% select_def() %>% filter(R==2 & sev==12))
ggsave("baseline_R2_sev12.png", bg="white", width=16, height=12)

plot_one_strat(betas %>% select_def(l=1095, scale_inf=2, scale_overcapacity=5) %>% filter(R==2.5 & sev==4), results %>% select_def(l=1095, scale_inf=2, scale_overcapacity=5) %>% filter(R==2.5 & sev==4))
ggsave("l1095_scale_inf2_overcap5_R25_sev4.png", bg="white", width=16, height=12)

plot_one_strat(betas %>% select_def(l=1095, scale_inf=2, scale_overcapacity=5) %>% filter(R==2.5 & sev==4), results %>% select_def(l=1095, scale_inf=2, scale_overcapacity=5) %>% filter(R==2.5 & sev==8))
ggsave("l1095_scale_inf2_overcap5_R25_sev8.png", bg="white", width=16, height=12)

plot_one_strat(betas %>% select_def(include_tisk=0) %>% filter(R==1.2 & sev==5), results %>% select_def(include_tisk=0) %>% filter(R==1.2 & sev==5))
ggsave("no_ttiq_12_5.png", bg="white", width=16, height=12)

ggplot(betas  %>% select_def() %>% filter(which=="best" & cost_beta_change==30 & R %in% c(1.4, 1.8, 2.2, 3.0) & sev %in% c(1,6,12, 18))%>% mutate(sev=factor(sev, levels=c(18, 12, 6, 1)))) + geom_line(aes(x=time, y=1-betas, group=1), size=1, span=0.1) + facet_grid(sev~R) + theme + scale_x_continuous("Time(days)", minor_breaks=c(0, 100, 200)) + ylab("Goverment Interventions (G(t))") + theme(text = element_text(size=)) + labs(color="Severity")
ggsave("fig_beta_grid.png", width=12, height=12)
ggsave("fig_beta_grid.eps", width=12, height=12)


# Plot histogram figures
plot_desc(res %>% select_def())
ggsave("fig_desc_strats_baseline1.png", width=25, height=20)

plot_desc(res %>% select_def(scale_beta=0.9))
ggsave("fig_desc_strats_generic_reduction2.png", width=25, height=20)

plot_desc(res %>% select_def(import=210))
ggsave("fig_desc_strats_generic_increased_import3.png", width=25, height=20)

plot_desc(res %>% select_def(vaccination_strat=2))
ggsave("fig_desc_strats_generic_weak_vaccine4.png", width=25, height=20)
plot_desc(res %>% select_def(cap=450))
ggsave("fig_desc_strats_generic_increased_icu_cap5.png", width=25, height=20)

plot_desc(res %>% select_def(include_tisk=0))
ggsave("fig_desc_strats_no_ttiq6.png", width=25, height=20)

plot_desc(res %>% select_def(voluntary_mode=0))
ggsave("fig_desc_strats_no_vsd7.png", width=25, height=20)

plot_desc(res %>% select_def(voluntary_mode=3))
ggsave("fig_desc_strats_increased_atruism8.png", width=25, height=20)

plot_desc(res %>% select_def(l=256, scale_ini=5))
ggsave("fig_desc_strats_unc_14d_suppr9.png", width=25, height=20)

plot_desc(res %>% select_def(l=256, scale_ini="R5"))
ggsave("fig_desc_strats_unc_14d_no_measure10.png", width=25, height=20)

plot_desc(res %>% select_def(l=550))
ggsave("fig_desc_strats_longer11.png", width=25, height=20)

plot_desc(res %>% select_def(scale_inf=0.5))
ggsave("fig_desc_strats_reduced_qaly12.png", width=25, height=20)

plot_desc(res %>% select_def(scale_inf=1.5, l=550))
ggsave("fig_desc_strats_long_increased_qaly14.png", width=25, height=20)

plot_desc(res %>% select_def(l=1095, scale_inf=1, scale_overcapacity=1))
ggsave("fig_desc_strats_mitigationI15.png", width=25, height=20)

plot_desc(res %>% select_def(l=1095, scale_inf=2, scale_overcapacity=1))
ggsave("fig_desc_strats_mitigationII16.png", width=25, height=20)

plot_desc(res %>% select_def(l=1095, scale_inf=2, scale_overcapacity=2))
ggsave("fig_desc_strats_mitigationIII17.png", width=25, height=20)

plot_desc(res %>% select_def(l=1095, scale_inf=2, scale_overcapacity=5))
ggsave("fig_desc_strats_mitigationIV18.png", width=25, height=20)



# Main figure

res$av_int <- 1 - res$beta_mean
manual <- list(list(R=c(1.45, 1.45), sev=c(0.5, 1.5)), list(R=c(1.55, 1.55), sev=c(1.5, 4.5)))
q1 <- add_contours(plot_strat(res %>% select_def(), pad=5), res=res %>% select_def(), manual=manual) + theme(legend.position="bottom", text=element_text(size=24),legend.title=element_blank(), legend.margin = margin(t = 0, l = 0, b = 0, r = 0))
q2 <- plot_var("av_int", res%>% select_def()%>% filter(which=="best"), 0,steps=c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9), color_low="#f1eef6", color_high="#034e7b") + geom_text(aes(x=R, y=sev, label=round(av_int,2)),size=8, data=res%>% select_def()%>% filter(which=="best" & sev %% 3==0 & (R*10-11) %% 3<1e-10  )) + theme(text=element_text(size=24))
q3 <- plot_var("frac_inf_rel", res %>%select_def() %>% filter(which=="best"), 0, steps=c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9), color_low="#e5f5f9", color_high="#005824")  + geom_text(aes(x=R, y=sev, label=round(frac_inf_rel,6)),size=8, data=res%>% select_def()%>% filter(which=="best" & sev %% 3==0 & (R*10-11) %% 3<1e-10  ))+ theme(text=element_text(size=24))

q4 <- plot_var("total", res %>% select_def() %>% filter(which=="best"),steps=c(20, 50, 100, 150, 200, 250, 300),label="bnNOK" ,barheight=20)  + theme( text=element_text(size=24)) + scale_fill_steps(breaks=c(20, 50, 100, 150, 200, 250, 300), low = "yellow", high = "red") + geom_text(aes(x=R, y=sev, label=round(total)),size=8, data=res%>% select_def()%>% filter(which=="best" & sev %% 3==0 & (R*10-11) %% 3<1e-10  )) + theme(text=element_text(size=28))
main <- cowplot::plot_grid(q1, q2, q3,q4, ncol=2, rel_heights=c(1,1,1,1), labels=c("A) Optimal Strategy", "B) Average intervention","C) Relative fraction infected", "D) Total Cost of Optimal Strategy"), label_size=25, vjust=1.0, hjust=0,  scale=0.92, align="hv")
ggsave("fig_main_strat.png", width=20, height=20, bg="white")
ggsave("fig_main_strat.eps", width=20, height=20)




q1 <- add_contours(plot_strat(res %>% select_def(scale_inf=0.5)) + theme(legend.position="bottom", text=element_text(size=34)), res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q2 <- add_contours(plot_strat(res %>% select_def(scale_inf=1.5)) + theme(legend.position="bottom", text=element_text(size=34),legend.title=element_blank()),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q3 <- add_contours(plot_strat(res %>% select_def(l=550)) + theme(legend.position="bottom", text=element_text(size=34),legend.title=element_blank()),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q4 <- add_contours(plot_strat(res %>% select_def(scale_inf=1.5, l=550)) + theme(legend.position="bottom", text=element_text(size=34),legend.title=element_blank()),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")

legend_b <- cowplot::get_plot_component(q4  +
    theme(legend.position = "bottom", legend.title=element_blank(), legend.spacing.x=unit(1.3, "cm")), "guide-box", return_all=TRUE)[[3]]



main <- cowplot::plot_grid(q1 + theme(legend.position="none"), q2  + theme(legend.position="none"), q3 + theme(legend.position="none"), q4 + theme(legend.position="none"), ncol=2, labels=c("A) 50% decreased cost of QALY", "B) 50% increased cost of QALY","C) 18 Months", "D) 18 Months and 50% increased cost of QALY"), label_size = 26, vjust=1.0, hjust=0, scale=0.95)
cowplot::plot_grid(main, legend_b ,ncol=1, rel_heights = c(1, 0.1))

ggsave("fig_strat_external.png", width=20, height=20, bg="white")
ggsave("fig_strat_external.eps", width=20, height=20)





q3 <- add_contours(plot_strat(res %>% select_def(beta_fact=0.9)) + theme(legend.position="bottom", text=element_text(size=25),legend.title=element_blank()),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q4 <- plot_diff( res %>% select_def(beta_fact=1.0) %>% filter(which=="best"),res %>% select_def(beta_fact=0.9) %>% filter(which=="best"), steps=c(-50, -40, -30, -20, -10, -5, 0), ,barheight=20, label="bnNOK", color_high="white", color_mid="blue", min_diff=-500)+ theme( text=element_text(size=25))
q5 <- add_contours(plot_strat(res %>% select_def(import=210)) + theme(legend.position="bottom", text=element_text(size=25),legend.title=element_blank()),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q6 <- plot_diff(res %>% select_def() %>% filter(which=="best"), res %>% select_def(import=210) %>% filter(which=="best"), steps=c(0,5, 10, 25, 50, 75, 100, 150, 200), label="bnNOK" ,barheight=20) + theme( text=element_text(size=25))
q7 <- add_contours(plot_strat(res %>% select_def(cap=450)) + theme(legend.position="bottom", text=element_text(size=25),legend.title=element_blank()),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q8 <- plot_diff(res %>% select_def() %>% filter(which=="best"),res %>% select_def(cap=450) %>% filter(which=="best"),  steps=-rev(c(0,2, 5, 10,15, 20,30)), label="bnNOK" ,barheight=20, color_high="white", color_mid="blue", min_diff=-500)  + theme( text=element_text(size=25))
q9 <- add_contours(plot_strat(res %>% select_def(vaccination_strat=2)) + theme(legend.position="bottom", text=element_text(size=25),legend.title=element_blank(), legend.spacing.x=unit(1, "cm")),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q10 <- plot_diff(res %>% select_def() %>% filter(which=="best"),res %>% select_def(vaccination_strat=2) %>% filter(which=="best"),  min_diff=-500, steps=-c(0,10, 20, 30, 50, 70, 90, 110), label="bnNOK" ,barheight=20, color_mid="blue", color_high="white") + theme( text=element_text(size=25))
cowplot::plot_grid(q3 + theme(legend.position="none"),q4,q5 + theme(legend.position="none"),q6,q7 + theme(legend.position="none"),q8,q9 ,q10, ncol=2, labels=c("A) Generic reduction", "B)", "C) Increased import", "D)", "E) Increased ICU capacity", "F)", "G)Vaccine with weak effect", "H)", "I) ICU prevalence < 250 and 18 months", "J)"), label_size=20, vjust=1, hjust=0, scale=0.95, align="hv")
ggsave("fig_strat_int.png", width=20, heigh=26)
ggsave("fig_strat_int.eps", width=20, heigh=26)





tmp  <-res %>% select_def(l=550, icu_prev_lim=250) %>% filter(which=="best")
max_resp <- results %>% select_def(l=550, icu_prev_lim=250) %>% filter(which=="best") %>% group_by(R, sev) %>% filter(time < 150) %>% summarise(max_resp=max(resp)) 
comb  <- tmp %>% left_join(max_resp, by=c("R", "sev"))
manual550 <- list(list(R=c(1.45, 1.45), sev=c(1.5, 14.5)), list(R=c(1.65, 1.65), sev=c(19.5, 20.5)))
q11 <- add_contours(plot_strat(res %>% select_def(l=550, icu_prev_lim=250), key="strategy4_plot"),res %>%select_def(l=550), manual=manual550) + theme(legend.position="bottom", text=element_text(size=22), legend.spacing.x=unit(10.0, "cm"),legend.title=element_blank()) 
q12 <- plot_diff(res %>% select_def(l=550) %>% filter(which=="best"), res %>% select_def(l=550, icu_prev_lim=250) %>% filter(which=="best"),  steps=c(-0, 100, 200, 300, 400), diverging=TRUE, label="bNOK" ,barheigh=25, min_diff=-200) + theme( text=element_text(size=25))#, legend.position="bottom")
cowplot::plot_grid(q11, q12, ncol=2, labels=c("A)", "B)" ), label_size=20, vjust=1, hjust=0, scale=0.95, align="hv")
ggsave("fig_strat_flatten.png", width=20, heigh=10)
ggsave("fig_strat_flatten.eps", width=20, heigh=10)


q1 <- add_contours(plot_strat(res %>% select_def(l=3*365)) + theme(legend.position="bottom", text=element_text(size=34)),res %>%select_def(), manual=manual)
q2 <- add_contours(plot_strat(res %>% select_def(l=3*365, scale_overcapacity=1, scale_inf=2)) + theme(legend.position="bottom", text=element_text(size=34)),res %>%select_def(), manual=manual)
q3 <- add_contours(plot_strat(res %>% select_def(l=3*365, scale_overcapacity=2, scale_inf=2)) + theme(legend.position="bottom", text=element_text(size=34)),res %>%select_def(), manual=manual)
q4 <- add_contours(plot_strat(res %>% select_def(l=3*365, scale_overcapacity=5, scale_inf=2)) + theme(legend.position="bottom", text=element_text(size=34)),res %>%select_def(), manual=manual)
legend_b <- cowplot::get_plot_component(
  q3+ 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.title=element_blank(), legend.spacing.x=unit(1.3, "cm")), "guide-box", return_all=TRUE
)[[3]]

main <- cowplot::plot_grid(q1 + theme(legend.position="none"), q2  + theme(legend.position="none"), q3 + theme(legend.position="none"), q4 + theme(legend.position="none"), ncol=2, labels=c("A) Baseline for 3 years", "B) 2X QALY Cost", "C) 2X QALY cost & additional 2X overcapacity", "D) 2X QALY cost & additional 5X overcapacity"), label_size = 26, vjust=1.0, hjust=0, scale=0.93)
cowplot::plot_grid(main, legend_b, ncol=1, rel_heights = c(1, 0.1))
ggsave("fig_strat_external_flatten.png", width=20, height=20, bg="white")




# VSD
q1 <- add_contours(plot_strat(res %>% select_def(voluntary_mode=1)) + theme(legend.position="bottom", text=element_text(size=34),legend.title=element_blank()),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q2 <- plot_diff(res %>% select_def() %>% filter(which=="best"), res %>% select_def(voluntary_mode=1) %>% filter(which=="best"), steps=c(-8, -5,-2, 0, 2,3,5, 8), label="bnNOK" ,barheight=20, min_diff=-200, diverging=TRUE, variable="total") + theme( text=element_text(size=30))
q3 <- add_contours(plot_strat(res %>% select_def(voluntary_mode=3)) + theme(legend.position="bottom", text=element_text(size=34),legend.title=element_blank()),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q4 <- plot_diff(res %>% select_def() %>% filter(which=="best"), res %>% select_def(voluntary_mode=3) %>% filter(which=="best"), steps=c(-30, 30),diverging=TRUE, label="bnNOK" ,barheight=20, min_diff=-40) + theme( text=element_text(size=30))
q5 <- add_contours(plot_strat(res %>% select_def(voluntary_mode=0)) + theme(legend.position="bottom", text=element_text(size=34),legend.title=element_blank(), legend.spacing.x=unit(1, "cm")),res %>%select_def(), manual=manual) #+ scale_fill_brewer("Dark")
q6 <- plot_diff(res %>% select_def() %>% filter(which=="best"), res %>% select_def(voluntary_mode=0) %>% filter(which=="best"),  steps=c( -20,-10,0,10, 20, 30, 50, 70, 90),color_low="blue", color_mid="white", diverging=TRUE, label="bnNOK" ,barheight=20, min_diff=-200) + theme( text=element_text(size=30))


legend_b <- cowplot::get_plot_component(
  q1 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"), "guide", return_all=TRUE
)[[3]]

cowplot::plot_grid(q3 + theme(legend.position="none"), q4 , q5 , q6, ncol=2,  labels=c("A) Increased Altruism", "B)", "C) No VSD", "D)"), label_size = 26, vjust=1.0, hjust=0, align="hv", scale=0.95)


ggsave("fig_strat_vsd.png", width=20, height=20)
ggsave("fig_strat_vsd.eps", width=20, height=20)



ggplot(results %>% select_def() %>% filter(which=="best" & time < 220 & include_tisk==1 & cost_beta_change==30) )+ geom_line(aes(x=time, y=incidence/5.37e6*10000, group=run_name), size=1, span=0.1) + facet_grid(strategy4~.) + theme + xlab("Time(days)") + ylab("Incidence per 10 000")
ggsave("fig_inc_all.png", height=20, width=17)
x_scale_ini <- 5
plot_unc(x_scale_ini=5)


res$R_desc <- paste0("R=", res$R)
q1 <- ggplot(res %>% filter(label=="cont_si" &  which=="best")) + geom_line(aes(x=scale_inf*1.4, y=av_rel_intervention, color=factor(sev)), size=2) + facet_wrap(.~R_desc) + theme + xlab("Cost of QALY (mNOK)") + ylab("Average relative intervention") + theme(text = element_text(size=30)) + labs(color="Severity    ") + scale_color_brewer(palette="Dark2")
q2 <- ggplot(res %>% filter(label=="cont_si" &  which=="best")) + geom_line(aes(x=scale_inf*1.4, y=frac_inf_rel, color=factor(sev)), size=2) + facet_wrap(.~R_desc) + theme + xlab("Cost of QALY (mNOK)") + ylab("Relative fraction infected") + theme(text = element_text(size=30)) + labs(color="Severity    ")+ scale_color_brewer(palette="Dark2")
cowplot::plot_grid(q1 + theme(legend.position="none"), q2 + theme(legend.position="bottom"), nrow=2, ncol=1,label_size = 36)
ggsave("fig_cont_scale_inf.png", width=20, height=16, bg="white")

q1 <- ggplot(betas %>% filter(label=="cont_si" & which=="best" & R==2 &  scale_inf %in% c(2/49*16, 2/49*32) & sev==7)) + geom_line(aes(x=time, y=1- betas, color=factor(round(scale_inf*1.4))), size=2) + theme + ylab("Goverment Intervention") + xlab("Time")+ theme(text = element_text(size=30)) 
q2 <- ggplot(results %>% filter(label=="cont_si" & which=="best" & R==2 &  scale_inf %in% c(2/49*16, 2/49*32) & sev==7)) + geom_line(aes(x=time, y=incidence, color=factor(round(scale_inf*1.4))), size=2)+ theme + ylab("Incidence") + xlab("Time")+ theme(text = element_text(size=30)) 
q3 <- ggplot(results %>% filter(label=="cont_si" & which=="best" & R==2 &  scale_inf %in% c(2/49*16, 2/49*32) & sev==7)) + geom_line(aes(x=time, y=tot_infected, color=factor(round(scale_inf*1.4))), size=2)+ theme + ylab("Cumulative incidence") + xlab("Time")+ theme(text = element_text(size=30))  + labs(color="QALY cost    ")
cowplot::plot_grid(q1 + theme(legend.position="none"), q2 + theme(legend.position="none"), q3+ theme(legend.position="bottom"), nrow=2, ncol=2)

ggsave("fig_cont_scale_inf_chosen.png", width=16, height=16)






q1  <- add_contours(plot_strat(res %>% select_def(include_tisk=0, l=270)), res=res %>% select_def(), manual=manual) + theme(legend.position="bottom", legend.title=element_blank())
ggsave("fig_no_tisk.png" ,width=18, height=15)




res$over_cap_base <- res$scale_overcapacity*2
q1 <- ggplot(res %>% filter(which=="best" & label=="test_flatten")) + geom_line(aes(x=sev, y=av_rel_intervention, color=factor(over_cap_base)), size=2) + theme_minimal() + xlab("Severity") + ylab("Average Relative Intervention") + scale_color_brewer("Increased Cost for Overcapacity         ", palette="Dark2") + theme(text = element_text(size=36), legend.position="bottom")
q2 <- ggplot(res %>% filter(which=="best"& label=="test_flatten")) + geom_line(aes(x=sev, y=frac_inf_rel, color=factor(over_cap_base)), size=2) + theme_minimal() + xlab("Severity") + ylab("Relative Fraction Infected") + scale_color_brewer("Cost overcapacity", palette="Dark2") + theme(text = element_text(size=36), legend.position="bottom")
q3 <- ggplot(res %>% filter(which=="best"& label=="test_flatten")) + geom_line(aes(x=sev, y=frac_max_inc, color=factor(over_cap_base)), size=2) + theme_minimal() + xlab("Severity") + ylab("Relative max incidence") + scale_color_brewer("Overcapacity", palette="Dark2") + theme(text = element_text(size=36), legend.position="bottom")

mr <- results %>%  filter(label=="test_flatten") %>% group_by(sev, scale_overcapacity, which) %>% summarise(max_resp=max(resp))
mr$over_cap_base <- mr$scale_overcapacity*2

a <- res %>% filter(which=="best") %>% left_join(mr %>% filter(which=="best"))

q4 <- ggplot(a) + geom_line(aes(x=sev, y=max_resp, color=factor(over_cap_base)), size=2) + theme_minimal() + xlab("Severity") + ylab("Peak ICU prevalence") + scale_color_brewer("Overcapacity", palette="Dark2") + theme(text = element_text(size=36), legend.position="bottom")

legend_b <- cowplot::get_plot_component(
  q1 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"), "guide", return_all=TRUE
)[[3]]



main <- cowplot::plot_grid(q1 + theme(legend.position="none"),q2+theme(legend.position="none"),q3+ theme(legend.position="none"), q4 + theme(legend.position="none"), scale=0.95 ,labels=c("A)", "B)", "C)", "D)"))
cowplot::plot_grid(main, legend_b, rel_heights=c(1, 0.1), ncol=1, nrow=2)

ggsave("fig_overcap.png", width=20, height=20, bg="white")

