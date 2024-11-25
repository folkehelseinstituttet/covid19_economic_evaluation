library(dplyr)
library(metapop)
library(metapopnorge)
library(ggplot2)
source("run_functions.R")
source("QALY.R")
source("parameters.r")
library(data.table)
param_file <- "parameter_files/parameters_vaccination.xlsx"
library(futile.logger)
flog.threshold(ERROR) 

run_vals <- function(x){
  cost(beta=x$beta, R=x$R, sev=x$sev, run_params=run_params, cost_params=cost_params, 
  ret_both=T,include_tisk=1, scale_ini=x$scale_ini, l=x$l, voluntary_mode=1) %>% mutate(beta=x$beta, R=x$R, severity=x$sev, scale_ini=x$scale_ini, lockdown_beta=1 - (1 - 1/x$R) / 0.8, l=x$l)

}

paramsets <- list()
for(R in c(1.5, 2)){
  for(sev in c(1,7, 15)){
    for(scale_ini in c(5)){
      for(l in c(270)){
        for(beta in seq(0.3, 1.0, by=0.025)){
        paramsets[[length(paramsets) +1]] <- list(R=R/0.93, sev=sev, beta=beta, scale_ini=scale_ini, l=l)
        }
      }
    }
  }
}


df <- parallel::mclapply(paramsets, run_vals, mc.cores=30)

df <- rbindlist(df)
plot_df <- df %>% mutate(hosp=hosp_loss+icu_loss) %>% select(beta, severity,R, lost_qaly, hosp, prod_loss_illness, loss_tisk, loss_interventions) %>% tidyr::pivot_longer(cols=c("lost_qaly", "hosp", "prod_loss_illness", "loss_tisk", "loss_interventions")) %>%mutate(type=recode(name, lost_qaly="Direct QALY loss", hosp="Hospital Capacity", prod_loss_illness="Production Loss Illness", loss_tisk="TTIQ", loss_interventions="Contact Reduction", inf="Total HC", beta_cost="Contact Reduction"))
plot_df$intervention <- 1 - plot_df$beta
ggplot(plot_df %>% filter(beta >0.5, R==1.5/0.93)) %>% + geom_line(aes(x=intervention, y=value, color=type), size=3)+ theme_minimal()+ xlab("Goverment Interventions(G)") + ylab("Costs (bnNok)") + theme(text = element_text(size=44)) + scale_color_brewer(palette="Dark2") + labs(color="Type") + facet_grid(.~severity, scales="free") 
ggsave("individual_costs.png", width=22, height=15)




#Plot age dists
scaling <- 0.92
x_R <- 2.4/scaling
x_sev <- 11
beta_from_scaling <- (1 - scaling)/run_params$combine_param_1
beta_ld <- 1- (1/(x_R) -1 + run_params$combine_param_1*beta_from_scaling)/(run_params$combine_param_2*beta_from_scaling - run_params$combine_param_1)

beta_opt <- betas %>% select_def() %>% filter(which=="best" & abs(R-2.4)<0.01 & sev==x_sev) %>% pull(betas)

cost_nothing <-  cost(beta=1, R=x_R, sev=x_sev, run_params=run_params, cost_params=cost_params, ret_all=T, scale_ini=5, l=270)
cost_supr <- cost(beta=beta_ld, R=x_R, sev=x_sev, run_params=run_params, cost_params=cost_params, ret_all=T, scale_ini=5, l=270)
cost_cut <- cost(beta=beta_opt, R=x_R, sev=x_sev, run_params=run_params, cost_params=cost_params, ret_all=T, scale_ini=5, n_cps=269, l=270)

cost_supr$cost
cost_cut$cost


age <- data.frame(age=rep(c("<10", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79","80+"), 3), cost=c(cost_supr$age_breakdown,
                                                                                                              cost_cut$age_breakdown,
                                                                                                              cost_nothing$age_breakdown),
                                                                                                              cost_per_person=c(cost_supr$age_breakdown/get_age_groups(),
                                                                                                              cost_cut$age_breakdown/get_age_groups(),
                                                                                                              cost_nothing$age_breakdown/get_age_groups())*1e6,
                  strategy=c(rep(c("Suppression   ", "Cut Overshoot   ", "No Interventions   "), each=9)))
age$strategy <- factor(age$strategy, levels=c( "Suppression   ","No Interventions   ",  "Cut Overshoot   "))
q1 <- ggplot(age) + geom_col(aes(x=age, y=cost, fill=strategy), position="dodge") + xlab("Age (years)") + scale_y_continuous("Cost (bnNOK)") +  theme_minimal() + theme(text = element_text(size=22)) + scale_fill_brewer( palette="Set1") 
q2 <- ggplot(age) + geom_col(aes(x=age, y=cost_per_person, fill=strategy), position="dodge") + xlab("Age (years)") + scale_y_continuous("Cost per person (1000 NOK)")+  theme_minimal() + theme(text = element_text(size=22)) + scale_fill_brewer( palette="Set1") 

cowplot::plot_grid(q1 + theme(legend.position="none"), q2+ theme(legend.position="bottom", legend.spacing.x=unit(1.3, "cm")), ncol=1, labels = c("A)", "B)"), label_size = 22, scale=0.9)

ggsave("age_dist.png", height=12, width=18, bg="white")


## VSD

add_theme <- function(q){
  q + theme_minimal() + theme(text = element_text(size=28), panel.spacing = unit(2, "lines")) + xlab("Time(days))") + scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0)))
    }

all_res <- list()
res <- list()

R_global_scaling <- 0.9381012
print(2/R_global_scaling)
for(R in c(1.05)){#c(1.5, 2, 3)){
  for( severity in c(1,9,15)){#,3,5,9)){
    for(x_vm in c(0,1,2,3,4,5,6)){
    r <- cost(1, R=R/R_global_scaling, sev=severity,include_tisk=1, vaccination_strat=1, run_params=run_params, cost_params=cost_params, scale_ini=5, l=270, n_cps=0, ret_all=1, voluntary_mode=x_vm, cost_beta_change=3)
    all_res[[length(all_res)+1]] <- r$results %>%mutate(vm=x_vm)
    res[[length(res)+1]] <- r$cost %>%mutate(vm=x_vm, sev=severity)
    }
  }
}
all_res <- data.table::rbindlist(all_res)
res <- data.table::rbindlist(res)

cc <- all_res %>% select(time, severity,vm, starts_with("contact_change")) %>% tidyr::pivot_longer(starts_with("contact_change"))

cc$name <- rep(c("<10", "10-19", "20-29", "30-39","40-49", "50-59", "60-69", "70-79", "80+"),  length(cc$name)/9)

cc <- cc %>% mutate(names=recode(vm, "0"="No VSD", "1"="Baseline", "2"="Double QALY cost", "3"="Altruism", "4"="Baseline2", "5"="4/5 Altruism", "6"="Avg", "7"="New"))

cc$names <- factor(cc$names, levels=c("No VSD", "Baseline", "No Altruism", "Altruism", "4/5 Own", "Increased Altruism", "Avg", "New"))
cc$severity <- factor(cc$severity, levels=c(20, 15, 9, 1))

add_theme(ggplot(cc %>% filter(time < 270 & vm %in% c(0, 1,3))) + geom_line(aes(x=time, y=value, color=name),size=1.5) + facet_grid(severity~names)) + ylab("Contact reduction") + xlab("Time(days)") + labs(color="Age group") + scale_color_brewer(palette="Spectral")
ggsave("vm_strats.png", width=18, height=10, bg="white")

all_res <-all_res  %>% mutate(names=factor(recode(vm, "0"="No VSD", "1"="Baseline", "2"="Double QALY cost", "3"="Altruism", "4"="Baseline", "5"="4/5 Altruism", "6"="Avg", "7"="New"),  levels=c("No VSD", "Baseline", "No Altruism", "Altruism", "4/5 Own", "Increased Altruism", "Avg", "New")))
all_res$severity <- factor(all_res$severity, levels=c(20, 15, 9, 1))
add_theme(ggplot(all_res %>% filter(vm %in% c(0,1,3) & time < 270 )) + geom_line(aes(x=time, y=incidence), size=1.5) + facet_grid(severity~names))
ggsave("vm_strats_inc.png", width=18, height=10, bg="white")
add_theme(ggplot(all_res%>% filter(vm %in% c(0,1,3) & time < 270)) + geom_line(aes(x=time, y=hosp_incidence), size=1.5) + facet_grid(severity~names, scale="free_y") + ylab("Hospital Admissions"))
ggsave("vm_strats_hosp_inc.png", width=18, height=10, bg="white")
