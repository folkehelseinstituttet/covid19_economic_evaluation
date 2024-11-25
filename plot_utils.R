
select_def <- function(data, ...){
  defaults=list(import=1, cap=350, vaccination_strat=1, scale_inf=1.0,
          voluntary_mode=1, scale_ini=5, l=270, include_tisk=1, cost_beta_change=30, beta_fact=1.0, scale_overcapacity=1, icu_prev_lim=0, label="default")
  in_vars <- as.list(match.call())
  vars <- list()
  for(key in names(defaults)){
    if(key %in% names(in_vars)){
      if(!is.null(in_vars[[key]]))
      vars[[key]] <- in_vars[[key]]
    }else{
      vars[[key]] <- defaults[[key]]
    }

  }
  filter_list <-list()
  for(key in names(vars)){
    var = sym(key)
    cond = vars[[key]]
    filter_list[[length(filter_list) +1 ]] <- expr(!!var == !!cond)

  }
 
  data %>% filter(!!!filter_list)
}

theme <- theme_minimal()  + theme(text = element_text(size=44))# + 
#color_sel <- scale_color_brewer("Dark") + scale_fill_brewer("Dark")

plot_one_strat <- function(betas, results){
  betas <- betas %>% mutate(which_n=recode(which, best="Optimal", lockdown="Force R=1", nothing="No Measures"))
  results <- results %>% mutate(which_n=recode(which, best="Optimal    ", lockdown="Force R=1    ", nothing="No Measures    "), which_n=factor(which_n, levels=c("No Measures    ","Force R=1    ","Optimal    ")))

  plot <- betas %>% filter(which=="best") %>% select(time, betas, which_n)
  plot <- rbind(plot, data.frame(time=plot$time, betas=1, which_n="No measures"),data.frame(time=plot$time, betas=betas$beta_mean_supr[1], which_n="Force R=1"))
  q1 <- ggplot(plot) + geom_line(aes(x=time, y=1-betas, color=which_n), size=1.5) + theme + xlab("Time(days)") + ylab("Goverment Interventions (G(t))") + theme(text = element_text(size=22)) + labs(color="Strategy      ") +  theme(legend.position="none") + scale_color_brewer(palette="Dark2")
  q2 <- ggplot(results %>% filter(which %in% c("best", "lockdown", "nothing"))) + geom_line(aes(x=time, y=incidence, group=run_name, color=factor(which_n)), size=1.5) + theme + xlab("Time(days)") + ylab("Incidence") + theme(text = element_text(size=22)) + labs(color="Strategy      ")+ theme(legend.position="bottom") +  scale_color_brewer(palette="Dark2") + geom_vline(xintercept=betas$l[1], alpha=0.5)
  q3 <- ggplot(results %>% filter(which %in% c("best", "lockdown", "nothing"))) + geom_line(aes(x=time, y=resp, group=run_name, color=factor(which_n)), size=1.5) + theme + xlab("Time(days)") + ylab("Prevalence ICU") + theme(text = element_text(size=22)) + labs(color="Strategy") + theme(legend.position="none")  + scale_color_brewer(palette="Dark2")+ geom_vline(xintercept=betas$l[1], alpha=0.5)
  q4 <- ggplot(results %>% filter(which %in% c("best", "lockdown", "nothing")) )+ geom_line(aes(x=time, y=tot_infected/5.37e6, group=run_name, color=factor(which_n)), size=1.5) + theme + xlab("Time(days)") + ylab("Cumulative Incidence ") + theme(text = element_text(size=22)) + labs(color="Strategy") + theme(legend.position="none") + scale_y_continuous(labels=scales::percent_format()) + scale_color_brewer(palette="Dark2")+ geom_vline(xintercept=betas$l[1], alpha=0.5)

  main <- cowplot::plot_grid(q1, q2 + theme(legend.position="none"),q3,q4, ncol=2)
legend_b <- cowplot::get_plot_component(q2  +
    theme(legend.position = "bottom", legend.spacing.x=unit(1.3, "cm")), "guide-box", return_all=TRUE)[[3]]

cowplot::plot_grid(main ,legend_b, ncol=1, rel_heights=c(1, 0.1))
#ggplot(results  %>% select_def(beta_fact=0.9) %>% filter(cost_beta_change==15 & which %in% c("nothing", "best") & R==1.3 & sev==1)) + geom_line(aes(x=time, y=incidence, color=factor(which)), size=1, span=0.1) + theme + xlab("Time(days)") + ylab("Incidence") + theme(text = element_text(size=36)) + labs(color="Severity")



}

add_contours <- function(q1, res, manual=FALSE){
  Rs <- list(c(), c(), c())
  sevs <- list(c(), c(), c())
  ignore="No Interventions"
  unique_R <- sort(unique(res$R))
  R_grid_size <- unique_R[2] - unique_R[1]
  unique_sev <- sort(unique(res$sev))
  sev_grid_size <- unique_sev[2] - unique_sev[1]
  i <- 1
  for(changes in list(c("Cut Overshoot", "Suppression"), c( "No Interventions", "Suppression"), c( "No Interventions","Cut Overshoot"))){
   previous <- ""#changes[1]
   first_time <- TRUE
    for(x_R in unique_R){
    
       for(x_sev in unique_sev){
   #   print(x_R)
         current <- as.character(res %>% filter(R==x_R & sev==x_sev &which=="best") %>% pull(strategy4))[1]

        if(!is.na(current) & !is.na(previous) &current!=previous & previous==changes[1] & current==changes[2] & first_time ){

             
              Rs[[i]] <- c(Rs[[i]], c(x_R - R_grid_size/2))
              sevs[[i]] <- c(sevs[[i]], x_sev - sev_grid_size/2)
              first_time <- FALSE
        }
      #       print(glue::glue("({x_R}, {x_sev}), Current: {current}, Previous: {previous}"))
      if(!is.na(current) & !is.na(previous) & current!= previous & current==changes[2] & previous==changes[1]){
    
   #      Rs[[i]] <- c(Rs[[i]], rep(x_R - R_grid_size/2, 2))
        Rs[[i]] <- c(Rs[[i]], c(x_R - R_grid_size/2,x_R + R_grid_size/2))
    #    sevs[[i]] <- c(sevs[[i]], x_sev - sev_grid_size/2, x_sev +sev_grid_size/2)
       sevs[[i]] <- c(sevs[[i]], rep(x_sev - sev_grid_size/2, 2))
      }

      # if(x_R==max(res$R)){
      #   if(max(Rs[[i]])!= x_R + R_grid_size/2){
      #    Rs[[i]] <- c(Rs[[i]],x_R + R_grid_size/2)
      #   sevs[[i]] <- c(sevs[[i]], x_sev - sev_grid_size/2)
      #   }

     #   }
     
      previous <- current
    }
  }
 # if(all(changes != c("suppression", "No interventions"))){
 # Rs[[i]] <- c(Rs[[i]],max(res$R) + R_grid_size/2)
 # sevs[[i]] <- c(sevs[[i]], sevs[[i]][length(sevs[[i]])])
 # }
  i <- i +1
  }

  j <- 1
  df <- data.frame(R=unlist(Rs), sev=unlist(sevs), type=c(rep(1, length(Rs[[1]])), rep(2, length(Rs[[2]])), rep(3, length(Rs[[3]]))))
  if(class(manual)=="list"){
    if(is.null(names(manual))){
      for( l in manual){
        df <- rbind(df, data.frame(R=l$R, sev=l$sev, type=3 + j))
        j <- j+1
      }
    }else{
     df <- rbind(df, data.frame(R=manual$R, sev=manual$sev, type=4))
    }
  }

 # print(df)
  return(q1 + geom_path(aes(x=R, y=sev, group=type), data=df))
}


prepare_result <- function(raw_results){
    
    
    fails <- c()
    res <- lapply(raw_results, function(x){tryCatch(x$costs, error=function(cond){ print(x); fails <- c(fails, x)})})
    mask <- unlist(lapply(res, function(x) {is.data.frame(x)}))
    res <- rbindlist(res[mask])
    if(!"max_n" %in% colnames(res)){
        res[, max_n:=1000]

    }

   # remove_dups <- duplicated(res[which=="best", .(R, sev, import, cap, vaccination_strat, scale_inf, include_tisk,l,scale_ini,voluntary_mode, cost_beta_change, beta_fact)])
    raw_results_nd <- raw_results[mask]#!remove_dups][mask]
    res <- rbindlist(lapply(raw_results_nd, function(x){x$costs}))
    res <- res %>% distinct()
    setDT(res)
    betas<- rbindlist(lapply(raw_results_nd, function(x) x$betas %>% mutate(period=1:n())), fill=TRUE)
    betas <- betas %>% distinct()
    setDT(betas)
    results <- rbindlist(lapply(raw_results_nd, function(x) x$results), fill=TRUE)
    results <- results %>% select(time, R, sev, cap,import, vaccination_strat, scale_inf, include_tisk, l, voluntary_mode, scale_ini, cost_beta_change, beta_fact, scale_overcapacity, icu_prev_lim, label, which, n_cps,start_point, local_algo, n_local, n_global,fit_cps,
                                    incidence, hosp, resp, D, tot_infected, tot_hosp, tot_resp, starts_with("tot_infected_age"))
    setDT(results)
    if(!"scale_overcapacity" %in% colnames(res)) {
      res <- res %>% mutate(scale_overcapacity=1)
      results <- results %>% mutate(scale_overcapacity=1)
      betas <- betas %>% mutate(scale_overcapacity=1)
    }
    if(!"icu_prev_lim" %in% colnames(res)) {
      res <- res %>% mutate(icu_prev_lim=0)
      results <- results %>% mutate(icu_prev_lim=0)
      betas <- betas %>% mutate(icu_prev_lim=0)

    } 
    
    sum_betas <- betas %>% filter(time<=l)  %>%
        group_by( R, sev, cap,import, vaccination_strat, scale_inf, include_tisk, l, voluntary_mode, scale_ini, cost_beta_change, beta_fact, scale_overcapacity, icu_prev_lim, label, which, n_cps,start_point, local_algo, n_local, n_global,fit_cps) %>% 
        mutate(lagged=betas - lag(betas)) %>% mutate(lags=tidyr::replace_na(lagged,0)) %>%
        summarise(min_beta=min(betas), days_open=sum(betas > 0.99), days_only_tisk = sum(betas >= 0.75 & betas < 0.99), 
        max_beta=max(betas), av_beta=mean(betas), med_beta=median(betas), sd_beta=sd(betas), mean_2_weeks=mean(ifelse(time<=14, betas, NA ),na.rm=T), mean_30_days=mean(ifelse(time<=30, betas, NA ), na.rm=T), mean_50_days=mean(ifelse(time<=50, betas, NA ), na.rm=T),
         days_open_or_tisk=sum(betas >= 0.745), days_closed=sum(betas <=0.05), change_score=sum(lags^2), max_jump=max(abs(lags), na.rm=T))
    sum_betas2 <- betas %>% filter(time<=l)  %>%
        group_by( R, sev, cap,import, vaccination_strat, scale_inf, include_tisk, l, voluntary_mode, scale_ini, cost_beta_change, beta_fact, scale_overcapacity, icu_prev_lim, label, which, n_cps,start_point, local_algo, n_local, n_global,fit_cps) %>% 
        filter(betas < 0.95) %>% summarise(first_intervention=min(time))
    
    sum_betas <- left_join(sum_betas, sum_betas2) %>% mutate(first_intervention=tidyr::replace_na(first_intervention, -1))
    sum_results <- results %>% filter(time<=l)  %>%
        group_by( R, sev, cap,import, vaccination_strat, scale_inf, include_tisk, l, voluntary_mode, scale_ini, cost_beta_change, beta_fact, scale_overcapacity, icu_prev_lim, label, which, n_cps,start_point, local_algo, n_local, n_global,fit_cps) %>% 
        filter(incidence == max(incidence)) %>% summarise(peak_time=min(time))
    res <- res %>% left_join(sum_betas,by=c("R", "sev", "import", "cap", "vaccination_strat", "scale_inf", "include_tisk", "l", "scale_ini", "voluntary_mode", "cost_beta_change", "beta_fact", "scale_overcapacity", "icu_prev_lim", "label", "which","start_point", "local_algo", "n_cps", "n_local", "n_global", "fit_cps"))
    res <- res %>% left_join(sum_results,by=c("R", "sev", "import", "cap", "vaccination_strat", "scale_inf", "include_tisk", "l", "scale_ini", "voluntary_mode", "cost_beta_change", "beta_fact", "scale_overcapacity", "icu_prev_lim", "label", "which","start_point", "local_algo", "n_cps", "n_local", "n_global", "fit_cps"))
    res_nothing <- res[which=="nothing", .(R, sev, import, cap, vaccination_strat, scale_inf, include_tisk,l,scale_ini,voluntary_mode, cost_beta_change,beta_fact,scale_overcapacity,icu_prev_lim,label, infected_no_measures=frac_infected, max_inc_no_measures=max_inc, cost_no_measures=total, peak_time_no_measures=peak_time)]
    res_nothing <- res_nothing[!duplicated(res_nothing)]
    res_ld <- res[which=="lockdown", .(R, sev, import, cap, vaccination_strat, scale_inf, include_tisk,l,scale_ini,voluntary_mode, beta_fact,cost_beta_change,scale_overcapacity,icu_prev_lim,label, beta_mean_supr=beta_mean, beta_mean_supr_last_two=beta_mean_remove_last_2, frac_infected_lockdown=frac_infected)]
    res_ld <- res_ld[!duplicated(res_ld)]    
    
    res <- res %>% left_join(res_nothing, by=c("R", "sev", "import", "cap", "vaccination_strat", "scale_inf", "include_tisk", "l", "scale_ini", "voluntary_mode", "cost_beta_change", "beta_fact", "scale_overcapacity", "icu_prev_lim", "label"))
    res <- res %>% left_join(res_ld, by=c("R", "sev", "import", "cap", "vaccination_strat", "scale_inf", "include_tisk", "l", "scale_ini","voluntary_mode", "cost_beta_change","beta_fact", "scale_overcapacity", "icu_prev_lim", "label"))


    res[, run_name:=paste(R, sev, import, cap, vaccination_strat, scale_inf, scale_ini, include_tisk, l, scale_ini,beta_fact, voluntary_mode, cost_beta_change, which, start_point, local_algo,n_cps, fit_cps, icu_prev_lim, label, fit_cps)]
    res[, max_rel_intervention:=(1-min_beta)/(1-beta_mean_supr)]
    res[, av_rel_intervention:=(1-av_beta)/(1-beta_mean_supr)]

    res[, med_rel_intervention:=(1-med_beta)/(1-beta_mean_supr)]
    res[,frac_inf_rel:=round((frac_infected)/infected_no_measures,2)]
    res[,frac_max_inc:=round(max_inc/max_inc_no_measures,2)]
    res[,beta_mean_rel:=round((1-beta_mean)/(1-beta_mean_supr),2)]
 #   res[,beta_mean_rel_4m:=round((1-beta_mean_4m)/(1-beta_mean_supr_4m),2)]
 #   res[beta_mean_rel < 0, beta_mean_rel:=0]
    res[, rel_min_beta:=1-((1-min_beta)/(1-beta_mean_supr))]




    res[, av_R0:=R*(1-0.8*((1-mean_2_weeks) + 0.0774) + 0.7*(1-mean_2_weeks)*0.0774) *beta_fact*reduction_in_beta_from_S/0.92]
    res[, av_R0_30:=R*(1-0.8*((1-mean_30_days) + 0.0774) + 0.7*(1-mean_30_days)*0.0774) *beta_fact*reduction_in_beta_from_S/0.92]
    res[, av_R0_50:=R*(1-0.8*((1-mean_50_days) + 0.0774) + 0.7*(1-mean_50_days)*0.0774) *beta_fact*reduction_in_beta_from_S/0.92]
    res[, av_R0_use:=av_R0_30]
    res[scale_ini=="R5", av_R0_use:=av_R0]

    res[, shape_d:="Uncertain"]
    res[av_rel_intervention>0.5, shape_d:="Closed"]
    res[ av_rel_intervention <0.15 , shape_d:="Open"]
    res[ av_rel_intervention <0.65 & av_rel_intervention > 0.15 & change_score < 0.3  , shape_d:="Medium - Low CS"]
    res[change_score > 0.25 & av_rel_intervention < 0.6, shape_d:="Cut Overshoot"]
    res[, shape:=shape_d]
    res[shape=="Closed - High CS", shape:="Closed"]
    res[shape=="Open - TISK", shape:="Open"]
    
    res[av_rel_intervention>0.5, strategy1:="Suppression"]
    res[av_rel_intervention<0.10, strategy1:="No Interventions"]
    res[av_rel_intervention>=0.15 & av_rel_intervention<=0.65, strategy1:="Mitigation"]
    
    res[,strategy2:="NULL"]
    res[av_rel_intervention>0.7, strategy2:="Suppression"]
    res[av_rel_intervention<0.10, strategy2:="No Interventions"]
    res[av_rel_intervention>=0.01& av_rel_intervention<=0.7 & change_score>1, strategy2:="Cut Overshoot"]
    res[av_rel_intervention>=0.1 & av_rel_intervention<=0.7 & change_score<1, strategy2:="Mitigation"]
#    res[abs(total - cost_no_measures) < 0.5, strategy2:="No Interventions"]
    #res[(av_rel_intervention>=0.1 | min_beta < 0.1) & av_rel_intervention<=0.65 & l>1000, strategy2:="Cut Overshoot"]
    res[, strategy3:=strategy2]
    res[,strategy4:="NULL"]
        
    res[peak_time >= 15 & (first_intervention > peak_time -5), strategy4:="Cut Overshoot"]
    res[strategy4!="Cut Overshoot" & av_rel_intervention>0.6, strategy4:="Suppression"]
    res[strategy4!="Cut Overshoot" & av_rel_intervention<0.01, strategy4:="No Interventions"]
    #res[av_rel_intervention>=0.15& av_rel_intervention<=0.7 & change_score>0.75, strategy4:="Cut Overshoot"]
    res[strategy4!="Cut Overshoot" & av_rel_intervention>=0.01 & av_rel_intervention<=0.6 , strategy4:="Mitigation"]
   #res[peak_time > 15 & (first_intervention > peak_time -5), strategy4:="Cut Overshoot"]
    
    res[,strategy5:="NULL"]
    res[peak_time > 15 & (first_intervention > peak_time_no_measures -5), strategy5:="Cut Overshoot"]
    res[strategy5!="Cut Overshoot" & av_rel_intervention>0.6, strategy5:="Suppression"]
    res[strategy5!="Cut Overshoot" & av_rel_intervention<0.05, strategy5:="No Interventions"]
    #res[av_rel_intervention>=0.15& av_rel_intervention<=0.7 & change_score>0.75, strategy5:="Cut Overshoot"]
    res[strategy5!="Cut Overshoot" & av_rel_intervention>=0.05 & av_rel_intervention<=0.6 , strategy5:="Mitigation"]
   # res[strategy4=="NULL" & change_score> 0.75, strategy4:="Cut Overshoot"]
    #res[,strategy4:=strategy5]

    #res[strategy2=="Suppression" & av_R0_use > 1.2 & av_rel_intervention<0.95, strategy3:="Mitigation"]
    #res[av_rel_intervention > 0.5 & av_rel_intervention<0.75, strategy3:="Mitigation"]
  

    res <- res %>% mutate(av_rel_2_weeks=(1-mean_2_weeks)/(1-beta_mean_supr))

    setDT(res)

 
    
    
    
    
    betas <- betas %>% left_join(res,  by=c("R"="R","sev", "import", "cap", "vaccination_strat", "scale_inf", "include_tisk", "l", "scale_ini", "scale_overcapacity", "voluntary_mode", "cost_beta_change", "which", "beta_fact", "start_point","n_cps", "local_algo", "fit_cps", "icu_prev_lim", "label", "n_local", "n_global"))
    
    #results <- results[!duplicated(results[, .(R, sev, import, cap, vaccination_strat, scale_inf,beta_fact, include_tisk,l,scale_ini,voluntary_mode, scale_overcapacity, cost_beta_change, time, which, start_point, fit_cps, icu_prev_lim, label, n_cps, local_algo)])]
    results <-results %>% left_join(res,  by=c("R"="R","sev","which", "start_point", "fit_cps","import", "cap", "vaccination_strat", "scale_inf", "include_tisk", "l", "scale_ini", "voluntary_mode", "cost_beta_change", "beta_fact", "scale_overcapacity", "n_cps", "local_algo","icu_prev_lim", "label", "n_local", "n_global"))
    res$shape <- factor(res$shape, levels=c("Open",  "Closed", "Cut Overshoot","Medium - Low CS", "Uncertain"))
    res$strategy1 <- factor(res$strategy1, levels=c("No Interventions",  "Suppression", "Mitigation","Medium - Low CS", "Uncertain"))
    res$strategy2 <- factor(res$strategy2, levels=c("No Interventions",  "Suppression", "Cut Overshoot","Medium - Low CS", "Uncertain"))
    res$strategy3 <- factor(res$strategy3, levels=c( "Suppression","Cut Overshoot", "No Interventions",  "Mitigation", "Medium - Low CS", "Uncertain"))

    res$strategy3_plot <- paste0(res$strategy3, "    ")
    res$strategy3_plot <- factor(res$strategy3_plot, levels=c( "Suppression    ","Cut Overshoot    ", "No Interventions    ",  "Mitigation    ", "Medium - Low CS", "Uncertain"))

    res$strategy4 <- factor(res$strategy4, levels=c( "Suppression","Cut Overshoot", "No Interventions",  "Mitigation", "Medium - Low CS", "Uncertain"))

    res$strategy4_plot <- paste0(res$strategy4, "    ")
    res$strategy4_plot <- factor(res$strategy4_plot, levels=c( "Suppression    ","Cut Overshoot    ", "No Interventions    ",  "Mitigation    ", "Medium - Low CS", "Uncertain"))


    return(list(results=results, betas=betas, res=res))


}

make_row <- function(res, results, x_R, x_sev,x_which, def, label, pop=5.37e6, icu_max_prev=350){

  res_i <- head(do.call(select_def, list(res, def)) %>% filter(R==x_R, sev==x_sev & which==x_which),1)
  results_i <- do.call(select_def, list(results, def)) %>% filter(R==x_R, sev==x_sev & which==x_which)


  return(data.frame(R=x_R, sev=x_sev,strategy=x_which, total_loss=res_i$total, QALY_loss=res_i$lost_qaly, HC_loss=res_i$icu_loss+res_i$hosp_loss, beta_change_loss=res_i$beta_change_cost, prod_loss_illness=res_i$prod_loss_illness, intervention_loss=res_i$beta_cost, TTIQ_loss=res_i$loss_tisk, average_intervention=1-res_i$beta_mean, fraction_infected=res_i$frac_infected*100,
        hospitalisations=tail(results_i,1)$tot_hosp/pop*100000, ICUs=tail(results_i,1)$tot_resp/pop*100000, deaths=tail(results_i,1)$D/pop*100000, max_icu_prevalence=max(results_i$resp)/icu_max_prev
  
   ))

}


make_table <- function(res, results, Rs, sevs, whichs,defs, labels){

  rows <- list()
  for(i in 1:length(Rs)){
    rows[[i]] <- make_row(res, results, Rs[i], sevs[i],whichs[i], defs[[i]], labels[i])
  }
  return(rbindlist(rows))
}

read_results <- function(filename=NULL, folder_path=NULL, max_n=100){
  if(!is.null(folder_path)){
    files <- list.files(folder_path)
    print(files)
    if(length(files)> max_n) files <- files[1:max_n]
    raw_results <- parallel::mclapply(files, function(file){
      if(grepl(filename, file)){
        print(file)
        tmp <- readRDS(fs::path(folder_path, file))


        print(paste(file, "-read"))
        nulls <- unlist(lapply(tmp, function(x) is.null(x)))  
        if(sum(nulls)==length(nulls)){
          prepd <- data.table()
        }else{               
          prepd <- prepare_result(tmp)
        }
        print(paste(file, "-done"))
        return(prepd) 

      }}, mc.cores=30)
      
  }else{
    raw_results <- readRDS(filename)
  }
  res <- rbindlist(lapply(raw_results, function(x) x$res))
  betas <- rbindlist(lapply(raw_results, function(x) x$betas))
  results <- rbindlist(lapply(raw_results, function(x) x$results), fill=TRUE)
  return(list(res=res, betas=betas, results=results))
}

plot_desc <- function(res_i){
  q1 <- ggplot(res_i%>% filter(which=="best"  )) + geom_histogram(aes(x=av_rel_intervention, y=stat(count/sum(count)), fill=strategy4_plot), bins=50) +theme + xlab("Average relative intervention") + ylab("Fraction of cases") +scale_y_continuous(labels = scales::percent_format())  + geom_segment(aes(x=0.01, xend=0.01, y=0, yend=0.4))  + geom_segment(aes(x=0.6, xend=0.6, y=0, yend=0.4))  + labs(fill="Strategy") + gg_fill + theme(legend.position="bottom", legend.title=element_blank())
  q2 <- ggplot(res_i %>% filter(which=="best")) + geom_point(aes(x=av_rel_intervention, y=frac_inf_rel, color=strategy4),size=4) + ylab("Relative fraction infected") + xlab("Average relative intervention") + theme + labs(color="Strategy") + gg_color 
  q3 <- ggplot(res_i %>% filter(which=="best")) + geom_point(aes(x=frac_max_inc, y=frac_inf_rel, color=strategy4),size=4) + ylab("Relative fraction infected") + xlab("Relative max incidence") + theme + labs(color="Strategy") + gg_color  #+ scale_fill_brewer("Dark")

  bottom_row <- cowplot::plot_grid(q2 + theme(legend.position="none"), q3 + theme(legend.position="none"), ncol=2, labels = c("B", "C"), label_size = 36)


  cowplot::plot_grid(q1, bottom_row, nrow=2, ncol=1, rel_heights=c(0.7,1), labels=c("A", ""), label_size = 36, scale=0.95)#, labels=c("Relativ fraction of infections", "Relative peak infections", "Relative mean beta", "Mean infection time (days)"), label_size=25)
}


plot_desc_comp <- function(res_i){
  q1 <- ggplot(res_i %>% filter(which=="best")) + geom_point(aes(x=frac_max_inc, y=frac_inf_rel, color=strategy3),size=4) + ylab("Relative fraction infected") + xlab("Relative max incidence") + theme + labs(color="Strategy") + gg_color  #+ scale_fill_brewer("Dark")
  q2 <- ggplot(res_i %>% filter(which=="best")) + geom_point(aes(x=frac_max_inc, y=frac_inf_rel, color=strategy4),size=4) + ylab("Relative fraction infected") + xlab("Relative max incidence") + theme + labs(color="Strategy") + gg_color  #+ scale_fill_brewer("Dark")

  bottom_row <- cowplot::plot_grid(q1 + theme(legend.position="none"), q2 + theme(legend.position="none"), ncol=2, labels = c("Old", "New"), label_size = 36)
}


cost <- function(res, x_import, capacity, x_vaccination_strat=1, x_scale_inf=1.0, suf="", steps=c(0,20,50,100), title=""){
  f <- res %>% filter( import==x_import, cap==capacity, vaccination_strat==x_vaccination_strat, scale_inf==x_scale_inf) %>% group_by(R, sev) %>%
    mutate(cost=min(total))
  q <- ggplot(f) +  geom_tile(aes(x=R,y=sev,fill=cost)) + theme_minimal()  + theme(text = element_text(size=44)) + xlab("Initial Reproduction Number")  +ylab("Effetive Severity")+ labs(fill="Cost(bNOK)") + scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0))) + scale_fill_steps(breaks=steps, low = "yellow", high = "red") + guides(fill = guide_colourbar(barwidth = 3, barheight = 25))# + geom_hline(yintercept=4)
  if(title!="")
    q <- q + ggtitle(title)
  ggsave(glue::glue("cost_{x_import}_{capacity}_{x_vaccination_strat}_{x_scale_inf}_{suf}.png"), width=20, height=15)
  q
}

plot_strat <- function(res, pad=FALSE, key="strategy4_plot"){

  
  ggplot(res %>% filter(which=="best")) + geom_raster(aes(x=R, y=sev, fill=get(key))) + theme + scale_fill_brewer(palette="Set1") + ylab("Effective Severity") + xlab("Initial Reproduction Number")+ scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0))) 
}

plot_var <- function(var, res, suf="", steps=c(0,20,50,100), title="", label="", barheight=25, color_low="yellow", color_high="red"){
#  f <- res %>% filter( import==x_import, cap==capacity, vaccination_strat==x_vaccination_strat, scale_inf==x_scale_inf, scale_ini==x_scale_ini, l==x_l , voluntary_mode==x_voluntary_mode, include_tisk==x_include_tisk) %>% group_by(R, sev) %>%filter(total==min(total))
  f <- res %>% group_by(R, sev) %>%filter(total==min(total))
  
  q <- ggplot(f) +  geom_tile(aes(x=R,y=sev,fill=get(var))) + theme_minimal() + theme(text = element_text(size=35)) + xlab("Initial Reproduction Number")  +ylab("Effective Severity")+ labs(fill=label) + scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0))) + scale_fill_steps(breaks=steps, low = color_low, high = color_high) + guides(fill = guide_colourbar(barwidth = 3, barheight = barheight))# + geom_hline(yintercept=4)
#  ggsave(glue::glue("frac_inf_{x_import}_{capacity}_{x_vaccination_strat}_{x_scale_inf}_{suf}.png"), width=20, height=15)
  q
}

plot_var_fact <- function(var, res, x_import, capacity, x_vaccination_strat=1, x_scale_inf=1.0,x_l=270,x_voluntary_mode=1, x_scale_ini=5, suf="", steps=c(0,20,50,100), title=""){
  f <- res %>% filter( import==x_import, cap==capacity, vaccination_strat==x_vaccination_strat, scale_inf==x_scale_inf, voluntary_mode==x_voluntary_mode, scale_ini==x_scale_ini, l==x_l)
  print(f)
  f <- f %>% group_by(R, sev) %>%filter(total==min(total))
                                                                                                                                       
  q <- ggplot(f) +  geom_tile(aes(x=R,y=sev,fill=get(var))) + theme_minimal()  + theme(text = element_text(size=44)) + xlab("Initial Reproduction Number")  +ylab("Effective Severity")+ labs(fill="Variable") + scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0))) + guides(fill = guide_colourbar(barwidth = 3, barheight = 25)) + ggtitle(title)# + geom_hline(yintercept=4)
#  ggsave(glue::glue("frac_inf_{x_import}_{capacity}_{x_vaccination_stra}_{x_scale_inf}_{suf}.png"), width=20, height=15)
  q
}


plot_diff <- function(res, res2, steps=c(0,1,2,3,4,5,6), trans="lin", min_diff=1, title="", label="Reduction(bNOK)", filename=NULL, barheight=5, diverging=FALSE, color_mid="white", color_low="blue", color_high="red", variable="total", direction="vertical"){
  f <- res %>% group_by(R, sev) %>%
    summarize(cost=min(get(variable)))
  f2 <- res2 %>% group_by(R, sev) %>% 
    summarize(cost2=min(get(variable)))
  f3 <- merge(f,f2, by=c("R","sev"))

  f3$diff <- f3$cost2 - f3$cost
  f3$diff[f3$diff < min_diff] <- 0
  print(f3)
  q <- ggplot(f3) +  geom_raster(aes(x=R,y=sev,fill=diff)) + theme_minimal() + theme(text = element_text(size=44)) + xlab("Initial Reproduction Number")  +ylab("Effective Severity")+ labs(fill=label) + scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0))) +  guides(fill = guide_colourbar(barwidth = 3, barheight = barheight, direction=direction))# + geom_hline(yintercept=4)
  if(diverging){
    #q <- q +  scale_fill_steps2(breaks=steps, low = color_low,mid=color_mid, high = color_high, limits=c(min(f3$diff), max(f3$diff))) 
    q <- q +  scale_fill_gradient2(low=color_low, mid=color_mid, high=color_high, midpoint=0)#,  limits=c(-max(abs(f3$diff)), max(abs(f3$diff))))

  }else{
    q <- q + scale_fill_steps(breaks=steps, low = color_mid, high = color_high) 
  }
  if(title!="")
    q <- q + ggtitle(title)

  if(!is.null(filename)){
    ggsave(filename, width=25, height=15)
  }
  q
}



plot_unc <- function(x_scale_ini=5){
  print(x_scale_ini)
  extra_lockdown <- res %>% select_def(scale_ini=x_scale_ini, l=256) %>% filter(which=="lockdown") %>% group_by(R, sev) %>% summarise(V1 = min(beta_cost)/ 256*14) %>% pull(V1)
  f <- res %>% select_def(scale_ini=x_scale_ini) %>% filter(which=="best")
  f <- f[!duplicated(f)]
  f <- f %>% group_by(R, sev) %>% summarise(strategy1=first(strategy4), total=min(total))
  f2 <-  res %>% select_def(scale_ini=as.character(glue::glue("R{x_scale_ini}")), l=256) %>% filter(which=="best") %>% group_by(R, sev) %>% summarise(tot2=min(total), strat2=first(strategy4))
  f3 <- res %>% select_def(scale_ini=x_scale_ini, l=256) %>% filter(which=="best") %>% group_by(R, sev, strategy4) %>% summarise(total=min(total)) %>% ungroup() %>% mutate(tot3=total + extra_lockdown) %>% select(R, sev, tot3, strat3=strategy4)

  plot <- f %>% left_join(f2) %>% left_join(f3)
  plot <- plot %>% mutate(diff=ifelse(strategy1=="Suppression", tot2 - total, tot3 - total))

  steps <- c(0,5, 10, 15,17, 20,30)
  barheight <- 20
  q1 <- ggplot(plot) +  geom_raster(aes(x=R,y=sev,fill=diff)) + theme_minimal() + theme(text = element_text(size=44)) + xlab("Initial Reproduction Number (R)")  +ylab("Effective Severity")+ labs(fill="bnNOK") + scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, 0))) + scale_fill_steps(breaks=steps, low = "white", high = "red") + guides(fill = guide_colourbar(barwidth = 3, barheight = barheight))# + geom_hline(yintercept=4)

  q2 <- add_contours(plot_strat(res %>% select_def(l=256, scale_ini=as.character(glue::glue("R{x_scale_ini}")))), res %>% select_def(scale_ini=x_scale_ini), manual=manual)
  q3 <- add_contours(plot_strat(res %>% select_def(l=256, scale_ini=x_scale_ini)), res %>% select_def(scale_ini=x_scale_ini), manual=manual)



  bottom_row <- cowplot::plot_grid(q2 + theme(legend.position="none"), q3 + theme(legend.position="none"), ncol=2, labels = c("B) 2 Weeks of No interventions", "C) 2 Weeks of Suppression"), vjust=1.0, hjust=0, label_size=36, scale=0.9)


  legend_b <- cowplot::get_plot_component(
    q2 + 
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom", legend.title=element_blank(), legend.spacing.x=unit(1.3, "cm")), "guide-box", return_all=TRUE
  )[[3]]

  cowplot::plot_grid(q1, bottom_row,legend_b, nrow=3, ncol=1, rel_heights=c(1.0,1, 0.1), label_size = 36, labels=c("A)", ""), scale=c(0.9,1))
  ggsave(glue::glue("fig_uncertainty_2_weeks_{x_scale_ini}.png"), width=20, height=16)

}