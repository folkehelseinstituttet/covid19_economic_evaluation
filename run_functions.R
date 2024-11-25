
cost_lin <- function(X, start=20, end=80, tot_cost_day=25/30){
  res <- tot_cost_day/(end - start)*(X-start)
  res[res<0] <- 0
  return(res)
}
cost_quad_3o4 <- function(X, cost_params){
  c1 <- cost_params$quad_cost_npi_1*100/75
  c2 <- cost_params$quad_cost_npi_2*(100/75)^2
  res <- (c1*X + c2*X^2)/30
  res[res<0] <- 0
  return(res)
}
cost_quad_100 <- function(X, cost_params){
  res <- (cost_params$quad_cost_npi_1*X + cost_params$quad_cost_npi_2*X^2)/30
  res[res<0] <- 0
  return(res)
}




get_beta <- function(beta, n, l){
  # Assume first n values are change points and next n+1 values are beta values
  if(n==0) return(rep(beta, l))
  cps <- sort(round(beta[1:n]))
  out <- rep(beta[n+1], cps[1])
  if(n > 1){
    for(i in 2:n){
      if(cps[i] != cps[i-1]){
        out <- c(out, rep(beta[n+i], cps[i] - cps[i-1]))
      }
    }
  }
  out <- c(out, rep(beta[length(beta)], l - cps[n]))
  return(out)
}




quality_red_func <-  function(x, max_loss=0.3, max_loss_frac=3, age_dist=c()){
  return(0)
}

quality_red_func_icu <-  function(x, max_loss=0.75, normal_risk=0.2,
                                  qaly_lost=10.4, increase_per_percent=0.007, age_dist=c()){
  
  max_loss_frac <- (max_loss/normal_risk -1 )/(increase_per_percent)/100 +1
  if(x < 1) return(0)
  if(x > max_loss_frac) return((max_loss-normal_risk)*qaly_lost)
  return( qaly_lost*increase_per_percent*((x-1)*100)*normal_risk)
}

update_severity <- function(params, severity, change_icu_prob=3.8, change_los_hosp=2.3, base_sev=15){
    params$hosp_prob <- params$hosp_prob*severity
    params$icu_prob <- params$icu_prob*min((1 + (change_icu_prob-1)*(severity-1)/(base_sev - 1)), change_icu_prob)
    params$length_hosp <- params$length_hosp*min((1 + (change_los_hosp-1)*(severity-1)/(base_sev - 1)), change_los_hosp)
    params$prob_death_hosp <- params$prob_death_icu <- params$prob_death_non_hosp <- params$prob_death_non_hosp*severity
    return(params)

}

get_basic_params <- function(param_file, daily_import=0,n_vac=1,
                               vac_pars=list(rr_inf = c(1,0.7),
                                             rr_hosp = c(1, 0.6),
                                             rr_death = c(1, 0.6),
                                             rr_icu = c(1, 0.6),
                                             rr_los_hosp = c(1, 0.85),
                                             rr_inf_asymp = c(1,0.8),
                                             rr_trans = c(1,1)), L=500) {

  params <- read_param_file_simplified(param_file)
  n_strain <- 1
  N <- 9
  params <- c(params,
              list(
                N_steps=L,
                n_vac=n_vac,
                n_strain=n_strain,
                dt=0.5,
                T_waning=array(1e10, dim=c(N, n_vac)),
                vaccinations=array(0,dim=c(L, N, n_vac)),
                beta_day=matrix(0.2, ncol=N, nrow=L),
                beta_strain=rep(1, n_strain),
                cross_protection=matrix(1, ncol=n_strain, nrow=n_strain),
                n=9,
                S_ini=matrix(c(get_age_groups(), get_age_groups()*0)[1:(9*n_vac)], nrow=N, ncol=n_vac),
                import_vec=array(0, dim=c(L,N, n_vac, n_strain)),
                I_ini=array(20, dim=c(N,n_vac,n_strain)),
                I_imp_ini=array(0, dim=c(N,n_vac,n_strain)),
                Ea_ini=array(0, dim=c(N,n_vac,n_strain)),
                Es_ini=array(0, dim=c(N,n_vac,n_strain)),
                A_ini=array(0, dim=c(N,n_vac,n_strain)),
                P_ini=array(0, dim=c(N,n_vac,n_strain)),
                H_ini=array(0, dim=c(N,n_vac,n_strain)),
                ICU_H_ini=array(0, dim=c(N,n_vac,n_strain)),
                ICU_R_ini=array(0, dim=c(N,n_vac,n_strain)),
                ICU_P_ini=array(0, dim=c(N,n_vac,n_strain)),
                B_D_ini=array(0, dim=c(N,n_vac,n_strain)),
                B_D_H_ini=array(0, dim=c(N,n_vac,n_strain)),
                B_D_ICU_ini=array(0, dim=c(N,n_vac,n_strain)),
                R_ini=array(0, dim=c(N,n_vac,n_strain)),
                D_ini=array(0, dim=c(N,n_vac,n_strain)),
                tot_infected_ini=array(0, dim=c(N,n_vac,n_strain)),
                tot_hosp_ini=array(0, dim=c(N,n_vac,n_strain)),
                tot_resp_ini=array(0, dim=c(N,n_vac,n_strain)),
                tot_vac_ini=array(0, dim=c(N, n_vac)),
                tot_vac_adm_ini=array(0, dim=c(N, n_vac)),
                beta_norm=get_age_groups(),
                reg_pop_long=get_age_groups(),
                N_regions=1,
                waning_immunity_vax = array(100000, dim=c(N,n_vac,n_strain)),
                waning_inf = 100000,
                age_groups=N,
                beta_mode=1
              )
              )
  if(n_vac==2){
    params$I_ini[,2,] <- 0
  }
  params$import_vec[1:L, 1:7, 1, 1] <- round(daily_import/7)
  basic_params <- fix_params(params, N, n_vac, n_strain, vac_pars)
}


beta_to_daily <- function(beta, l , n_cps, equaly_spaced_cps=TRUE){
  if(!equaly_spaced_cps){
    new_beta <- get_beta(beta,n_cps, l)
  }else{
     if(n_cps == 0) {
      new_beta <- rep(beta, l)
    }else{
    cps <- ( l*(1:(n_cps)))%/%(n_cps+1)
    lagged <- lag(cps)
    lagged[1] <- 0

    spacings <- c(cps, l ) - c(lagged,max(cps))

    new_beta <- rep(beta, spacings)
     
  }
  }
  return(new_beta)
}

prepare_params <- function(beta, R, severity, import_x,vaccination_strat, run_params,cost_params, 
              n_cps, l,dt,scale_ini=1, cool_down=0, cool_down_R=0, voluntary_mode=0, beta_type=NULL, thin=FALSE, equaly_spaced_cps=NULL, beta_fact=beta_fact, fixed_cps=c(), ihr=NULL, ifr=NULL){
  orig_beta <- c(fixed_cps,beta)
  beta <- beta_to_daily(orig_beta, l - cool_down, n_cps, equaly_spaced_cps=equaly_spaced_cps)
 
  if(cool_down > 0){
      beta[(l-cool_down+1):l] <-  (1 - (1 - 1/R) / run_params$combine_param_1)*cool_down_R
  }
  if(sum(is.na(beta)> 0)){
    print("ERROR")
    print(beta)
    print(fixed_cps)
    print(l)
    print(l - cool_down)
    print((l - cool_down)/(n_cps+1))

  }
  
   import <- import_x
  L <- l/dt 
  daily <- rep(beta, each=1/dt)
  n_vac <- 1
  if(vaccination_strat>1){
    n_vac <- 2
  }
  new_params <- get_basic_params(param_file, daily_import=import, L=L, n_vac=n_vac)
  if(!is.null(ihr) & n_vac==1){
    new_params$hosp_prob[,,1] <- ihr    
    new_params$prob_death_non_hosp[,,1] <- ifr 
    new_params$prob_death_hosp[,,1] <- ifr 
    new_params$prob_death_icu[,,1] <- ifr 
  }
  new_params$interpolation_time <- 1:L
  if(n_cps > 0 & orig_beta[1]>1 & !all(round(orig_beta[1:n_cps]) == orig_beta[1:n_cps])){
 
     spots <- round(round(orig_beta[1:n_cps]))
    new_params$interpolation_time[spots] <- orig_beta[1:n_cps]
  
  }
  new_params$include_waning <- 0
  new_params$daily <- daily
  if(thin){
    new_params$daily <- beta
  }
  
  
  
  if(scale_ini %in% c("R5", "R10")){
    growth_rates <- readRDS("parameter_files/growth_rates.RDS")
    fact <- 10
    use_R <- R
    if("global_R_scaling" %in% names(run_params)){
      use_R <- R*run_params$global_R_scaling
    }
    if(scale_ini =="R5") fact <- 5
    scale_ini <- fact*exp(growth_rates[[as.character(use_R)]]*14)
  }
  new_params$I_ini <- new_params$I_ini*scale_ini

  new_params$dt <- dt
  new_params$S_div_vac <- array(0, dim=c(9,n_vac))
  beta_ini <- new_params$beta_day[1,]
  beta_1 <- fix_beta_large(new_params, new_params$S_ini, new_params$I_ini, 1, beta=beta_ini, use_eig=TRUE)
  beta_0 <- 0.2*beta_1*R
  new_params$beta_day <- matrix(beta_0*daily*beta_fact, ncol=9, nrow=L)

  new_params <- update_severity(new_params, severity, change_icu_prob=run_params$change_icu_prob, change_los_hosp=run_params$change_los_hosp)
  if(vaccination_strat == 2){
    vac_uptake <- run_params$vac_uptake
    tot <- new_params$S_ini[,1]
    new_params$S_ini[,1] <- tot*(1-vac_uptake)
    new_params$S_ini[,2] <- tot*vac_uptake
  }
  beta_after_vac <- fix_beta_large(new_params, new_params$S_ini, new_params$I_ini, 1, beta=beta_ini*beta_fact, use_eig=TRUE)
  
  new_params$reduction_in_beta_from_S <- beta_1/beta_after_vac
  
  new_params$beta_mode <- 4
  if(voluntary_mode==0){
    exp_health_loss <- -expected_qaly_loss(new_params, severity, qaly_cost=cost_params$qaly_cost)*5*1e9*100
  }else if(voluntary_mode== 1 ){
    exp_health_loss <- expected_qaly_loss(new_params, severity, qaly_cost=cost_params$qaly_cost)*1e9*run_params$combine_param_1
  }else if(voluntary_mode==2){
    exp_health_loss <- expected_qaly_loss(new_params, severity, qaly_cost=cost_params$qaly_cost)*1e9*2*run_params$combine_param_1
  } else if(voluntary_mode==3){
    exp_health_loss_orig <- expected_qaly_loss(new_params, severity, qaly_cost=cost_params$qaly_cost)*1e9*1*run_params$combine_param_1
    #mean_exp <- mean(exp_health_loss_orig)
    #exp_health_loss <- (exp_health_loss_orig+mean_exp)/2
    exp_health_loss <- (1*exp_health_loss_orig + 1*(new_params$mixing_matrix[, ] /rowSums(new_params$mixing_matrix)) %*% exp_health_loss_orig)/1
 } else if(voluntary_mode==4){
    exp_health_loss_orig <- expected_qaly_loss(new_params, severity, qaly_cost=cost_params$qaly_cost)*1e9*1*run_params$combine_param_1#*0.5
    #mean_exp <- mean(exp_health_loss_orig)
    #exp_health_loss <- (exp_health_loss_orig+mean_exp)/2
    exp_health_loss <- (4/5*exp_health_loss_orig + 1/5*(new_params$mixing_matrix[, ] /rowSums(new_params$mixing_matrix)) %*% exp_health_loss_orig)/1
  }else if(voluntary_mode==5){
    exp_health_loss_orig <- expected_qaly_loss(new_params, severity, qaly_cost=cost_params$qaly_cost)*1e9*1*run_params$combine_param_1#*0.5
    #mean_exp <- mean(exp_health_loss_orig)
    #exp_health_loss <- (exp_health_loss_orig+mean_exp)/2
    exp_health_loss <- (2/3*exp_health_loss_orig + 1/3*(new_params$mixing_matrix[, ] /rowSums(new_params$mixing_matrix)) %*% exp_health_loss_orig)/1
  }else if(voluntary_mode==6){
    exp_health_loss_orig <- expected_qaly_loss(new_params, severity, qaly_cost=cost_params$qaly_cost)*1e9*1*run_params$combine_param_1
    #mean_exp <- mean(exp_health_loss_orig)
    #exp_health_loss <- (exp_health_loss_orig+mean_exp)/2
    exp_health_loss <- (0*exp_health_loss_orig + sum(get_age_groups()* exp_health_loss_orig)/sum(get_age_groups()))
  }

  new_params$expected_health_loss <- array(exp_health_loss, dim=c(new_params$n,new_params$n_vac))
  T <- cost_params$time_voluntary
  cost_idle <- cost_params$cost_idle / 31
  a <- 4/(0.9*cost_idle*T)
  if(voluntary_mode==0){
    a <- 100
  }
  new_params$spont_behav_change_params <- c(beta_0*beta_fact,T, cost_idle,a , 
       run_params$combine_param_1, run_params$combine_param_2,cost_params$max_cost_idle*T/31)
  new_params$spont_behav_mode <- 2

   steps_per_day <- 1/dt
  imp_vec <- array(0, dim=c(l,9,n_vac,1))
  imp_vec[,,1,] <- import_x/7
  imp_vec <- abind::abind(imp_vec, array(0, dim=c(1, dim(imp_vec)[2:4])), along=1)
  rows <- c()
  for(i in 1:(dim(imp_vec)[1]-1)){
    rows <- c(rows, c(i, rep(dim(imp_vec)[1], steps_per_day-1)))
  }
  a <- imp_vec[rows,,,]
  if(n_vac ==1){
    dim(a) <- c(dim(a), 1,1)
  }else{
    dim(a) <- c(dim(a), 1)
  }
  new_params$import_vec <-a

  return(new_params) 
}

calc_loss_prod_loss_and_tisk <- function(res, new_params, cost_params, dt,l, include_tisk=1, cool_down=0, reduced_cost_from_NPI=0 ){
  setDT(res)			     
  max_beta_tisk <- cost_params$max_beta_tisk
  daily <- new_params$daily[1:((l-cool_down)/dt)]
  add_tisk <- numeric(length(daily))
  add_tisk[daily < 1.0 & daily > max_beta_tisk] <- (1 - daily[daily < 1.0 & daily > max_beta_tisk])/(1 - max_beta_tisk)
  add_tisk[daily <= max_beta_tisk] <- 1
  add_tisk <- add_tisk*include_tisk
  prod_loss <- calc_prod_loss(res[1:((l-cool_down)/dt),], new_params, tisk_on=add_tisk, reduced_cost_from_NPI=reduced_cost_from_NPI)
  is_tisk_on = add_tisk > 0
  fixed_costs <- sum(cost_params$tisk_per_day * is_tisk_on)*dt
  variable_costs <- sum(res[time <=(l-cool_down) & time %% 1 ==0, mean(incidence), by=t]$V1*.colMeans(add_tisk, 1/dt, l-cool_down))*cost_params$cost_per_tisk
  return(list(prod_loss=prod_loss, direct_tisk=fixed_costs + variable_costs))
}

calc_loss_hosp_capacity <- function(res, capacity, cost_params, dt, l, scaling){
  
  hosp_loss <- c()
  icu_loss <- c()
  for(i in unique(res$sim)){
    hosp_loss <- c(hosp_loss, calc_qual_post_loss( res %>% dplyr::filter(sim==i, time%%1==0) %>% pull(hosp_incidence),
    res %>% dplyr::filter(sim==i, time%%1==0) %>% pull(hosp),
      capacity[1],cost_params$acute_hosp, cost_params$elective_hosp, cost_params$los_hosp_acute,
      cost_params$los_hosp_elective, cost_params$qaly_loss_per_day_hosp,quality_red_func,
                                                  dt=1, ret_full=F,loss_per_QALY=cost_params$qaly_cost*scaling))
    res_tmp <- res %>% filter(sim==i)
    icu_quality_red_func <- function(x, age_dist=c()) quality_red_func_icu(x, max_loss=cost_params$icu_max_loss,
                             normal_risk=cost_params$icu_normal_risk,
                                  qaly_lost=cost_params$icu_qaly_lost,
                                    increase_per_percent=cost_params$icu_risk_increase_per_percent, age_dist=age_dist)
    resp_incidence <- c(res_tmp[2:(l/dt),] %>% dplyr::pull(tot_resp) - res_tmp[1:(l/dt-1),] %>% dplyr::pull( tot_resp),0)
    resp_incidence_day <- .colSums(resp_incidence, 1/dt, l)
    icu_loss <- c(icu_loss, calc_qual_post_loss(resp_incidence_day,res %>% dplyr::filter(sim==i, time%%1==0) %>% pull(resp),
               capacity[2],cost_params$acute_icu, cost_params$elective_icu, cost_params$los_icu_acute, 
               cost_params$los_icu_elective, cost_params$qaly_loss_per_day_icu, icu_quality_red_func, 
               dt=1, ret_full=F,  loss_per_QALY=cost_params$qaly_cost*scaling))
   }
  return(list(hosp_loss=hosp_loss, icu_loss=icu_loss))
}

calc_loss_NPI <- function(res, new_params, run_params, cost_params,dt, l ,include_tisk=1, cool_down=0){
  setDT(res)
  z <- as.matrix(res[1:((l-cool_down)/dt), .(mean(get("contact_change[1]")),
                     mean(get("contact_change[2]")),
                     mean(get("contact_change[3]")),
                     mean(get("contact_change[4]")),
                     mean(get("contact_change[5]")),
                     mean(get("contact_change[6]")),
                     mean(get("contact_change[7]")),
                     mean(get("contact_change[8]")),
                     mean(get("contact_change[9]"))

                     ), by=t] %>% select(-t))
  
  r <- 1 - new_params$daily[1:((l-cool_down)/dt)]
  cost_z <- cost_quad_100(100*z, cost_params)*dt
  if(include_tisk){
    cost_p <- cost_quad_3o4(100* pmax(r - (1 - cost_params$max_beta_tisk), 0),cost_params)*dt
  }else{
    cost_p <- cost_quad_100(100* r, cost_params)*dt
  }
 
  max_intervention <- (2*run_params$combine_param_1 - run_params$combine_param_2)
  total_costs <- 1/max_intervention*(
     run_params$combine_param_1*(cost_p + cost_z) - run_params$combine_param_2/(cost_params$max_npi_cost/30)*(cost_p*cost_z))
  return(total_costs)
}
 
 

calc_beta_change_cost <- function(beta, cost_param){
   if(length(beta)==1) return(0)
  abs_diffs <- abs(beta[2:length(beta)] - beta[1:(length(beta)-1)])
  abs_diffs[beta[2:length(beta)]==1] <- 0
    return(sum((1/(1 + exp(-abs_diffs*5)) -0.5 )*2)*cost_param)
}


#' Main function to run the model and calculate the costs
cost <- function(beta=1, R=1.2, severity=1, import_x=0, vaccination_strat=1, capacity=c(11000, 350),scale_inf=1,
                 n_cps=0, run_params=NULL,cost_params=NULL, include_tisk=1, l=168, do_print=FALSE, cost_beta_change=1,
                 ret_both=F, ret_res=F, dt=1/50, cool_down=150, cool_down_R=0.8,voluntary_mode=1,ihr=NULL, ifr=NULL,
                 n_parts=10, n_threads=1, deterministic_stoch_model=FALSE, scale_ini=1, beta_type=NULL, thin=TRUE, beta_fact=1,scale_overcapacity=1, fixed_cps=c(),
                 equaly_spaced_cps=TRUE, use_determinsitic_model=TRUE, ret_all=FALSE, refine_type="minimal", icu_prev_lim=NULL){
 

  if(do_print){
    print(glue::glue("Called with"))
    print(beta)
  }
  if(use_determinsitic_model){
    dt <- 1
  }
  l <- l+cool_down
  thinned_dt <- dt
  thin_before_refine <- NULL
  if(thin){
    thinned_dt <- 1
    thin_before_refine <- 1
  }

  new_params <- prepare_params(beta, R, severity, import_x,vaccination_strat, 
         run_params, cost_params, n_cps, l,dt,scale_ini=scale_ini, cool_down=cool_down, cool_down_R=cool_down_R, voluntary_mode = voluntary_mode, beta_type=beta_type, thin=thin, equaly_spaced_cps = equaly_spaced_cps, beta_fact=beta_fact, fixed_cps=fixed_cps, ihr=ihr, ifr=ifr)
  
   
  for(j in 1:5){
			  cost_tmp <- tryCatchLog::tryCatchLog({     
			    res <- metapop::run_params(new_params, L=l, n_parts, n_threads, deterministic=deterministic_stoch_model, thin_before_refine = thin_before_refine, use_determinsitic_model=use_determinsitic_model,
          refine_type=refine_type)
			    setDT(res)
          setkey(res, time)
	        res <- as.data.frame(res)
          res <- res %>% mutate( severity =severity, R=R,scenario=1, vac_scenario=vaccination_strat)
          hosp_icu_loss <- calc_loss_hosp_capacity(res, capacity, cost_params, thinned_dt,l, scaling=scale_inf*scale_overcapacity)
          hosp_loss <- hosp_icu_loss$hosp_loss
          icu_loss <- hosp_icu_loss$icu_loss
          qaly_loss <- calc_qaly(res, qaly_cost=cost_params$qaly_cost)
          tot_qaly_loss <- sum(qaly_loss$value_lost_QALYs_death_and_disease)
      
          beta_change_cost <- calc_beta_change_cost(c(1, new_params$daily[1:(l-cool_down)]), cost_params$cost_change_beta*cost_beta_change)

          total_cost_NPI <- calc_loss_NPI(res, new_params, run_params, cost_params, thinned_dt, l, include_tisk=include_tisk, cool_down=cool_down)
    
          cost_NPI <- sum(total_cost_NPI*new_params$beta_norm)/sum(new_params$beta_norm)


          prod_loss <- calc_loss_prod_loss_and_tisk(res, new_params, cost_params, thinned_dt,l,include_tisk=include_tisk, cool_down=cool_down, reduced_cost_from_NPI=total_cost_NPI)
        
          icu <- as.numeric(res %>% filter(t==max(t)) %>% select(starts_with("tot_resp_age")))
          hosp <- as.numeric(res %>% filter(t==max(t)) %>% select(starts_with("tot_hosp_age")))
      
          loss_by_age <- (hosp_loss*hosp/sum(hosp) + icu_loss*icu/sum(icu) + qaly_loss$value_lost_QALYs_death_and_disease)*scale_inf + colSums(total_cost_NPI*new_params$beta_norm)/sum(new_params$beta_norm) + prod_loss$prod_loss_illness + prod_loss$prod_loss_isolation + prod_loss$prod_loss_quarantine + prod_loss$direct_tisk
          costs_from_infections <- tot_qaly_loss*scale_inf +
                          sum(prod_loss$prod_loss_illness) + mean(hosp_loss)*scale_inf + mean(icu_loss)*scale_inf
      
          loss_tisk <- sum(prod_loss$prod_loss$prod_loss_isolation + prod_loss$prod_loss$prod_loss_quarantine)+ 
                                          prod_loss$direct_tisk
          data.frame(inf=costs_from_infections,
                              beta_cost=cost_NPI+
                                              loss_tisk + beta_change_cost,
                              lost_qaly=tot_qaly_loss*scale_inf,
                              hosp_loss=mean(hosp_loss)*scale_inf,
                          #   cost_no_self_reg=sum(cost_NPI$cost_no_self_reg*new_params$dt),
                              icu_loss=mean(icu_loss)*scale_inf,
                              prod_loss_illness=sum(prod_loss$prod_loss$prod_loss_illness),
                              beta_change_cost=beta_change_cost,
                              mean_inf_time=sum(res[, "incidence"]*res[, time])/sum(res[, "incidence"]),
                              max_inc=mean(res[, max(incidence), by=sim]$V1),
                              loss_interventions=cost_NPI,
                              frac_infected=mean(res[time==l-cool_down, tot_infected])/sum(new_params$beta_norm),
                              frac_infected_tot=mean(res[time==l, tot_infected])/sum(new_params$beta_norm),
                              loss_tisk=loss_tisk,
                              reduction_in_beta_from_S=new_params$reduction_in_beta_from_S
                                        
                              )
        },
           error=function(cond){
            print("Internal Error")
            print(cond)
          #  print(head(res))
            print("Internal Error")
				    return(NULL)  
           }
    , silent.warnings=TRUE)  
    if(!is.null(cost_tmp)) break;
  }
  cost_tmp$total <- cost_tmp$inf + cost_tmp$beta_cost
  if(!is.null(icu_prev_lim)){
    max_resp <- max(res$resp)
    if(max_resp>icu_prev_lim){
      cost_tmp$total<- cost_tmp$total +  1e5*(max_resp - icu_prev_lim )  
    }
    

  }
  
  
  if(do_print){
    print(glue::glue("Total cost {cost_tmp$total}"))
  }
  if(ret_all){
    res <- res %>% mutate(across(.cols=starts_with("beta_spont"), .fns= function(x) x/new_params$spont_behav_change_params[1], .names= "{.col}_norm"))
    setDT(res)
    return(list(cost=cost_tmp, results=res, age_breakdown=loss_by_age))
  }
  if(ret_both){
    return(cost_tmp)
  }
  if(ret_res){
    return(res)
  }
   return(-cost_tmp$total)
}




