create_cut_overshoot_change_points <- function(n, l, inc, peak_frac=0.8, early_frac=0.1, late_frac=0.1){
  if(n < 5){stop("Need more than five change points")}
  inc <- inc[1:l]
  max_t <- which(inc==max(inc))
  init <- min(which(inc > max(inc)*0.2))
  end <- max(which(inc > max(inc)*0.01))
  before_peak <- round(seq(init, max_t+1, length.out=ceiling(peak_frac*n*0.3)))
  n_after <- ceiling(peak_frac*n) - length(before_peak)
  after_peak <- round(seq(max_t + (end-max_t)/n_after, end, length.out=n_after))
  
  n_late =max(1, floor(late_frac*n))
  n_early = max(1, n - n_late - length(after_peak)-length(before_peak))
   early_cps <- round(seq(1, min(before_peak), length.out=n_early+2)[2:(n_early+1)])
  late_cps <- round(seq(max(after_peak), l, length.out=n_late+2)[2:(n_late+1)])
  return(list(cps=c(early_cps, before_peak,after_peak, late_cps), start=c(rep(1, n_early+length(before_peak)), rep(0, 3), rep(0.75, n_late+length(after_peak)-2) )))
}

create_cut_overshoot_cps <- function(l, inc, n_cps){
  if(n_cps < 2){
     stop("Need more than 1 cp")
  }

  inc <- inc[1:l]
  max_t <- min(which(inc>max(inc)*0.95))
  end <- max(which(inc > max(inc)*0.25))
  
  if(n_cps %in% c(2,3,4) & (end- max_t > 30)) end <- max_t + 18
  if(n_cps > (end - max_t +1)){

    print(max_t)
    frac <- (max_t-1) / (l-end + max_t -1)
    print(frac)
    n_before <- round(frac*(n_cps-end + max_t))
    n_after <- n_cps - end + max_t - n_before -1
    print(n_before)
    print(n_after)
    if(n_before == 0){
      before <- c()
    }else{
     before <- round(seq((max_t-1)/n_before, max_t-1, length.out=n_before))
    }
    cps <- max_t:end
    after <- round(seq(end + 1, l, length.out=n_after))
    cps <- c(before, cps, after)
    start <- c(rep(1, n_before+1), rep(0, length(max_t:end)), rep(0.75, n_after))

  }else{
    cps <- round(seq(max_t, end, length.out=n_cps))
    start <- c(1, rep(0, length(cps) -1), 0.75)
  }
  return(list(cps= cps, start=start))
}

create_flatten_cps <- function(l, inc, n_cps){

  return(list(cps=c(90, 120), start=c(1, 0.75, 1)))
}


#' Main function to optimise beta
optim_betas <-  function(R=3, sev=1,import=0, vaccination_strat=1, capacity=c(11000, 350), scale_inf=1,scale_ini=5,l=270, run_params=NULL, voluntary_mode=1, cost_beta_change=1, beta_fact=1,scale_overcapacity=1,
                 cost_params=NULL, n_cps, include_tisk=1, opt_control, n_parallel_inner=1, flatten_the_curve_thresholds=c(), icu_prev_lim=NULL, label=NULL, global_R_scaling=NULL){
    strat_columns <- list(R=R, sev=sev, cap=capacity, import=import, vaccination_strat=vaccination_strat, scale_inf=scale_inf, include_tisk=include_tisk, l=l, scale_ini=scale_ini, voluntary_mode=voluntary_mode, cost_beta_change=cost_beta_change, beta_fact=beta_fact, scale_overcapacity=scale_overcapacity, icu_prev_lim=icu_prev_lim)
 
    if(!is.null(global_R_scaling)){
        no_measures <- do.call(cost, c(strat_columns, list(ret_all=T, thin=TRUE, run_params=run_params, cost_params=cost_params)))
        scaling <- no_measures$results$`beta_spont_behaviour[1]_norm`[1]
        run_params$global_R_scaling <- scaling
        strat_columns$R <- R/scaling

    }else{
        run_params$global_R_scaling<- 1
    }
    run_list <- lapply(opt_control, function(x) list(opt_control=x, strat_columns=strat_columns ))




    results <- lapply(run_list, function(x) optim_beta_inner(x, run_params, cost_params))
    costs <- rbindlist(lapply(results, function(x) x$cost)) %>% mutate(which=ifelse(total==min(total), "best", "opt_res"))
    betas <- rbindlist(lapply(results, function(x) x$betas))
    betas$which = rep(costs$which,each= l)
    res <- rbindlist(lapply(results, function(x) x$results))
    res$which=rep(costs$which, each=l+150)
    no_measures <- do.call(cost, c(strat_columns, list(ret_all=T, thin=TRUE, run_params=run_params, cost_params=cost_params)))

    ls <- list()
    other_res <- list() 

    ls[[length(ls) +1]] <- no_measures$cost%>% mutate(which="nothing", n_cps=NA, beta_mean=1)
    other_res[[length(other_res) + 1]] <- no_measures$results %>% mutate(which="nothing",  n=NA)
    #Lockdown

    beta_from_scaling <- (1 - run_params$global_R_scaling)/run_params$combine_param_1
    beta_lockdown <- 1- (1/(R/run_params$global_R_scaling*no_measures$cost$reduction_in_beta_from_S) -1 + run_params$combine_param_1*beta_from_scaling)/(run_params$combine_param_2*beta_from_scaling - run_params$combine_param_1)

    if(beta_lockdown < 0) beta_lockdown <- 0
    if(beta_lockdown > 1) beta_lockdown <- 1
    ld <- results[[1]]$partial_cost_function(rep(beta_lockdown, l), n_cps=l-1, equaly_spaced_cps=TRUE, ret_all=TRUE)
    ls[[length(ls) +1]] <- ld$cost %>% mutate(which="lockdown", n_cps=NA, beta_mean=beta_lockdown)
    other_res[[length(other_res) + 1]] <- ld$results %>% mutate(which="lockdown",  n_cps=NA)
    #Flatten the curve
    continue <- length(flatten_the_curve_thresholds)>0
    beta <- beta_lockdown + 0.05
    i <- 1
    while(i <= length(flatten_the_curve_thresholds) & beta < 1){
        current_threshold <- flatten_the_curve_thresholds[i]
        max_hosp <- 0
        while(max_hosp < current_threshold & beta<=1){
        tmp_res <- results[[1]]$partial_cost_function(rep(beta, l), n_cps=l-1, equaly_spaced_cps=TRUE, ret_all=TRUE)$results
        max_hosp <- max(tmp_res[, mean(hosp), by=.(t)]$V1)
        beta <- beta + 0.01
        }
        threshold_res <- results[[1]]$partial_cost_function(rep(beta-0.01, l), n_cps=l-1, equaly_spaced_cps=TRUE, ret_all=TRUE)
        ls[[length(ls) +1]] <- threshold_res$cost %>% mutate(which=as.character(glue::glue("flatten_{current_threshold}")), n_cps=NA, beta_mean=beta)
        other_res[[length(other_res) + 1]] <- threshold_res$results %>% mutate(which=as.character(glue::glue("flatten_{current_threshold}")), n_cps=NA)
    i <- i + 1
    }
    strat_columns$cap <- strat_columns$cap[2]
    if(is.null(icu_prev_lim))strat_columns$icu_prev_lim <- 0

    strat_columns$R <- R
    new_columns <- as.data.frame(strat_columns)
    if(is.null(label)){
        x_label="default"
    }else{
        x_label <- label
    }
    ls <- rbind(costs, rbindlist(ls, fill=TRUE), fill=TRUE) %>% mutate(new_columns) %>% mutate(label=x_label)
    results <- rbind(res, rbindlist(other_res, fill=TRUE), fill=TRUE) %>% mutate(new_columns) %>% mutate(label=x_label)
    betas <- betas %>% mutate(new_columns)%>% mutate(label=x_label)


    return(list(costs=ls, 
            betas=betas,
            results = results))
}  

optim_beta_inner <-function(run, run_params, cost_params){
  results <-list()
  ls <- list()
  betas <- list()
  no_measures <- do.call(cost, c(run$strat_columns, list(ret_all=T, thin=TRUE, run_params=run_params, cost_params=cost_params)))
  inc_no_measures <- no_measures$results$inc
  n_cps <- run$opt_control$n_cps
  start_point <- run$opt_control$start_point
  start_point <- run$opt_control$start_point
  n_local <- run$opt_control$n_local
  local_algo <- run$opt_control$local_algo
  n_global <- run$opt_control$n_cps
  fit_cps <- run$opt_control$fit_cps

  global_algo <- "direct"

  lwr <- rep(0, n_cps+1)
  upr <- rep(1.0, n_cps+1)

   beta_from_scaling <- (1 - run_params$global_R_scaling)/run_params$combine_param_1
   beta_lockdown <- 1- (1/(run$strat_columns$R*no_measures$cost$reduction_in_beta_from_S) -1 + run_params$combine_param_1*beta_from_scaling)/(run_params$combine_param_2*beta_from_scaling - run_params$combine_param_1)

   fixed_cps <- c()
    equaly_spaced_cps <- TRUE
    if(!fit_cps){
        if(start_point=="nothing"){
          start <- rep(1, n_cps+1)
        }else if(start_point=="lockdown"){
          
          start <- rep(beta_lockdown,n_cps+1)
        }else if(start_point=="cut"){
          tmp <- create_cut_overshoot_cps(run$strat_columns$l, inc_no_measures, n_cps )
         
          start <- tmp$start
          fixed_cps <- tmp$cps
          equaly_spaced_cps <- FALSE
      }else if(start_point=="flatten"){
          tmp <- create_flatten_cps(run$strat_columns$l, inc_no_measures, n_cps )
          
          start <- tmp$start
          fixed_cps <- tmp$cps
          equaly_spaced_cps <- FALSE
      }
    }else{
      lwr <- c(rep(3, n_cps), lwr)
      upr <- c(rep(run$strat_columns$l, n_cps), upr)
      if(start_point=="nothing"){
        start <- c(seq(10, run$strat_columns$l-10, length.out=n_cps), rep(1, n_cps+1))
        }else if(start_point=="lockdown"){

          start <- c(seq(10, run$strat_columns$l-10, length.out=n_cps), rep(beta_lockdown, n_cps+1))
          
        }else if(start_point=="cut"){
          tmp <- create_cut_overshoot_cps(run$strat_columns$l, inc_no_measures, n_cps )
   
          start <- c(tmp$cps, tmp$start)
          if(start[1]<3) start[1] <- 3
          if(start[2]<3) start[2] <- 3

        }else if(start_point=="flatten"){
          tmp <- create_flatten_cps(run$strat_columns$l, inc_no_measures, n_cps )
          
          start <- c(tmp$cps, tmp$start)
         
         
        }
         equaly_spaced_cps <- FALSE
      equaly_spaced_cps <- FALSE
    }
   partial_cost_function <-  do.call(purrr::partial, c(list(cost), run$strat_columns, list(run_params=run_params, cost_params=cost_params)))
   cost_function <-purrr::partial(partial_cost_function, equaly_spaced_cps=equaly_spaced_cps,n_cps=n_cps, fixed_cps=fixed_cps)
 

    if(start_point=="global"){
      if(global_algo=="direct"){
         sol <- nloptr::direct(function(...) -cost_function(...), lower=lwr, upper=upr, control=list(maxeval=n_global))
         start <- sol$par
      }
    }
   if(local_algo == "optim"){
        op <- optim(start, function(...) -cost_function(...), gr=NULL, lower=lwr, upper=upr, method="Brent", control=list(maxit=n_local, trace=0))
        solution <- op$par
      }else if(local_algo=="GA"){
        GA <- ga(type="real-value", fitness=cost_function, lower=lwr, upper=upr, maxiter = floor(n_local)/25,
                popSize=25, parallel = 1)
        solution = GA@solution
    
     }else if(local_algo=="sbplx"){
     
      solution <- nloptr::sbplx(start, function(...) -cost_function(...), lower=lwr, upper=upr, control=list(maxeval=n_local, trace=FALSE, tol=0.1), nl.info=FALSE)$par
     }else if(local_algo=="bobyqa"){
        solution <- nloptr::bobyqa(start, function(...) -cost_function(...), lower=lwr, upper=upr, control=list(maxeval=n_local, trace=FALSE, tol=0.1), nl.info=FALSE)$par

     }
    solution <- beta_to_daily(c(fixed_cps,solution), run$strat_columns$l, n_cps , equaly_spaced_cps=equaly_spaced_cps)
    
      beta_mean = mean(solution)
      beta_mean_remove_last_2 = mean(solution[1:(run$strat_columns$l-2*30)])
 
      best_cost <- partial_cost_function(solution, n_cps=run$strat_columns$l-1, equaly_spaced_cps=TRUE, ret_all=TRUE)
      cost <- best_cost$cost  %>% mutate(beta_mean=beta_mean, beta_mean_remove_last_2=beta_mean_remove_last_2) %>% mutate(as.data.frame(run$opt_control))
      betas <- data.frame(betas=solution, time=1:length(solution)) %>% mutate(as.data.frame(run$opt_control))
      results <- best_cost$results %>% mutate(as.data.frame(run$opt_control))
  return(list(solution=solution, partial_cost_function=partial_cost_function,cost=cost, betas=betas, results=results))
}



