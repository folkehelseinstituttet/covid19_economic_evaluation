library(dplyr)
library(data.table)
library(metapop)
library(metapopnorge)
source("parameters.r")
source("run_functions.R")
source("QALY.R")
param_file <- "parameter_files/parameters_vaccination.xlsx"
#set# PARAMETERS
RhpcBLASctl::omp_set_num_threads(1)
setDTthreads(1)
library(futile.logger)
flog.threshold(ERROR) 

#global_R_scaling <- 0.9381012


paramsets <- list()
for(R in seq(1.2, 3, by=0.1)){
  for(sev in 1:20){
    for(beta_cost in c(30)){#c(0,5,15,30)){
      paramsets[[length(paramsets) +1]] <- list(R=R,
                                                 sev=sev,
                                                 cap=c(11000, 350),
                                                 import=1,
                                                 vaccination_strat=1,
                                                 scale_inf=1,
                                                 scale_ini=5,
                                                 l=270,
                                                 cost_beta_change=beta_cost,
 						                                      voluntary_mode=1,
                                                 include_tisk=1,
                                                 scale_overcapacity=1,
                                                  beta_fact=1.0)
     }
     for(cap in list(c(11000, 450))){
       paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=cap,
                                                  import=1,
                                                  vaccination_strat=1,
                                                  scale_inf=1,
                                                  scale_ini=5,
                                                  l=270,
                                                  cost_beta_change=30,
  						                                     voluntary_mode=1,
                                                  include_tisk=1,

                                                 scale_overcapacity=1,
                                                  beta_fact=1.0)
      }
  
       for(include_tisk in c(0)){
         paramsets[[length(paramsets) +1]] <- list(R=R,
                                                   sev=sev,
                                                   cap=c(11000, 350),
                                                   import=1,
                                                   vaccination_strat=1,
                                                   scale_inf=1,
                                                   scale_ini=5,
                                                   l=270,
                                                   cost_beta_change=30,
                                                   include_tisk=include_tisk,
                                                   scale_overcapacity=1,
                                                   voluntary_mode=1,
                                                  beta_fact=1.0)

      }
      for(import in c(210)){
        paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=import,
                                                  vaccination_strat=1,
                                                  scale_inf=1,
                                                   scale_ini=5,
                                                   cost_beta_change=30,
                                                  l=270,
                                                  include_tisk=1,
                                                  voluntary_mode=1,

                                                 scale_overcapacity=1,
                                                  beta_fact=1.0)
      }
    for (vs in c(2)){
      paramsets[[length(paramsets) +1]] <- list(R=R,
                                                sev=sev,
                                               cap=c(11000, 350),
                                               import=1,
                                               vaccination_strat=vs,
                                               scale_inf=1,
                                                scale_ini=5,
                                                cost_beta_change=30,
                                               l=270,
                                               include_tisk=1,
                                               voluntary_mode=1,

                                                 scale_overcapacity=1,
                                                  beta_fact=1.0)
   }
      for (si in c(0.5, 1.5)){
        paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=1,
                                                  vaccination_strat=1,
                                                  scale_ini=5,
                                                  cost_beta_change=30,
                                                  l=270,
                                                  scale_inf=si,
                                                  include_tisk=1,
                                                  voluntary_mode=1,

                                                 scale_overcapacity=1,
                                                  beta_fact=1.0) 
       paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=1,
                                                  vaccination_strat=1,
                                                  scale_ini=5,
                                                  cost_beta_change=30,
                                                  l=550,
                                                  scale_inf=si,
                                                  include_tisk=1,
                                                  voluntary_mode=1,

                                                 scale_overcapacity=1,
                                                  beta_fact=1.0)
      }
      for (l in c(550)){
        paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=1,
                                                  vaccination_strat=1,
                                                  cost_beta_change=30,
                                                  scale_inf=1,
                                                  scale_ini=5,
                                                  l=l,
                                                  include_tisk=1,
                                                  voluntary_mode=1,

                                                 scale_overcapacity=1,
                                                  beta_fact=1.0)
     }
    #  for(scale_ini in c(1,10)){
    #     paramsets[[length(paramsets) +1]] <- list(R=R,
    #                                               sev=sev,
    #                                               cap=c(11000, 350),
    #                                               import=1,
    #                                               vaccination_strat=1,
    #                                               scale_inf=1,
    #                                               scale_ini=scale_ini,
    #                                               cost_beta_change=30,
    #                                               l=270,
    #                                               include_tisk=1,
    #                                               voluntary_mode=1,
    #                                               beta_fact=1.0)
    # 						  }
    for(vm in c(0,3,4)){
        paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=1,
                                                  vaccination_strat=1,
                                                  scale_inf=1,
                                                  scale_ini=5,
                                                  l=270,
                                                  cost_beta_change=30,
                                                  include_tisk=1,
                                                  voluntary_mode=vm,

                                                 scale_overcapacity=1,
                                                  beta_fact=1.0)
                                                 }
    
    for(beta_fact in c(0.9)){
        paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=1,
                                                  vaccination_strat=1,
                                                  scale_inf=1,
                                                  scale_ini=5,
                                                  l=270,
                                                  cost_beta_change=30,
                                                  include_tisk=1,
                                                  voluntary_mode=1,

                                                 scale_overcapacity=1,
                                                  beta_fact=beta_fact)
    }

   
     
      paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=1,
                                                  vaccination_strat=1,
                                                  scale_inf=1,
                                                  scale_ini=5,
                                                  l=550,
                                                  cost_beta_change=30,
                                                  include_tisk=1,
                                                  voluntary_mode=1,
                                                  scale_overcapacity=1,
                                                  beta_fact=1, 
                                                  icu_prev_lim=250)

     paramsets[[length(paramsets) +1]] <- list(R=R,
                                                 sev=sev,
                                                 cap=c(11000, 350),
                                                 import=1,
                                                 vaccination_strat=1,
                                                 scale_inf=1,
                                                 scale_ini="R5",
                                                 l=270-14,
                                                 cost_beta_change=30,
                                                 scale_overcapacity=1,
 						                                     voluntary_mode=1,
                                                 include_tisk=1,
                                                beta_fact=1.0)
    
  paramsets[[length(paramsets) +1]] <- list(R=R,
                                                 sev=sev,
                                                 cap=c(11000, 350),
                                                 import=1,
                                                 vaccination_strat=1,
                                                 scale_inf=1,
                                                 scale_ini=5,
                                                 l=270-14,
                                                 cost_beta_change=30,
                                                 scale_overcapacity=1,
 						                                    voluntary_mode=1,
                                                 include_tisk=1,
                                                beta_fact=1.0)

   paramsets[[length(paramsets) +1]] <- list(R=R,
                                                 sev=sev,
                                                 cap=c(11000, 350),
                                                 import=1,
                                                 vaccination_strat=1,
                                                 scale_inf=2,
                                                 scale_ini=5,
                                                 l=3*365,
                                                 cost_beta_change=beta_cost,
 						                                    voluntary_mode=1,
                                                include_tisk=1,
                                                
                                                scale_overcapacity=1,
                                                beta_fact=1.0)
    paramsets[[length(paramsets) +1]] <- list(R=R,
                                                 sev=sev,
                                                 cap=c(11000, 350),
                                                 import=1,
                                                 vaccination_strat=1,
                                                 scale_inf=2,
                                                 scale_ini=5,
                                                 l=3*365,
                                                 cost_beta_change=beta_cost,
 						                                    voluntary_mode=1,
                                                include_tisk=1,
                                                
                                                scale_overcapacity=5,
                                                beta_fact=1.0)
   paramsets[[length(paramsets) +1]] <- list(R=R,
                                                 sev=sev,
                                                 cap=c(11000, 350),
                                                 import=1,
                                                 vaccination_strat=1,
                                                 scale_inf=2,
                                                 scale_ini=5,
                                                 l=3*365,
                                                 cost_beta_change=beta_cost,
 						                                    voluntary_mode=1,
                                                include_tisk=1,
                                                
                                                scale_overcapacity=2,
                                                beta_fact=1.0)
  paramsets[[length(paramsets) +1]] <- list(R=R,
                                                 sev=sev,
                                                 cap=c(11000, 350),
                                                 import=1,
                                                 vaccination_strat=1,
                                                 scale_inf=1,
                                                 scale_ini=5,
                                                 l=3*365,
                                                 cost_beta_change=beta_cost,
 						                        voluntary_mode=1,
                                                include_tisk=1,
                                                
                                                scale_overcapacity=1,
                                                beta_fact=1.0)

  }
}

#paramsets <- list()



for(R in c(2, 2.5)){
  for(sev in c(3, 7, 13)){
    for (si in seq(0, 2, length.out=50)){
        paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=1,
                                                  vaccination_strat=1,
                                                  scale_ini=5,
                                                  cost_beta_change=30,
                                                  l=270,
                                                  scale_inf=si,
                                                  include_tisk=1,
                                                  voluntary_mode=1,
                                                  label="cont_si",
                                                 scale_overcapacity=1,
                                                  beta_fact=1.0) 
    
  }
}
}


for(R in c(2)){
  for(sev in seq(1,18.9, by=0.1)){   
    for(oc in c(1,2,3,4,5)){
        paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=1,
                                                  vaccination_strat=1,
                                                  scale_ini=5,
                                                  cost_beta_change=30,
                                                  l=3*365,
                                                  scale_inf=2,
                                                  include_tisk=1,
                                                  voluntary_mode=1,
                                                  label="test_flatten",
                                                 scale_overcapacity=oc,
                                                  beta_fact=1.0) 
    
  }
}
}

for(R in c(1.05)){
  for(sev in seq(1,3.9, by=0.1)){   

        paramsets[[length(paramsets) +1]] <- list(R=R,
                                                  sev=sev,
                                                  cap=c(11000, 350),
                                                  import=1,
                                                  vaccination_strat=1,
                                                  scale_ini=5,
                                                  cost_beta_change=30,
                                                  l=270,
                                                  scale_inf=1,
                                                  include_tisk=1,
                                                  voluntary_mode=1,
                                                  label="test_postpone",
                                                 scale_overcapacity=1,
                                                  beta_fact=1.0) 
    
  }
}





#x <- paramsets[[7]]

#cost(c(16, 140, 0.75, 0.75, 1),1.4, 1, 0, 1, c(11000, 350), 1,2, deterministic=TRUE, n_parts=1, do_print=TRUE)
N_threads_external <- 115
n_cps <- 4

retry_n <- 2
n <- 50000




opt_control <- list(
  list(n_cps=15, n_local=n, local_algo="bobyqa",start_point="lockdown", n_global=10, fit_cps=FALSE),
  list(n_cps=15, n_local=n, local_algo="bobyqa",start_point="nothing", n_global=10, fit_cps=FALSE),
  list(n_cps=15, n_local=n, local_algo="bobyqa",start_point="global", n_global=n, fit_cps=FALSE),
  list(n_cps=25, n_local=n, local_algo="bobyqa",start_point="global", n_global=n, fit_cps=FALSE),
  list(n_cps=15, n_local=n, local_algo="bobyqa",start_point="cut", n_global=10, fit_cps=FALSE),
  list(n_cps=4, n_local=n, local_algo="sbplx",start_point="cut", n_global=10, fit_cps=TRUE),
  list(n_cps=2, n_local=n, local_algo="sbplx",start_point="flatten", n_global=10, fit_cps=TRUE),
  list(n_cps=2, n_local=n, local_algo="bobyqa",start_point="flatten", n_global=10, fit_cps=TRUE),
  list(n_cps=3, n_local=n, local_algo="sbplx",start_point="cut", n_global=10, fit_cps=TRUE),#,
  list(n_cps=3, n_local=n, local_algo="bobyqa",start_point="cut", n_global=10, fit_cps=TRUE),#,
  list(n_cps=2, n_local=n, local_algo="sbplx",start_point="cut", n_global=10, fit_cps=TRUE)#,
  #list(n_cps=2, n_local=n, local_algo="sbplx",start_point="cut", n_global=10, fit_cps=TRUE)

)

run_id <- function(id){
  print(paste("Running ", id))
for(j in 1:retry_n){
  out <- tryCatchLog::tryCatchLog({
   latest_res <- data.frame()
   x <- paramsets[[as.numeric(id)]]
   
   optim_betas(R=x$R, sev=x$sev, import=x$import, vaccination_strat=x$vaccination_strat, capacity=x$cap,scale_inf=x$scale_inf,l=x$l, run_params=run_params, cost_params=cost_params, opt_control=opt_control, include_tisk=x$include_tisk,voluntary_mode=x$voluntary_mode, scale_ini=x$scale_ini, beta_fact=x$beta_fact, cost_beta_change=x$cost_beta_change,scale_overcapacity=x$scale_overcapacity,icu_prev_lim=x$icu_prev_lim, label=x$label, global_R_scaling=TRUE)},
                         #  optim_betas(x$R, x$sev, x$import, x$vaccination_strat, x$cap,x$scale_inf,x$l, run_params, cost_params,ns=c(0,-5), include_tisk=x$include_tisk, n_parallel=FALSE,voluntary_mode=x$voluntary_mode, max_n=max_n, type=c("sbplx"), scale_ini=x$scale_ini,flatten_the_curve_thresholds = c(1000,2000), cost_beta_change=x$cost_beta_change)},
   error=function(cond){
   print("Error")
	    print(latest_res)
    print(cond)
   print("Error")
	return(NULL)
 }, silent.warnigns=TRUE)
 if(!is.null(out)) break;
}
return(out)

}


N_nodes <- 30
N_per_node <- length(paramsets)/N_nodes
cores_per_node <- 40

id <- Sys.getenv("SLURM_ARRAY_TASK_ID")

if(id=="") id <- 1
id <- as.numeric(id)

ids <- 1:(N_per_node) + (id-1)*N_per_node
print(ids)

results <- parallel::mclapply(
      ids, run_id, mc.cores=cores_per_node
)




saveRDS(results, glue::glue("/projects/ec194/economic_evaluation_results/raw_results/raw_results_long{id}.RDS"))
## results <- lapply(
##                        paramsets[1:10],
##                        function (x) optim_betas(x$R, x$sev, x$import, x$vaccination_strat, x$cap, 0))



#res[, strat:="Supression"]
#res[inf/total > 0.1 & inf/total < 0.9, strat:="Mitigation"]
#res[inf/total > 0.9, strat:="No measures"]




#saveRDS(res, "main_results.RDS")
