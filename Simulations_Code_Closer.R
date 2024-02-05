###################################################################
###################################################################
###################################################################
## Simulations Code: 
###################################################################
###################################################################
###################################################################
.libPaths("/panfs/roc/groups/12/ssafo/haine108/R/x86_64-pc-linux-gnu-library/4.1")
.libPaths("/panfs/jay/groups/17/ssafo/haine108/R/x86_64-pc-linux-gnu-library/4.1")
library(Rcpp); library(RcppNumerical); library(RcppEigen);  library(MASS); 
library(boot); library(dplyr); library(tidyr); library(parallel)

i = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

## get arguments and values under study: 
args=(commandArgs(T))

print(args[1])
print(args[2])

method = gsub("\\method=", "", args[1])
scenario = gsub("scenario=", "", args[2], fixed = T)
#print(n1)
int_perc = seq(0.2, 0.55, 0.025)
int_val= int_perc[i]

if(scenario == "ideal_full_overlap"){
  conf = 0; overlap_deg = 1
}
if(scenario == "ideal_50"){
  conf = 0; overlap_deg = 0.50

  
}
if(scenario == "ideal_no_overlap"){
  conf = 0; overlap_deg = 0

  
}


nreps = 5000; 

no_cores = 100
print(no_cores)
print(int_val)
print(overlap_deg)
print(scenario)
#print(cur_mix)

codewd <- "/home/murra484/haine108/Futility/rcode/"
savewd <- "/home/murra484/haine108/Futility/Results/"

# codewd <- "~/Desktop/Dissertation/Paper3/Futility /"
# savewd <- "~/Desktop/Dissertation/Paper3/Futility /"

source(paste0(codewd, "functions.R"))
Rcpp::sourceCpp(paste0(codewd, "MEM_cpp.cpp") )
Rcpp::sourceCpp(paste0(codewd, "SS_MIX_cpp.cpp") )


### run all reps: 

#### read in all alternate trial data:
alt_list = readRDS(paste0(savewd, "/data/", "overlap_", overlap_deg, "_conf_", conf, "_alt",  ".RDS"))
#### read in all null trial data:
null_list = readRDS(paste0(savewd, "/data/", "overlap_", overlap_deg, "_conf_", conf, "_null",  ".RDS"))

int_list_alt =   lapply(1:nreps, function(x) get_interim_dat(alt_list[[x]], int_perc = int_val) )
int_list_null =  lapply(1:nreps, function(x) get_interim_dat(null_list[[x]], int_perc = int_val) )

#### run all sims:
start_time <- Sys.time()
int_to = mclapply(1:nreps, function(x) run_trial(int_val, theta = 0.975,
                                                 nRWD = 1000, nRCT = 400, A = 1,
                                                 nX = 5, res_conf = conf, overlap_deg = overlap_deg,
                                                 seed = x, no_cores,
                                                 data_alt =int_list_alt[[x]],
                                                 data_null = int_list_null[[x]],
                                                 int_flag = T, method = "TO", alt_ss = NULL), mc.cores = no_cores)


end_time <- Sys.time()
end_time-start_time
# 
# saveRDS(int_to, paste0(savewd, "/both_version/", "overlap_", overlap_deg, "_n", int_val, "_conf_", conf, "_int_to.RDS"))
# 
start_time <- Sys.time()
final = mclapply(1:nreps, function(x) run_trial(int_val, theta = 0.975,
                                                nRWD = 1000, nRCT = 400, A = 1,
                                                nX = 5, res_conf = conf, overlap_deg = overlap_deg,
                                                seed = x, no_cores,
                                                data_alt = alt_list[[x]],
                                                data_null = null_list[[x]], int_flag = F, method = "TO"), mc.cores = no_cores)


end_time <- Sys.time()
end_time-start_time
# 
# saveRDS(final, paste0(savewd, "/both_version/", "overlap_", overlap_deg, "_n", int_val, "_conf_", conf, "_final.RDS"))
# 
start_time <- Sys.time()
int_ss_alt = mclapply(1:nreps, function(x) run_trial(int_val, theta = 0.975,
                                                     nRWD = 1000, nRCT = 400, A = 1,
                                                     nX = 5, res_conf = conf, overlap_deg = overlap_deg,
                                                     seed = x, no_cores-25,
                                                     data_alt = int_list_alt[[x]],
                                                     data_null = int_list_null[[x]],
                                                     int_flag = T, method = "SS", alt_ss = T), mc.cores = 10)


end_time <- Sys.time()
end_time-start_time

start_time <- Sys.time()
int_ss_null = mclapply(1:nreps, function(x) run_trial(int_val, theta = 0.975,
                                                      nRWD = 1000, nRCT = 400, A = 1,
                                                      nX = 5, res_conf = conf, overlap_deg = overlap_deg,
                                                      seed = x, no_cores-25,
                                                      data_alt = int_list_alt[[x]],
                                                      data_null = int_list_null[[x]],
                                                      int_flag = T, method = "SS", alt_ss = F), mc.cores = 10)


end_time <- Sys.time()
end_time-start_time
# 
# saveRDS(int_ss, paste0(savewd, "/both_version/", "overlap_", overlap_deg, "_n", int_val, "_conf_", conf, "_int_ss.RDS"))


all_trials = list(ss_alt = int_ss_alt,
                  ss_null = int_ss_null,
                  ref = int_to, 
                  final = final)



saveRDS(all_trials, paste0(savewd, "/both_version/", "overlap_", overlap_deg, "_n", int_val, "_conf_", conf, "_closer", ".RDS"))
