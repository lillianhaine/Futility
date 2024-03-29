###################################################################
###################################################################
###################################################################
## Example: Simulations Code to run 1000 simulations :
## we assume the following scenario: 
## method is the SS-MIX-MEM
## confounding = 0 - so no unobserved confounding
## overlap = 0.50, so we're in a mixed population 
## note the null and alternative data needs to be saved already to be read in
###################################################################
###################################################################
###################################################################

## load libraries: 
library(Rcpp); library(RcppNumerical); library(RcppEigen);  library(MASS); 
library(boot); library(dplyr); library(tidyr); library(parallel)

## set where your Rcode is housed: 
codewd <- "/home/murra484/haine108/Futility/rcode/"

## set where you want to save results: 
savewd <- "/home/murra484/haine108/Futility/Results/"

#### read in all alternate trial data:
alt_list = readRDS(paste0(savewd, "/data/", "overlap_", overlap_deg, "_conf_", conf, "_alt",  ".RDS"))
#### read in all null trial data:
null_list = readRDS(paste0(savewd, "/data/", "overlap_", overlap_deg, "_conf_", conf, "_null",  ".RDS"))

## scenario is either ideal_full_overlap, ideal_50, ideal_no_overlap - varying the mixing parameter between 0, 1, and 0.50
## we adding confounding as well, but don't offer that up in our scenarios defined here
scenario = "ideal_50"
#print(n1)
#int_perc = seq(0.10, 0.75, 0.05)
int_val= 0.30 
conf = 0; overlap_deg = 0.50

## set the number of sims
nreps = 1000; 

## set the number of cores
no_cores = 4

## source relevant functions: 
source(paste0(codewd, "functions.R"))
Rcpp::sourceCpp(paste0(codewd, "MEM_cpp.cpp") )
Rcpp::sourceCpp(paste0(codewd, "SS_MIX_cpp.cpp") )


### get interim data
int_list_alt =   lapply(1:nreps, function(x) get_interim_dat(alt_list[[x]], int_perc = int_val) )
int_list_null =  lapply(1:nreps, function(x) get_interim_dat(null_list[[x]], int_perc = int_val) )

#### run all sims and all methods (TO and SS) using the function specified:
start_time <- Sys.time()
int_to = mclapply(1:nreps, function(x) run_trial(int_val, theta = 0.975,
                                                 nRWD = 1000, nRCT = 400, A = 1,
                                                 nX = 3, res_conf = conf, overlap_deg = overlap_deg,
                                                 seed = x, no_cores,
                                                 data_alt =int_list_alt[[x]],
                                                 data_null = int_list_null[[x]],
                                                 int_flag = T, method = "TO", alt_ss = NULL), mc.cores = no_cores)


end_time <- Sys.time()
end_time-start_time


# to_fn = paste0(savewd, "/both_version/", "overlap_", overlap_deg, "_n", int_val, "_conf_", conf, "_int_to.RDS")
# saveRDS(int_to, to_fn)
# rm(int_to)



start_time <- Sys.time()
final = mclapply(1:nreps, function(x) run_trial(int_val, theta = 0.975,
                                              nRWD = 1000, nRCT = 400, A = 1,
                                              nX = 3, res_conf = conf, overlap_deg = overlap_deg,
                                              seed = x, no_cores,
                                              data_alt = alt_list[[x]],
                                              data_null = null_list[[x]], int_flag = F, method = "TO"), mc.cores = no_cores)


end_time <- Sys.time()
end_time-start_time
 
# saveRDS(final, paste0(savewd, "/both_version/", "overlap_", overlap_deg, "_n", int_val, "_conf_", conf, "_final.RDS"))
# rm(final)

start_time <- Sys.time()
int_ss_alt = mclapply(1:nreps, function(x) run_trial(int_val, theta = 0.975,
                                                 nRWD = 1000, nRCT = 400, A = 1,
                                                 nX = 3, res_conf = conf, overlap_deg = overlap_deg,
                                                 seed = x, no_cores,
                                                 data_alt = int_list_alt[[x]],
                                                 data_null = int_list_null[[x]],
                                                 int_flag = T, method = "SS", alt_ss = T), mc.cores = no_cores)


end_time <- Sys.time()
end_time-start_time

# saveRDS(int_ss_alt, paste0(savewd, "/both_version/", "overlap_", overlap_deg, "_n", int_val, "_conf_", conf, "_ss_alt.RDS"))
# rm(int_ss_alt)

start_time <- Sys.time()
int_ss_null = mclapply(1:nreps, function(x) run_trial(int_val, theta = 0.975,
                                                     nRWD = 1000, nRCT = 400, A = 1,
                                                     nX = 5, res_conf = conf, overlap_deg = overlap_deg,
                                                     seed = x, no_cores,
                                                     data_alt = int_list_alt[[x]],
                                                     data_null = int_list_null[[x]],
                                                     int_flag = T, method = "SS", alt_ss = F), mc.cores = no_cores)


end_time <- Sys.time()
end_time-start_time
 
#saveRDS(int_ss, paste0(savewd, "/both_version/", "overlap_", overlap_deg, "_n", int_val, "_conf_", conf, "_int_ss.RDS"))


## get everything                     
all_trials = list(ss_alt = int_ss_alt,
                  ss_null = int_ss_null,
                  ref = int_to,
                  final = final)
# 
# unlink(to_fn, ss_alt_fn, ss_ref_fn, final_fn)


## save everything: 
saveRDS(all_trials, paste0(savewd, "/both_version/", "overlap_", overlap_deg, "_n", int_val, "_conf_", conf, ".RDS"))
