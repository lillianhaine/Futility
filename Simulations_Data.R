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

conf_vect = seq(-2, 2, 0.5)

if(scenario == "full_overlap"){
  overlap_deg = 1
}
if(scenario == "50"){
  overlap_deg = 0.50
}
if(scenario == "no_overlap"){
  overlap_deg = 0
}

conf = conf_vect[i]

nreps = 2000; 

no_cores = 24
print(no_cores)
#print(int_val)
print(overlap_deg)
print(scenario)
#print(cur_mix)

codewd <- "/home/murra484/haine108/Futility/rcode/"
savewd <- "/home/murra484/haine108/Futility/Results/"

# codewd <- "~/Desktop/Dissertation/Paper3/Futility /"
# savewd <- "~/Desktop/Dissertation/Paper3/Futility /"

source(paste0(codewd, "functions.R"))

### run and save all data: 
save_dat(nRWD = 1000, nRCT = 400, A = 1, 
           nX = 5, res_conf = conf, overlap_deg = overlap_deg,
           nreps = nreps, no_cores = no_cores)
