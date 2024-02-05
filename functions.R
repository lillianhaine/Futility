########################################################## 
########################################################## 
## R code to run ss-mix and rct only
########################################################## 
########################################################## 

########################################################## 
########################################################## 
## posterior predictive function: 

prob_pt_lt_pc = function(ac, bc, at, bt){
  f = function(p) {pbeta(p, at, bt)*(dbeta(p, ac, bc) )}
  
  #print(k)
  
  ret.val =  integrate(f, 0, 1)$value
  
  #print(ret.val)
  
  return( ret.val)
  
}

########################################################## 
########################################################## 
## rct only: 
draws_reference = function(dat, nMCMC, N, theta, no_cores){
  #### filter data to be only RCT: 
  rct_dat = dplyr::filter(dat, S == 1)
  ctrl = dplyr::filter(rct_dat, arm == 0); trt = dplyr::filter(rct_dat, arm == 1)
  print(nrow(rct_dat))
  ### summary statistics: 
  rc = sum(ctrl$Y); nc = nrow(ctrl)
  rt = sum(trt$Y); nt = nrow(trt)
  
  #if(N > nc){print("INTERIM")}
  #if(N == nc){print("FINAL")}
  
  #### get draws from pc and pt posteriors: 
  ### pc posterior: ac = rc+1, bc = 1+nc-rc
  ac_post = 1+rc; bc_post = 1+nc-rc
  pc = rbeta(nMCMC, ac_post, bc_post)
  
  ### pt posterior: ac = rc+1, bc = 1+nc-rc
  at_post = 1+rt; bt_post = 1+nt-rt
  pt = rbeta(nMCMC, at_post, bt_post)
  
  #### diff: 
  diff = pc-pt
  #### for each pc, pt draw get response draws: 
  rc_star = sapply(1:nMCMC, function(x) rbinom(1, N-nc, prob = pc[x]) )
  rt_star = sapply(1:nMCMC, function(x) rbinom(1, N-nt, prob = pt[x]) )
  
          
  #### get Pr(pt < pc | D, D*) 
  ac_vect = 1+rc_star+rc; at_vect = 1+rt_star+rt
  bc_vect = 1+N-(rc_star+rc); bt_vect = 1+N-(rt_star+rt)
  pr_draws = sapply(1:nMCMC, function(x) 
    prob_pt_lt_pc(ac_vect[x], bc_vect[x], at_vect[x], bt_vect[x]) )
  
  #### get Pr( Pr(pt < pc | D, D*) > theta | D)
  PPP = sum(pr_draws > theta)/length(pr_draws)
  if(N == nc){
    PPP = ifelse(prob_pt_lt_pc(ac_post, bc_post, at_post, bt_post) > theta, 1, 0)
  }

  ret.list = list(pc = list(mean = mean(pc), sd = sd(pc)), 
                  pt = list(mean = mean(pt), sd = sd(pt)), 
                  diff = list(mean = mean(diff), sd = sd(diff)), 
                  pr_draws = list(mean = mean(pr_draws), sd = sd(pr_draws) ), 
                  PPP = PPP)
  ret.val = PPP
  return(ret.val)
}

########################################################## 
########################################################## 
## ss-mix-mem: 
draws_ss_mix_mem = function(dat, nMCMC1, nMCMC2, 
                            nburnin1, nburnin2, N, theta, no_cores){
  dat_ps = mutate(dat, logitPS = logit(P_RCT))
  tot_MCMC1 = nMCMC1+nburnin1
  trial = dplyr::filter(dat_ps, S == 1) ; rct_logitps = trial$logitPS; 
  #print(nrow(trial))
  rwd = dplyr::filter(dat_ps, S == 0); rwd_logitps = rwd$logitPS; yrwd = rwd$Y
  
  ### step 1: run the ss-mix step: 
  ss_mix_draws = SS_MIX(ps_trial = rct_logitps, ps_ext = rwd_logitps, Mall = tot_MCMC1,  
                        y_ext = yrwd,
                        sd_muct = 0.00001, a_tauct = .001, b_tauct = .001, 
                        sd_muu =  0.00001, a_tauu = .001,  b_tauu = .001, 
                        a_w = 1, b_w = 1)
  
  ### set up sufficient statistics for MEM: 
  r_ext = ss_mix_draws$yobsIN[-c(1:nburnin1)]; n_ext = ss_mix_draws$nobsIN[-c(1:nburnin1)]
  rm(ss_mix_draws)
  ctrl = filter(trial, arm == 0); trt = filter(trial, arm == 1)
  rc = sum(ctrl$Y); nc = nrow(ctrl); rt = sum(trt$Y); nt = nrow(trt)
  ### step 2: run the mem step: 
  nstar = N-nc
  mem_draws = mclapply(1:nMCMC1, function(x) MEM_binary(rc, nc, rt, nt, 
                                                      r_ext[x], n_ext[x], 0, 0, 
                                                      nburnin2, nMCMC2, 0.50, 0, 1, 1, 1, 1, 1, 1, 1, 1, nstar), 
                       mc.cores = no_cores-5)
  
  #pc = do.call(rbind, lapply(1:nMCMC1, function(x) mem_draws[[x]]$pc[-c(1:nburnin2)] ) )
  #pt = do.call(rbind, lapply(1:nMCMC1, function(x) mem_draws[[x]]$pt[-c(1:nburnin2)] ) )
  #rc_star = do.call(rbind, lapply(1:nMCMC1, function(x) mem_draws[[x]]$rcstar[-c(1:nburnin2)] ) )
  #rt_star = do.call(rbind, lapply(1:nMCMC1, function(x) mem_draws[[x]]$rtstar[-c(1:nburnin2)] ) )
  
  
  #totMCMC = nMCMC1*nMCMC2
  #rc_star = sapply(1:totMCMC, function(x) rbinom(1, N-nc, prob = pc[x]) )
  #rt_star = sapply(1:totMCMC, function(x) rbinom(1, N-nt, prob = pt[x]) )
  
  
  #### get Pr(pt < pc | D, D*) 
  #ac_vect = 1+rc_star+rc; at_vect = 1+rt_star+rt
  #bc_vect = 1+N-(rc_star+rc); bt_vect = 1+N-(rt_star+rt)
  #pr_draws_noC = do.call(rbind, lapply(1:totMCMC, function(x) prob_pt_lt_pc(ac_vect[x], bc_vect[x], at_vect[x], bt_vect[x])) )
  
  pr_draws = do.call(rbind, lapply(1:nMCMC1, function(x) mem_draws[[x]]$PPP[-c(1:nburnin2)] ) )
  rm(mem_draws)
  #pc_draws = do.call(rbind, lapply(1:nMCMC1, function(x) mem_draws[[x]]$pc[-c(1:nburnin2)] ) )
  #pt_draws = do.call(rbind, lapply(1:nMCMC1, function(x) mem_draws[[x]]$pt[-c(1:nburnin2)] ) )

  #### get Pr( Pr(pt < pc | D, D*) > theta | D)
  #PPP = sum(pr_draws > theta, 1, 0)/totMCMC
  
  PPP = sum(pr_draws > theta)/length(pr_draws)
  
  
  # ret.val = list(PPP = PPP, 
  #                pc_vals = list(mean = mean(pc_draws), 
  #                               sd = sd(pc_draws)), 
  #                pr_vals = list(mean = mean(pr_draws), 
  #                               sd = sd(pr_draws)), 
  #                pt_vals = list(mean = mean(pt_draws), 
  #                               sd = sd(pt_draws)), 
  #                n_in = list(mean = mean(n_ext), 
  #                            sd = sd(n_ext))
  # )
  return(PPP)
  
  
}


########################################################## 
########################################################## 
## data generation
data_generation = function(nRCT, nRWD, A, nX = 5,
                           res_conf, cur_beta, mixing_param,
                           b0, tau, b0_pt, seed){
  #print(seed)
  ### from previous run for getting b0 and tau:
  set.seed(seed)
  tot_RCT = (A+1)*nRCT
  nTOT = tot_RCT+nRWD
  arms_inds = seq(0, A)
  
  logit = function(p){
    ret.val = log(p/(1-p))
    return(ret.val)
  }
  
  inv_logit = function(p){
    ret.val = exp(p)/(1+exp(p))
    return(ret.val)
  }
  
  ###### generative mixture components: 
  # components <- sample(1:2,prob=c(mixing_param,
  #                                 1-mixing_param),size = 10000,replace=TRUE)
  # ext_dat_mean = c(0, -1.5)
  #components_all <- c( rep(1, 10000), components) 
  # X_all <- do.call(rbind, lapply(components_all,
  #                                function(x){
  #                                  mvrnorm(n=1,
  #                                          mu=rep(ext_dat_mean[x], nX),
  #                                          Sigma = diag(1, nX, nX)) } ) )
  X_all <- mvrnorm(n = 100000, 
                   mu = rep(0, nX), 
                   Sigma = diag(1, nX, nX) )
  # print(mixing_param)
  # print( sum(components == 1)/length(components) )
  # print( sum(components_all == 1)/length(components_all) )
  # X_ext = data.frame(X_ext)
  # X_ext$trueZ = components
  # 
  # ##### generate mvn X matrices for rct participants
  # X_rct = mvrnorm(n = 10000, 
  #                 mu = rep(0, nX), 
  #                 Sigma = diag(10, nX, nX))
  # X_rct = data.frame(X_rct)
  # X_rct$trueZ = 1
  
  ##### now put them all together: 
  # X_all <- rbind(X_ext, X_rct)
  
  #print(dim(X_all))
  
  dat_X = data.frame(X_all)
  ##### Pr(Trial | X) (get PS for trial membership)
  dat_ps = mutate(dat_X, P_RCT = inv.logit(cur_beta*X1+cur_beta*X2+cur_beta*X3) )
  dat_ps$Z = rbinom(nrow(dat_ps), 1, prob = dat_ps$P_RCT)
  #### sample RCT based on Pr(Trial | X)
  rct_dat = slice_sample(filter(dat_ps, Z == 1), n = tot_RCT)
  rwd_dat = slice_sample(filter(dat_ps, Z == 0), n = nRWD)%>% 
    mutate(S = 0, arm = 0)
  
  rct_dat = mutate(rct_dat, S = 1) 
  rct_dat$arm = as.vector(sapply(1:(nrow(rct_dat)/(A+1)), function(p) sample(c(arms_inds)) ) )
  
  #### Pr(Y = 1 | X) (generate event probabilities for outcomes) and Y variable
  dat = data.frame(rbind(rct_dat, rwd_dat)) %>% 
    mutate(y_prob = inv.logit(b0+tau*arm+(1-S)*res_conf+X1+X2+X3))
  dat$Y = rbinom(nTOT, 1, prob = dat$y_prob)
  
  
  ret.list = list(dat = dat, 
                  ctrl = mean(filter(dat, S == 1, arm == 0)$Y), 
                  trt = mean(filter(dat, S == 1, arm == 1)$Y))
  
  return(dat)
}

########################################################## 
########################################################## 
## interim data:
### get the specified % of total RCT for interim analysis: 
get_interim_dat = function(all_dat, int_perc){
  
  rwd = dplyr::filter(all_dat, S == 0)
  ### filter ctrl and treated: 
  ctrl = dplyr::filter(all_dat, S == 1, arm == 0)
  trt = dplyr::filter(all_dat, S == 1, arm == 1)
  
  ### now get the interim data: 
  n_int = ceiling((int_perc/2)*nrow(ctrl))
  indices = sample(1:nrow(ctrl), size = n_int, replace = F)
  #print(paste0("N at Interim:", n_int) )
  int_ctrl = ctrl[indices, ]
  int_trt = trt[indices, ]
  
  ### all together: 
  int_dat = data.frame(rbind(int_ctrl, int_trt, rwd))
  return(int_dat)
}

save_dat = function(nRCT, nRWD, A, nX = 5, 
                    res_conf, overlap_deg, nreps, no_cores){
    ### step 1: generate trial data: 
    if(overlap_deg == 1){
      b0 = -1.265625; tau = -0.790625; cur_beta = 0; b0_pt = 0; mixing_param = 1
    }
    if(overlap_deg == 0.50){
      b0 = -2.03125; tau =  -0.728125; cur_beta = 0.75; b0_pt = 0; mixing_param = 0.50
    }
    if(overlap_deg == 0){
      b0 = -2.296875; tau =  -0.68125; cur_beta = 2; b0_pt = 0; mixing_param = 0
    }
    
    all_data_alt = mclapply(1:nreps, function(x) data_generation(nRCT = nRCT, nRWD = nRWD, A = A, nX = 5, 
                                                                 res_conf = res_conf, cur_beta = cur_beta, mixing_param = mixing_param,
                                   b0 = b0, tau = tau, b0_pt = b0_pt, x), mc.cores = no_cores )
    all_data_null = mclapply(1:nreps, function(x) data_generation(nRCT = nRCT, nRWD = nRWD, A = A, nX = 5, 
                                                                  res_conf = res_conf, cur_beta = cur_beta, mixing_param = mixing_param,
                                    b0 = b0, tau = 0, b0_pt = b0_pt, x), mc.cores = no_cores )
    
    saveRDS(all_data_alt, paste0(savewd, "/data/", "overlap_", overlap_deg, "_conf_", conf, "_alt",  ".RDS"))
    
    saveRDS(all_data_null, paste0(savewd, "/data/", "overlap_", overlap_deg, "_conf_", conf, "_null",  ".RDS"))
    
    
  
}

########################################################## 
########################################################## 
##  run a trial:
run_trial = function(int_perc, theta, nRCT, nRWD, A, nX = 5, 
                     res_conf, overlap_deg, seed, num_cores, 
                     data_alt, data_null, int_flag, method, alt_ss = T){
  set.seed(seed)
  print(paste0("start iter:", seed))
  ### step 2: get the interim data for treated, control: 
  #print(paste0("N data RCT interim:", nrow(int_dat_alt)-nRWD) )
  #print(paste0("percent at interim: ", int_perc*nRCT))
  
   if(int_flag == T){
     if(method == "TO"){
  
      int_alt_ref = draws_reference(data_alt,10000, nRCT, theta )

      int_null_ref = draws_reference(data_null,10000, nRCT, theta )
      
      ret.list = list(alt = int_alt_ref, 
                      null = int_null_ref)
      
      print(paste0("done with ref"))


      }
    if(method == "SS"){
      if(alt_ss == T){
        int_alt_ss = draws_ss_mix_mem(data_alt, nMCMC1 = 200, nMCMC2 = 200,
                                                           nburnin1 = 1000, nburnin2 = 1000, N = nRCT, theta, num_cores ) 
        ret.list = int_alt_ss
        print(paste0("done with ss-mix-mem alt"))
        }

      if(alt_ss == F){
        int_null_ss = draws_ss_mix_mem(data_null, nMCMC1 = 200, nMCMC2 = 200,
                                       nburnin1 = 1000, nburnin2 = 1000, N = nRCT, theta, num_cores )
        ret.list = int_null_ss
        print(paste0("done with ss-mix-mem null"))
      }

      

      }

  
  }
  # rm(all_data_null)
  # rm(all_data_alt)
  # 
  # ret.list$ref_int = list(alt = int_alt_ref, null = int_null_ref)
  # ret.list$final = list(alt = final_alt, null = final_null)
  # rm(int_alt_ref)
  # rm(int_null_ref)
  if(int_flag == F){
    final_alt = draws_reference(data_alt, 10000, nRCT, theta)
    final_null = draws_reference(data_null, 10000, nRCT, theta)
    ret.list = list(alt = final_alt, 
                 null = final_null)
  }
  print(paste0("done with iter:", seed))
  
  
  #}
  # if(method == "reference"){
  #   final_alt = draws_reference(all_data_alt, 10000, nRCT, theta)
  #   
  #   final_null = draws_reference(all_data_null, 10000, nRCT, theta)
  #   
  #   int_alt = draws_reference(int_dat_alt,10000, nRCT, theta )
  #   int_null = draws_reference(int_dat_null,10000, nRCT, theta )
  #   
  #   ret.list = list(int = list( alt = int_alt, 
  #                               null = int_null), 
  #                   final = list(alt = final_alt, 
  #                                null = final_null) )
  #   print(paste0("done with ref"))
  #   
  # }
  # if(method == "ss_mix_mem"){
  #   final_alt = draws_reference(all_data_alt, 10000, nRCT, theta)
  #   
  #   final_null = draws_reference(all_data_null, 10000, nRCT, theta)
  #   int_alt = draws_ss_mix_mem(int_dat_alt,nMCMC1 = 1000, nMCMC2 = 1000,
  #                              nburnin1 = 1000, nburnin2 = 1000, N = nRCT, theta, no_cores )
  #   int_null = draws_ss_mix_mem(int_dat_null,nMCMC1 = 1000, nMCMC2 = 1000,
  #                               nburnin1 = 1000, nburnin2 = 1000, N = nRCT, theta, no_cores )
  #   ret.list = list(int = list( alt = int_alt, 
  #                               null = int_null), 
  #                   final = list(alt = final_alt, 
  #                                null = final_null) )
  #   print(paste0("done with ss-mix-mem"))
  # }
  return(ret.list)
}
