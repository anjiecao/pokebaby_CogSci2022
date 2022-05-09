# ------------------ tidy up raw simulation results  ------------------ #
tidy_up_raw_simulation_results <- function(df, im_type){
  
  if(im_type == "basic"){
    tidy_df <- df %>% 
      rowwise() %>% 
      mutate(
        log_surprisal = log(surprisal), 
        log_kl = log(kl),
        complexity = case_when(
          grepl("nf_6_of_1",stim_info) ~ "simple", 
          grepl("nf_6_of_3", stim_info) ~ "complex"
        ), 
        sequence_scheme_print = case_when(
          sequence_scheme == "BBBBBB" ~ "No Deviant", 
          sequence_scheme == "BDBBBB" ~ "Deviant at 2nd Trial", 
          sequence_scheme == "BBBDBB" ~ "Deviant at 4th Trial", 
          sequence_scheme == "BBBBBD" ~ "Deviant at Last Trial")
      ) %>% 
      ungroup() %>% 
      select(params_info, complexity, trial_number, sequence_scheme, sequence_scheme_print,surprisal, log_surprisal, kl, log_kl)
    
    
  }else if (im_type == "random"){
  
    tidy_df <- df %>% 
      mutate(
        trial_number = stimulus_idx, 
        log_sample_n = log(sample_n),
        sequence_scheme =  sub(".*_", "", stim_info),
        stimuli_rep = str_replace(stim_info, paste0("_", "ss_", sequence_scheme), ""), 
        complexity = case_when(
          stimuli_rep == "nf_6_of_1" ~ "simple", 
          stimuli_rep == "nf_6_of_3" ~ "complex"
        ), 
        sequence_scheme_print = case_when(
          sequence_scheme == "BBBBBB" ~ "No Deviant", 
          sequence_scheme == "BDBBBB" ~ "Deviant at 2nd Trial", 
          sequence_scheme == "BBBDBB" ~ "Deviant at 4th Trial", 
          sequence_scheme == "BBBBBD" ~ "Deviant at Last Trial")
      ) %>% 
      ungroup() %>% 
      mutate(params_info = NA) %>% 
      select(params_info, complexity, trial_number, sequence_scheme, sequence_scheme_print, sample_n, log_sample_n) %>% 
      mutate(im_type = im_type)
  
  }else{
    tidy_df <- df %>% 
      mutate(
        trial_number = stimulus_idx, 
        log_sample_n = log(sample_n),
        sequence_scheme =  sub(".*_", "", stim_info),
        stimuli_rep = str_replace(stim_info, paste0("_", "ss_", sequence_scheme), ""), 
        complexity = case_when(
          stimuli_rep == "nf_6_of_1" ~ "simple", 
          stimuli_rep == "nf_6_of_3" ~ "complex"
        ), 
        sequence_scheme_print = case_when(
          sequence_scheme == "BBBBBB" ~ "No Deviant", 
          sequence_scheme == "BDBBBB" ~ "Deviant at 2nd Trial", 
          sequence_scheme == "BBBDBB" ~ "Deviant at 4th Trial", 
          sequence_scheme == "BBBBBD" ~ "Deviant at Last Trial")
      ) %>% 
      ungroup() %>% 
      select(params_info, complexity, trial_number, sequence_scheme, sequence_scheme_print, sample_n, log_sample_n) %>% 
      mutate(im_type = im_type)

  }
   return (tidy_df)
}


# ------------------ tidy up a summary results  ------------------ #
summarise_tidy_sim_res <- function(tidy_df){
  tidy_sim_res_sum <- tidy_df %>% 
    group_by(im_type, params_info, complexity, trial_number, sequence_scheme, sequence_scheme_print) %>% 
    summarise(
      trial_sim_res = mean(sample_n), 
      log_trial_sim_res = mean(log_sample_n)
    ) %>% 
    filter(!is.na(complexity)) %>% 
    select(im_type, params_info, complexity, trial_number, sequence_scheme, sequence_scheme_print, trial_sim_res, log_trial_sim_res)
  
  return(tidy_sim_res_sum)
}

# ------------------ tidy up a  scaled summary results  ------------------ #

summarise_scaled_tidy_sim_res <- function(scaled_tidy_df){
  
  scaled_tidy_sim_res_sum <- scaled_tidy_df %>% 
    group_by(im_type, params_info, complexity, trial_number, sequence_scheme, sequence_scheme_print) %>% 
    summarise(
      mean_scaled_log_sample = mean(scaled_log_sample_n)
    ) %>% 
    filter(!is.na(complexity)) %>% 
    select(im_type, params_info, complexity, trial_number, sequence_scheme, sequence_scheme_print, mean_scaled_log_sample)
  return(scaled_tidy_sim_res_sum)

  }

# ------------------ scaling the model results to match mean and sd of the behavioral results according to the best fit parameters   ------------------  #

scale_model_results <- function(behavioral_df, sim_df, best_params){
  
  if(is.na(best_params)){
    scaled_sim_df <- behavioral_df %>% 
      left_join(sim_df, by = c("trial_number", "complexity", "sequence_scheme", "sequence_scheme_print")) %>% 
      group_by(params_info) %>% 
      summarise(
        mean_log_trial_sim_res = mean(log_sample_n),  # m1
        sd_log_trial_sim_res = sd(log_sample_n), # s1
        mean_log_trial_looking_time = mean(log(trial_looking_time)), # m2
        sd_log_trial_looking_time = sd(log(trial_looking_time)) # s2
      ) %>% 
      left_join(sim_df, by = "params_info") %>% 
      mutate(
        multiply_const = sd_log_trial_looking_time / sd_log_trial_sim_res # s2/s1
      ) %>% 
      mutate(
        scaled_log_sample_n = mean_log_trial_looking_time + (log_sample_n - mean_log_trial_sim_res) * multiply_const
      ) %>% 
      filter(!is.na(complexity)) %>% 
      select(im_type, params_info, complexity, trial_number, sequence_scheme, sequence_scheme_print, sample_n, log_sample_n, scaled_log_sample_n)
    
    
  }else{
    
    scaled_sim_df <- behavioral_df %>% 
      left_join(sim_df %>% filter(params_info == best_params), by = c("trial_number", "complexity", "sequence_scheme", "sequence_scheme_print")) %>% 
      group_by(params_info) %>% 
      summarise(
        mean_log_trial_sim_res = mean(log_sample_n),  # m1
        sd_log_trial_sim_res = sd(log_sample_n), # s1
        mean_log_trial_looking_time = mean(log(trial_looking_time)), # m2
        sd_log_trial_looking_time = sd(log(trial_looking_time)) # s2
      ) %>% 
      left_join(sim_df, by = "params_info") %>% 
      mutate(
        multiply_const = sd_log_trial_looking_time / sd_log_trial_sim_res # s2/s1
      ) %>% 
      mutate(
        scaled_log_sample_n = mean_log_trial_looking_time + (log_sample_n - mean_log_trial_sim_res) * multiply_const
      ) %>% 
      filter(!is.na(complexity)) %>% 
      select(im_type, params_info, complexity, trial_number, sequence_scheme, sequence_scheme_print, sample_n, log_sample_n, scaled_log_sample_n)
    
    
  }
  
  
  return(scaled_sim_df)
}

# a slightly different version of scaling for the basic model; 
# since it doesn't have trial level variation, scaling based on all the mean and sd of the summary as well 
scale_basic_model <- function(behavioral_df, sim_df, best_params){
  scaled_sim_df <- sim_df %>% 
    filter(params_info == best_params) %>% 
    left_join(behavioral_df, by = c("trial_number", "complexity", "sequence_scheme", "sequence_scheme_print")) %>% 
    group_by(params_info) %>% 
    summarise(
      
      mean_log_surprisal = mean(log_surprisal), 
      sd_log_surprisal = sd(log_surprisal), 
      mean_log_kl = mean(log_kl), 
      sd_log_kl = sd(log_kl),
      
      mean_log_trial_looking_time = mean(log(trial_looking_time)), 
      sd_log_trial_looking_time = sd(log(trial_looking_time))
    ) %>% 
      left_join(sim_df, by = "params_info") %>% 
    mutate(
      scaled_log_surprisal = mean_log_trial_looking_time + (log_surprisal - mean_log_surprisal) * (sd_log_trial_looking_time / sd_log_surprisal), 
      scaled_log_kl =  mean_log_trial_looking_time + (log_kl - mean_log_kl) * (sd_log_trial_looking_time / sd_log_kl)
    ) %>% 
    left_join(behavioral_df,  by = c("trial_number", "complexity", "sequence_scheme", "sequence_scheme_print"))
  
  return(scaled_sim_df)
}



# ------------------ calculate correlation between behavioral and simulation results  ------------------  #
calculate_correlation <- function(sim_b_df){
  d <- sim_b_df %>% 
    ungroup() %>% 
    nest_by(params_info) 
  
  d$r <- unlist(map(d$data, function(x){
    r <- cor(x$trial_sim_res, x$trial_looking_time, method = "pearson")
  }))
  d$r_in_log <- unlist(map(d$data, function(x){
    r_in_log <- cor(x$log_trial_sim_res, x$log_trial_looking_time, method = "pearson")
  }))
  
  d$rmse <- unlist(map(d$data, function(x){
    rmse_log <- Metrics::rmse(x$trial_looking_time, x$trial_sim_res)
  }))
  
  d$rmse_in_log <- unlist(map(d$data, function(x){
    rmse_log <- Metrics::rmse(x$log_trial_looking_time, x$log_trial_sim_res)
  }))
  
  d %>% 
    arrange(-r_in_log)
}

# ------------------ calculate correlation between behavioral and simulation results after scaling for the best parameter ------------------  #

calculate_scaled_correlation <- function(sim_b_df){
  
  rmse_interval <- function(rmse, deg_free, p_lower = 0.025, p_upper = 0.975){
    tibble(.pred_lower = sqrt(deg_free / qchisq(p_upper, df = deg_free)) * rmse,
           .pred_upper = sqrt(deg_free / qchisq(p_lower, df = deg_free)) * rmse)
  }
  
  
  d <- sim_b_df %>% 
    ungroup() %>% 
    nest_by(params_info) 
  
  
  d$r_in_log <- unlist(map(d$data, function(x){
    r_in_log <- cor(x$mean_scaled_log_sample, x$log_trial_looking_time, method = "pearson")
  }))
  
  d$r_in_log_ci <- unlist(map(d$data, function(x){
    r_conf <- confintr::ci_cor(x$mean_scaled_log_sample, x$log_trial_looking_time, method = "pearson", type = "bootstrap")["interval"][[1]]
    r_conf_print <- paste0("[", round(r_conf[1],2), ", ", round(r_conf[2],2), "]")
    return(r_conf_print)
  }))
 
  d$rmse_in_log <- unlist(map(d$data, function(x){
    rmse_log <- Metrics::rmse(x$log_trial_looking_time, x$mean_scaled_log_sample)
  }))
  
  d$rmse_in_log_ci <- unlist(map(d$data, function(x){
    rmse_log <- Metrics::rmse(x$log_trial_looking_time, x$mean_scaled_log_sample)
    rmse_log_interval <- rmse_interval(rmse_log, nrow(x))
    rmse_conf_print <- paste0("[", round(rmse_log_interval$.pred_lower,2), ", ", 
                              round(rmse_log_interval$.pred_upper,2), "]")
    return (rmse_conf_print)
  }))
  
  d %>% 
    arrange(-r_in_log)
}


calculate_scaled_correlation_for_sa <- function(sim_b_df){
  d <- sim_b_df %>% 
    ungroup() %>% 
    nest_by(params_info) 
  
  
  d$r_in_log <- unlist(map(d$data, function(x){
    r_in_log <- cor(x$log_trial_sim_res, x$log_trial_looking_time, method = "pearson")
  }))
  
  
  d$rmse_in_log <- unlist(map(d$data, function(x){
    rmse_log <- Metrics::rmse(x$log_trial_looking_time, x$log_trial_sim_res)
  }))
  
  d %>% 
    arrange(-r_in_log)
}


# ------------------ calculate correlation for the basic results  ------------------  #
calculate_correlation_for_basic_model <- function(basic_sim_d_df, best_for){
  d <- basic_sim_d_df %>% 
    nest_by(params_info)
  
  
  d$kl_r <- unlist(map(d$data, function(x){
    r <- cor(x$trial_looking_time, x$kl, method = "pearson")
  }))
  
  d$kl_r_log <- unlist(map(d$data, function(x){
    r <- cor(x$log_trial_looking_time, log(x$kl), method = "pearson")
  }))
  
  
  d$surprisal_r <- unlist(map(d$data, function(x){
    r <- cor(x$trial_looking_time, x$surprisal, method = "pearson")
  }))
  
  d$surprisal_r_log <- unlist(map(d$data, function(x){
    r <- cor(x$log_trial_looking_time, log(x$surprisal), method = "pearson")
  }))
  
  if(best_for == "KL"){
    return (
      d %>% 
        arrange(-kl_r_log)
    )
  }else if(best_for == "surprisal"){
    return (
      d %>% 
        arrange(-surprisal_r_log)
    )
  }
  
}

calculate_correlation_for_scaled_basic_model <- function(basic_sim_d_df){
  d <- basic_sim_d_df %>% 
    nest_by(params_info)
  
  
 
  d$kl_r_log <- unlist(map(d$data, function(x){
    r <- cor(x$log_trial_looking_time, log(x$scaled_log_kl), method = "pearson")
  }))
  
  
  d$surprisal_r_log <- unlist(map(d$data, function(x){
    r <- cor(x$log_trial_looking_time, x$scaled_log_surprisal, method = "pearson")
  }))
  
  d$kl_rmse_log <- unlist(map(d$data, function(x){
    rmse_log <- Metrics::rmse(x$log_trial_looking_time, x$scaled_log_kl)
  }))
  
  d$surprisal_rmse_log <- unlist(map(d$data, function(x){
    rmse_log <- Metrics::rmse(x$log_trial_looking_time, x$scaled_log_surprisal)
  }))
  
  
  
  return (d)
  
}


