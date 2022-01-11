# ------------------ tidy up raw simulation results  ------------------ #
tidy_up_raw_simulation_results <- function(df, im_type){
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



# ------------------ scaling the model results to match mean and sd of the behavioral results according to the best fit parameters   ------------------  #
scale_model_results <- function(behavioral_df, sim_df, best_params){
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
  
  return(scaled_sim_df)
}




# ------------------ calculate correlation between behavioral and simulation results  ------------------  #
calculate_correlation <- function(sim_b_df){
  d <- sim_b_df %>% 
    ungroup() %>% 
    nest_by(params_info) 
  
  d$perason_corr <- unlist(map(d$data, function(x){
    pearson_r <- cor(x$trial_sim_res, x$trial_looking_time, method = "pearson")
  }))
  d$perason_corr_in_log <- unlist(map(d$data, function(x){
    pearson_r_in_log <- cor(x$log_trial_sim_res, x$log_trial_looking_time, method = "pearson")
  }))
  
  
  d %>% 
    arrange(-perason_corr_in_log)
}