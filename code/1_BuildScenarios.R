params_3_VOC <- list()
for(i in seq_len(nrow(model_selected_ur))){
  params_3_VOC[[i]] <- gen_country_basics(country = model_selected_ur$country_name[i],
                                          waning_nat = 52*7*3,
                                          R0_assumed  = model_selected_ur$r[i],
                                          date_start = as.character(ymd("2019-12-01") + model_selected_ur$t[i]),
                                          date_end = "2022-12-31",
                                          processes = burden_processes,
                                          deterministic = TRUE)
  
  add_vac <- function(params){
    res <- params %>%
      update_vac_char(.,
                      ve_i   = ve$ve_i_o[1],  # infection blocking VE post 1 dose
                      v2e_i  = ve$ve_i_o[2],  # infection blocking VE post 2 doses
                      ve_d   = ve$ve_d[1],    # clinical fraction among breakthrough post 1 dose
                      v2e_d  = ve$ve_d[2],    # clinical fraction among breakthrough post 2 doses
                      wv = 1/360) %>%
      update_vac_char(.,
                      ve_i   = ve$ve_i_o[1],  # infection blocking VE post 1 dose
                      v2e_i  = ve$ve_i_o[2],  # infection blocking VE post 2 doses
                      ve_d   = ve$ve_d[1],    # clinical fraction among breakthrough post 1 dose
                      v2e_d  = ve$ve_d[2],    # clinical fraction among breakthrough post 2 doses
                      wv = 1/360) %>% # 1/ waning duration
      change_VOC(.,
                 date_switch = c("2021-04-15", "2021-12-15"),
                 rc_severity = c(1.5, 1.1),
                 rc_transmissibility = c(1.5, 1.5),
                 rc_ve = c(0.5, 0.5))
    return(res);rm(res)
  }
  
  params_3_VOC[[i]] %<>% add_vac(.)

}

add_vp <- function(params){
  res <- vac_policy(params,
                    # these two parameters define the supply conditions
                    milestone_date = c("2021-03-01", # start from 0
                                       "2021-06-30", # 0.03
                                       "2021-12-31", # all population; 0.2
                                       "2022-12-31"), # 0.6
                    milestone_cov = c(0,
                                      0.03,
                                      0.15,
                                      0.3),
                    # prioritisation, assume 60+  all prioritised
                    priority = c(NA, NA, NA, NA,
                                 2,  2,  2,  2,
                                 2,  2,  2,  2,
                                 1,  1,  1,  1),
                    # maximum feasible uptakes
                    cov_max = c(rep(0,2),
                                rep(0.7, 10),
                                rep(0.9, 4)),
                    # the proportion of age groups at which we will compare
                    # rolling out dose 2s or more dose 1s
                    # p_change = 0.7,
                    supply_delay = 24, # unit = weeks
                    dose_interval = 4,
                    switch_date = c("2022-01-15", 
                                    "2022-03-15"))
  return(res);rm(res)
}

params_3_VOC_vp <- list()
for(i in 3:length(params_3_VOC)){
  params_3_VOC_vp[[i]] <- add_vp(params_3_VOC[[i]])
  
  params_3_VOC_vp[[i]]$scenarios %>% 
    lapply(., "[[", "daily_vac_scenarios") %>% 
    map(., function(x){
      x %>% 
        dplyr::select(date, starts_with("Y", ignore.case = F)) %>% 
        mutate_at(vars(starts_with("Y", ignore.case = F)), cumsum) %>% 
        pivot_longer(starts_with("Y", ignore.case = F)) %>% 
        separate(name, into = c("ag", "dose")) %>% 
        mutate(ag = parse_number(ag)) %>% 
        filter(ag > 4) %>% 
        ggplot(., aes(x = date, y = value, group = dose, color = dose)) +
        geom_line() +
        facet_wrap(~ag, scales = "free")
    }) -> p_list
  
  tmp_dir <- paste0("figs/intermediate/allocations/",
                    model_selected_ur$country_name[i])
  
  if(!file.exists(tmp_dir)){
    dir.create(path = tmp_dir)
  }
  
  lapply(1:length(p_list), function(x){
    ggsave(filename = paste0(tmp_dir,"/scenario_",x,".png"),plot = p_list[[x]])
  })
         
  
  print(round(i*100/nrow(model_selected_ur),0))
  
  qs::qsave(params_3_VOC_vp,
       file = "data/intermediate/params_3_vp_18.qs")

  }
