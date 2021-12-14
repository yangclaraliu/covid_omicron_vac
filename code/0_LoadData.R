pacman::p_load(
  tidyverse, sf, countrycode, rnaturalearth, magrittr, data.table,
  ggsflabel, mgcv, pspline, viridis, ggsci, mgcv, imputeTS, ggpattern,
  ggpubr, gridExtra, grid,   tidyverse, sf, countrycode, rnaturalearth, 
  magrittr, data.table, ggsflabel, mgcv, pspline, viridis, ggsci, mgcv, 
  imputeTS, cowplot, qs, testthat
)

path_euro <- "C:/Users/eideyliu/Documents/GitHub/COVID_Vac_Delay/"
path_dropbox <- "C:/Users/eideyliu/Dropbox/Github_Data/COVID-Vac_Delay/"

cm_path <- paste0(path_euro, "code/covidm_for_fitting/")
cm_force_rebuild <- F
cm_build_verbose <- T
cm_version <- 2
source(paste0(cm_path, "/R/covidm.R"))
source(paste0("code/util_functions.R"))

model_selected_ur <- read_rds(paste0(path_euro,
                                     "data/intermediate/DEoptim3_selected.rds"))
country_dictionary <- read_rds(paste0(path_euro, 
                                      "data/country_dictionary.rds"))
members_all <- model_selected_ur$iso3c
members_remove <- read_rds(paste0(path_dropbox, 
                                  "/intermediate/members_remove.rds"))
euro_lmic <- c("ALB","ARM","AZE","BLR","BIH","BGR","GEO","KAZ","XKX","KGZ",
               "MDA","MNE","MKD","RUS","SRB","TJK","TUR","TKM","UKR","UZB")
euro_inuse <- setdiff(euro_lmic, members_remove)

#### contact data ####
# updated contact matrices
load(paste0(path_dropbox, "contact_all.rdata"))
load(paste0(path_dropbox, "contact_work.rdata"))
load(paste0(path_dropbox, "contact_home.rdata"))
load(paste0(path_dropbox, "contact_school.rdata"))
load(paste0(path_dropbox, "contact_others.rdata"))

model_selected_ur %<>% 
  mutate(to_replace = iso3c %in% names(contact_all)) %>% 
  left_join(country_dictionary, c("iso3c" = "iso"))

##### age group labels #####
tmp <- cm_parameters_SEI3R("Thailand")
ag_labels <- tmp$pop[[1]]$group_names; rm(tmp)

##### assign Prem 2021 contact matrices #####
for(i in 1:nrow(model_selected_ur)){
  if(model_selected_ur$to_replace[i]){
    cm_matrices[[model_selected_ur$country_name[i]]]$home <-
      as.matrix(contact_home[[model_selected_ur$iso3c[i]]]) %>% 
      set_colnames(ag_labels) %>% 
      set_rownames(ag_labels)
    
    cm_matrices[[model_selected_ur$country_name[i]]]$work <-
      as.matrix(contact_work[[model_selected_ur$iso3c[i]]]) %>% 
      set_colnames(ag_labels) %>% 
      set_rownames(ag_labels)
    
    cm_matrices[[model_selected_ur$country_name[i]]]$school <-
      as.matrix(contact_school[[model_selected_ur$iso3c[i]]]) %>% 
      set_colnames(ag_labels) %>% 
      set_rownames(ag_labels)
    
    cm_matrices[[model_selected_ur$country_name[i]]]$other <-
      as.matrix(contact_others[[model_selected_ur$iso3c[i]]]) %>% 
      set_colnames(ag_labels) %>% 
      set_rownames(ag_labels)
  }
}

# fix discrepancies among country names
names(cm_matrices)[names(cm_matrices) == "TFYR of Macedonia"] <- 
  "North Macedonia"

# proxy for contact matrices (currently unavailable ?)
cm_matrices[["Republic of Moldova"]] <- cm_matrices$Romania
cm_matrices[["Turkmenistan"]] <- cm_matrices$Uzbekistan


#### load epidemic parameters ####
#####  Clinical Fraction #####
# (based on Davies et al, Nature paper) 
cf <- c(
  0.2904047, 0.2904047, 0.2070468, 0.2070468, 0.2676134,
  0.2676134, 0.3284704, 0.3284704, 0.3979398, 0.3979398,
  0.4863355, 0.4863355, 0.6306967, 0.6306967, 0.6906705, 0.6906705
)

##### susceptibility #####
# (based on Davies et al, Nature paper)
sus <- c(
  0.3956736, 0.3956736, 0.3815349, 0.3815349, 0.7859512,
  0.7859512, 0.8585759, 0.8585759, 0.7981468, 0.7981468,
  0.8166960, 0.8166960, 0.8784811, 0.8784811, 0.7383189, 0.7383189
)

#### Google mobility ####
gm_type <- c("retail", "grocery", "parks", "transit", "work", "residential")
gm <- read_rds(paste0(path_dropbox, "gm.rds"))

##### mobility scalers ######
curves <- data.table(
  work_scaler = c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0.008, 0.021, 0.033, 0.046, 0.058, 0.071, 0.083, 0.096, 0.108, 0.121, 0.133,
    0.146, 0.158, 0.171, 0.183, 0.196, 0.208, 0.221, 0.233, 0.246, 0.258, 0.271,
    0.283, 0.296, 0.308, 0.321, 0.334, 0.346, 0.359, 0.371, 0.384, 0.397, 0.41,
    0.422, 0.435, 0.448, 0.461, 0.474, 0.487, 0.5, 0.513, 0.526, 0.539, 0.552,
    0.566, 0.579, 0.592, 0.606, 0.619, 0.633, 0.646, 0.66, 0.674, 0.687, 0.701,
    0.715, 0.729, 0.743, 0.757, 0.771, 0.785, 0.799, 0.813, 0.828, 0.842, 0.856,
    0.87, 0.885, 0.899, 0.914, 0.928, 0.942, 0.957, 0.971, 0.986, 1, 1.014, 1.029,
    1.043, 1.058, 1.072, 1.087, 1.101, 1.115, 1.13, 1.144, 1.159, 1.173, 1.188,
    1.202, 1.216, 1.231, 1.245, 1.26, 1.274, 1.289, 1.303, 1.317, 1.332, 1.346, 1.361
  ),
  other_scaler = c(
    0.064, 0.066, 0.067, 0.068, 0.069, 0.071, 0.072, 0.073, 0.075, 0.076, 0.077, 0.078,
    0.08, 0.081, 0.082, 0.084, 0.085, 0.086, 0.087, 0.089, 0.09, 0.091, 0.092, 0.094,
    0.095, 0.096, 0.098, 0.099, 0.1, 0.101, 0.103, 0.104, 0.105, 0.106, 0.108, 0.109,
    0.11, 0.112, 0.113, 0.114, 0.116, 0.118, 0.119, 0.121, 0.123, 0.125, 0.128, 0.13,
    0.132, 0.135, 0.137, 0.14, 0.143, 0.146, 0.15, 0.154, 0.159, 0.164, 0.169, 0.175,
    0.182, 0.19, 0.198, 0.207, 0.217, 0.228, 0.24, 0.252, 0.266, 0.28, 0.295, 0.31,
    0.327, 0.344, 0.361, 0.379, 0.398, 0.418, 0.438, 0.459, 0.48, 0.502, 0.525, 0.549,
    0.572, 0.597, 0.621, 0.647, 0.672, 0.698, 0.725, 0.751, 0.778, 0.805, 0.833, 0.86,
    0.888, 0.916, 0.944, 0.972, 1, 1.028, 1.056, 1.084, 1.112, 1.14, 1.168, 1.196, 1.224,
    1.252, 1.28, 1.308, 1.337, 1.365, 1.393, 1.421, 1.449, 1.477, 1.505, 1.533, 1.561,
    1.589, 1.617, 1.645, 1.673, 1.701
  ),
  perc = round(seq(0, 1.25, 0.01), 2)
)

si <- qread(paste0(path_dropbox, "si.qs"))
oxcgrt <- qread(paste0(path_dropbox,"oxcgrt.qs"))
print(paste0("The Stringency Index and School Closure Data in used was downloaded on ", 
             file.info(paste0(path_dropbox, "si.qs"))$mtime, "."))

si %>%
  filter(!is.na(StringencyIndex)) %>%
  group_by(date) %>%
  tally() %>%
  mutate(
    n_max = max(n),
    missing = (n_max - n) / n_max
  ) %>%
  filter(missing > 0.1,
         date > "2021-01-01") %>%
  pull(date) %>%
  min() -> si_stopdate

owid_epi <- qread(paste0(path_euro, "data/epi.qs")) 
owid_vac <- qread(paste0(path_dropbox, "owid_vac.qs"))
schedule_raw <- read_rds(paste0(path_dropbox,"schedule_raw.rds"))

# schedule_raw %>% 
#   ggplot(., aes(x = date, y = work, color = status)) +
#   geom_point() +
#   facet_wrap(~wb)

#### vaccine efficacy tested ####
data.table(ve_i_o = c(0.67, 0.68),
           ve_d_o = c(0.67, 0.78)) %>% 
  mutate(
    ve_d = exp_ve(ve_d_o, ve_i_o),
    ve_h = c(0.845, 0.9),
    ve_mort = c(0.845, 0.95)) -> ve

#### burden processes ####
critical2 <- 0
picu_cocin_func <- function(age) {
  x <- c(-0.1309118, 0, 17.2398874, 65.7016492, 100)
  y <- c(-2.1825091, -2.1407043, -1.3993552, -1.2344361, -8.8191062)
  p <- splinefun(x, y)(age)
  exp(p) / (1 + exp(p))
}
picu_cocin <- picu_cocin_func(0:85)

# Infection fatality rate (derived from Levin et al., preprint)
ifr_levin <- 100 * exp(-7.56 + 0.121 * 0:85) / (100 + exp(-7.56 + 0.121 * 0:85)) / 100
# Infection hospitalisation rate (derived from Salje et al., Science)
ihr_salje <- exp(-7.37 + 0.068 * 0:85) / (1 + exp(-7.37 + 0.068 * 0:85))
# Amalgamate probabilities
probabilities <- data.table(age = 0:85, ihr = ihr_salje, ifr = ifr_levin, picu = picu_cocin)
probabilities[, age_group := pmin(15, age %/% 5)]
probabilities <- probabilities[, lapply(.SD, mean), by = age_group, .SDcols = 2:4]

# Create model burden processes
P.critical <- probabilities[, ihr * picu]
P.severe <- probabilities[, ihr * (1 - picu)]
P.death <- probabilities[, ifr]
P.hosp <- P.critical + P.severe

delay_2death <- cm_delay_gamma(26, 5, 60, 0.25)$p
delay_2severe <- cm_delay_gamma(8.5, 5, 60, 0.25)$p
delay_2hosp <- cm_delay_gamma(14.6, 5, 60, 0.25)$p

burden_processes <- list(
  cm_multinom_process("E",       
                      data.frame(death = P.death),                   
                      delays = data.frame(death = delay_2death), report = "o"),
  cm_multinom_process("Ev",      
                      data.frame(death = P.death*(1-ve$ve_mort[1])), 
                      delays = data.frame(death = delay_2death), report = "o"),
  cm_multinom_process("Ev2",     
                      data.frame(death = P.death*(1-ve$ve_mort[2])), 
                      delays = data.frame(death = delay_2death), report = "o"),
  
  
  cm_multinom_process("E",       
                      data.frame(to_hosp = P.hosp),                  
                      delays = data.frame(to_hosp = delay_2severe)),
  cm_multinom_process("Ev",      
                      data.frame(to_hosp = P.hosp*(1-ve$ve_h[1])),   
                      delays = data.frame(to_hosp = delay_2severe)),
  cm_multinom_process("Ev2",     
                      data.frame(to_hosp = P.hosp*(1-ve$ve_h[2])),   
                      delays = data.frame(to_hosp = delay_2severe)),
  
  cm_multinom_process("to_hosp", 
                      data.frame(hosp = rep(1,16)),                  
                      delays = data.frame(hosp = delay_2hosp),   
                      report = "ip")
)
