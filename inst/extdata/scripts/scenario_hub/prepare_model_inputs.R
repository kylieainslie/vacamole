# Script for preparing and specifying model inputs for running scenarios for the 
# European Scenario Hub
# URL: https://github.com/covid19-forecast-hub-europe/covid19-scenario-hub-europe#readme

# preamble ---------------------------------------------------------------------
# This script will specify all of the model inputs and (if necessary)
# write them as RDS files, which will be called in main_script.R
# ------------------------------------------------------------------------------

# osiris case data --------------------------------------------------
# this must be run on the RIVM R servers
path <- "/rivm/r/COVID-19/Surveillance/Data/OSIRIS/Geschoond/"
file <- list.files(path, pattern = ".rds")
if (identical(file, character(0))) {
  path <- paste0(path,"Previous/")
  file <- list.files(path, pattern = ".rds") %>%
    max()
}

osiris <- readRDS(paste0(path,file)) # read in file from path

osiris_tally <- osiris %>%           # aggregate for number of cases per day
  # this removes any identifiable data
  select(OSIRISNR, INFECTIEZIEKTE, ZIE1eZiekteDt, Land) %>%
  filter(Land == "Nederland",
         INFECTIEZIEKTE %in% c("NCOV", "Weak Positive", 
                               "Antitgen Pos. + Symptoms", 
                               "PCR Positief", "Antigen Positief")) %>%
  select(-Land) %>%
  rename(date = ZIE1eZiekteDt) %>%
  group_by(date) %>%
  summarise(inc = n()) %>%
  filter(!is.na(date)) %>%
  complete(date = seq.Date(min(date), max(date), by="day"), fill = list(inc = 0))

cutoff_date <- as.Date("2022-03-12")

osiris1 <- osiris_tally %>%
  filter(date <= cutoff_date)

# vaccinations params ----------------------------------------------------------
# delay to protection ----------------------------------------------------------
delays <- list(
  pfizer = c(14, 7, 7, 7, 7),
  moderna = c(14, 7, 7, 7, 7), 
  astrazeneca = c(14, 7),
  jansen = c(14, 7)
)

# VE against infection by dose and vaccine -------------------------------------
ve_inf_list <- list(
  wildtype = list(
    pfizer = c(0.926, 0.948, c(rep(0.948, 3))), # from clinical trial
    moderna = c(0.896, 0.941, c(rep(0.948, 3))), # from clinical trial
    astrazeneca = c(0.583, 0.621), # from clinical trial
    jansen = c(0.661, 0.661) # from clinical trial
  ),
  alpha = list(
    pfizer = c(0.66, 0.8, c(rep(0.8, 3))), # from Pritchard et al. 2021 Nature
    moderna = c(0.66, 0.8, c(rep(0.8, 3))), # assumed to be the same as pfizer
    astrazeneca = c(0.61, 0.79), # from Pritchard et al. 2021 Nature 
    jansen = c(0.767, 0.767) # from Corchado-Garcia et al. 2021 medRxiv (need to check if this is against alpha!)
  ),
  delta = list( # from Dutch data sources
    pfizer = c(0.57, 0.69, 0.93, 0.93, 0.93), 
    moderna = c(0.66, 0.82, 0.93, 0.93, 0.93), 
    astrazeneca = c(0.41, 0.54), 
    jansen = c(0.5, 0.5)
  ),
  omicron = list( # from https://www.medrxiv.org/content/10.1101/2022.02.06.22270457v2
    pfizer = c(0, 0.33, 0.68, 0.65, 0.65),  # Grewal et al. 2022 (4th dose) - assuming same for 5th dose
    moderna = c(0, 0.33, 0.68, 0.65, 0.65), # Grewal et al. 2022 (4th dose) - assuming same for 5th dose
    astrazeneca = c(0, 0.33), 
    jansen = c(0, 0.33)
  )
)

# VE against transmission (for vaccine failures) by dose and vaccine -----------
ve_trans_list <- list(
  wildtype = list( # same as alpha !!!
    pfizer = c(0.26, 0.70, c(rep(0.70, 3))),      
    moderna = c(0.51, 0.88, c(rep(0.88, 3))),     
    astrazeneca = c(0.15, 0.58),
    jansen = c(0.77, 0.77)
  ),
  alpha = list(  # de Gier et al.
    pfizer = c(0.26, 0.70, c(rep(0.70, 3))),      
    moderna = c(0.51, 0.88, c(rep(0.88, 3))),     
    astrazeneca = c(0.15, 0.58),
    jansen = c(0.77, 0.77)
  ),
  delta = list( # de Gier et al. (updated)
    pfizer = c(0.46, 0.52, c(rep(0.52, 3))),      
    moderna = c(0.66, 0.24, c(rep(0.24, 3))),     
    astrazeneca = c(0, 0.25),
    jansen = c(0.42, 0.42) 
  ),
  omicron = list( # MADE-UP VALUES
    pfizer = c(0.25, 0.33, 0.4, 0.4, 0.4),      
    moderna = c(0.25, 0.33, 0.4, 0.4, 0.4),     
    astrazeneca = c(0, 0.25),
    jansen = c(0.25, 0.33) 
  )
)

# VE against hospitalisation (for vaccine failures) by dose and vaccine --------
ve_hosp_list <- list(
  wildtype = list( # same as alpha variant!!!
    pfizer = c(0.81, 0.95, 0.95, 0.95, 0.95),      # Dutch data
    moderna = c(0.81, 0.95, 0.95, 0.95, 0.95),     # assumed same as pfizer because Dutch estimates were weird
    astrazeneca = c(0.83, 0.95), # Dutch data
    jansen = c(0.85, 0.85)             # from RIVM website: https://www.rivm.nl/en/covid-19-vaccination/vaccines/efficacy-and-protection
  ),
  alpha = list(
    pfizer = c(0.81, 0.95, 0.95, 0.95, 0.95),      # Dutch data
    moderna = c(0.81, 0.95, 0.95, 0.95, 0.95),     # assumed same as pfizer because Dutch estimates were weird
    astrazeneca = c(0.83, 0.95), # Dutch data
    jansen = c(0.85, 0.85)             # from RIVM website: https://www.rivm.nl/en/covid-19-vaccination/vaccines/efficacy-and-protection
  ),
  delta = list(
    pfizer = c(0.89, 0.96, 0.98, 0.98, 0.98),      # from Brechje (pre-print)
    moderna = c(0.95, 0.85, 0.98, 0.98, 0.98),     # from Brechje (pre-print)
    astrazeneca = c(0.88, 0.94),       # from Brechje (pre-print)
    jansen = c(0.92, 0.92)                   # from Brechje (pre-print)
  ),
  omicron = list( # from https://www.rivm.nl/documenten/effectiviteit-van-covid-19-vaccinatie-tegen-ziekenhuis-en-intensive-care-opname-in-8
    pfizer = c(0, 0.56, 0.88, 0.92, 0.92),   # Grewal et al. 2022 (4th dose) - assuming same for 5th dose
    moderna = c(0, 0.49, 0.82, 0.92, 0.92),  # Grewal et al. 2022 (4th dose) - assuming same for 5th dose
    astrazeneca = c(0, 0.39),
    jansen = c(0.73, 0.73)
  )
)

# hospitalisations multiplier
# calculated as (1-ve_hosp)/(1-ve_inf)
calc_mult_fun <- function(ve_hosp,ve_inf){
  h_multiplier <- list()
  
  for (i in 1:length(ve_hosp)){
    h_multiplier[[i]] <- list(
      pfizer = (1-ve_hosp[[i]]$pfizer)/(1-ve_inf[[i]]$pfizer),
      moderna = (1-ve_hosp[[i]]$moderna)/(1-ve_inf[[i]]$moderna),
      astrazeneca = (1-ve_hosp[[i]]$astrazeneca)/(1-ve_inf[[i]]$astrazeneca),
      jansen = (1-ve_hosp[[i]]$jansen)/(1-ve_inf[[i]]$jansen)
    )
  }
  names(h_multiplier) <- names(ve_hosp)
  return(h_multiplier)
}

h_mult_list <- calc_mult_fun(ve_hosp = ve_hosp_list, ve_inf = ve_inf_list)

# create list of vaccine params and save as RDS file ---------------------------
ve_list <- list(delays = delays,
                ve_inf = ve_inf_list, 
                ve_hosp = h_mult_list,
                ve_trans = ve_trans_list)
saveRDS(ve_list, "inst/extdata/inputs/ve_params.rds")

# Vaccination schedule ---------------------------------------------------------
# Read in file and change column names for booster doses
# vac_path <- "C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/manuscripts/impact_vac/data/vaccination_scenarios/"
vac_path <- "/rivm/s/ainsliek/data/"

vac_sched <- read_csv(paste0(vac_path,"Cum_upt20220503.csv")) %>%
  rename_with(~ gsub("B1", "d3", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("B2", "d4", .x, fixed = TRUE)) %>%
  select(-starts_with("X")) %>%
  mutate(date = as.Date(date, format = "%m/%d/%Y"))

# add columns for 5th dose of pfizer and moderna vaccines
new_columns <- c(paste0("pf_d5_", 1:10), paste0("mo_d5_", 1:10))
vac_sched[,new_columns] <- 0

# rearrange columns to preserve correct order (also exclude Novovax doses)
vac_sched1 <- vac_sched %>%
  select(date, pf_d1_1:pf_d4_10, pf_d5_1:pf_d5_10,
         mo_d1_1:mo_d4_10, mo_d5_1:mo_d5_10,
         az_d1_1:az_d2_10, ja_d1_1:ja_d2_10)

# create "empty" values from 1/1/2020 until start of vac sched (1/4/2021)
n_cols <- dim(vac_sched1)[2]-1 # exclude date column
n_rows <- vac_sched1$date[1] - osiris1$date[1]
empty_mat <- matrix(rep(0, n_cols * n_rows), nrow = n_rows)
dates <- seq.Date(osiris1$date[1], vac_sched1$date[1]-1, by = "day")
my_df <- data.frame(date = dates, empty_mat)
names(my_df) <- names(vac_sched1)
vac_schedule <- bind_rows(my_df, vac_sched1)

# write out to directory
saveRDS(vac_schedule,"inst/extdata/inputs/vac_schedule_real_w_4th_and_5th_dose.csv")

