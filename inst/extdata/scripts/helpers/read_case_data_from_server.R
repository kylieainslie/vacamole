# ----------------------------------------------------
# read in updated OSIRIS data from RIVM server
# ----------------------------------------------------

# get file path --------------------------------------
path <- "/rivm/r/COVID-19/Surveillance/Data/OSIRIS/Geschoond/"
file <- list.files(path, pattern = ".rds")
if (identical(file, character(0))) {
  path <- paste0(path,"Previous/")
  file <- list.files(path, pattern = ".rds") %>%
    max()
}

# read in data ---------------------------------------
osiris <- readRDS(paste0(path,file))

# aggregate for number of cases per day --------------
# this removes any identifiable data

osiris_tally <- osiris %>%
  select(OSIRISNR, INFECTIEZIEKTE, ZIE1eZiekteDt, Land) %>%
  filter(Land == "Nederland",
         INFECTIEZIEKTE %in% c("NCOV", "Weak Positive", "Antitgen Pos. + Symptoms", "PCR Positief", "Antigen Positief")) %>%
  select(-Land) %>%
  rename(date = ZIE1eZiekteDt) %>%
  group_by(date) %>%
  summarise(inc = n()) %>%
  filter(!is.na(date)) %>%
  complete(date = seq.Date(min(date), max(date), by="day"), fill = list(inc = 0))

# calculate rolling average and remove last 3 days ---
osiris1 <- osiris_tally %>%
  mutate(roll_avg = zoo::rollmean(inc, k = 7, fill = 0)) %>%
  filter(date < max(date)-2) # remove last 3 days due to reporting delay

last_date_in_osiris <- tail(osiris1$date, 1)
# save output ----------------------------------------
saveRDS(osiris1, file = paste0("inst/extdata/data/case_data_upto_", last_date_in_osiris, ".rds"))


# aggregate for number of cases per day by age group -
osiris_tally_by_age <- osiris %>%
  select(OSIRISNR, INFECTIEZIEKTE, ZIE1eZiekteDt, Land, Leeftijdsgroep10) %>%
  filter(Land == "Nederland",
         INFECTIEZIEKTE %in% c("NCOV", "Weak Positive", "Antitgen Pos. + Symptoms", "PCR Positief", "Antigen Positief")) %>%
  select(-Land) %>%
  rename(date = ZIE1eZiekteDt,
         age_group = Leeftijdsgroep10) %>%
  mutate(age_group = factor(age_group, levels = c("0-9", "10-19", "20-29", "30-39", "40-49",
                                                  "50-59", "60-69", "70-79", "80+"))) %>%
  group_by(date, age_group) %>%
  summarise(inc = n()) %>%
  filter(!is.na(date)) %>%
  complete(date = seq.Date(min(date), max(date), by="day"), fill = list(inc = 0))

last_date_in_osiris <- tail(osiris_tally_by_age$date, 1)
# save output ----------------------------------------
saveRDS(osiris_tally_by_age, file = paste0("../case_data_by_age_group_upto_", last_date_in_osiris, ".rds"))
