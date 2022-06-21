# Process data sources for input into model

# Load required packages
library(dplyr)
library(readr)
library(here)

### Read in data files
file_path <- "C:/Users/ainsliek/Dropbox/Kylie/Projects/RIVM/scenario_hub_europe/round0/data/"
cases <- read_csv(paste0(file_path, "truth_JHU-Incident Cases.csv"))

## vaccination
# Get 32 Euro hub locations
hub <- read_csv("https://raw.githubusercontent.com/covid19-forecast-hub-europe/covid19-scenario-hub-europe/main/data-locations/locations_eu.csv") %>%
  select("location_name", "location", "population")

# Source: Our World In Data
# Metadata: https://github.com/owid/covid-19-data/tree/master/public/data

# Get total population vaccination data
owid_vax <- read_csv("https://github.com/owid/covid-19-data/blob/master/public/data/vaccinations/vaccinations.csv?raw=true") %>%
  filter(location %in% hub$location_name) %>%
  rename("location_name" = "location") %>%
  select(-"iso_code") %>%
  left_join(hub, by = "location_name") %>%
  select("location", "location_name", "population", "date", everything())

# Get vaccination data by age
# - Note: age data not available for Denmark (DE), Greece (GR), Netherlands (NL), United Kingdom (GB)
owid_vax_age <- read_csv("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/vaccinations/vaccinations-by-age-group.csv") %>%
  filter(location %in% hub$location_name) %>%
  rename("location_name" = "location") %>%
  left_join(hub, by = "location_name") %>%
  select("location", "location_name", "population", "date", everything())

# contact matrices (NL only)
list_contact_files <- list.files(path = paste0(file_path, "contact_matrices/NL/"), pattern="*.csv")
files <- paste0(file_path, "contact_matrices/NL/",list_contact_files)
my_contact_files <- files %>%
  map(read_csv)
names(my_contact_files) <- gsub('.csv', '', list_contact_files)

# subset cases and vaccination by country
# we'll focus on NL for now
cases_nl <- cases %>%
  filter(location == "NL")

vac_nl <- vac_by_age_group %>%
  filter(location == "Netherlands")
