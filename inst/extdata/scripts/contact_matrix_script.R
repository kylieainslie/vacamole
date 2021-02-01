# Create symmetrical contact matrix for model input

# read in data
contact_matrix_wo_hh <- readRDS("/s-schijf/ainsliek/vac-a-mole/inst/extdata/data/Contactmatrix_withoutHH_2020-12-04.rds")
populationNL_2020 <- readRDS("/s-schijf/ainsliek/vac-a-mole/inst/extdata/data/populationNL_2020.rds")

# summarise pop data into age groups
pop_NL_summary <- populationNL_2020 %>%
  mutate(age_group = case_when(
    age < 10 ~ 
  ))