getPopulationData <- function(dir = "./", year) {
  
  tmp <- str_detect(list.files(dir), "population.csv")
  if(any(tmp)) {
    # Read-in downloaded CBS data ----
    pop.data <- read_csv(file = paste0(dir, "/", list.files(dir)[tmp]), col_types = "ifii")
  } else {
    # Load packages and download CBS data ----
    library(cbsodataR)
    library(stringi)
    
    # # 1. See what is available
    # cbs.list <- cbs_get_toc(Language = "nl") %>%
    #   as.data.frame #%>% View
    # 
    # # 2. Search for tables with bevolking, leeftijd, geslacht
    # cbs.list %>%
    #   select(Identifier, Title) %>%
    #   filter(
    #     Title %>% str_detect(pattern = ".evolking") &
    #       Title %>% str_detect(pattern = ".eeftijd") &
    #       Title %>% str_detect(pattern = ".eslacht"))
    # 
    # # 3. Choose table and see structure of metadata (the 'codebook')
    # #    From this we can see what to download
    # cbs_get_meta(id = "7461bev") %>% str
    
    # Download NL age- and sex-specific population data 2017
    pop.data <- cbs_get_data(
      id = "7461bev",
      select = c("Perioden", "Leeftijd", "Geslacht", "BurgerlijkeStaat", "Bevolking_1"),
      #Perioden = paste0(2007:2017,"JJ00"),
      Leeftijd = c(10010, seq(from = 10100, to = 19800, by = 100), 22100) %>% as.character,
      Geslacht = c("3000   ", "4000   "),
      BurgerlijkeStaat = "T001019") %>% 
      transmute(
        year = Perioden %>% str_replace_all(pattern = "JJ00", replacement = "") %>% as.integer,
        gender = Geslacht %>% as.integer %>% factor(levels = c(3000, 4000), labels = c("Male", "Female")),
        age = Leeftijd %>% factor(levels = c(10010, seq(from = 10100, to = 19800, by = 100), 22100), labels = 0:99) %>% as.character %>% as.integer,
        population = as.integer(Bevolking_1))
    
    # write_csv(pop.data, path = "data/NL_population.csv")
  }
  
  yr <- year
  return(pop.data %>% filter(year == yr))
}
