## Download and read the IntCal13, Marine13, and SHCal13 calibration curves
library(tidyverse)

# IntCal13
download.file("http://www.radiocarbon.org/IntCal13%20files/intcal13.14c", "data-raw/intcal13.14c")
intcal13 <- readr::read_csv("data-raw/intcal13.14c", skip = 9)[-1, ] %>%
  dplyr::rename(`CAL BP` = `# CAL BP`) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `CAL BP`,
      `14C age`,
      Error
    ),
    as.integer
  ) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `Delta 14C`,
      Sigma
    ),
    as.numeric
  )

usethis::use_data(intcal13,
  overwrite = T
)

# Marine13
download.file("http://www.radiocarbon.org/IntCal13%20files/marine13.14c", "data-raw/marine13.14c")
marine13 <- readr::read_csv("data-raw/marine13.14c", skip = 9)[-1, ] %>%
  dplyr::rename(`CAL BP` = `# CAL BP`) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `CAL BP`,
      `14C age`,
      Error
    ),
    as.integer
  ) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `Delta 14C`,
      Sigma
    ),
    as.numeric
  )

usethis::use_data(marine13,
  overwrite = T
)

# SHCal13
download.file("http://www.radiocarbon.org/IntCal13%20files/shcal13.14c", "data-raw/shcal13.14c")
shcal13 <- readr::read_csv("data-raw/shcal13.14c", skip = 9)[-1, ] %>%
  dplyr::rename(`CAL BP` = `# CAL BP`) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `CAL BP`,
      `14C age`,
      Error
    ),
    as.integer
  ) %>%
  dplyr::mutate_at(
    .vars = dplyr::vars(
      `Delta 14C`,
      Sigma
    ),
    as.numeric
  )

usethis::use_data(shcal13,
  overwrite = T
)

# Raw data for IntCal13
# You must have access to the SQL server where the data are hosted; see http://intcal.qub.ac.uk/shcal13/query/Rhelp.html
library(DBI)
library(RMySQL)
library(writexl)

con <- dbConnect(MySQL(),
                 user = "customer",
                 dbname = "intcalx",
                 host = "intcal.qub.ac.uk")

intcal13_raw <- dbListTables(con) %>%
  magrittr::set_names(.,.) %>%
  purrr::map(function(x){
    con %>%
      dbGetQuery(paste0("select * from ", x)) %>%
      tibble::as_tibble()
  }) %T>%
  writexl::write_xlsx("data-raw/intcal13_raw.xlsx")

usethis::use_data(intcal13_raw,
                  overwrite = T
)

