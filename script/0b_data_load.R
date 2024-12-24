
# packages ----------------------------------------------------------------

library(tidyverse)
library(openxlsx)

# load data ---------------------------------------------------------------

datafile_cont_BA2 <- read.xlsx('./data/Omicron_BA2_contact.xlsx')
datafile_chains_BA2 <- read.xlsx('./data/Omicron_BA2_chains.xlsx') |> filter(!is.na(from))
datafile_info_BA2 <- read.xlsx('./data/Omicron_BA2.xlsx')

datafile_cont_BA1 <- read.xlsx('./data/Omicron_BA1_contact.xlsx')
datafile_chains_BA1 <- read.xlsx('./data/Omicron_BA1_chains.xlsx') |> filter(!is.na(from))
datafile_info_BA1 <- read.xlsx('./data/Omicron_BA1.xlsx')

datafile_cont_Delta <- read.xlsx('./data/Delta_contact.xlsx')
datafile_chains_Delta <- read.xlsx('./data/Delta_chains.xlsx') |> filter(!is.na(from))
datafile_info_Delta <- read.xlsx('./data/Delta.xlsx')

# convert date -------------------------------------------------------------

datafile_cont_BA2 <- datafile_cont_BA2 |> 
  mutate_at(vars(contains('date')), convertToDate) |> 
  mutate(vaccine = if_else(vaccine > 0 & 
                             vaccinelastdate >= lastexposedate - 14 &
                             !is.na(vaccinelastdate) &
                             !is.na(lastexposedate),
                           vaccine - 1,
                           vaccine))
datafile_info_BA2 <- datafile_info_BA2 |> 
  mutate_at(vars(contains('date')), convertToDate)

datafile_cont_BA1 <- datafile_cont_BA1 |> 
  mutate_at(vars(contains('date')), convertToDate)
datafile_info_BA1 <- datafile_info_BA1 |> 
  mutate_at(vars(contains('date')), convertToDate)

datafile_cont_Delta <- datafile_cont_Delta |> 
  mutate_at(vars(contains('date')), convertToDate) |> 
  mutate(vaccine = if_else(vaccine > 0 & 
                             vaccinelastdate >= lastexposedate - 14 &
                             !is.na(vaccinelastdate) &
                             !is.na(lastexposedate)&
                             drop == 'N',
                           vaccine - 1,
                           vaccine)) |> 
  select(-drop)

datafile_info_Delta <- datafile_info_Delta |> 
  mutate_at(vars(contains('date')), convertToDate)


save.image('./data/sars_2_cov.Rdata')

datafile_cont_Delta$gender[!datafile_cont_Delta$gender %in% c('Male', 'Female')] <- ''
