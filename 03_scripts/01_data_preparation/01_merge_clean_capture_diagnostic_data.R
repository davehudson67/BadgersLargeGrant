library(tidyverse)
library(lubridate)

rm(list=ls())

## load capture data ----------------------------------------------------------------------------
captures <- read_csv("Data/tblCaptures.csv")

captures <- captures %>%
  select(tattoo, date, sett, pm, socg)

captures$date <- dmy(captures$date)

#captures <- captures %>%
#  mutate(capture_yr = year(date))

#m change 040l to 040Lin capture
captures <- as.data.frame(captures)
captures[captures$tattoo=="040l"] <- "040L"

## Load badger data------------------------------------------------------------------------------
badgers <- read_csv("Data/tblBadger.csv")

badgers <- badgers %>%
  select(tattoo, age_fc, year_fc, sex)

badgers$sex <- as.factor(badgers$sex)
summary(badgers)

## change NAs to UNKNOWN for sex
badgers$sex[is.na(badgers$sex)] <- "UNKNOWN"

## check for multiple sex entries
badgers %>%
  group_by(tattoo, sex) %>%
  count(sort = TRUE)

## check
summary(badgers)

## combine datasets
bc <- left_join(captures, badgers, by = "tattoo")

## adjust pm variable
bc$pm[bc$pm == "YES"] <- "Yes"
bc$pm[bc$pm == "NO"] <- "No"

bc$age_fc[bc$tattoo == "067B"]  <- "ADULT"

## checks
stopifnot(any(!is.na(bc$age_fc)))
stopifnot(any(!is.na(bc$year_fc)))

rm(badgers)
rm(captures)
#bc <- bc %>%
#  group_by(tattoo) %>%
#  mutate(year_fc = min(capture_yr)) %>%
#  ungroup() %>%
#  mutate(age_fc_yrs = ifelse(age_fc == "CUB", 0, 
#                             ifelse(age_fc == "YEARLING", 1, NA))) %>%
#  filter(!is.na(age_fc_yrs)) %>%
#  mutate(age_yrs = capture_yr - year_fc)

## load infection data CULTURES---------------------------------------------------------------------
culture <- read_csv("Data/tblCulture.csv")

## set date
culture$date <- dmy(culture$date)

## adjust pm variable
culture$pm[culture$pm == "YES"] <- "Yes"
culture$pm[culture$pm == "NO"] <- "No"

#change levels for cult$results
culture <- culture %>%
  mutate(result2 = ifelse(result == "M.BOVIS", 1, 
                          ifelse(result == "NEGATIVE", 0, NA))) %>%
  mutate(sample2 = ifelse(sample == "OA" | sample == "TA"| sample == "URINE" | sample=="FAECES", sample, "OTHER"))

## remove NAs
culture <- culture[!is.na(culture$result2),]

#generate columns for overall result and breakdown according to sample type (standard 4 clinical plus 'other')
culture <- culture %>%
  group_by(tattoo, date) %>%
  mutate(Cult_Urine = ifelse(sample2 == "URINE", result2, 0)) %>%
  mutate(Cult_Urine = max(Cult_Urine)) %>%
  mutate(Cult_Faeces = ifelse(sample2 == "FAECES", result2, 0)) %>%
  mutate(Cult_Faeces = max(Cult_Faeces)) %>%
  mutate(Cult_OA = ifelse(sample2 == "OA", result2, 0)) %>%
  mutate(Cult_OA = max(Cult_OA)) %>%
  mutate(Cult_TA = ifelse(sample2 == "TA", result2, 0)) %>%
  mutate(Cult_TA = max(Cult_TA)) %>%
  mutate(Cult_Other = ifelse(sample2 == "OTHER", result2, 0)) %>%
  mutate(Cult_Other = max(Cult_Other)) %>%
  ungroup() %>%
  distinct(tattoo, date, .keep_all = TRUE)

culture <- culture %>%
  group_by(tattoo, date) %>%
  mutate(Cult_Sum = sum(c(Cult_Urine, Cult_Faeces, Cult_OA, Cult_TA, Cult_Other))) %>%
  ungroup() %>%
  select(tattoo, date, pm, Cult_Urine, Cult_Faeces, Cult_OA, Cult_TA, Cult_Other, Cult_Sum)

## Load DPP data---------------------------------------------------------------------------------------------------
DPP <- read_csv("Data/tblDPPTest.csv")

## adjust pm variable
DPP$PM[DPP$PM == "YES"] <- "Yes"
DPP$PM[DPP$PM == "NO"] <- "No"

DPP$VisualLine1 <- as.factor(DPP$VisualLine1)
summary(DPP$VisualLine1)

DPP$VisualLine2 <- as.factor(DPP$VisualLine2)
summary(DPP$VisualLine2)

DPP <- DPP %>%
  mutate(date = dmy(CaptureDate)) %>%
  mutate(tattoo = BadgerID) %>%
  mutate(pm = PM) %>%
  mutate(v1 = if_else(VisualLine1 == "N" | VisualLine1 == "Nx", 0, 1)) %>%
  mutate(v2 = if_else(VisualLine2 == "N" | VisualLine2 == "Nx", 0, 1))

## some NAs in v2
#DPP$v2[is.na(DPP$v2)] <- -1

DPP <- DPP %>%
  group_by(tattoo, date) %>%
  mutate(v1s = max(v1, na.rm = TRUE)) %>%
  mutate(v2s = max(v2, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(tattoo, date, .keep_all = TRUE) %>%
  select(tattoo, date, pm, v1s, v2s)


## Load IFNGamma data-------------------------------------------------------------------------------------------------
IFN <- read_csv("Data/tblIFN_Gamma_ELISA.csv")

IFN$date <- dmy(IFN$CaptureDate)
IFN$tattoo <- IFN$BadgerId
IFN$GAMMA <- IFN$Result

IFN <- IFN %>%
  select(tattoo, pm, date, GAMMA)

## check number of records per date
IFN %>%
  group_by(tattoo, date, GAMMA) %>%
  count(sort = TRUE)

IFN <- IFN %>%
  mutate(GAMMA = ifelse(IFN$GAMMA == "Positive", 1, 0))

#IFN$GAMMA[is.na(IFN$GAMMA)] <- -1

IFN <- IFN %>%
  group_by(tattoo, date) %>%
  mutate(GAMMA = max(GAMMA, na.rm = TRUE))

IFN$GAMMA[IFN$GAMMA == -Inf] <- NA

IFN <- IFN %>%
  arrange(tattoo, date)

IFN <- IFN %>%
  distinct(tattoo, date, GAMMA, .keep_all = TRUE)

## Combine data------------------------------------------------------------------------------------------------------

bc1 <- full_join(culture, IFN, by = c("tattoo", "date", "pm"))
bc1 <- full_join(bc1, DPP, by = c("tattoo", "date", "pm"))
bc <- full_join(bc, bc1, by = c("tattoo", "date", "pm"))

## Load previous data
pd <- read_csv("Data/all.diag.results.csv")

## keep only BROCK and statpak test info
pd <- pd %>%
  #  filter(!is.na(brock) | !is.na(statpak)) %>%
  select(tattoo, date, pm, brock, statpak)
#  filter(!is.na(age_years))

#colnames(pd)[which(names(pd) == "age_years")] <- "age_yrs"
pd$date <- dmy(pd$date)

bc1 <- full_join(bc, pd, by = c("tattoo", "date", "pm"))

bc1$age_fc[is.na(bc1$age_fc)] <- 0 # Make NAs = 0
bc1$age_fc <- as.factor(bc1$age_fc)
summary(bc1$age_fc)
bc1$age_fc <- as.numeric(bc1$age_fc)

bc1 <- bc1 %>%
  mutate(capture_yr = year(date)) %>%
  group_by(tattoo) %>%
  mutate(Kage = ifelse(any(age_fc == 3 | age_fc == 4), 1, 0)) %>% # 3 = cub, 4 = yearling
  ungroup() %>%
  filter(Kage == 1) # keep only known age badgers

## some badgers have more than one year_fc
errors <- bc1 %>%
  group_by(tattoo) %>%
  distinct(year_fc) %>%
  count(sort = TRUE)
errors <- errors$tattoo[errors$n == 2]

## adjust NAs to zero
bc1$year_fc[is.na(bc1$year_fc)] <- 0

## check that it was only due to NAs...
bc1 %>%
  filter(year_fc > 0) %>%
  group_by(tattoo) %>%
  distinct(year_fc) %>%
  count(sort = TRUE)

## yes - so can adjust
summary(bc1)
bc1 <- bc1 %>%
  group_by(tattoo) %>%
  mutate(year_fc = max(year_fc)) %>%
  ungroup() %>%
  filter(year_fc > 0) %>%
  mutate(age_fc_yrs = ifelse(age_fc == 3, 0, 1)) %>%  # age first caught is 3 which = cub so age in years is 0, age first caught is not 3 = yearling so age in years is 1 
  mutate(age = capture_yr - year_fc + age_fc_yrs)

summary(bc1)

bc1 <- arrange(bc1, date)

## create a trap month variable
bc1 <- bc1 %>%
  mutate(trap_month = as.numeric(month(date))) %>%
  mutate(capture_yr = as.numeric(year(date)))

## split trap month into 4 different 'seasons', year split into 4 quarters
bc1$trap_season <- NA
for (i in 1:nrow(bc1)){
  if(bc1$trap_month[i] == 1 | bc1$trap_month[i] == 2 | bc1$trap_month[i] == 3){
    bc1$trap_season[i] <- 1
  } else if (bc1$trap_month[i] == 4 | bc1$trap_month[i] == 5 | bc1$trap_month[i] == 6){
    bc1$trap_season[i] <- 2
  } else if (bc1$trap_month[i] == 7 | bc1$trap_month[i] == 8 | bc1$trap_month[i] == 9){
    bc1$trap_season[i] <- 3
  } else if (bc1$trap_month[i] == 10 | bc1$trap_month[i] == 11 | bc1$trap_month[i] == 12){
    bc1$trap_season[i] <- 4
  }
}

## create occasion variable by combining capture year and trap season
bc1 <- unite(bc1, "occasion", capture_yr, trap_season, sep = ".", remove = FALSE)

bc1 <- arrange(bc1, date) 

## check there is atrap season for every possibility
unique(bc1$occasion)

## missing 1977.4, 1979.1, 1981.1, 

## create occ
bc1$occ <- bc1 %>%
  group_by(occasion) %>%
  group_indices()

## make adjustments for missing dates
bc1$occ[bc1$occ > 3] <- bc1$occ[bc1$occ > 3] + 1
bc1$occ[bc1$occ > 7] <- bc1$occ[bc1$occ > 7] + 1
bc1$occ[bc1$occ > 15] <- bc1$occ[bc1$occ > 15] + 1

## add 5 on to all occ to allow births to take place before the study began.
bc1$occ <- bc1$occ + 5
max(bc1$occ)

## Sort out different entries of No and NO for pm
bc1$pm[bc1$pm == "no"] <- "No"
bc1$pm <- as.factor(bc1$pm) 
summary(bc1$pm)

## Sort out missing Sex
bc1 %>%
  group_by(tattoo) %>%
  distinct(sex) %>%
  count(sort = TRUE)

bc1$sex <- as.numeric(bc1$sex)
bc1$sex[is.na(bc1$sex)] <- 0

## check that each individual only has one defined sex and NAs not mixed sex
bc1 %>%
  group_by(tattoo) %>%
  filter(sex > 0) %>%
  distinct(sex) %>%
  count(sort  = TRUE)

bc1 <- bc1 %>%
  group_by(tattoo) %>%
  mutate(sex = max(sex))
  
bc1$sex <- as.factor(bc1$sex)
levels(bc1$sex) <- c("FEMALE", "MALE", "UNKNOWN")
summary(bc1$sex)

saveRDS(bc1, "BadgersSG2024.rds")
