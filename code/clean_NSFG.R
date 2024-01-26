
##==============================================================================
##
##  Project title: Modeling of childlessness
##
##  This code loads, cleans up, and tidy the NSFG data used to estimate
##  the proportion of childless women.
##
## 
##  Author: Benjamin Schlüter
##  Date: December 2023
##==============================================================================
##
##  Notes
## ------
## 1. The NSFG is designed to be nationally representative of women and men 
## ages 15 to 49 years in the household-based population of the United States. 
## Prior to 2002, the NSFG included only women ages 15 to 44 years. 
## Men ages 15 to 44 years were included in the survey in 2002; 
## the age range was expanded to 15 to 49 years in 2015.
## 2. For 2015-2017, interviews were conducted with 5,554 women and 4,540 men, 
## for a total sample size of 10,094. The response rate was 66.7% for female 
## respondents and 63.6% for male respondents. The overall response rate was 
## 65.3%.
##
##==============================================================================

rm(list = ls())


## LOAD PACKAGES ===============================================================

# Install/load packages
packages <- c("tidyverse", "here", "survey",
              "SAScii", "readr")
for(p in packages){
    if(!require(p,character.only = TRUE)) install.packages(p)
    library(p,character.only = TRUE)
}



## FUNCTIONS ===================================================================



## LOAD DATA ===================================================================

# Code taken from
# https://asdfree.com/national-survey-of-family-growth-nsfg.html

# Need to add period 2006_2010
dwld.period <- c("2006_2010", "2011_2013", "2013_2015", "2015_2017", "2017_2019")

list.nsfg.df <- lapply(dwld.period,
                       function(p) {
                           
                           # File name changed after 2006_2010
                           if (p != "2006_2010") {
                               file_name <- "_FemRespData.dat"
                           } else { 
                               file_name <- "_FemResp.dat"
                           }
                           # Url to download data
                           dat_url <- paste0(
                               "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NSFG/",
                               p,
                               file_name
                           )
                           
                           sas_url <-
                               file.path(dirname(dat_url) , 
                                         paste0("sas/",
                                                p,
                                                "_FemRespSetup.sas") )
                           
                           sas_positions <-
                               parse.SAScii( sas_url )
                           
                           sas_positions[ , 'varname' ] <-
                               tolower( sas_positions[ , 'varname' ] )
                           
                           sas_positions[ , 'column_types' ] <-
                               ifelse( sas_positions[ , 'char' ] , "c" , "d" )
                           
                           # Load data
                           nsfg_tbl <-
                               read_fwf(
                                   dat_url ,
                                   fwf_widths( 
                                       abs( sas_positions[ , 'width' ] ) , 
                                       col_names = sas_positions[ , 'varname' ] 
                                   ) ,
                                   col_types = paste0( sas_positions[ , 'column_types' ] , collapse = "" ) ,
                                   na = c( "" , "." )
                               )
                           
                           # Weights have different names in 2006_2010
                           if (p != "2006_2010") {
                               wgt_var <- paste0("wgt", p)
                           } else { 
                               wgt_var <- "wgtq1q16"
                           }
                           # Select variables of interest
                           nsfg_df <- data.frame( nsfg_tbl ) |> 
                               dplyr::select(
                                   age_r, hisprace2, parity, secu, sest, !!sym(wgt_var)
                               ) |> 
                               rename(
                                   # Harmonize weight var name to row bind later
                                   "wgt" = !!sym(wgt_var)
                               ) |> 
                               mutate(
                                   # Var for period
                                   period = p
                               )
                           }
                       )
nsfg.df <- do.call("rbind", list.nsfg.df)

# Store data
saveRDS(nsfg.df, 
        here(
            "data_private",
            "2006_2019_NSFG_FemRespData.rds"
        ), 
        compress = FALSE )



## TIDY DATA ===================================================================

## Filter data of interest
df <- nsfg.df |> 
    mutate(
        ## Create 5 variables for parity (to be used with svyby())
        p0 = ifelse(parity == 0, 1, 0),
        p1 = ifelse(parity == 1, 1, 0),
        p2 = ifelse(parity == 2, 1, 0),
        p3 = ifelse(parity == 3, 1, 0),
        p4 = ifelse(parity >= 4, 1, 0),
        ## Name abbridged race
        race_eth = case_when(
            hisprace2 == 1 ~ "Hispanic",
            hisprace2 == 2 ~ "NH-White",
            hisprace2 == 3 ~ "NH-Black",
            hisprace2 == 4 ~ "NH-Other"
        )
    )
## Complex survey desing
nsfg_design <- 
    svydesign( 
        id = ~ secu , 
        strata = ~ sest , 
        data = df , 
        weights = ~ wgt, 
        nest = TRUE 
    )

## Check similarity with workshop slides (sl. 24)
# https://www.bgsu.edu/content/dam/BGSU/college-of-arts-and-sciences/center-for-family-and-demographic-research/documents/Workshops/2022-nsfg-workshop.pdf
# svyby(~ parity, ~ race_eth + period, nsfg_design, svymean)
## Weighted parity by age

df.p.nsfg <- tibble(
    svyby( ~ p0 + p1 + p2 + p3 + p4 , ~ age_r + period, nsfg_design , svymean )
) |> 
    dplyr::select(
        period, age_r, p0, p1, p2, p3, p4
        ) |> 
    pivot_longer(
        p0:p4, 
        names_to = "parity", 
        values_to = "p"
        )

## Weighted parity by age and race
df.p.race.nsfg <- tibble(
    svyby( ~ p0 + p1 + p2 + p3 + p4 , ~ age_r + race_eth + period, nsfg_design , svymean )
) |> 
    dplyr::select(
        period, age_r, race_eth, p0, p1, p2, p3, p4
        ) |> 
    pivot_longer(
        p0:p4, 
        names_to = "parity",
        values_to = "p"
        )



## VISUALIZATION ===============================================================

## NSFG by race and % parent
df.p.race.nsfg |> 
    filter(
        race_eth != "NH-Other",
        parity == "p0"
    ) |> 
    mutate(
        parent = 1 - p
    ) |> 
    ggplot(aes(x = age_r, 
               y = parent,
               group = period,
               col = period)) + 
    facet_wrap( ~ race_eth) +
    geom_line(linewidth = 1) +
    geom_point() +
    theme_bw()

                           
