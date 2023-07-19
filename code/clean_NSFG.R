
rm(list=ls())

library(tidyverse)
library(here)
library(HMDHFDplus)
library(survey)
source(here("code", "secrets.R"))


## install_github( "ajdamico/lodown" , dependencies = TRUE )
library(lodown)



## Load data -------------------------------------------------------------------

## NSFG

# # Examine all available NSFG microdata files
# nsfg_cat <-
#     get_catalog( "nsfg" ,
#                  output_dir = here("data_raw") )
# 
# # 2017-2019 only
# nsfg_cat <- subset( nsfg_cat , grepl( "2017_2019" , full_url ) )
# # Download the microdata to your local computer
# nsfg_cat <- lodown( "nsfg" , nsfg_cat )
## Load data
nsfg_df <- readRDS( here("data_raw" , "2017_2019_FemRespData.rds" ) )

## HFD

exppa <- readHFDweb(CNTRY = "USA",
                    item = "exposRRpa",
                    username = usr,
                    password = pwd)

## Tidy data -------------------------------------------------------------------

## Filter data of interest
df <- nsfg_df |> 
    dplyr::select(age_r, hisprace2, parity, secu, sest, wgt2017_2019) |> 
    mutate(
        ## Create 5 variables for parity (to be used with svyby())
        p0 = ifelse(parity == 0, 1, 0),
        p1 = ifelse(parity == 1, 1, 0),
        p2 = ifelse(parity == 2, 1, 0),
        p3 = ifelse(parity == 3, 1, 0),
        p4 = ifelse(parity >= 4, 1, 0),
        ## Name abbridged race
        abr_race = case_when(
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
        weights = ~ wgt2017_2019 , 
        nest = TRUE 
    )

## Check similarity with workshop slides
# https://www.bgsu.edu/content/dam/BGSU/college-of-arts-and-sciences/center-for-family-and-demographic-research/documents/Workshops/2022-nsfg-workshop.pdf
# svyby(~ parity, ~ abr_race, nsfg_design, svymean)
## Weighted parity by age

df.p.nsfg <- tibble(
    svyby( ~ p0 + p1 + p2 + p3 + p4 , ~ age_r , nsfg_design , svymean )
) |> 
    dplyr::select(age_r, p0, p1, p2, p3, p4) |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p")

## Weighted parity by age and race
df.p.race.nsfg <- tibble(
    svyby( ~ p0 + p1 + p2 + p3 + p4 , ~ age_r + abr_race , nsfg_design , svymean )
) |> 
    dplyr::select(age_r, abr_race, p0, p1, p2, p3, p4) |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p")

## HFD
df.p.hfd <- exppa |> 
    ## Sum over the same period as NSFG
    filter(Year %in% 2017:2019) |> 
    group_by(Age) |> 
    summarise(across(E0x:E4px, ~sum(.x))) |> 
    ungroup() |> 
    mutate(Ex = rowSums(across(c(E0x, E1x, E2x, E3x, E4px))),
           across(E0x:E4px, ~.x/Ex)) |> 
    dplyr::select(-Ex) |> 
    `colnames<-`(c("age", "p0", "p1", "p2", "p3", "p4")) |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p")

## Combine and merge NSFG & HFD 
combined.df <- bind_rows(
    df.p.nsfg |> 
        rename("age" = age_r) |> 
        mutate(source = "NSFG"),
    df.p.hfd |> 
        mutate(source = "HFD")
    
)



## Visu ------------------------------------------------------------------------

## NSFG by age
df.p.nsfg |> 
    ggplot(aes(x = age_r, y = p, 
               col = parity)) + 
    geom_line(linewidth = 1) +
    theme_bw()

## NSFG vs HFD
combined.df |> 
    ggplot(aes(x = age, y = p, group = interaction(parity, source))) +
    geom_line(aes(col = parity, linetype = source),
              linewidth = 0.8) +
    theme_bw()

## NSFG by race
df.p.race.nsfg |> 
    ggplot(aes(x = age_r, y = p, 
               col = parity)) + 
    facet_wrap( ~ abr_race) +
    geom_line(linewidth = 1) +
    theme_bw()