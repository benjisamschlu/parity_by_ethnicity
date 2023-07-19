
rm(list=ls())

library(tidyverse)
library(here)
library(HMDHFDplus)

source(here("code", "secrets.R"))



## Load data -------------------------------------------------------------------

## Adapt code below to work on frever
## Careful about age of mother
# https://cps.ipums.org/cps-action/variables/FREVER#universe_section

## read_csv allows to open compressed .gz 
data <- read_csv(here("data_private", "cps_00005.csv.gz"))



exppa <- readHFDweb(CNTRY = "USA",
                    item = "exposRRpa",
                    username = usr,
                    password = pwd)



## Tidy data ---------------------------------------------------------------

## HFD
df.p.hfd <- exppa |> 
    mutate(Ex = rowSums(across(c(E0x, E1x, E2x, E3x, E4px))),
           across(E0x:E4px, ~.x/Ex))

## Lowers col names
names(data) <- tolower(names(data))



## Number of minor per person --------------------------------------------------

## Consider individuals aged >= 10 yo (to be filtered later)

## What about weigths?

## Minimum nber of col to optimize code
df <- data |> 
    ## Only keep year where frever (nber of children
    ## ever born) was recorded
    filter(!is.na(frever)) |> 
    dplyr::select(
        year, month, serial,pernum, ## identification
        momloc, momloc2, poploc, poploc2, frever,## parenthood
        wtfinl, ## weights
        age, sex, race, hispan) |> ## demographic
    mutate(nminor = 0)

## Keep only women since 1980, 
## not having 999 for frever
df.p <- df |> 
    filter(year >= 1990,
           sex ==2,
           frever != 999)

## Parity per age and sex
df.p <- df.p |> 
    mutate(
        ## Rescale weights
        w = wtfinl/1000,
        ## Create 4+ nber of minor
        nchild = ifelse(frever >= 4, 4, frever)
    ) |> 
    group_by(age, sex, year,
             .drop = F) |> 
    summarise(p0 = sum(w[nchild == 0])/sum(w),
              p1 = sum(w[nchild == 1])/sum(w),
              p2 = sum(w[nchild == 2])/sum(w),
              p3 = sum(w[nchild == 3])/sum(w),
              p4 = sum(w[nchild == 4])/sum(w)) 

## Combine and merge CPS & HFD data
combined.df <- bind_rows(
    df.p |> 
        pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
        mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
        filter(sex == "Female") |> 
        mutate(source = "CPS"),
    df.p.hfd |> 
        dplyr::select(-c(OpenInterval, Ex)) |> 
        `colnames<-`(c("year", "age", "p0", "p1", "p2", "p3", "p4")) |> 
        pivot_longer(p0:p4, names_to = "parity", values_to = "p") |>
        mutate(sex = "Female",
               source = "HFD")
    
) 



## Visu ------------------------------------------------------------------------

## CPS Fertility data 2019
df.p |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
    mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
    filter(year == 2020) |> 
    ggplot(aes(x = age, y = p, 
               group = parity,
               col = parity)) + 
    geom_line(linewidth = 1) +
    theme_bw()


## HFD vs CPS
combined.df |> 
    filter(year == 2020) |> 
    ggplot(aes(x = age, y = p, group = interaction(parity, source))) +
    geom_line(aes(col = parity, linetype = source),
              linewidth = 0.8) +
    theme_bw()
## Pros:
# Both gender
# Important sample size

## Cons:
# Survey is HH-based
# As parent age, less child in HH
# Divorce and (non) shared custody
# Constraining child to be <18 yo


