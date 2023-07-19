
rm(list=ls())

library(tidyverse)
library(here)
library(HMDHFDplus)

source(here("code", "secrets.R"))



## Load data -------------------------------------------------------------------

exppa <- readHFDweb(CNTRY = "USA",
                    item = "exposRRpa",
                    username = usr,
                    password = pwd)
## Mortality data
lt.US <- readRDS("../parental_loss_by_race/data/lt_US.rda") |> 
    ## Compute mx by causes
    mutate(mx_other = (n_deaths - (n_drug + n_firearm))/pop,
           mx_drug = n_drug/pop,
           mx_firearm = n_firearm/pop)



## Tidy data -------------------------------------------------------------------

## HFD
df.p.hfd <- exppa |> 
    filter(Year %in% 1999:2020) |> 
    mutate(Ex = rowSums(across(c(E0x, E1x, E2x, E3x, E4px))),
           across(E0x:E4px, ~.x/Ex)) |> 
    dplyr::select(-c(OpenInterval, Ex)) |> 
    `colnames<-`(c("year", "age", "p0", "p1", "p2", "p3", "p4")) 
## Extrapol oldest age (>55)
extrapol <- function(d) {
    
    A <- max(d$age)
    to.extrapol <- d[d$age == A,] |> 
        dplyr::select(-age)
    temp <- rbind(to.extrapol, to.extrapol[rep(1, 84-A), ]) |> 
        mutate(age = (A+1):85)
    out <- bind_rows(d,
                     temp)
    return(out)
}
## Extend years in HFD
df.p.hfd <- df.p.hfd |> 
    group_by(year) |> 
    group_modify(~ extrapol(.x))

## Combine data frame
combined_df <- lt.US |> 
    filter(year >= 1999,
           race_eth == "total",
           age >= 12) |> 
    dplyr::select(year, sex, age, pop, mx_drug, mx_firearm) |> 
    ## Use female distribution for both sexes
    left_join(df.p.hfd,
              by = c("year", "age")) |> 
    mutate(
        ## parity multiplier
        mult = 0*p0 + 1*p1 + 2*p2 + 3*p3 + 4*p4,
        ## calculation of children affected
        nchild_drug = pop * mx_drug * mult,
        nchild_firearm = pop * mx_firearm * mult,
    )



## Visu ------------------------------------------------------------------------

combined_df |>
    group_by(year, sex) |> 
    summarise(drug = sum(nchild_drug),
              firearm = sum(nchild_firearm)) |> 
    pivot_longer(drug:firearm, 
                 names_to = "cause", values_to ="N") |> 
    ggplot(aes(x = year, y = N/1000, group = cause)) + 
    facet_grid(cause ~ sex) +
    geom_line() +
    theme_bw() +
    labs(y = "Number of children impacted by 
        \nparental loss (in thousands)")
## Significantly higher 
## * MPM stops at age 18 yo of children (indirectly constraints age of parents)
## * Fathers have their child later in MPM (shift), reduce even more the windows
## while currently, male mortality at young age might have more weight
## which might explain the bigger difference for firearms 
## -> Sensitivity check with male = female fertility

## Cumulative number    
combined_df |>
    group_by(year, sex) |> 
    summarise(drug = sum(nchild_drug),
              firearm = sum(nchild_firearm)) |> 
    ungroup() |> 
    group_by(sex) |> 
    mutate(drug = cumsum(drug),
           firearm = cumsum(firearm)) |> 
    pivot_longer(drug:firearm, 
                 names_to = "cause", values_to ="cumN") |> 
    ggplot(aes(x = year, y = cumN/1000, group = cause)) + 
    facet_grid(cause ~ sex) +
    geom_line() +
    theme_bw() +
    labs(y = "Cum. number of children impacted by 
        \nparental loss (in thousands)")

