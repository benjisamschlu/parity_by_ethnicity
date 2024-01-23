
##==============================================================================
##
##  Project title: Modeling of childlessness
##
##  This code does some exploratory data analysis
##
## 
##  Author: Benjamin Schlüter
##  Date: December 2023
##==============================================================================
##
##  Notes
## ------
## 1. 
## 
##==============================================================================



## LOAD PACKAGES ===============================================================

# Install/load packages
packages <- c("tidyverse", "here", "ggplot2")
for(p in packages){
    if(!require(p,character.only = TRUE)) install.packages(p)
    library(p,character.only = TRUE)
}



## FUNCTIONS ===================================================================



## LOAD DATA ===================================================================

# Childless corrected with hh info
# for the years 2008 and 2010
df <- readRDS(
    here(
        "data", 
        "df_childness.rds"
    )
) 

# Childlessness from NSFG data
df.p.nsfg <- readRDS(
    here(
        "data",
        "df_childness_us_nsfg.rds"
    )
)

df.p.race.nsfg <- readRDS(
    here(
        "data",
        "df_childness_race_nsfg.rds"
    )
)



## VISUALIZATIONS ==============================================================

# Proportion childless over the ages at the national level
df |> 
    summarise(
        # National level by year and 5y age group
        .by = c(year, age),
        
        across(y:n, ~ sum(.x)),
        # Compute prop childless
        p = (y / n) * 100
    ) |>
    ggplot(aes(x = age,
               y = p)) +
    facet_wrap(~ year) +
    geom_point() +
    geom_line() +
    theme_bw()

# Proportion childless by race over the ages at the national level
df |> 
    summarise(
        # National level by year and 5y age group
        .by = c(year, age, race_eth),
        
        across(y:n, ~ sum(.x)),
        # Compute prop childless
        p = (y / n) * 100
    ) |>
    ggplot(aes(x = age,
               y = p,
               group = race_eth,
               col = race_eth)) +
    facet_wrap(~ year) +
    geom_point() +
    geom_line() +
    theme_bw()

# Proportion childless by states over the ages at the national level
df |> 
    summarise(
        # National level by year and 5y age group
        .by = c(year, age, name),
        
        across(y:n, ~ sum(.x)),
        # Compute prop childless
        p = (y / n) * 100
    ) |>
    ggplot(aes(x = age,
               y = p,
               group = name,
               col = name)) +
    facet_wrap(~ year) +
    geom_point() +
    geom_line() +
    theme_bw()

# Comparison of childlessness 
# from CPS and NSFG by race 
bind_rows(
    df |> 
        summarise(
            # National level by year and 5y age group
            .by = c(year, age, race_eth),
            
            across(y:n, ~ sum(.x))
        ) |> 
        filter(year == 2018) |> 
        mutate(
            p = y / n,
            p = 1 - p,
            data = "CPS"
        ) |> 
        dplyr::select(age, race_eth, p, data),
    
    df.p.race.nsfg |> 
        rename(age = "age_r") |> 
        mutate(
            p = 1 - p,
            data = "NSFG"
        )
) |> 
    ggplot(aes(x = age,
               y = p,
               group = data,
               col = data)) +
    facet_wrap(~ race_eth) +
    geom_point() +
    geom_line() +
    theme_bw()
    

