
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



## ===== FUNCTIONS =============================================================



## ===== LOAD DATA =============================================================

# Prop of mother at US level
df.us <- readRDS(
    here(
        "data", 
        "df_us.rds"
    )
) 

# Prop of mother at RACE level
df.race <- readRDS(
    here(
        "data",
        "df_race.rds"
    )
)



## ===== VISUALIZATION =========================================================

# Prop of mother over ages by race/ethnicity and survey
df.race |> 
    ggplot(aes(x = age,
               y = p,
               ymin = l95,
               ymax = u95,
               group = year,
               col = year)) +
    facet_grid(survey ~ race_eth) +
    geom_line() +
    geom_point() +
    geom_pointrange() +
    theme_bw()

# Check if weighted prop. agrees with counts in CPS
df.race |> 
    filter(
        survey == "CPS"
        ) |> 
    ggplot(aes(x = age,
               y = p)) +
    facet_grid(race_eth ~ year) +
    geom_line(aes(y = p, col = "Prop."),
              alpha = .2) +
    geom_point(aes(y = p, col = "Prop."),
               alpha = .2) +
    geom_line(aes(y = y / n, col = "Count"),
              alpha = .2) +
    geom_point(aes(y = y / n, col = "Count"),
               alpha = .2) +
    theme_bw() +
    scale_color_manual(values = c("Prop." = "red4",
                                  "Count" = "skyblue3"))
# The two corresponds!


# Compare data from NSFG and CPS
df.race |> 
    ggplot(aes(x = age,
               y = y / n,
               group = survey,
               col = survey)) +
    facet_grid(race_eth ~ year) +
    geom_line() +
    geom_point() +
    theme_bw() +
    labs(
        y = "Prop. Mother",
        x = "Age"
    )






