
##==============================================================================
##
##  Project title: Modeling of childlessness
##
##  This code compares proportion of childless women with the ones
##  obtained from the US Census Bureau's report.
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

# Careful: no observation for all years (~ every 2 y)

# Childless by race and state
df <- readRDS(
    here(
        "data", 
        "df_childness.rds"
        )
    ) 
# Childless corrected with hh info
# for the years 2008 and 2010
df_cor <- readRDS(
    here(
        "data", 
        "df_childness_cor_2008_10.rds"
    )
) 



## VALIDATION ==================================================================

# Reproduce plot in US Census Bureau report
cols <- c("30" = "blue",
          "35" = "red",
          "40" = "green")
# Fig 1
df |> 
    mutate(
        data = "initial"
    ) |> 
    bind_rows(
        df_cor |> 
            mutate(
                data = "corrected"
            )
    ) |> 
    mutate(
        # Create 5y age groups
        age5 = cut(age,
                   breaks = seq(15, 45, 5),
                   labels = seq(15, 40, 5),
                   right = FALSE)
    ) |> 
    summarise(
        # National level by year and 5y age group
        .by = c(year, age5, data),
        
        across(y:n, ~ sum(.x)),
        # Compute prop childless
        p = (y / n) * 100
    ) |> 
    filter(
        # Focus on years plotted on report's Fig 1
        year >= 2000,
        year <= 2014,
        # Focus on 3 age groups
        age5 %in% seq(30, 40, 5),
        # Check data without correction
        data == "initial"
    ) |> 
    ggplot(aes(x = year, y = p, 
               group = age5, col = age5)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_y_continuous(breaks = seq(14, 32, 2),
                       limits = c(14, 32)) +
    scale_color_manual(values = cols)
# Same values

lines <- c("initial" = "solid",
           "corrected" = "dashed")

# Fig 2
df |> 
    mutate(
        data = "initial"
    ) |> 
    bind_rows(
        df_cor |> 
            mutate(
                data = "corrected"
            )
    ) |> 
    mutate(
        # Create 5y age groups
        age5 = cut(age,
                   breaks = seq(15, 45, 5),
                   labels = seq(15, 40, 5),
                   right = FALSE)
    ) |> 
    summarise(
        # National level by year and 5y age group
        .by = c(year, age5, data),
        
        across(y:n, ~ sum(.x)),
        # Compute prop childless
        p = (y / n) * 100
    ) |> 
    filter(
        # Focus on years plotted on report's Fig 1
        year >= 2006,
        year <= 2014,
        # Focus on 3 age groups
        age5 %in% seq(30, 40, 5)
    ) |> 
    ggplot(aes(x = year, 
               y = p, 
               group = interaction(age5, data), 
               col = age5,
               linetype = data)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_y_continuous(breaks = seq(14, 32, 2),
                       limits = c(14, 32)) +
    scale_color_manual(values = cols) +
    scale_linetype_manual(values = lines)
# Same values except a really small diff for 35



