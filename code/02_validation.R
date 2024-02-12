
##==============================================================================
##
##  Project title: Modeling of childlessness
##
##  This code compares proportion of childless women after
##  correction of 2008 and 2010 data with the ones
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

# Motherhood at US level
df.us <- readRDS(
    here(
        "data", 
        "df_us.rds"
        )
    ) 



## VALIDATION ==================================================================

# Fig 1
cols <- c("30" = "blue",
          "35" = "red",
          "40" = "green")

df.us |> 
    filter(
        survey == "CPS"
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
        .by = c(year, age5),
        
        across(y:n, ~ sum(.x)),
        # Compute prop childless
        p = (1 - (y / n)) * 100
    ) |> 
    filter(
        # Focus on years plotted on report's Fig 1
        year >= 2000,
        year <= 2014,
        # Focus on 3 age groups
        age5 %in% seq(30, 40, 5)
    ) |> 
    ggplot(aes(x = year, y = p, 
               group = age5, col = age5)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    scale_y_continuous(breaks = seq(14, 32, 2),
                       limits = c(14, 32)) +
    scale_color_manual(values = cols)
    

