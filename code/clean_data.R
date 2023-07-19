
rm(list=ls())

library(tidyverse)
library(here)
library(HMDHFDplus)

source(here("code", "secrets.R"))



## Load data -------------------------------------------------------------------

## read_csv allows to open compressed .gz 
data <- read_csv(here("data_private", "cps_00004.csv.gz"))

exppa <- readHFDweb(CNTRY = "USA",
                    item = "exposRRpa",
                    username = usr,
                    password = pwd)



## Tidy data ---------------------------------------------------------------

## HFD
df.p.hfd <- exppa |> 
    filter(Year == 2019) |> 
    mutate(Ex = rowSums(across(c(E0x, E1x, E2x, E3x, E4px))),
           across(E0x:E4px, ~.x/Ex))

## Lowers col names
names(data) <- tolower(names(data))



## Number of minor per person --------------------------------------------------

## Consider individuals aged >= 10 yo (to be filtered later)

## What about weigths?

## Minimum nber of col to optimize code
df <- data |>  
    dplyr::select(
        year, month, serial,pernum, ## identification
        momloc, momloc2, poploc, poploc2, ## parenthood
        wtfinl, ## weights
        age, sex, race, hispan) |> ## demographic
    mutate(nminor = 0)

## Nber of minor per individual
## Function
count.minor <- function(d) {
    
    ## Who is parent of minor
    find <- (d$momloc != 0 | d$momloc2 != 0 | d$poploc != 0 | d$poploc2 != 0) & d$age < 18
    ## How many time are they listed parent
    n <- table(c(d$momloc[find], d$momloc2[find], d$poploc[find], d$poploc2[find]))
    ## Do not account for 0 (construction of "find" will always contain a 0)
    n <- tail(n, -1)
    id <- as.numeric(names(n))
    ## Store
    d$nminor[d$pernum %in% id] <- n
    
    return(d)
}

## Apply function and get counts
df.minor <- df |> 
    group_by(year, month, serial) |> 
    group_modify(~ count.minor(.x))

## Check warning was printed after function
# saveRDS(df.minor,
#         here("data", "df_minor.rda"))
# df.minor <- readRDS(here("data", "df_minor.rda"))

## Parity per age and sex
df.p <- df.minor |> 
    ## Consider ind from 10 yo 
    ## (interested in parents)
    filter(age >= 10) |> 
    mutate(
        ## Rescale weights
        w = wtfinl/1000,
        ## Create 4+ nber of minor
        nminor = ifelse(nminor >= 4, 4, nminor)
    ) |> 
    group_by(age, sex,
             .drop = F) |> 
    summarise(p0 = sum(w[nminor == 0])/sum(w),
              p1 = sum(w[nminor == 1])/sum(w),
              p2 = sum(w[nminor == 2])/sum(w),
              p3 = sum(w[nminor == 3])/sum(w),
              p4 = sum(w[nminor == 4])/sum(w)) 

## Combine and merge CPS & HFD data
combined.df <- bind_rows(
    df.p |> 
        pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
        mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
        filter(sex == "Female") |> 
        mutate(source = "CPS"),
    df.p.hfd |> 
        dplyr::select(-c(Year, OpenInterval, Ex)) |> 
        `colnames<-`(c("age", "p0", "p1", "p2", "p3", "p4")) |> 
        pivot_longer(p0:p4, names_to = "parity", values_to = "p") |>
        mutate(sex = "Female",
               source = "HFD")
    
)



## Visu ------------------------------------------------------------------------

## CPS Fertility data 2019
df.p |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
    mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
    ggplot(aes(x = age, y = p, 
               group = interaction(parity, sex),
               col = parity,
               linetype = sex)) + 
    geom_line(linewidth = 1) +
    theme_bw()


## HFD vs CPS
combined.df |> 
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


