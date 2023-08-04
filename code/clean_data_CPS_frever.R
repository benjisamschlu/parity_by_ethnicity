
rm(list=ls())

library(tidyverse)
library(here)
library(HMDHFDplus)

source(here("code", "secrets.R"))



## Load data -------------------------------------------------------------------


## Careful about age of mother: not the same boundaries over the years
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

## Tidy
df <- data |> 
    ## Only keep year where frever (nber of children
    ## ever born) was recorded.
    ## Recorded only for females.
    filter(!is.na(frever)) |> 
    dplyr::select(
        year, month, serial,pernum, ## identification
        momloc, momloc2, poploc, poploc2, frever,## parenthood
        wtfinl, ## weights,
        statefip, ## states
        age, sex, race, hispan) |> ## demographic
    mutate(
        ## Create abridged race (!! NEED TO BE CHECKED !!)
        race_eth = case_when(race == 100 & hispan == 0 ~ "NH-White",
                             race == 200 & hispan == 0 ~ "NH-Black",
                             !(race %in% c(100, 200)) & hispan == 0 ~ "NH-Other",
                             hispan != 0 ~ "Hispanic")
    ) |> 
    ## Keep only women since 1990, 
    ## aged 15-45 not having 999
    ## only in latest CPS that itwees
    ## are until 50 yo
    ## !! Men don't have info on frever !!
    filter(year >= 1990,
           age >= 15,
           age < 45,
           frever != 999)

## Childlessness
df.chldness <- df |>  
    mutate(
        ## Rescale weights
        w = wtfinl/1000,
        ## Create 4+ nber of minor
        childless = ifelse(frever == 0, 1, 0)
    ) |> 
    group_by(age, year, 
             .drop = F) |> 
    summarise(y = sum(w[childless == 1]),
              n = sum(w)) 

## Store data for modeling
saveRDS(df.chldness,
        here("data", "df_chldness.rda"))

## Parity per age 
df.p <- df |> 
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

## Store data for modeling
saveRDS(df.p,
        here("data", "df_p.rda"))

## Parity per age, and race
df.p.race <- df |> 
    mutate(
        ## Rescale weights
        w = wtfinl/1000,
        ## Create 4+ nber of minor
        nchild = ifelse(frever >= 4, 4, frever)
    ) |> 
    group_by(age, sex, year, race_eth, 
             .drop = F) |> 
    summarise(p0 = sum(w[nchild == 0])/sum(w),
              p1 = sum(w[nchild == 1])/sum(w),
              p2 = sum(w[nchild == 2])/sum(w),
              p3 = sum(w[nchild == 3])/sum(w),
              p4 = sum(w[nchild == 4])/sum(w)) 

## Store data for modeling
saveRDS(df.p.race,
        here("data", "df_p_race.rda"))

## Parity per age, and states
df.p.states <- df |> 
    mutate(
        ## Rescale weights
        w = wtfinl/1000,
        ## Create 4+ nber of minor
        nchild = ifelse(frever >= 4, 4, frever)
    ) |> 
    group_by(age, sex, year, statefip, 
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
        mutate(source = "CPS"),
    df.p.hfd |> 
        dplyr::select(-c(OpenInterval, Ex)) |> 
        `colnames<-`(c("year", "age", "p0", "p1", "p2", "p3", "p4")) |> 
        pivot_longer(p0:p4, names_to = "parity", values_to = "p") |>
        mutate(sex = "Female",
               source = "HFD")
    
) 



## Visu ------------------------------------------------------------------------

## CPS Fertility data 2020
df.p |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
    mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
    filter(year == 2020) |> 
    ggplot(aes(x = age, y = p, 
               group = parity,
               col = parity)) + 
    facet_wrap(~ sex) +
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
# Important sample size

## CPS Fertility data 2020 by race
df.p.race |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
    mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
    filter(year == 2020) |> 
    ggplot(aes(x = age, y = p, 
               group = parity,
               col = parity)) + 
    facet_wrap(~ race_eth) +
    geom_line(linewidth = 1) +
    theme_bw()

## CPS Fertility data 2020 by states
df.p.states |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
    mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
    filter(year == 2020) |> 
    ggplot(aes(x = age, y = p, 
               group = parity,
               col = parity)) + 
    facet_wrap(~ statefip) +
    geom_line(linewidth = 1) +
    theme_bw()



## Parametric function to model p ----------------------------------------------

## Logistic 
## Janoscheck is easy to fit and has flex inflexion point
## Richards (gen logistic) is hard to fit but super flex

## Focusing on parity 0
logistic_fct <- function(x, x0=30, U=1, k) {
    
    out <- U / (1 + exp( k*(x - x0) ))
    return(out)
}

gen_log_fct <- function(x, x0=30, L , U=1, k) {
    
    out <- L + ((U - L) / (1 + exp( k*(x - x0) )))
    return(out)
}

janoscheck_fct <- function(x, x0=30, L, U=1, k) {
    
    out <- U - (U - L)*exp(k*x^x0)
    return(out)
}

richards_fct <- function(x, x0=30, L, U=1, d, k) {
    
    out <- L + (U - L)/(1 + (d - 1)*exp(k * (x - x0) ))^(1 - d)
    return(out)
}

df.p |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
    mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
    filter(year == 2020,
           parity == "p0") |> 
    mutate(fit_p0 = logistic_fct(age, k = 0.24),
           gen_fit_p0 = gen_log_fct(age, x0= 27, L=0.15, k = 0.24),
           jano_fit_p0 = janoscheck_fct(age, L= 0.2, k = 0.24),
           rich_fit_p0 = richards_fct(age, L=0.2, k=0.24, d=2)) |> 
    ggplot(aes(x = age, y = p)) + 
    geom_point() +
    geom_line(aes(y = fit_p0)) +
    geom_line(aes(y = gen_fit_p0), col = "red3") +
    # geom_line(aes(y = rich_fit_p0), col = "pink") +
    # geom_line(aes(y = jano_fit_p0), col = "green3") +
    theme_bw()



## Trying to use optim() -------------------------------------------------------------------------

df.rss <- df.p |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
    mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
    filter(year == 2020,
           parity == "p0") |> 
    ungroup() |> 
    dplyr::select(age, p)

min.RSS <- function(data, par) {
    with(data, sum( ((par[1] + ((1-par[1])/(1+exp(par[2]*(age-par[3]))))) -  p)^2 ))
}

pars <- optim(par = c(0.15, 0.24, 27), 
              fn = min.RSS, 
              data = df.rss,
              method = "L-BFGS-B",
              lower = c(0, 0, 0),
              upper = c(0.5, 0.5, 40))$par


df.p |> 
    pivot_longer(p0:p4, names_to = "parity", values_to = "p") |> 
    mutate(sex = ifelse(sex == 2, "Female", "Male")) |> 
    filter(year == 2020,
           parity == "p0") |> 
    mutate(gen_fit_p0 = gen_log_fct(age, x0= pars[3], L=pars[1], k = pars[2])) |> 
    ggplot(aes(x = age, y = p)) + 
    geom_point() +
    geom_line(aes(y = gen_fit_p0), col = "red3") +
    # geom_line(aes(y = rich_fit_p0), col = "pink") +
    # geom_line(aes(y = jano_fit_p0), col = "green3") +
    theme_bw()
