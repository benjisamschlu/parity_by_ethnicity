
## EDA LOGISTIC PARAMETERS -----------------------------------------------------

rm(list=ls())

library(tidyverse)
library(here)
library(rstan)
library(splines)



## Load data -------------------------------------------------------------------

## Careful: no observation for all years (~ every 2 y)

## Childless by race and state
df <- readRDS(here("data", "df_chldness.rda"))



## Key dimensions in code ------------------------------------------------------

ages <- unique(df$age)
n.ages <- length(ages)
years <- unique(df$year)
n.years <- length(years)
states <-  unique(df$name)
n.states <- length(states)
races <- unique(df$race_eth)
n.races <- length(races)



## EDA - LOGISTIC PARS FOR US ==================================================

# Childless 
df.stan.us <- df |>
    ungroup() |> 
    # create factors to keep
    # all combinations with
    # group_by
    mutate(
        age = factor(age,
                     levels = ages,
                     labels = ages),
        year = factor(year,
                      levels = years,
                      labels = years),
        name = factor(name, 
                      levels = name,
                      labels = name),
        race_eth = factor(race_eth,
                          levels = races,
                          labels = races)
    ) |> 
    summarise(
        .by = c(year, age),
        across(y:n, sum)
    ) |> 
    ungroup() |> 
    # Arrange for STAN
    arrange(year, age) 

# STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.stan.us)[1],
    Y_obs = n.years,
    A = n.ages,
    age = ages,
    G = 1,
    # Careful !!
    # Structured as vector but order super important !!!
    # Needs to coincide with logmx_ord in STAN (multiple to_vector() )
    y = round(df.stan.us$y, 0), 
    n = round(df.stan.us$n, 0)
)

# STAN set-up
niter = 2000
nwarmup = niter/2
nchains = 4
nsim = (niter-nwarmup)*nchains

options(mc.cores = parallel::detectCores()-1)

# STAN fit
pars <- c("x0", "k", "U","p")

fit = stan(here("code", "stan", "p0_logistic_pars_eda.stan"),
           pars  = pars,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data
)

# Store and tidy logistic pars 
pars.eda <- lapply(pars[pars != "p"], function(x) {
    
    # Store in arrays
    par.draw <- array(as.matrix(fit, x),
                      c(nsim, n.years),
                      dimnames = list(1:nsim,
                                      years))
    # Get CI
    par.ci <- apply(par.draw, 
                    2, 
                    quantile, probs = c(0.975, 0.5, 0.025))
    # Store in df
    out <- as.data.frame.table(par.ci) %>%
        rename("year" = Var2) %>%
        mutate(Var1 = case_when( 
            Var1 == "97.5%" ~ "upper95",
            Var1 == "50%" ~ "median",
            Var1 == "2.5%" ~ "lower95")) |> 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        mutate(year = as.character(year) %>% as.numeric,
               par = x)
})

# Combine list items into one data frame
df.pars.eda.us <- do.call("rbind", pars.eda)

saveRDS(df.pars.eda.us,
        here("data", "df_logistic_pars_eda_us.rda"))

# df.pars.eda.us <- readRDS(here("data", "df_logistic_pars_eda_us.rda"))

## Store logistic fit 
p.fit <- array(as.matrix(fit, "p"),
               c(nsim, n.ages, n.years),
               dimnames = list(1:nsim,
                               ages,
                               years
               ))
## Quantile of p
p.ci <- apply(p.fit, c(2:3), quantile, probs = c(0.975, 0.5, 0.025))

df.fit.eda.us <- as.data.frame.table(p.ci) %>%
    rename("age" = Var2,
           "year" = Var3) %>%
    mutate(Var1 = case_when( 
        Var1 == "97.5%" ~ "upper95",
        Var1 == "50%" ~ "median",
        Var1 == "2.5%" ~ "lower95")) |> 
    pivot_wider(names_from = Var1, values_from = Freq) %>% 
    mutate(age = as.character(age) %>% as.numeric,
           year = as.character(year) %>% as.numeric)

## Store estimates
saveRDS(df.fit.eda.us,
        here("data", "df_logistic_fit_eda_us.rda"))

## Load estimates
# df.p.us <- readRDS(here("data", "df_logistic_fit_eda_us.rda"))


## EDA - LOGISTIC PARS FOR RACES ===============================================

# Childless by race_eth
df.stan.race <- df |>
    ungroup() |> 
    # create factors to keep
    # all combinations with
    # group_by
    mutate(
        age = factor(age,
                     levels = ages,
                     labels = ages),
        year = factor(year,
                      levels = years,
                      labels = years),
        name = factor(name, 
                      levels = name,
                      labels = name),
        race_eth = factor(race_eth,
                          levels = races,
                          labels = races)
    ) |> 
    summarise(
        .by = c(year, age, race_eth),
        across(y:n, sum)
    ) |> 
    ungroup() |> 
    # Arrange for STAN
    arrange(race_eth, year, age) 

# STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.stan.race)[1],
    Y_obs = n.years,
    A = n.ages,
    age = ages,
    G = n.races,
    # Careful !!
    # Structured as vector but order super important !!!
    # Needs to coincide with logmx_ord in STAN (multiple to_vector() )
    y = round(df.stan.race$y, 0), 
    n = round(df.stan.race$n, 0)
)

# STAN set-up
niter = 2000
nwarmup = niter/2
nchains = 4
nsim = (niter-nwarmup)*nchains

options(mc.cores = parallel::detectCores()-1)

# STAN fit
pars <- c("x0", "k", "U","p")

fit = stan(here("code", "stan", "p0_logistic_pars_eda.stan"),
           pars  = pars,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data
)

# Store and tidy logistic pars 
pars.eda <- lapply(pars[pars != "p"], function(x) {
    
    # Store in arrays
    par.draw <- array(as.matrix(fit, x),
                      c(nsim, n.years, n.races),
                      dimnames = list(1:nsim,
                                      years,
                                      races))
    # Get CI
    par.ci <- apply(par.draw, 
                    2:3, 
                    quantile, probs = c(0.975, 0.5, 0.025))
    # Store in df
    out <- as.data.frame.table(par.ci) %>%
        rename("year" = Var2,
               "race_eth" = Var3) %>%
        mutate(Var1 = case_when( 
            Var1 == "97.5%" ~ "upper95",
            Var1 == "50%" ~ "median",
            Var1 == "2.5%" ~ "lower95")) |> 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        mutate(year = as.character(year) %>% as.numeric,
               par = x)
})

# Combine list items into one data frame
df.pars.eda.race <- do.call("rbind", pars.eda)

saveRDS(df.pars.eda.race,
        here("data", "df_logistic_pars_eda_race.rda"))

# df.pars.eda.race <- readRDS(here("data", "df_logistic_pars_eda_race.rda"))

## Store logistic fit 
p.fit <- array(as.matrix(fit, "p"),
               c(nsim, n.races, n.ages, n.years),
               dimnames = list(1:nsim,
                               races,
                               ages,
                               years
               ))
## Quantile of p
p.ci <- apply(p.fit, c(2:4), quantile, probs = c(0.975, 0.5, 0.025))

df.fit.eda.race <- as.data.frame.table(p.ci) %>%
    rename("race_eth" = Var2,
           "age" = Var3,
           "year" = Var4) %>%
    mutate(Var1 = case_when( 
        Var1 == "97.5%" ~ "upper95",
        Var1 == "50%" ~ "median",
        Var1 == "2.5%" ~ "lower95")) |> 
    pivot_wider(names_from = Var1, values_from = Freq) %>% 
    mutate(age = as.character(age) %>% as.numeric,
           year = as.character(year) %>% as.numeric)

## Store estimates
saveRDS(df.fit.eda.race,
        here("data", "df_logistic_fit_eda_race.rda"))

## Load estimates
# df.p.race <- readRDS(here("data", "df_logistic_fit_eda_race.rda"))


## EDA - LOGISTIC PARS FOR STATES ==============================================

# Childless by race_eth
df.stan.state <- df |>
    ungroup() |> 
    # create factors to keep
    # all combinations with
    # group_by
    mutate(
        age = factor(age,
                     levels = ages,
                     labels = ages),
        year = factor(year,
                      levels = years,
                      labels = years),
        name = factor(name, 
                      levels = name,
                      labels = name),
        race_eth = factor(race_eth,
                          levels = races,
                          labels = races)
    ) |> 
    group_by(
        year, age, name,
        .drop = FALSE
    ) |> 
    summarise(
        across(y:n, sum)
    ) |> 
    ungroup() |> 
    # Arrange for STAN
    arrange(name, year, age) 

# STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.stan.state)[1],
    Y_obs = n.years,
    A = n.ages,
    age = ages,
    G = n.states,
    # Careful !!
    # Structured as vector but order super important !!!
    # Needs to coincide with logmx_ord in STAN (multiple to_vector() )
    y = round(df.stan.state$y, 0), 
    n = round(df.stan.state$n, 0)
)

# STAN set-up
niter = 2000
nwarmup = niter/2
nchains = 4
nsim = (niter-nwarmup)*nchains

options(mc.cores = parallel::detectCores()-1)

# STAN fit
pars <- c("x0", "k", "U","p")

fit = stan(here("code", "stan", "p0_logistic_pars_eda.stan"),
           pars  = pars,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data
)

# Store and tidy logistic pars 
pars.eda <- lapply(pars[pars != "p"], function(x) {
    
    # Store in arrays
    par.draw <- array(as.matrix(fit, x),
                      c(nsim, n.years, n.states),
                      dimnames = list(1:nsim,
                                      years,
                                      states))
    # Get CI
    par.ci <- apply(par.draw, 
                    2:3, 
                    quantile, probs = c(0.975, 0.5, 0.025))
    # Store in df
    out <- as.data.frame.table(par.ci) %>%
        rename("year" = Var2,
               "name" = Var3) %>%
        mutate(Var1 = case_when( 
            Var1 == "97.5%" ~ "upper95",
            Var1 == "50%" ~ "median",
            Var1 == "2.5%" ~ "lower95")) |> 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        mutate(year = as.character(year) %>% as.numeric,
               par = x)
})

# Combine list items into one data frame
df.pars.eda.state <- do.call("rbind", pars.eda)

saveRDS(df.pars.eda.state,
        here("data", "df_logistic_pars_eda_state.rda"))

# df.pars.eda.state <- readRDS(here("data", "df_logistic_pars_eda_state.rda"))

## Store logistic fit 
p.fit <- array(as.matrix(fit, "p"),
               c(nsim, n.states, n.ages, n.years),
               dimnames = list(1:nsim,
                               states,
                               ages,
                               years
               ))
## Quantile of p
p.ci <- apply(p.fit, c(2:4), quantile, probs = c(0.975, 0.5, 0.025))

df.fit.eda.state <- as.data.frame.table(p.ci) %>%
    rename("name" = Var2,
           "age" = Var3,
           "year" = Var4) %>%
    mutate(Var1 = case_when( 
        Var1 == "97.5%" ~ "upper95",
        Var1 == "50%" ~ "median",
        Var1 == "2.5%" ~ "lower95")) |> 
    pivot_wider(names_from = Var1, values_from = Freq) %>% 
    mutate(age = as.character(age) %>% as.numeric,
           year = as.character(year) %>% as.numeric)

## Store estimates
saveRDS(df.fit.eda.state,
        here("data", "df_logistic_fit_eda_state.rda"))

## Load estimates
# df.p.state <- readRDS(here("data", "df_logistic_fit_eda_state.rda"))


## VISUALIZATIONS ==============================================================

# US pars
df.pars.eda.us |> 
    ggplot(aes(x = year, y = median,
               ymin = lower95, ymax = upper95)) +
    facet_wrap(~ par,
               scales = "free_y") +
    geom_point() +
    geom_errorbar() +
    geom_vline(xintercept = 2012, col = "red4") +
    theme_bw()

# US fit
df.p.us |> 
    left_join(df.stan.us |> 
                  mutate(
                      p = y/n,
                      year = as.character(year) |> as.numeric(),
                      age = as.character(age) |> as.numeric()
                  ),
              by = c("year", "age")) |> 
    ggplot(aes(x = age, y = median,
               ymin = lower95, ymax = upper95)) +
    facet_wrap(~ year) +
    geom_point(aes(y = p), size = 1, alpha = .5) +
    geom_line(linewidth = 1, col = "red3") +
    geom_ribbon(col = NA, alpha = .3) +
    theme_bw()

# Races pars
df.pars.eda.race |> 
    filter(
        race_eth != "NH-Other"
    ) |> 
    ggplot(aes(x = year, y = median,
               ymin = lower95, ymax = upper95,
               col = race_eth)) +
    facet_wrap(~ par,
               scales = "free_y") +
    geom_point() +
    geom_errorbar() +
    geom_vline(xintercept = 2012) +
    theme_bw()

# Races fit
df.p.race |>
    left_join(df.stan.race |> 
                  mutate(
                      p = y/n,
                      year = as.character(year) |> as.numeric(),
                      age = as.character(age) |> as.numeric()
                  ),
              by = c("year", "age", "race_eth")) |> 
    filter(year %in% c(1995, 2004, 2012, 2020),
           race_eth != "NH-Other") |> 
    ggplot(aes(x = age, y = median,
               ymin = lower95, ymax = upper95,
               group = race_eth, col = race_eth)) +
    facet_wrap( ~ year) +
    geom_point(aes(y = p), size = 1, alpha = .5) +
    geom_line(linewidth = 1) +
    geom_ribbon(col = NA, alpha = .3) +
    theme_bw()

# States pars
df.pars.eda.state |> 
    ggplot(aes(x = year, y = median,
               ymin = lower95, ymax = upper95,
               col = name)) +
    facet_wrap(~ par,
               scales = "free_y") +
    geom_line() +
    # geom_errorbar() +
    geom_vline(xintercept = 2012) +
    theme_bw() +
    theme(legend.position = "none")

# States fit
