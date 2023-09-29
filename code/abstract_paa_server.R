
rm(list=ls())

library(tidyverse)
library(here)
library(rstan)
library(splines)



## Load data -------------------------------------------------------------------

## Careful: no observation for all years (~ every 2 y)

## Childless by race and state
df <- readRDS(here("data", "df_chldness.rda")) |> 
    filter(
        ## Focus on 3 main bridged race_eth
        race_eth != "NH-Other"
    )



## Key dimensions in code ------------------------------------------------------

ages <- unique(df$age)
n.ages <- length(ages)
years <- unique(df$year)
n.years <- length(years)
states <-  unique(df$name)
n.states <- length(states)
races <- unique(df$race_eth)
n.races <- length(races)



## Model states with full model ----------------------------------------

## Childless by race_eth
df.stan <- df |>
    ## create factors to keep
    ## all combinations with
    ## group_by
    mutate(
        age = factor(age,
                     levels = ages,
                     labels = ages),
        year = factor(year,
                      levels = years,
                      labels = years),
        name = factor(name, 
                          levels = name,
                          labels = name)
    ) |> 
    group_by(
        age, year, name,
        .drop = FALSE
    ) |> 
    summarise(
        across(y:n, sum)
    ) |> 
    ungroup() |> 
    ## Arrange for STAN
    arrange(name, year, age) 

## Create linear splines basis
B <- bs(ages, knots=seq(20, 40, 10), degree=3)
nalpha = ncol(B)

## STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.stan)[1],
    Y_obs = n.years,
    A = n.ages,
    age = ages,
    G = n.states,
    B = B,
    K = nalpha,
    # Careful !!
    # Structured as vector but order super important !!!
    # Needs to coincide with logmx_ord in STAN (multiple to_vector() )
    y = round(df.stan$y, 0), 
    n = round(df.stan$n, 0)
)

## STAN set-up
niter = 2000
nwarmup = niter/2
nchains = 4
nsim = (niter-nwarmup)*nchains

options(mc.cores = parallel::detectCores()-1)

## STAN fit
pars <- c("p")

fit = stan(here("code", "stan", "p0_logistic_splines.stan"),
           pars  = pars,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data
)

## Store p 
p.fit <- array(as.matrix(fit, "p"),
               c(nsim, n.states, n.ages, n.years),
               dimnames = list(1:nsim,
                               states,
                               ages,
                               years
               ))
## Quantile of p
p.ci <- apply(p.fit, c(2:4), quantile, probs = c(0.975, 0.5, 0.025))

df.p <- as.data.frame.table(p.ci) %>%
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
saveRDS(df.p,
        here("data", "df_p0_state_logistic_splines_fit.rda"))



