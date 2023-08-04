
rm(list=ls())

library(tidyverse)
library(here)
library(rstan)



## Load data -------------------------------------------------------------------

## Careful: no observation for all years (~ every 2 y)
df <- readRDS(here("data", "df_chldness.rda"))



## EDA -------------------------------------------------------------------------

df |> 
    mutate(p = y/n) |> 
    ggplot(aes(x = age, y = p)) + 
    facet_wrap(~ year) +
    geom_point() +
    theme_bw()



## STAN section ----------------------------------------------------------------

## STAN data
df.stan <- df |> 
    ## Arrange for STAN
    arrange(year, age) 

years <- unique(df.stan$year)
n.years <- length(years)
ages <- unique(df.stan$age)
n.ages <- length(ages)
years.all <- 1990:2020
n.years.all <- length(years.all)
year.i <- which(years.all %in% years)

stan_data <- list(# Dimensions
    N = dim(df.stan)[1],
    Y_obs = n.years,
    A = n.ages,
    age = ages,
    Y = n.years.all,
    t_i = year.i,
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
fit = stan(here("code", "stan", "p0_logistic_fct.stan"),
           pars  = c("p"),
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data
)

## Estimates
fit

## Store
p.fit <- array(as.matrix(fit, "p"),
               c(nsim, n.ages, n.years.all),
               dimnames = list(1:nsim,
                               ages,
                               years.all))
## Quantile of p
p.ci <- apply(p.fit, c(2,3), quantile, probs = c(0.975, 0.5, 0.025))

df.p <- as.data.frame.table(p.ci) %>%
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
saveRDS(df.p,
        here("data", "df_p0_fitted.rda"))


## Visu ------------------------------------------------------------------------

## Fit to points
df.p |> 
    left_join(
        df.stan,
        by = c("age", "year")
    ) |> 
    mutate(
        p = y/n
    ) |>
    ggplot(aes(x = age)) +
    facet_wrap(~ year) +
    geom_point(aes(y = p)) +
    geom_line(aes(y = median), col = "red4") +
    geom_ribbon(aes(ymin = lower95, ymax = upper95)) +
    theme_bw()

## Change over time
df.p |> 
    left_join(
        df.stan,
        by = c("age", "year")
    ) |> 
    mutate(
        p = y/n
    ) |>
    filter(year %in% seq(1990, 2020, 5)) |> 
    ggplot(aes(x = age, group = year)) +
    geom_line(aes(y = median, col = year)) +
    geom_point(aes(y = p, col = year)) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = year), alpha = 0.3) +
    theme_bw()
    
