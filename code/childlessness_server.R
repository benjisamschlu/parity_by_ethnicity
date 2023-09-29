
rm(list=ls())

library(tidyverse)
library(here)
library(rstan)



## Load data -------------------------------------------------------------------

## Careful: no observation for all years (~ every 2 y)

## Childless by race and state
df <- readRDS(here("data", "df_chldness.rda"))



## STAN section state level -----------------------------------------------------

ages <- unique(df$age)
n.ages <- length(ages)
years <- unique(df$year)
n.years <- length(years)
states <-  unique(df$name)
n.states <- length(states)

## STAN data
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
                      levels = states,
                      labels = states)
    ) |> 
    group_by(age, year, name,
             .drop = FALSE) |> 
    summarise(across(y:n, sum)) |> 
    ## Arrange for STAN
    arrange(name, year, age) 

## EXPLORATORY WORK BY SETTING
## PARS ~ N(0,1)

stan_data <- list(# Dimensions
    N = dim(df.stan)[1],
    Y_obs = n.years,
    A = n.ages,
    age = ages,
    G = n.states,
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

## STAN fit for EDA

# pars <- c("L", "k", "x0")
# 
# fit = stan(here("code", "stan", "p0_group_logistic_eda.stan"),
#            pars  = pars,
#            include = TRUE,
#            iter = niter,
#            warmup = nwarmup,
#            chains = 4,
#            data = stan_data
# )
# 
# ## Store ad tidy estimates 
# pars.fit <- lapply(as.character(pars), function(x) {
#     
#     ## Store in arrays
#     par.draw <- array(as.matrix(fit, x),
#                       c(nsim, n.years, n.states),
#                       dimnames = list(1:nsim,
#                                       years,
#                                       states))
#     ## Get CI
#     par.ci <- apply(par.draw, 
#                     c(2,3), 
#                     quantile, probs = c(0.975, 0.5, 0.025))
#     ## Store in df
#     out <- as.data.frame.table(par.ci) %>%
#         rename("year" = Var2,
#                "state" = Var3) %>%
#         mutate(Var1 = case_when( 
#             Var1 == "97.5%" ~ "upper95",
#             Var1 == "50%" ~ "median",
#             Var1 == "2.5%" ~ "lower95")) |> 
#         pivot_wider(names_from = Var1, values_from = Freq) %>% 
#         mutate(year = as.character(year) %>% as.numeric,
#                par = x)
#     
#     return(out)
# })
# 
# df.pars.eda <- do.call("rbind", pars.fit)
# 
# saveRDS(df.pars.eda,
#         here("data", "df_p0_state_eda.rda"))
# 


## STAN FIT RW

## STAN fit
pars <- c("L", "k", "x0", "p")

fit = stan(here("code", "stan", "p0_group_logistic.stan"),
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
    rename("state" = Var2,
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
        here("data", "df_p0_state_rw_fit.rda"))

# df.p <- readRDS(here("data", "df_p0_state_rw_fit.rda"))


## Store and tidy logistic estimates 
pars <- pars[pars != "p"]
pars.fit <- lapply(as.character(pars), function(x) {
    
    ## Store in arrays
    par.draw <- array(as.matrix(fit, x),
                      c(nsim, n.years, n.states),
                      dimnames = list(1:nsim,
                                      years,
                                      states))
    ## Get CI
    par.ci <- apply(par.draw, 
                    c(2,3), 
                    quantile, probs = c(0.975, 0.5, 0.025))
    ## Store in df
    out <- as.data.frame.table(par.ci) %>%
        rename("year" = Var2,
               "state" = Var3) %>%
        mutate(Var1 = case_when( 
            Var1 == "97.5%" ~ "upper95",
            Var1 == "50%" ~ "median",
            Var1 == "2.5%" ~ "lower95")) |> 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        mutate(year = as.character(year) %>% as.numeric,
               par = x)
    
    return(out)
})

df.pars.state.fit <- do.call("rbind", pars.fit)

saveRDS(df.pars.state.fit,
        here("data", "df_p0_state_rw_pars.rda"))

# df.pars.state.fit <- readRDS(here("data", "df_p0_state_rw_pars.rda"))

