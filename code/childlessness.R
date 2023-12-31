
rm(list=ls())

library(tidyverse)
library(here)
library(rstan)



## Load data -------------------------------------------------------------------

## Careful: no observation for all years (~ every 2 y)

## Childless by race and state
df <- readRDS(here("data", "df_chldness.rda"))



## EDA -------------------------------------------------------------------------

## Shape of childlessness over time
df |> 
    group_by(age, year) |> 
    summarise(across(y:n, sum)) |> 
    mutate(p = y/n) |> 
    ggplot(aes(x = age, y = p)) + 
    facet_wrap(~ year) +
    geom_point() +
    theme_bw()

## Temporal variation
df |> 
    group_by(age, year) |> 
    summarise(across(y:n, sum)) |> 
    mutate(p = y/n) |> 
    ggplot(aes(x = age, y = p, group = year, col = year)) + 
    geom_line() +
    theme_bw()

## Temporal variation by states
states.names <- unique(df$name)
states.drawn <- sample(states.names, 10)
df |> 
    #filter(name %in% states.drawn) |> 
    group_by(age, year, name) |> 
    summarise(across(y:n, sum),
              p = y/n) |> 
    ungroup() |> 
    ggplot(aes(x = age, y = p, group = year, col = year)) + 
    facet_wrap(~ name,
               ncol = 10,
               nrow = 6) +
    geom_line() +
    theme_bw()

## Temporal variation by race
df |> 
    group_by(age, year, race_eth) |> 
    summarise(across(y:n, sum),
              p = y/n) |> 
    ungroup() |> 
    ggplot(aes(x = age, y = p, group = year, col = year)) + 
    facet_grid(~ race_eth) +
    geom_point() +
    theme_bw()

## Temporal variation by state & race
pop_states <- c("California", "Texas", "Florida", "New York", "Pennsylvania",
                "Illinois", "Ohio", "Georgia", "North Carolina", "Michigan")
df |> 
    filter(name %in% pop_states,
           race_eth %in% c("Hispanic", "NH-Black", "NH-White")) |> 
    group_by(age, year, name, race_eth) |> 
    summarise(across(y:n, sum),
              p = y/n) |> 
    ungroup() |> 
    ggplot(aes(x = age, y = p, group = year, col = year)) + 
    facet_grid(race_eth ~ name) +
    geom_line() +
    theme_bw()



## STAN section national level -------------------------------------------------

## STAN data
df.stan <- df |> 
    group_by(age, year) |> 
    summarise(across(y:n, sum)) |> 
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

# df.p <- readRDS(here("data", "df_p0_fitted.rda"))



## STAN section race level -----------------------------------------------------

## STAN data
df.stan <- df |> 
    group_by(age, year, race_eth) |> 
    summarise(across(y:n, sum)) |> 
    mutate(race_eth = factor(race_eth)) |> 
    ## Arrange for STAN
    arrange(race_eth, year, age) 

## EXPLORATORY WORK BY SETTING
## PARS ~ N(0,1)

ages <- unique(df.stan$age)
n.ages <- length(ages)
years <- unique(df.stan$year)
n.years <- length(years)
races <-  unique(df.stan$race_eth)
n.races <- length(races)

stan_data <- list(# Dimensions
    N = dim(df.stan)[1],
    Y_obs = n.years,
    A = n.ages,
    age = ages,
    G = n.races,
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
pars <- c("L", "k", "x0")

fit = stan(here("code", "stan", "p0_group_logistic_eda.stan"),
           pars  = pars,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data
)

## Store ad tidy estimates 
pars.fit <- lapply(as.character(pars), function(x) {
    
    ## Store in arrays
    par.draw <- array(as.matrix(fit, x),
                 c(nsim, n.years, n.races),
                 dimnames = list(1:nsim,
                                 years,
                                 races))
    ## Get CI
    par.ci <- apply(par.draw, 
                 c(2,3), 
                 quantile, probs = c(0.975, 0.5, 0.025))
    ## Store in df
    out <- as.data.frame.table(par.ci) %>%
        rename("year" = Var2,
               "race" = Var3) %>%
        mutate(Var1 = case_when( 
            Var1 == "97.5%" ~ "upper95",
            Var1 == "50%" ~ "median",
            Var1 == "2.5%" ~ "lower95")) |> 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        mutate(year = as.character(year) %>% as.numeric,
               par = x)
    
    return(out)
})

df.pars.race.eda <- do.call("rbind", pars.fit)

saveRDS(df.pars.race.eda,
        here("data", "df_p0_race_eda.rda"))

# df.pars.race.eda <- readRDS(here("data", "df_p0_race_eda.rda"))
        


## BUILD MODEL FROM 
## EXPLORATORY WORK

options(mc.cores = parallel::detectCores()-1)

## STAN fit
pars <- c("L", "k", "x0")

fit = stan(here("code", "stan", "p0_race_logistic.stan"),
           pars  = pars,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data
)

## Store
p.fit <- array(as.matrix(fit, "p"),
               c(nsim, n.races, n.ages, n.years),
               dimnames = list(1:nsim,
                               races,
                               ages,
                               years
                               ))
## Quantile of p
p.ci <- apply(p.fit, c(2:4), quantile, probs = c(0.975, 0.5, 0.025))

df.p <- as.data.frame.table(p.ci) %>%
    rename("race" = Var2,
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
        here("data", "df_p0_race_fit.rda"))

# df.pars.race.eda <- readRDS(here("data", "df_p0_race_fit.rda"))


## BUILD MODEL FROM 
## DISCUSSION WITH MONICA
## -> RW2 and RW1

options(mc.cores = parallel::detectCores()-1)

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
               c(nsim, n.races, n.ages, n.years),
               dimnames = list(1:nsim,
                               races,
                               ages,
                               years
               ))
## Quantile of p
p.ci <- apply(p.fit, c(2:4), quantile, probs = c(0.975, 0.5, 0.025))

df.p <- as.data.frame.table(p.ci) %>%
    rename("race" = Var2,
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
        here("data", "df_p0_race_rw_fit.rda"))

# df.pars.race.fit <- readRDS(here("data", "df_p0_race_fit.rda"))

## Store and tidy logistic estimates 
pars <- pars[pars != "p"]
pars.fit <- lapply(as.character(pars), function(x) {
    
    ## Store in arrays
    par.draw <- array(as.matrix(fit, x),
                 c(nsim, n.years, n.races),
                 dimnames = list(1:nsim,
                                 years,
                                 races))
    ## Get CI
    par.ci <- apply(par.draw, 
                 c(2,3), 
                 quantile, probs = c(0.975, 0.5, 0.025))
    ## Store in df
    out <- as.data.frame.table(par.ci) %>%
        rename("year" = Var2,
               "race" = Var3) %>%
        mutate(Var1 = case_when( 
            Var1 == "97.5%" ~ "upper95",
            Var1 == "50%" ~ "median",
            Var1 == "2.5%" ~ "lower95")) |> 
        pivot_wider(names_from = Var1, values_from = Freq) %>% 
        mutate(year = as.character(year) %>% as.numeric,
               par = x)
    
    return(out)
})

df.pars.race.fit <- do.call("rbind", pars.fit)

saveRDS(df.pars.race.fit,
        here("data", "df_p0_race_rw_pars.rda"))

# df.pars.race.eda <- readRDS(here("data", "df_p0_race_eda.rda"))
        




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

## STAN fit
pars <- c("L", "k", "x0")

fit = stan(here("code", "stan", "p0_group_logistic_eda.stan"),
           pars  = pars,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = 4,
           data = stan_data
)

## Store ad tidy estimates 
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

df.pars.fit <- do.call("rbind", pars.fit)

saveRDS(df.pars.fit,
        here("data", "df_p0_state_eda.rda"))

# df.pars.states.eda <- readRDS(here("data", "df_p0_state_eda.rda"))


## BUILD MODEL FROM 
## EXPLORATORY WORK
options(mc.cores = parallel::detectCores()-1)

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




## Visu ------------------------------------------------------------------------


## NATIONAL 

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
    geom_point(aes(y = p), col = "lightblue") +
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



## BY RACE

## EDA setting logistic par priors to N(0,1)
df.pars.race.fit |> 
    filter(race != "NH-Other") |> 
    ggplot(aes(x = year, y = median, group = race, col = race)) +
    facet_wrap(~ par,
               scales = "free_y") +
    geom_pointrange(aes(ymin = lower95, ymax = upper95)) +
    geom_line() +
    theme_bw()
## Temporal trend in x0: race around a global mean, mean ~ RW
## Might assume that U is 1
## k noisy and seems to be race specific without temporal trend over time, hiear
## L RW(1) by race
    
## Plot of fit vs data points
## All races and years
df.p |> 
    filter(race != "NH-Other") |> 
    left_join(
        df |> 
            group_by(age, year, race_eth) |> 
            summarise(across(y:n, sum),
                      p = y/n) |> 
            ungroup() |> 
            dplyr::select(age, year, race = race_eth, p),
        by = c("age", "year", "race")
    ) |> 
    ggplot(aes(x = age, group = race, col = race, fill = race)) +
    facet_wrap(~ year,
               scales = "free_y") +
    geom_point(aes(y = p)) +
    geom_line(aes(y = median), linewidth = 1) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95),
                alpha = .3, col = NA) +
    theme_bw()

## NH-Black over years
df.p |> 
    left_join(
        df |> 
            group_by(age, year, race_eth) |> 
            summarise(across(y:n, sum),
                      p = y/n) |> 
            ungroup() |> 
            dplyr::select(age, year, race = race_eth, p),
        by = c("age", "year", "race")
    ) |> 
    filter(race == "NH-Black") |> 
    ggplot(aes(x = age, group = race, col = race)) +
    facet_wrap(~ year,
               scales = "free_y") +
    geom_point(aes(y = p)) +
    geom_line(aes(y = median), linewidth = 1) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95)) +
    theme_bw()


## BY STATES
# EDA setting logistic par priors to N(0,1)
ex.states <- sample(unique(df.pars.states.eda$state), 10)
df.pars.states.eda |> 
    filter(state %in% ex.states) |> 
    ggplot(aes(x = year, y = median, group = state, col = state)) +
    facet_wrap(~ par,
               scales = "free_y") +
    #geom_line() +
    geom_pointrange(aes(ymin = lower95, ymax = upper95)) +
    theme_bw() +
    theme(legend.position = "none")

## Look at fit for each state
df.p.states <- df.p |> 
    left_join(
        df |> 
            group_by(age, year, name) |> 
            summarise(across(y:n, sum),
                      p = y/n) |> 
            ungroup() |> 
            dplyr::select(age, year, state = name, p),
        by = c("age", "year", "state")
    )

fig.fit.state <- lapply(unique(df$name), function(x) {
    
    
    fig.out <- df.p.states |> 
        filter(state == x) |> 
        ggplot(aes(x = age)) +
        facet_wrap(~ year,
                   scales = "free_y") +
        geom_point(aes(y = p)) +
        geom_line(aes(y = median), linewidth = 1) +
        geom_ribbon(aes(ymin = lower95, ymax = upper95),
                    alpha = .3, col = NA) +
        theme_bw() +
        labs(title = x)
    
    return(fig.out)
})

## Save list of plots
ggsave(
    filename = here("plots", "p0_state_rw_fit.pdf"), 
    plot = gridExtra::marrangeGrob(fig.fit.state, nrow=1, ncol=1), 
    width = 15, height = 9
)

    
