
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
        race_eth != "NH-Other",
        # Year > 2000
        year >= 2000
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



## MODEL RACES WITH FULL MODEL =================================================

## Childless by race_eth
df.stan <- df |>
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
        race_eth = factor(race_eth)
    ) |> 
    group_by(
        age, year, race_eth,
        .drop = FALSE
    ) |> 
    summarise(
        across(y:n, sum)
    ) |> 
    ungroup() |> 
    ## Arrange for STAN
    arrange(race_eth, year, age) 

## Create Bsplines basis
B <- bs(ages, knots=30, degree=3)
nalpha = ncol(B)
# matplot(15:44, B, type="l")

## STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.stan)[1],
    Y_obs = n.years,
    A = n.ages,
    age = ages,
    G = n.races,
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

fit = stan(here("code", "stan", "p0_relational_logistic.stan"),
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
saveRDS(df.p,
        here("data", "df_p0_race_rel_logistic_fit.rda"))

## Load estimates
# df.p <- readRDS(
#     here("data", "df_p0_race_logisticrw2_splines_fit.rda")
# )



## MODEL SAMPLE STATES WITH FULL MODEL =========================================

# Sample States
sample.states <- c("Alaska",
                   "California",
                   "Florida",
                   "Georgia",
                   "New York",
                   "North Dakota",
                   "Oklahoma",
                   "Pennsylvania",
                   "Rhode Island",
                   "Tennessee",
                   "Texas",
                   "Washington",
                   "Wisconsin",
                   "Wyoming")

# Childless by sample states
df.stan <- df |>
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
        race_eth = factor(race_eth)
    ) |> 
    group_by(
        age, year, name,
        .drop = FALSE
    ) |> 
    summarise(
        across(y:n, sum)
    ) |> 
    ungroup() |> 
    filter(
        name %in% sample.states
    ) |> 
    ## Arrange for STAN
    arrange(name, year, age) 

## Create Bsplines basis
B <- bs(ages, knots=30, degree=3)
nalpha = ncol(B)
# matplot(15:44, B, type="l")

## STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.stan)[1],
    Y_obs = n.years,
    A = n.ages,
    age = ages,
    G = length(sample.states),
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

fit = stan(here("code", "stan", "p0_logisticrw2_splines.stan"),
           pars  = pars,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains,
           data = stan_data
)

saveRDS(fit,
        here("data", "fit_backup.rda"))

## Store p 
p.fit <- array(as.matrix(fit, "p"),
               c(nsim, length(sample.states), n.ages, n.years),
               dimnames = list(1:nsim,
                               sample.states,
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
        here("data", "df_p0_sample_state_logisticrw2_splines_fit.rda"))

## Load estimates
# df.p <- readRDS(
#     here("data", "df_p0_sample_state_logisticrw2_splines_fit.rda")
# )




## VISUALIZATION OF ESTIMATES --------------------------------------------------

# Race
df.p |> 
    left_join(
        df.stan |> 
            mutate(age = as.character(age) |> as.numeric(),
                   year = as.character(year) |> as.numeric()
            ),
        by = c("age", "year", "race_eth")
    ) |> 
    filter(
        year %in% c(2000, 2008, 2014, 2020)
    ) |> 
    ggplot(aes(x = age)) +
    facet_grid(race_eth ~ year) +
    geom_line(aes(y = median),
              col = "red4") +
    geom_point(aes(y = y/n), size = 1) +
    theme_bw()

df.p |> 
    left_join(
        df.stan |> 
            mutate(age = as.character(age) |> as.numeric(),
                   year = as.character(year) |> as.numeric()
            ),
        by = c("age", "year", "race_eth")
    ) |> 
    filter(year >= 2008) |> 
    ggplot(aes(x = age)) +
    facet_wrap( ~ race_eth) +
    geom_line(aes(y = median, col = year, group = year)) +
    theme_bw()




# Sample States
## Store state plots in list
df.p.state <-
    df.p |> 
    left_join(
        df.stan |> 
            mutate(age = as.character(age) |> as.numeric(),
                   year = as.character(year) |> as.numeric()
            ),
        by = c("age", "year", "name")
    ) 

fig.fit.state <- lapply(sample.states, function(x) {
    
    
    fig.out <- df.p.state |> 
        filter(name == x) |> 
        ggplot(aes(x = age)) +
        facet_wrap(~ year,
                   scales = "free_y") +
        geom_point(aes(y = y/n)) +
        geom_line(aes(y = median), linewidth = 1) +
        geom_ribbon(aes(ymin = lower95, ymax = upper95),
                    alpha = .3, col = NA) +
        theme_bw() +
        labs(title = x)
})

## Save list of plots
ggsave(
    filename = here("plots", "p0_sample_state_logisticrw2_Cspline_fit.pdf"), 
    plot = gridExtra::marrangeGrob(fig.fit.state, nrow=1, ncol=1), 
    width = 15, height = 9
)


