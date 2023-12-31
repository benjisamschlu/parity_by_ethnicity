
rm(list=ls())

library(tidyverse)
library(here)
library(rstan)
library(splines)



## Load data -------------------------------------------------------------------

## Careful: no observation for all years (~ every 2 y)

## Childless by race and state
df <- readRDS(here("data", "df_chldness.rda")) 



## STAN section ----------------------------------------------------------------

## Key dimensions
ages <- unique(df$age)
n.ages <- length(ages)
years <- unique(df$year)
n.years <- length(years)
states <-  unique(df$name)
n.states <- length(states)

## Create linear splines basis
B   <- bs(ages, knots=seq(20, 40, 10), degree=1)
nalpha = ncol(B)


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
    ungroup() |> 
    ## Take complement to model
    ## p = (1 - childlessness)
    mutate(
        y = n - y
    ) |> 
    ## Arrange for STAN
    arrange(name, year, age) 

stan_data <- list(# Dimensions
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
        here("data", "df_p0_state_logistic_splines_fit.rda"))

df.p <- readRDS(
    here("data", "df_p0_state_logistic_splines_fit.rda")
)


## Visualize outputs
## Look at fit for each state
df.p.states <- df.p |> 
    left_join(
        df.stan |> 
            mutate(
                p = y/n,
                age = as.character(age) |> as.numeric(),
                year = as.character(year) |> as.numeric()
            ) |> 
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
    filename = here("plots", "p0_state_logistic_spline_fit.pdf"), 
    plot = gridExtra::marrangeGrob(fig.fit.state, nrow=1, ncol=1), 
    width = 15, height = 9
)



