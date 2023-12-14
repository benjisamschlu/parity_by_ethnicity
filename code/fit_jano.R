
##==============================================================================
##
##  Project title: Modeling of childlessness
##
##  This code estimate proportion of parent (1-childless) with a Janoshek
##  curve.
##
## 
##  Author: Benjamin Schlüter
##  Date: December 2023
##==============================================================================
##
##  Notes
## ------
## 1. Only consider years from 2008 as previous years might over-estimate
##    childlessness (see US Census Bureau report)
## 
##==============================================================================



## ==== LOAD PACKAGES ==========================================================

# Install/load packages
packages <- c("tidyverse", "here", "rstan")
for(p in packages){
    if(!require(p,character.only = TRUE)) install.packages(p)
    library(p,character.only = TRUE)
}




## ==== LOAD DATA ==============================================================

# Childless corrected with hh info
# for the years 2008 and 2010
df_cor <- readRDS(
    here(
        "data", 
        "df_childness_cor_2008_10.rds"
    )
) 
# Childless by race and state
df <- readRDS(
    here(
        "data", 
        "df_childness.rds"
    )
) |> 
    filter(
        # Focus on years using coresident info
        year >= 2012
    )

# Bind corrected years with the rest
df <- bind_rows(
    df_cor, df
)

# Key dimensions in code 
ages <- unique(df$age)
n.ages <- length(ages)
years <- unique(df$year)
n.years <- length(years)



## ==== INITIAL PAR VALUES =====================================================

# Data of childlessness at US level
df.us <- df |>
    ungroup() |> 
    summarise(
        .by = c(year, age),
        
        across(y:n, sum),
        # Perspective of prop parent instead of childless (Janoshek fct)
        y = n - y,
        # Compute prop parent (Janoshek fct is increasing)
        p = y / n
    ) |> 
    ungroup() |> 
    # Arrange for STAN
    arrange(year, age)

# Visualize the data
df.us |> 
    ggplot(aes(x = age, y = p)) +
    facet_wrap(~ year) +
    geom_point() +
    theme_bw()

# Janoscheck growth model
get_jano <- function(U, L, k, c, x) {
    
    y <- U - (U - L) * exp(-k * x^c)
}

# Get plausible initial pars value by trial error
# for mid-period year: 2010
plot(x = 0:29, y = df.us$p[df.us$year == 2010])
lines(x = 0:29, y = get_jano(0.83, 0.01, 0.012, 1.8, 0:29), type = "l")
par.init <- c(0.83, 0.01, 0.012, 1.8)

## Data used in optim() to minimize RSS
df.rss <- 
    df.us |> 
    mutate(
        x = age - 15
    ) |> 
    dplyr::select(
        -c(y, n, age)
    ) 

# Function minimizing the squared differences between observed p
# and predicted p by janoscheck
min.RSS <- function(data, par) {
    with(data, 
         sum(((par[1] - (par[1] - par[2]) * exp(-par[3] * x^par[4])) - p)^2)
         )
}

# Minimize and extract estimated parameters
# over each year
pars <- sapply(years, function(y) {
    
    out <- optim(par = par.init, 
                  fn = min.RSS, 
                  data = df.rss[df.rss$year == y, ],
                  method = "L-BFGS-B",
                  lower = c(0.7, 0, 0, 1),
                  upper = c(1, 0.05, 1, 3)
                 )$par
    
}, 
simplify = FALSE,
USE.NAMES = TRUE)

pars <- do.call("rbind", pars)

# Exctract range and variance 
# to build initial par draws in STAN
range.pars <- apply(pars, 2, range)
var.pars <- apply(pars, 2, var)

# Check fit of obtained pars
for (y in years) {
    
    i <- which(years == y)
    
    plot(x = 0:29, y = df.us$p[df.us$year == y], main = y)
    lines(x = 0:29, y = get_jano(pars[i, 1], pars[i, 2], pars[i, 3], pars[i, 4], 0:29), type = "l")
}



## ==== FIT IN STAN ============================================================

# STAN set-up ------------------------------------------------------------------

niter = 2000
nwarmup = niter/2
nchains = 4
nsim = (niter-nwarmup)*nchains
options(mc.cores = parallel::detectCores()-1)
pars <- c("p")

# Inits -----------------------------------------------------------------------

# Consider US, so only one group
n.groups <- 1

stanInit = function() {
    range.U <- range.pars[, 1]
    range.L <- range.pars[, 2]
    range.k <- range.pars[, 3]
    range.c <- range.pars[, 4]
    
    var.U <- var.pars[1]
    var.L <- var.pars[2]
    var.k <- var.pars[3]
    var.c <- var.pars[4]
    
    list(
        U = matrix( runif(n.years*n.groups, range.pars[1, 1], range.pars[2, 1]), n.years, n.groups),
        L = matrix( runif(n.years*n.groups, range.pars[1, 2], range.pars[2, 2]), n.years, n.groups),
        k = matrix( runif(n.years*n.groups, range.pars[1, 3], range.pars[2, 3]), n.years, n.groups),
        c = matrix( runif(n.years*n.groups, range.pars[1, 4], range.pars[2, 4]), n.years, n.groups),
        sigma_U = runif(1, var.pars[1]/10, var.pars[1]*10),
        sigma_L = runif(1, var.pars[2]/10, var.pars[2]*10),
        sigma_k = runif(1, var.pars[3]/10, var.pars[3]*10),
        sigma_c = runif(1, var.pars[4]/10, var.pars[4]*10)
    )
}

## US-level --------------------------------------------------------------------

# STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.us)[1],
    Y_obs = n.years,
    A = n.ages,
    age = 0:29,
    G = n.groups,
    # Careful !!
    # Structured as vector but order super important !!!
    # Needs to coincide with logmx_ord in STAN (multiple to_vector() )
    y = round(df.us$y, 0), 
    n = round(df.us$n, 0)
)

fit = stan(here("code", "stan", "p0_jano.stan"),
           data = stan_data,
           pars  = pars,
           init = stanInit,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains
)

# Store p 
p.fit <- array(as.matrix(fit, "p"),
               c(nsim, n.ages, n.years),
               dimnames = list(1:nsim,
                               ages,
                               years
               ))
## Quantile of p
p.ci <- apply(p.fit, c(2:3), quantile, probs = c(0.975, 0.5, 0.025))

df.p.us <- as.data.frame.table(p.ci) %>%
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
saveRDS(df.p.us,
        here(
            "data", 
            "df_p0_us_jano_fit.rds"
            )
        )

## Load estimates
# df.p.us <- readRDS(
#     here(
#         "data", 
#         "df_p0_us_jano_fit.rds"
#         )
#     )

## Race/eth-level --------------------------------------------------------------

# Data of childlessness at racial/ethnic level
df.race <- df |>
    # Focus on most populous race/ethnicity
    filter(
        race_eth != "NH-Other"
        ) |> 
    ungroup() |> 
    summarise(
        .by = c(year, race_eth, age),
        
        across(y:n, sum)
    ) |> 
    ungroup() |> 
    mutate(
        # Janoscheck requires to define
        # the outcome in terms of % parent
        y = n - y,
        p = y / n
    ) |> 
    ## Arrange for STAN
    arrange(race_eth, year, age)

races <- unique(df.race$race_eth)
n.races <- length(races)

# Define grouping in STAN
n.groups <- n.races

# STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.race)[1],
    Y_obs = n.years,
    A = n.ages,
    age = 0:29,
    G = n.groups,
    # Careful !!
    # Structured as vector but order super important !!!
    # Needs to coincide with logmx_ord in STAN (multiple to_vector() )
    y = round(df.race$y, 0), 
    n = round(df.race$n, 0)
)

fit = stan(here("code", "stan", "p0_jano.stan"),
           data = stan_data,
           pars  = pars,
           init = stanInit,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains
)

# Store p 
p.fit <- array(as.matrix(fit, "p"),
               c(nsim, n.races, n.ages, n.years),
               dimnames = list(1:nsim,
                               races, 
                               ages,
                               years
               ))
# Quantile of p
p.ci <- apply(p.fit, c(2:4), quantile, probs = c(0.975, 0.5, 0.025))

df.p.race <- as.data.frame.table(p.ci) %>%
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

# Store estimates
saveRDS(df.p.race,
        here(
            "data", 
            "df_p0_race_jano_fit.rds"
            )
        )

# Load estimates
# df.p.race <- readRDS(
#     here(
#         "data", 
#         "df_p0_race_jano_fit.rds"
#         )
#     )

## State-level -----------------------------------------------------------------

# Sample of States
# Make computation quicker during the dvlpt of model
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

# Data of childlessness at state level
df.state <- df |>
    ungroup() |> 
    mutate(
        # Define factor to keep group 
        # of 0
        name = factor(name),
        year = factor(year),
        age = factor(age)
    ) |> 
    group_by(
        year, name, age,
        .drop = FALSE
    ) |> 
    summarise(
        across(y:n, sum)
    ) |> 
    mutate(
        # Janoscheck requires to define
        # the outcome in terms of % parent
        y = n - y,
        p = y / n
    ) |> 
    filter(
        name %in% sample.states
    ) |> 
    ## Arrange for STAN
    arrange(name, year, age)

states <- unique(df.state$name)
n.states <- length(states)

# Define grouping in STAN
n.groups <- n.states

## STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.state)[1],
    Y_obs = n.years,
    A = n.ages,
    age = 0:29,
    G = n.groups,
    # Careful !!
    # Structured as vector but order super important !!!
    # Needs to coincide with logmx_ord in STAN (multiple to_vector() )
    y = round(df.state$y, 0), 
    n = round(df.state$n, 0)
)

fit = stan(here("code", "stan", "p0_jano.stan"),
           data = stan_data,
           pars  = pars,
           init = stanInit,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains
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

df.p.state <- as.data.frame.table(p.ci) %>%
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
saveRDS(df.p.state,
        here(
            "data", 
            "df_p0_state_jano_fit.rds"
            )
        )

## Load estimates
# df.p.race <- readRDS(
#     here(
#         "data", 
#         "df_p0_state_jano_fit.rds"
#         )
#     )


## === VISUALIZATION OF FIT ====================================================

# US
df.p.us |> 
    left_join(
        df.us |> 
            mutate(age = as.character(age) |> as.numeric(),
                   year = as.character(year) |> as.numeric()
            ),
        by = c("age", "year")
    ) |> 
    ggplot(aes(x = age)) +
    facet_wrap( ~ year) +
    geom_ribbon(aes(ymin = lower95,
                    ymax = upper95),
                col = NA,
                fill = "red4",
                alpha = .3) + 
    geom_line(aes(y = median),
              col = "red4",
              linewidth = 1) +
    # geom_ribbon(col = NA) +
    geom_point(aes(y = y / n), 
               col = "skyblue3",
               alpha = .6) +
    theme_bw() +
    labs(x = "Age",
         y = "Proportion parent")

# Race over time
df.p.race |> 
    left_join(
        df.race |> 
            mutate(age = as.character(age) |> as.numeric(),
                   year = as.character(year) |> as.numeric()
            ),
        by = c("age", "year", "race_eth")
    ) |> 
    ggplot(aes(x = age,
               group = race_eth,
               col = race_eth)) +
    facet_wrap(~ year) +
    geom_line(aes(y = median)) +
    # geom_ribbon(col = NA) +
    geom_point(aes(y = y / n), 
               alpha = .6) +
    theme_bw() +
    theme(legend.position = c(0.75, 0.15)) +
    labs(x = "Age",
         y = "Proportion parent")

# Race 
df.p.race |> 
    left_join(
        df.race |> 
            mutate(age = as.character(age) |> as.numeric(),
                   year = as.character(year) |> as.numeric()
            ),
        by = c("age", "year", "race_eth")
    ) |> 
    ggplot(aes(x = age,
               col = race_eth)) +
    facet_grid(race_eth ~ year) +
    geom_line(aes(y = median)) +
    # geom_ribbon(col = NA) +
    geom_point(aes(y = y / n), 
               alpha = .6) +
    theme_bw() +
    theme(legend.position = "top",
          legend.title = element_blank()) +
    labs(x = "Age",
         y = "Proportion parent")

# State
## Join fit and data
df.p.state <- 
    df.p.state |> 
    left_join(
        df.state |> 
            mutate(
                age = as.character(age) |> as.numeric(),
                year = as.character(year) |> as.numeric()
            ) |> 
            dplyr::select(age, year, name, p),
        by = c("age", "year", "name")
    )

## Store state plots in list
fig.fit.state <- lapply(sample.states, function(x) {
    
    
    fig.out <- df.p.state |> 
        filter(name == x) |> 
        ggplot(aes(x = age)) +
        facet_wrap(~ year,
                   scales = "free_y") +
        geom_line(aes(y = median), 
                  linewidth = 1,
                  col = "red4") +
        geom_ribbon(aes(ymin = lower95, ymax = upper95),
                    alpha = .3, col = NA, fill = "red4") +
        geom_point(aes(y = p),
                   col = "skyblue3",
                   alpha = .6) +
        theme_bw() +
        labs(title = x)
})

## Save list of plots
ggsave(
    filename = here("plots", "p0_state_jano_fit.pdf"), 
    plot = gridExtra::marrangeGrob(fig.fit.state, nrow=1, ncol=1), 
    width = 15, height = 9
)



## === JANOSCHEK & FX ==========================================================

# Derivative of Janoschek fct
get_r <- function(U, L, k, c, x) {
    
    r <- k * c * (U - L) * x^(c - 1) * exp(-k * x^c)
}
fx <- get_r(pars[10,1], pars[10,2], pars[10,3], pars[10,4], 0:29)
plot(x = 15:44, y = fx, type = "l") 

# Normalise to sum to 1
fx <- fx/sum(fx)

# Multiply by TFR in 2010
fx <- fx * 1.73

# Not conclusive and TFR refers to all children
# Modify by:
# Using STAN pars instead of optim
# What to replace TFR with?
    

