
##==============================================================================
##
##  Project title: Modeling of childlessness
##
##  This code estimate proportion of parent (1-childless) with a Janoshek
##  curve using both CPS and NSFG data.
##
## 
##  Author: Benjamin Schlüter
##  Date: February 2024
##==============================================================================
##
##  Notes
## ------
## 1. 
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

# Prop of mother at US level
df.us <- readRDS(
    here(
        "data", 
        "df_us.rds"
    )
) 

# Prop of mother at RACE level
df.race <- readRDS(
    here(
        "data",
        "df_race.rds"
    )
)



## ===== KEY DIMENSIONS IN CODE ================================================

ages <- unique(df.us$age)
n.ages <- length(ages)
years <- unique(df.us$year)
years.all <- min(years):max(years)
n.years <- length(years)
n.years.all <- length(min(years):max(years))



## ===== INITIAL PARS VALUES FOR STAN ==========================================

# Janoscheck growth model
get_jano <- function(U, L, k, c, x) {
    
    y <- U - (U - L) * exp(-k * x^c)
}

# Get plausible initial pars value by trial error
# for mid-period year: 2016 (relative to 2014,
# it has the advantage to go until age 50 yo)
plot(x = 0:35, y = df.us$p[df.us$year == 2016 & df.us$survey=="CPS"])
lines(x = 0:35, y = get_jano(0.85, 0.01, 0.009, 1.85, 0:35), type = "l")

par.init <- c(0.85, 0.01, 0.009, 1.85)

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

# Check fit of obtained pars on CPS only
for (y in years) {
    
    if (y %in% c(2008, 2010)) {
        x <- 0:29
    } else {x <- 0:35}
    
    i <- which(years == y)
    
    plot(x = x, y = df.us$p[df.us$year == y & df.us$survey=="CPS"], main = y)
    lines(x = x, y = get_jano(pars[i, 1], pars[i, 2], pars[i, 3], pars[i, 4], x), type = "l")
}



## ===== FIT IN STAN ===========================================================

# STAN set-up ------------------------------------------------------------------

niter = 2000
nwarmup = niter/2
nchains = 4
nsim = (niter-nwarmup)*nchains
options(mc.cores = parallel::detectCores()-1)
pars <- c("mu", "U", "L", "k", "c")

# Inits ------------------------------------------------------------------------

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
        U = matrix( runif(n.years.all*n.groups, range.pars[1, 1], range.pars[2, 1]), n.years.all, n.groups),
        L = matrix( runif(n.years.all*n.groups, range.pars[1, 2], range.pars[2, 2]), n.years.all, n.groups),
        k = matrix( runif(n.years.all*n.groups, range.pars[1, 3], range.pars[2, 3]), n.years.all, n.groups),
        c = matrix( runif(n.years.all*n.groups, range.pars[1, 4], range.pars[2, 4]), n.years.all, n.groups),
        sigma_U = runif(1, var.pars[1]/10, var.pars[1]*10),
        sigma_L = runif(1, var.pars[2]/10, var.pars[2]*10),
        sigma_k = runif(1, var.pars[3]/10, var.pars[3]*10),
        sigma_c = runif(1, var.pars[4]/10, var.pars[4]*10)
    )
}



## ===== US-level ==============================================================

# Create indexing for observed data,
# order, and convert SE to the right scale (logit)
df.us <-
    df.us |> 
    mutate(
        # Create indexing var used in STAN to vectorize for 
        # survey, group (race or state), year, and age
        survey_i = as.factor(survey) |> as.numeric(),
        # At the US level, we only consider one group
        group_i = 1,
        # Year index starts at 1
        year_i = year - 2007,
        # Age index starts at 1
        age_i = age - 14,
        # Allows to index observed data 
        obs_i = 1,
        # Errors of survey estimate
        # need to be on the logit scale.
        # There are 6 observations where se == 0.
        # Use of delta method of logit transformation
        logit_se = ifelse( se != 0, se / (p * (1 - p)), 0)
    ) |> 
    arrange(
        # !!! Order super important for STAN !!!
        # Has to match with indexing vars
        survey, group_i, year, age
    )



## DEBUGGIN ====================================================================
df.us |> 
    ggplot(aes(x = age, 
               y = log((y/n) / (1 - (y/n))), 
               col = survey, 
               group = survey)) +
    facet_wrap(~ year) +
    geom_point() +
    geom_errorbar(aes(ymin = log((y/n) / (1 - (y/n)))-1.96*logit_se,
                      ymax = log((y/n) / (1 - (y/n)))+1.96*logit_se)) +
    theme_bw()

# Prob of 1 at age 45 or 50 in NSFG !!! Outliers
df.us |>
    filter(
        survey == "NSFG",
        year == 2016,
        age == 50
    )


df.us |> 
    summarise(
        .by = c(year, survey),
        
        N = sum(n)
    )
## =============================================================================
    

# Create indexing for estimated combination i (gp, year, age) 
df.i <- expand.grid(
    # !!!! Same order as order imposed above !!!!
    age_i = 1:n.ages,
    year_i = 1:n.years.all,
    group_i = 1:n.groups
    ) |> 
    mutate(
        # Create index for all combination of 
        # group, year, and t 
        i = row_number()
    )

# Index of each observed data in each survey
obs_i <- lapply(c("CPS", "NSFG"), function(s) {
    
    df.i |> 
        left_join(
            df.us |> 
                filter(
                    survey == s
                ) |> 
                dplyr::select(group_i, year_i, age_i, obs_i),
            by = c("group_i", "year_i", "age_i")
        ) |> 
        filter(
            !is.na(obs_i)
        ) |> 
        pull(i) 
})
obs_i <- unlist(obs_i)

# STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.us)[1],
    T_all = n.years.all,
    A_all = n.ages,
    age = 0:35,
    G = n.groups,
    # Indexing 
    i = obs_i,
    # Careful !!
    # Structured as vector but order super important !!!
    # Needs to coincide with logmx_ord in STAN (multiple to_vector() )
    y = round(df.us$y, 0), 
    n = round(df.us$n, 0),
    se = df.us$logit_se
)

fit = stan(here("code", "stan", "p0_jano_vec.stan"),
           data = stan_data,
           pars  = pars,
           init = stanInit,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains
)

# Store p 
p.fit <- array(as.matrix(fit, "mu"),
               c(nsim, n.ages, n.years.all),
               dimnames = list(1:nsim,
                               ages,
                               years.all
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

# Store estimates
saveRDS(df.p.us,
        here(
            "data", 
            "df_p0_us_jano_fit.rds"
        )
)

# # Load estimates
# df.p.us <- readRDS(
#     here(
#         "data",
#         "df_p0_us_jano_fit.rds"
#         )
#     )

# Store pars 
pars.us <- lapply(pars[pars != "mu"], function(x) {
    
    # Store in arrays
    par.draw <- array(as.matrix(fit, x),
                      c(nsim, n.years.all),
                      dimnames = list(1:nsim,
                                      years.all))
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
df.pars.us <- do.call("rbind", pars.us)

# Store estimates
saveRDS(df.pars.us,
        here(
            "data", 
            "df_pars_us_jano_fit.rds"
        )
)

# Load estimates
# df.pars.us <- readRDS(
#     here(
#         "data",
#         "df_pars_us_jano_fit.rds"
#     )
# )





## ===== Race-level ============================================================

# Create indexing for observed data,
# order, and convert SE to the right scale (logit)
df.race <-
    df.race |> 
    mutate(
        # Create indexing var used in STAN to vectorize for 
        # survey, group (race or state), year, and age
        survey_i = as.factor(survey) |> as.numeric(),
        # At the US level, we only consider one group
        group_i = as.factor(race_eth) |> as.numeric(),
        # Year index starts at 1
        year_i = year - 2007,
        # Age index starts at 1
        age_i = age - 14,
        # Allows to index observed data 
        obs_i = 1,
        # Errors of survey estimate
        # need to be on the logit scale.
        # There are 6 observations where se == 0.
        # Use of delta method of logit transformation
        
        
        # SHOULDN'T BE ANY SE ==0
        
        
        logit_se = ifelse( se != 0, se / (p * (1 - p)), 0)
    ) |> 
    arrange(
        # !!! Order super important for STAN !!!
        # Has to match with indexing vars
        survey, group_i, year, age
    )



## DEBUGGIN ====================================================================
df.race |> 
    ggplot(aes(x = age, 
               y = log((y/n) / (1 - (y/n))), 
               col = survey, 
               group = survey)) +
    facet_grid(race_eth ~ year) +
    geom_point() +
    geom_errorbar(aes(ymin = log((y/n) / (1 - (y/n)))-1.96*logit_se,
                      ymax = log((y/n) / (1 - (y/n)))+1.96*logit_se)) +
    theme_bw()

df.race |> 
    ggplot(aes(x = age, 
               y = y/n, 
               col = survey, 
               group = survey)) +
    facet_grid(race_eth ~ year) +
    geom_point() +
    theme_bw()

df.race |> filter(se == 0)
# se==0 would have a lot of influence in model !!!

## =============================================================================

n.groups <- length(unique(df.race$race_eth))
races <- unique(df.race$race_eth)

# Create indexing for estimated combination i (gp, year, age) 
df.i <- expand.grid(
    # !!!! Same order as order imposed above !!!!
    age_i = 1:n.ages,
    year_i = 1:n.years.all,
    group_i = 1:n.groups
) |> 
    mutate(
        # Create index for all combination of 
        # group, year, and t 
        i = row_number()
    )

# Index of each observed data in each survey
obs_i <- lapply(c("CPS", "NSFG"), function(s) {
    
    df.i |> 
        left_join(
            df.race |> 
                filter(
                    survey == s
                ) |> 
                dplyr::select(group_i, year_i, age_i, obs_i),
            by = c("group_i", "year_i", "age_i")
        ) |> 
        filter(
            !is.na(obs_i)
        ) |> 
        pull(i) 
})
obs_i <- unlist(obs_i)

# STAN data
stan_data <- list(
    # Dimensions
    N = dim(df.race)[1],
    T_all = n.years.all,
    A_all = n.ages,
    age = 0:35,
    G = n.groups,
    # Indexing 
    i = obs_i,
    # Careful !!
    # Structured as vector but order super important !!!
    # Needs to coincide with logmx_ord in STAN (multiple to_vector() )
    y = round(df.race$y, 0), 
    n = round(df.race$n, 0),
    se = df.race$logit_se
)

fit = stan(here("code", "stan", "p0_jano_vec.stan"),
           data = stan_data,
           pars  = pars,
           init = stanInit,
           include = TRUE,
           iter = niter,
           warmup = nwarmup,
           chains = nchains
)


# Store p 
p.fit <- array(as.matrix(fit, "mu"),
               c(nsim, n.groups, n.ages, n.years.all),
               dimnames = list(1:nsim,
                               races,
                               ages,
                               years.all
               ))
## Quantile of p
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

# # Load estimates
# df.p.us <- readRDS(
#     here(
#         "data",
#         "df_p0_us_jano_fit.rds"
#         )
#     )

# Store pars 
pars.us <- lapply(pars[pars != "mu"], function(x) {
    
    # Store in arrays
    par.draw <- array(as.matrix(fit, x),
                      c(nsim, n.years.all),
                      dimnames = list(1:nsim,
                                      years.all))
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
df.pars.us <- do.call("rbind", pars.us)

# Store estimates
saveRDS(df.pars.us,
        here(
            "data", 
            "df_pars_us_jano_fit.rds"
        )
)

# Load estimates
# df.pars.us <- readRDS(
#     here(
#         "data",
#         "df_pars_us_jano_fit.rds"
#     )
# )





## ===== VISUALIZATION =========================================================

# p - US
expand.grid(
    # !!!! Same order as order imposed above !!!!
    age = 15:49,
    year = 2008:2022,
    survey = c("CPS", "NSFG")
    ) |> 
    left_join(
        df.us,
        by = c("year", "age", "survey")
        
    ) |> 
    left_join(
        df.p.us,
        by = c("year", "age")
    ) |>
    ggplot(aes(x = age,
               ymin = lower95,
               ymax = upper95,
               group = survey,
               col = survey)) +
    facet_wrap( ~ year) +
    geom_ribbon(col = NA,
                fill = "red4",
                alpha = .3) + 
    geom_line(aes(y = median),
              col = "red4",
              linewidth = 1) +
    geom_point(aes(y = p), 
               alpha = .6) +
    theme_bw() +
    labs(x = "Age",
         y = "Proportion mothers")


# p - race
expand.grid(
    # !!!! Same order as order imposed above !!!!
    age = 15:49,
    year = 2008:2022,
    survey = c("CPS", "NSFG"),
    race_eth = races
) |> 
    left_join(
        df.race,
        by = c("year", "age", "survey", "race_eth")
        
    ) |> 
    left_join(
        df.p.race,
        by = c("year", "age", "race_eth")
    ) |>
    ggplot(aes(x = age,
               ymin = lower95,
               ymax = upper95,
               group = survey,
               col = survey)) +
    facet_grid(race_eth ~ year) +
    geom_ribbon(col = NA,
                fill = "red4",
                alpha = .3) + 
    geom_line(aes(y = median),
              col = "red4",
              linewidth = 1) +
    geom_point(aes(y = p), 
               alpha = .6) +
    theme_bw() +
    labs(x = "Age",
         y = "Proportion mothers")

# p - race: comparison of races
expand.grid(
    # !!!! Same order as order imposed above !!!!
    age = 15:49,
    year = 2008:2022,
    survey = c("CPS", "NSFG"),
    race_eth = races
) |> 
    left_join(
        df.race,
        by = c("year", "age", "survey", "race_eth")
        
    ) |> 
    left_join(
        df.p.race,
        by = c("year", "age", "race_eth")
    ) |>
    ggplot(aes(x = age,
               ymin = lower95,
               ymax = upper95,
               group = race_eth,
               col = race_eth)) +
    facet_wrap( ~ year) +
    geom_line(aes(y = median),
              linewidth = 1) +
    theme_bw() +
    labs(x = "Age",
         y = "Proportion mothers")

# p - race: race over time
expand.grid(
    # !!!! Same order as order imposed above !!!!
    age = 15:49,
    year = 2008:2022,
    survey = c("CPS", "NSFG"),
    race_eth = races
) |> 
    left_join(
        df.race,
        by = c("year", "age", "survey", "race_eth")
        
    ) |> 
    left_join(
        df.p.race,
        by = c("year", "age", "race_eth")
    ) |>
    ggplot(aes(x = age,
               ymin = lower95,
               ymax = upper95,
               group = year,
               col = year)) +
    facet_wrap( ~ race_eth) +
    geom_line(aes(y = median),
              linewidth = 1) +
    theme_bw() +
    labs(x = "Age",
         y = "Proportion mothers")

