
rm(list=ls())

library(tidyverse)
library(here)
library(rstan)
library(splines)



## Load data -------------------------------------------------------------------

## Childless by race and state
df <- readRDS(here("data", "df_chldness.rda"))



## Simulate data ---------------------------------------------------------------

## Age range
ages <- 15:45

## Logistic function
compute_logistic <- function(L, U, k, x0) {
    
    p = L + (U - L) / (1 + exp(k * (ages - x0)))
    
    return(p)
}

## Mean value
mu_L <- 0.175
mu_U <- 1
mu_k <- 0.4
mu_x0 <- 23

## Generate 4 groups (= races)
n.gp <- 4

df.pars.draws <- tibble(L = rnorm(n.gp, mu_L, 0.03),
                        U = 1.03,
                        k = rnorm(n.gp, mu_k, 0.05),
                        x0 = rnorm(n.gp, mu_x0, 1)) 

list.p.draws = pmap(df.pars.draws, compute_logistic)
df.p.draws = do.call("rbind", list.p.draws) |> 
    t() |> 
    as.data.frame() |> 
    mutate(
        age = ages
    ) |> 
    pivot_longer(V1:V4, names_to = "gp", values_to = "p")



## Simulate Hierarchical data -------------------------------------------------

## Logistic function with data frame as output
compute_logistic_df <- function(L, k, x0) {
    
    p = L + (1 - L) / (1 + exp(k * (ages - x0)))
    
    out <- tibble(
        p0 = p,
        age = ages, 
        L = L, 
        k = k, 
        x0 = x0
    )
    
    return(out)
}

## 31 ages and 31 years

## simulate a random walk for x0

    ## generate values around it means for each year

## k and L have each a mean around a global mean
    


## Assess impact of spline around logistic function ----------------------------

## Check rows where n=0

df |> 
    filter(
        n < 1
        )

## Logistic function
compute_logistic <- function(L, k, x0) {
    
    p = L + (1 - L) / (1 + exp(k * (ages - x0)))
    
    return(p)
}

## Logistic with values
theta <- compute_logistic(0.175, 0.4, 23)

## Logistic with splines
position.knots <- seq(20, 40, 10)
B   <- bs(ages, knots=position.knots, degree=1)
alpha <- c(1,1,-1.5,1)

omega <- qlogis(theta) + B %*% alpha

p <- 1 / (1 + exp(-omega))

df <- tibble(
    fit = c(theta, p),
    type = c(rep("logistic", length(ages)), 
             rep("logistic+splines", length(ages))),
    age = rep(ages, 2)
)

## Visu
df |> 
    ggplot(aes(x = age, y = fit, col = type, group = type)) +
    geom_line() +
    geom_vline(xintercept = position.knots,
               linetype = "dashed") +
    theme_bw()

## Visu splines
df.splines <- tibble(
    splines = B %*% alpha,
    age = ages
)
df.splines |> 
    ggplot(aes(x = age, y = splines)) +
    geom_line() +
    theme_bw()

## EDA -------------------------------------------------------------------------

df.p.draws |> 
    ggplot(aes(x = age, y = p, col = gp)) +
    geom_line() +
    theme_bw()
