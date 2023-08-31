
rm(list=ls())

library(tidyverse)
library(here)
library(rstan)



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
    




## EDA -------------------------------------------------------------------------

df.p.draws |> 
    ggplot(aes(x = age, y = p, col = gp)) +
    geom_line() +
    theme_bw()
