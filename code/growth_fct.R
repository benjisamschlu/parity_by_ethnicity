
##==============================================================================
##
##  Project title: Model childlessness
##
##  Explore different growth function
##
## 
##  Author: Benjamin Schlüter
##  Date: Oct. 2023
##==============================================================================
##
##  Notes
## ------
## 1. 
## 2.
## 
##==============================================================================

rm(list = ls())




## LOAD PACKAGES ---------------------------------------------------------------

# Install/load packages
packages <- c("tidyverse", "here")
for(p in packages){
    if(!require(p,character.only = TRUE)) install.packages(p)
    library(p,character.only = TRUE)
}



## FUNCTIONS -------------------------------------------------------------------


get_janoscheck <- function(U, L, k, p, x) {
    
    y <- U - (U - L) * exp(-k * x^p)
}

get_logistic <- function(U, L, k, x0, x) {
    
    y <- L + (U - L) / (1 + exp(k * (x - x0)))
}

get_logit <- function(p) {
    
    y <- log(p / (1 - p))
}

get_expit <- function(l) {
    
    y <- exp(l) / (1 + exp(l))
}

linear_transf <- function(a, b, x) {
    
    y = a + b * x
}



## DATA -------------------------------------------------------------------

ages <- 0:30

# Janoschek
plot(x = ages, y = get_janoscheck(1, 0.2, 1, 1.2, ages))

ages <- 15:45

# Logistic
fit <- get_logistic(1, 0.2, 0.8, 25, ages)

plot(x = ages, y = fit, type = "l")

logit_fit <- get_logit(fit)
a <- 0
b <- 0.5
trans_logit_fit <- linear_transf(a, b, logit_fit)
new_fit <- get_expit(trans_logit_fit)

plot(x = ages, y = new_fit, type = "l")



## ELEMENTS OF FUNCTION -----------------------------------------------------

# Janoscheck eq: U - (U - L) * exp(-k * t^p)
t <- 0:10
plot(x = x, y = exp(-t), type = "l",
     xlim = range(x),
     ylim = c(0,exp(max(-x))))

# think in terms of exp(-k * t^p) with t>0

# exp(-t) starts at 1 and goes to zero

# small k takes more time for exp(-k * t) to go towards zero
plot(x = x, y = exp(-0.5*t), type = "l",
     xlim = range(t),
     ylim = c(0,exp(max(-t))))
plot(x = x, y = exp(-2*t), type = "l",
     xlim = range(t),
     ylim = c(0,exp(max(-t))))
# p by making small number smaller and big bigger creates the sigmoid
plot(x = t, y = exp(-0.1*(t^1.5)), type = "l",
     xlim = range(t),
     ylim = c(0,exp(max(-t))))

t <- 0:30
plot(x = t, y = exp(-0.03*(t^1.5)), type = "l")
plot(x = t, y = get_janoscheck(0.98, 0.2, 0.01, 1.8, t), type = "l")



## JANOSCHECK ON DATA ------------------------------------------------------------------

# Load data
df <- readRDS(here("data", "df_chldness.rda")) |> 
    group_by(
        age, year,
        .drop = FALSE
    ) |> 
    summarise(
        across(y:n, sum)
    ) |> 
    ungroup() |> 
    filter(year == 2014) |> 
    mutate(p = 1 - y/n)

ages <- 0:29
plot(x = ages, y = df$p)
lines(x = ages, y = get_janoscheck(0.87, 0.005, 0.01, 1.8, ages), type = "l")

## SECTION 1 ===================================================================


## SECTION 2 ===================================================================



