
##==============================================================================
##
##  Project title: Model childlessness
##
##  Explore Janosheck function
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


get_janoscheck <- function(U, beta, b, c, x) {
    
    y <- U * (1 - beta * exp(-b * x^c))
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
plot(x = ages, y = get_janoscheck(1, 1, 0.8, 0.2, ages))

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



## KEY DIMENSIONS OF CODE -----------------------------------------------------



## CLEAN DATA ------------------------------------------------------------------



## SECTION 1 ===================================================================


## SECTION 2 ===================================================================



