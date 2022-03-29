
## Load libraries

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(bioinactivation)

## Functions for deSolve

two_compartment <- function(Time, state, pars, get_temp, get_k) {
    
    pars <- as.list(pars)
    state <- as.list(state)
    
    temp <- get_temp(Time)
    k <- get_k(pars, temp)
    
    dP <- k * (1 - (state$P/pars$Pmax)^pars$m)
    
    normal_D <- 10^( log10(pars$D_R) - (temp - pars$temp_ref)/pars$z )
    new_D <- 1 + pars$C * state$P
    
    dN <- -log(10)/normal_D/new_D * state$N
    
    return(list(c(dP, dN)))
    
}

two_compartment_tail <- function(Time, state, pars, get_temp, get_k) {
    
    pars <- as.list(pars)
    state <- as.list(state)
    
    temp <- get_temp(Time)
    k <- get_k(pars, temp)
    
    dP <- k * (1 - (state$P/pars$Pmax)^pars$m)
    
    normal_D <- 10^( log10(pars$D_R) - (temp - pars$temp_ref)/pars$z )
    new_D <- 1 + pars$C * state$P
    
    dN <- -log(10)/normal_D/new_D * state$N * (1 - pars$Nmin/state$N)
    
    return(list(c(dP, dN)))
    
}

## Models for k

my_k <- function(pars, temp) {
    
    pars <- as.list(pars)
    
    k <- ifelse(temp <= pars$temp_min, 
                0, 
                pars$A*exp(-pars$E/(temp - pars$temp_min))
                )
    
}

step_k <- function(pars, temp) {
    
    pars <- as.list(pars)
    
    k <- ifelse(temp <= pars$temp_min, 0, pars$kmax)
    
    k
}

linear_k <- function(pars, temp) {
    
    pars <- as.list(pars)
    
    k <- ifelse(temp <= pars$temp_min, 0, (temp - pars$temp_min)*pars$A)
    
    k
    
}

exp_k <- function(pars, temp) {
    
    pars <- as.list(pars)
    
    lnk <- pars$k0 + pars$A * temp
    
    exp(lnk)
    
}

# ## Functions for calculating predictions
# 
# make_prediction <- function(my_pars, yini, times, get_temp, get_k) {
#     
#     out   <- ode(yini, times, two_compartment, my_pars,
#                  get_temp = get_temp, get_k = get_k) %>%
#         as.data.frame()
#     
#     out
# }
# 
# calculate_residuals <- function(this_pars, this_data, known_pars,
#                                 yini, get_temp, get_k) {
#     
#     my_pars <- c(this_pars, known_pars)
#     times <- seq(0, max(this_data$time), length = 1000)
#     
#     my_prediction <- make_prediction(my_pars, yini, times, get_temp, get_k) %>%
#         mutate(logN = log10(N))
#     
#     my_cost <- modCost(model = select(my_prediction, time, logN),
#                        obs = select(this_data, time, logN))
#     
#     my_cost
# }



















