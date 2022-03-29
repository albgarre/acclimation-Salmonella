
library(tidyverse)
library(readxl)
library(cowplot)
library(FME)
library(deSolve)

library(wesanderson)

## Fit the isothermal data

source("isothermal_model.R")
source("acclimation_models.R")

## Load data

sheet_names <- excel_sheets("./data/Salmonella Enteritidis Heat resistance_normal pH.xlsx")

my_data <- sheet_names[!grepl("isotherm", sheet_names)] %>%
    set_names(., .) %>%
    map(., ~ read_excel("./data/Salmonella Enteritidis Heat resistance_normal pH.xlsx", sheet = .))

# sheet_names <- excel_sheets("./data/old/for_R.xlsx")
# 
# my_data <- sheet_names[grepl("Enteritidis", sheet_names)] %>%
#     set_names(., .) %>%
#     map(., ~ read_excel("./data/old/for_R.xlsx", sheet = .))

## Visualize the data

ggplot(my_data[[1]]) +
    geom_point(aes(time, log10(N)))

ggplot(my_data[[2]]) +
    geom_point(aes(time, log10(N)))

ggplot(my_data[[3]]) +
    geom_point(aes(time, log10(N)))

ggplot(my_data[[4]]) +
    geom_point(aes(time, log10(N)))

ggplot(my_data[[5]]) +
    geom_point(aes(time, log10(N)))

## Estimate initial counts

my_data <- my_data %>%
    map(.,
        ~ mutate(., logN = log10(N))
    ) %>%
    map(., 
        ~ mutate(., logN0 = mean(ifelse(time ==  0, logN, NA) , na.rm = TRUE))
    )

## Compare predictions and observations

logN0_data <- my_data %>%
    map(.,
        ~ filter(., time == 0)
    ) %>%
    map_dfr(.,
            ~ summarize(., logN0 = mean(logN))
    ) %>%
    pull(logN0)

dyna_predictions <- my_data %>%
    map2(., logN0_data,
        ~ predict_inactivation("Bigelow", seq(0, max(.x$time), length = 1000),
                               c(coef(iso_Ent$nls), temp_ref = iso_Ent$parameters$temp_ref, logN0 = .y),
                               .x
        ) 
    )


## Table 2 (I)

dyna_predictions %>%
    map(., ~ .$simulation) %>%
    map(., ~ select(., time, logN)) %>%
    map2(., my_data,
         ~ modCost(
             model = as.data.frame(.x),
             obs = as.data.frame(select(.y, time, logN))
         )
         ) %>%
    map(., ~.$residuals) %>%
    map(., ~ summarize(., 
                       RMSE = sqrt(mean(res^2)),
                       ME = mean(res)
                       )
        ) %>%
    imap_dfr(., ~ mutate(.x, profile = .y))

##

aa <- map2(dyna_predictions, my_data, 
     ~ plot(.x, plot_temp = TRUE, label_y1 = "Microbial count (log CFU/g)",
            label_y2 = "Temperature (ºC)") +
         geom_point(aes(x = time, y = logN), data = .y) +
         # coord_cartesian(ylim = c(min(.y$logN), max(.y$logN))) +
         theme_bw() +
         theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 16),
               strip.text = element_text(size = 16)) +
         xlab("Treatment time (min)")
)

my_names <- str_split(names(aa), "_") %>%
    map_chr(.,
            ~ paste0(.[2], "ºC/min")
    )

aa %>%
    set_names(., my_names) %>%
    imap(.,
         ~ .x + ggtitle(.y) + coord_cartesian(ylim = c(0, 6.5))
    ) %>%
    plot_grid(plotlist = .)
##

my_data[c(2, 3, 4)] %>%
    map(.,
        ~ filter(., temperature >= 57.5)
        ) %>%
    map(.,
        ~ ggplot(.) + geom_point(aes(x =time, y = logN))
        )

# aa <- my_data[c(2, 3, 4)] %>%
#     map(.,
#         ~ filter(., temperature >= 57.5)
#     ) %>%
#     map(.,
#         ~ mutate(., time = time - min(time))
#         ) %>%
#     map(.,
#         ~ fit_inactivation_MCMC(experiment_data = as.data.frame(.),
#                                 "Bigelow",
#                                 temp_profile = .,
#                                 c(D_R = 20, logN0 = 3, z = 7),
#                                 c(D_R = 100, logN0 = 5, z = 12),
#                                 c(D_R = .1, logN0 = 1, z = 3),
#                                 c(temp_ref = 60),
#                                 niter = 3000
#         )
#     )
# 
# 
# aa %>% map(plot)
# aa %>% map(., ~ unlist(.$best_prediction$model_parameters))
# aa %>% map(., ~ pairs(.$modMCMC))

## Acclimation model

make_prediction <- function(my_pars, yini, times, get_temp) {
    
    out   <- ode(yini, times, two_compartment, my_pars,
                 get_temp = get_temp, get_k = my_k) %>%
        as.data.frame()
    
    out
}

calculate_residuals <- function(this_pars, this_data, known_pars,
                                yini, times, get_temp,
                                NA_value = -5) {
    
    my_pars <- c(this_pars, known_pars)
    
    my_prediction <- make_prediction(my_pars, yini, times, get_temp) %>%
        mutate(logN = log10(N)) %>%
        mutate(logN = ifelse(is.na(logN), NA_value, logN),  # So it actually calculates residuals
               logN = ifelse(logN < NA_value, NA_value, logN)
               )
    
    my_cost <- modCost(model = select(my_prediction, time, logN),
                       obs = select(this_data, time, logN))
    
    my_cost
}

## Guess the initial values

yini <- c(P = 0, N = 1e6)

iso_pars <- coef(iso_Ent$nls)
# iso_pars <- c(D_R = 10^(-1.239), z = 3)

known_pars <- c(temp_min = 37, Pmax = 1,
                m = 1,
                Nmin = 0,
                temp_ref = 57.5)

# adapt_pars <- c(A = 1e-1, E = 1, C = 6)  # lower
adapt_pars <- c(A = 1e0, E = 50, C = 6)  # start
# adapt_pars <- c(A = 1, E = 100, C = 6)  # upper

my_data %>%
    map(.,
        ~ list(temp = approxfun(.$time, .$temperature, rule = 2),
               max_t = max(.$time)*1.1
               )
        ) %>%
    map(.,
        ~ make_prediction(c(iso_pars, known_pars, adapt_pars),
                          yini,
                          seq(0, .$max_t, length = 1000),
                          .$temp
                          )
        ) %>%
    map(.,
        ~ mutate(., logN = log10(N))
        ) %>%
    map(.,
        ~ ggplot(.) + geom_line(aes(x = time, y = logN))
        ) %>%
    imap(.,
         ~ .x + ggtitle(.y)
         ) %>%
    map2(., my_data,
         ~ .x + geom_point(aes(x = time, y = logN), data = .y) + coord_cartesian(ylim = c(0, 6))
         ) %>%
    plot_grid(plotlist = .)

my_data %>%
    map(.,
        ~ list(temp = approxfun(.$time, .$temperature, rule = 2),
               max_t = max(.$time)
        )
    ) %>%
    map(.,
        ~ make_prediction(c(iso_pars, known_pars, adapt_pars),
                          yini,
                          seq(0, .$max_t, length = 1000),
                          .$temp
        )
    ) %>%
    map(.,
        ~ ggplot(.) + geom_line(aes(x = time, y = P*6))
    ) %>%
    imap(.,
         ~ .x + ggtitle(.y)
    ) %>%
    map2(., my_data,
         ~ .x + geom_point(aes(x = time, y = logN), data = .y)
    ) %>%
    plot_grid(plotlist = .)

## Fit the model

set.seed(12341)

# d <- as.data.frame(my_data$S.enter._Dynamic_1C_min)
# yini["N"] <- 10^logN0_data[2]

d <- as.data.frame(my_data$S.enter._Dynamic_0.5C_min)
yini["N"] <- 10^logN0_data[1]
    
my_temp <- approxfun(d$time, d$temperature, rule = 2)

fit_acclimation <- modMCMC(calculate_residuals,
                           adapt_pars,
                           lower = c(A = 1e-1, E = 1e0, C = 3),
                           # upper = c(A = 1e2, E = 1e2, C = 15),
                           upper = c(A = 20, E = 1e2, C = 15),
                           this_data = d, 
                           known_pars = c(iso_pars, known_pars),
                           yini = yini, 
                           times = seq(0, max(d$time), length = 100),
                           get_temp = my_temp,
                           # niter = 1000
                           niter = 5000,  updatecov = 500,
                           burninlength = 1000
                           )

my_MCMC <- mcmc(fit_acclimation$pars)
summary(my_MCMC)
plot(my_MCMC)
pairs(fit_acclimation)
# 
# acfplot(my_MCMC)
# autocorr.plot(my_MCMC)
# rejectionRate(my_MCMC)
# geweke.diag(my_MCMC)
# pnorm(geweke.diag(my_MCMC)$z)
# geweke.plot(my_MCMC)
# heidel.diag(my_MCMC)

make_prediction(c(fit_acclimation$bestpar, known_pars, iso_pars), yini, 
                seq(0, max(d$time), length = 100),
                my_temp) %>%
    mutate(logN = log10(N)) %>%
    ggplot() +
    geom_line(aes(x = time, y = logN)) +
    geom_point(aes(x = time, y = logN), data = d)

## Compare against data

my_data %>%
    map(.,
        ~ list(temp = approxfun(.$time, .$temperature, rule = 2),
               max_t = max(.$time)*1.1
        )
    ) %>%
    map(.,
        ~ make_prediction(c(iso_pars, known_pars, fit_acclimation$bestpar),
                          yini,
                          seq(0, .$max_t, length = 1000),
                          .$temp
        )
    ) %>%
    map(.,
        ~ mutate(., logN = log10(N))
    ) %>%
    map(.,
        ~ ggplot(.) + geom_line(aes(x = time, y = logN))
    ) %>%
    imap(.,
         ~ .x + ggtitle(.y)
    ) %>%
    map2(., my_data,
         ~ .x + geom_point(aes(x = time, y = logN), data = .y) + coord_cartesian(ylim = c(0, 6))
    ) %>%
    plot_grid(plotlist = .)

my_data %>%
    map(.,
        ~ list(temp = approxfun(.$time, .$temperature, rule = 2),
               max_t = max(.$time)
        )
    ) %>%
    map(.,
        ~ make_prediction(c(iso_pars, known_pars, fit_acclimation$bestpar),
                          yini,
                          seq(0, .$max_t, length = 1000),
                          .$temp
        )
    ) %>%
    map(.,
        ~ ggplot(.) + geom_line(aes(x = time, y = P*6))
    ) %>%
    imap(.,
         ~ .x + ggtitle(.y)
    ) %>%
    map2(., my_data,
         ~ .x + geom_point(aes(x = time, y = logN), data = .y)
    ) %>%
    plot_grid(plotlist = .)

## Compare both model predictions

lines_1 <- dyna_predictions %>%
    map(.,
        ~ plot(., plot_temp = TRUE, 
               label_y1 = "Microbial count (log CFU/g)",
               label_y2 = "Temperature (ºC)")
    )


# my_points <- my_data %>%
#     map(.,
#         ~ geom_point(aes(x = time, y = logN), data = .)
#         )
# 
# map2(lines_1, my_points, 
#      ~ .x + .y + ylim(c(0, 7))
#      )

iso_plots<- map2(dyna_predictions, my_data, 
                  ~ plot(.x, plot_temp = TRUE, label_y1 = "Microbial count (log CFU/g)",
                         label_y2 = "Temperature (ºC)",
                         ylims = c(0, max(.y$logN)+.2)
                         ) +
                      geom_point(aes(x = time, y = logN), data = .y) +
                     # ylim(c(0, max(.y$logN)+.2)) +
                      # coord_cartesian(ylim = c(min(.y$logN)-.2, max(.y$logN)+.2)) +
                      # theme_bw() +
                      theme(axis.text = element_text(size = 14),
                            axis.title = element_text(size = 16),
                            strip.text = element_text(size = 16)) +
                      xlab("Treatment time (min)")
                  )

my_lines <- my_data %>%
    map(.,
        ~ list(temp = approxfun(.$time, .$temperature, rule = 2),
               max_t = max(.$time)*1.1
        )
    ) %>%
    map2(., logN0_data,
        ~ make_prediction(c(iso_pars, known_pars, fit_acclimation$bestpar),
                          c(P = 0, N = 10^.y),
                          seq(0, .x$max_t, length = 1000),
                          .$temp
        )
    ) %>%
    map(.,
        ~ mutate(., logN = log10(N))
    ) %>%
    map(.,
        ~  geom_line(aes(x = time, y = logN), 
                     linetype = 3,
                     data =.)
    )

out_plots <- map2(iso_plots, my_lines,
     ~ .x + .y + theme_cowplot()
     )
## Figure 4

# plot_grid(plotlist = out_plots[c(1, 2, 4, 6)],
#           labels = "AUTO")

# out_plots[[1]]
out_plots[[1]] 
## Figure 5

# plot_grid(plotlist = out_plots[c(2, 4, 3, 6)],
#           labels = "AUTO")

plot_grid(plotlist = out_plots[c(2, 3, 6)],
          labels = "AUTO", 
          nrow = 1)


## Figure 6

# my_data[c(1, 2, 3, 6)] %>%
#     map(.,
#         ~ list(temp = approxfun(.$time, .$temperature, rule = 2),
#                max_t = max(.$time)
#         )
#     ) %>%
#     map(.,
#         ~ make_prediction(c(iso_pars, known_pars, fit_acclimation$bestpar),
#                           yini,
#                           seq(0, .$max_t, length = 1000),
#                           .$temp
#         )
#     ) %>%
#     map(.,
#         ~ ggplot(.) + geom_line(aes(x = time, 
#                                     # y = P*fit_acclimation$bestpar["C"])
#                                     y = P*6)
#                                 )
#     ) %>%
#     map2(., my_data[c(1, 2, 3, 6)],
#          ~ .x + geom_point(aes(x = time, y = logN), data = .y) + 
#              scale_y_continuous(name = "Microbial count (log CFU/ml)",
#                                 sec.axis = sec_axis(~./6*100, name = "Degree of adaptation (%)")) +
#              theme_cowplot()
#                                 # sec.axis = sec_axis(~.*100, name = "Increase in D-value (%)"))
#     ) %>%
#     plot_grid(plotlist = ., labels = "AUTO")

my_data[c(1, 2, 3, 6)] %>%
    map(.,
        ~ list(temp = approxfun(.$time, .$temperature, rule = 2),
               max_t = max(.$time)
        )
    ) %>%
    map(.,
        ~ make_prediction(c(iso_pars, known_pars, fit_acclimation$bestpar),
                          yini,
                          seq(0, .$max_t, length = 1000),
                          .$temp
        )
    ) %>%
    map(.,
        ~ ggplot(.) + geom_line(aes(x = time, 
                                    # y = P*fit_acclimation$bestpar["C"])
                                    y = P*6)
        )
    ) %>%
    map2(., my_data[c(1, 2, 3, 6)],
         ~ .x + geom_point(aes(x = time, y = logN), data = .y) + 
             scale_y_continuous(name = "Microbial count (log CFU/ml)",
                                limits = c(0, 6.5),
                                sec.axis = sec_axis(~./6*100, name = "Degree of acclimation (%)")) +
             theme_cowplot()
         # sec.axis = sec_axis(~.*100, name = "Increase in D-value (%)"))
    ) %>%
    map2(., dyna_predictions[c(1, 2, 3, 6)],
         
         ~ .x + geom_line(aes(x = time, y = logN),
                          data = .y$simulation,
                          linetype = 2
                          ) 
         ) %>%
    plot_grid(plotlist = ., labels = "AUTO")


## Table 2 (II)

my_data[c(1, 2, 4, 3, 6)] %>%
    map(.,
        ~ list(temp = approxfun(.$time, .$temperature, rule = 2),
               max_t = max(.$time)
        )
    ) %>%
    map(.,
        ~ make_prediction(c(iso_pars, known_pars, fit_acclimation$bestpar),
                          yini,
                          seq(0, .$max_t, length = 1000),
                          .$temp
        )
    ) %>%
    map(., as.data.frame) %>%
    map(., ~ mutate(., logN = log10(N))) %>%
    map(., ~ select(., time, logN)) %>%
    map2(., my_data[c(1, 2, 4, 3, 6)],
         ~ modCost(
             model = .x,
             obs = as.data.frame(select(.y, time, logN))
             )
         ) %>%
    map(.,
        ~ .$residuals
        ) %>%
    map(.,
        ~ summarize(., ME = mean(res), RMSE = sqrt(mean(res^2)))
        ) %>%
    imap_dfr(., ~ mutate(.x, profile = .y))

my_data %>%
    map(.,
        ~ lm(temperature ~ time, data = .)
        )



