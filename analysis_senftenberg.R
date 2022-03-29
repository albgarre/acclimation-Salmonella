
library(tidyverse)
library(readxl)
library(cowplot)
library(stringr)
library(deSolve)
library(FME)
library(wesanderson)

## Fit the isothermal data

source("isothermal_model.R")

## Load data

sheet_names <- excel_sheets("./data/Salmonella SFT Heat resistance_normal pH.xlsx")

my_data <- sheet_names[!grepl("Isothermal", sheet_names)] %>%
    set_names(., .) %>%
    map(., ~ read_excel("./data/Salmonella SFT Heat resistance_normal pH.xlsx", sheet = .)) %>%
    map(.,
        ~ mutate(., temperature = approxfun(.$time, .$temperature, rule = 2)(.$time))
        )


# sheet_names <- excel_sheets("./data/for_R.xlsx")
# 
# my_data <- sheet_names[grepl("Senftenberg", sheet_names)] %>%
#     set_names(., .) %>%
#     map(., ~ read_excel("./data/for_R.xlsx", sheet = .))


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
        ) %>%
    map(.,
        ~ mutate(.,
                 slope = max(logN)/(max(temperature) - min(temperature)),
                 intercept =  (min(temperature) * (max(logN)-0)/(max(temperature) - min(temperature)) - 0),
                 fake_temp = temperature * slope - intercept
                 )
        ) 

## Predictions by Peleg

peleg_predictions <- my_data %>%
    map(., 
        ~ predict_inactivation("Peleg", seq(0, max(.$time), length = 1000),
                               c(coef(peleg_Senf$nls), temp_ref = 60, logN0 = .$logN0[1]),
                               .
                               ) 
        ) %>%
    map(.,
        ~ .$simulation
    ) %>%
    map(.,
        ~ select(., time, logS)
    ) %>%
    map(.,
        ~ mutate(., model = "Peleg")
        )

## Predictions by Mafart

mafart_predictions <- my_data %>%
    map(., 
        ~ predict_inactivation("Mafart", seq(0, max(.$time), length = 1000),
                               c(coef(mafart_Senf$nls), temp_ref = 60, logN0 = .$logN0[1]),
                               .
        ) 
    ) %>%
    map(.,
        ~ .$simulation
        ) %>%
    map(.,
        ~ select(., time, logS)
        ) %>%
    map(.,
        ~ mutate(., model = "Mafart")
        )

## Predictions by mod-Mafart

mod_Mafart <- function(Time, State, Parms, get_temp) {  # With the simplification
    
    with(as.list(c(State, Parms)), {
        
        temperature <- get_temp(Time)
        logdelta <- log10(delta_ref) - (temperature - ref_temp)/z
        delta <- 10^logdelta
        
        time_star <- (-logS*delta^p)^(1/p)
        dlogS <- -1/delta^p * p * time_star^(p - 1)
        
        return(list(dlogS))
    })
}

modMafart_predictions <- my_data %>%
    map(.,
        ~ ode(y = c(logS = 0-1e-6),
              times = seq(1e-6, max(.$time), length = 1000),
              func = mod_Mafart,
              parms = c(coef(mafart_Senf$nls), ref_temp = 60),
              get_temp = approxfun(.$time, .$temperature, rule = 2)
        )
        ) %>%
    # map(.,
    #     ~ ode(y = c(logS = 0-1e-6),
    #           times = seq(1e-6, 3, length = 1000),
    #           func = mod_Mafart,
    #           parms = c(coef(mafart_Senf$nls), ref_temp = 60),
    #           get_temp = .
    #           )
    #     ) %>%
    map(as.data.frame) %>%
    map(.,
        ~ mutate(., model = "mod-Mafart")
        )

## Figure 2


map2(mafart_predictions, modMafart_predictions,
     bind_rows
     ) %>%
    map2(., peleg_predictions,
         bind_rows
         ) %>%
    map2(., map(my_data, ~ head(., 1)$logN0),
         ~ mutate(.x, logN = logS + .y)
    ) %>%
    map2(., my_data,
        ~ ggplot(.x) + 
            geom_line(aes(x = time, y = logN, colour = model, 
                          linetype = model),
                      size = 1) +
            geom_line(aes(x = time, y = fake_temp),
                      data = .y, inherit.aes = FALSE) +
            geom_point(aes(x = time, y = logN), data = .y,
                       inherit.aes = FALSE) +
            coord_cartesian(ylim = c(0, 6)) +
            theme_bw(base_size = 14) +
            # theme(legend.position = "none") +
            xlab("Treatment time (min)") +
            ylab("Microbial concentration (log CFU/g)") +
            scale_linetype_manual(values = 2:4) +
            theme(legend.position = "none") +
            scale_color_manual(values = wes_palette("Darjeeling1", 3, type = "discrete"))
        ) %>%
    plot_grid(plotlist = ., # ncol=1,
              labels = "AUTO",
              ncol = 2)

## Look at the residuals

map2(mafart_predictions, map(my_data, ~ head(., 1)$logN0),
     ~ mutate(.x, logN = logS + .y)
     ) %>%
    map2(., my_data,
         ~ modCost(model = as.data.frame(select(.x, time, logN)),
                   obs = as.data.frame(select(.y, time, logN)))
         ) %>%
    map(., ~ .$residuals) %>%
    map(., ~ mutate(., res2 = res^2)) %>%
    map(., ~ summarize(., 
                       ME = mean(res),
                       MSE = mean(res2),
                       RMSE = sqrt(MSE))) %>%
    imap_dfr(., ~ mutate(.x, cond = .y))

map2(peleg_predictions, map(my_data, ~ head(., 1)$logN0),
     ~ mutate(.x, logN = logS + .y)
) %>%
    map2(., my_data,
         ~ modCost(model = as.data.frame(select(.x, time, logN)),
                   obs = as.data.frame(select(.y, time, logN)))
    ) %>%
    map(., ~ .$residuals) %>%
    map(., ~ mutate(., res2 = res^2)) %>%
    map(., ~ summarize(., 
                       ME = mean(res),
                       MSE = mean(res2),
                       RMSE = sqrt(MSE))) %>%
    imap_dfr(., ~ mutate(.x, cond = .y))

map2(modMafart_predictions, map(my_data, ~ head(., 1)$logN0),
     ~ mutate(.x, logN = logS + .y)
) %>%
    map2(., my_data,
         ~ modCost(model = as.data.frame(select(.x, time, logN)),
                   obs = as.data.frame(select(.y, time, logN)))
    ) %>%
    map(., ~ .$residuals) %>%
    map(., ~ mutate(., res2 = res^2)) %>%
    map(., ~ summarize(., 
                       ME = mean(res),
                       MSE = mean(res2),
                       RMSE = sqrt(MSE))) %>%
    imap_dfr(., ~ mutate(.x, cond = .y))


