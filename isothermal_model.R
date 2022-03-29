
library(tidyverse)
library(readxl)
library(bioinactivation)

## Load the data

data_Ent <- read_excel("./data/Salmonella Enteritidis Heat resistance_normal pH.xlsx", sheet = "S.enteritidis_isotherm")
data_Senf <- read_excel("./data/Salmonella SFT Heat resistance_normal pH.xlsx", sheet = "Senftenberg_Isothermal")

# data_Ent <- read_excel("./data/old/isothermal_forR.xlsx", sheet = "Enteritidis")
# data_Senf <- read_excel("./data/old/isothermal_forR.xlsx", sheet = "Senftenberg")

ggplot(data_Ent, aes(x = time, y = logN)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap("temperature", scales = "free")

ggplot(data_Senf, aes(x = time, y = logN)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap("temperature", scales = "free")

## Estimate initial counts

data_Ent <- data_Ent %>%
    group_by(temperature) %>%
    mutate(logN0 = mean(ifelse(time == 0, logN, NA) , na.rm = TRUE)) %>% 
    mutate(log_diff = logN - logN0) %>%
    rename(temp = temperature)

data_Ent %>%
    mutate(temp = paste0(temp, "ºC")) %>%
    ggplot(., aes(x = time, y = log_diff)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap("temp", scales = "free") +
    xlab("Treatment time (min)") +
    ylab("Microbial count (log CFU/g)") + 
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))


data_Senf <- data_Senf %>%
    group_by(time, temperature) %>%
    summarize(logN = mean(logN)) %>%
    ungroup() %>%
    group_by(temperature) %>%
    mutate(logN0 = mean(ifelse(time == 0, logN, NA) , na.rm = TRUE)) %>% 
    mutate(log_diff = logN - logN0) %>%
    rename(temp = temperature)

data_Senf %>%
    mutate(temp = paste0(temp, "ºC")) %>%
    ggplot(., aes(x = time, y = logN)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap("temp", scales = "free") +
    xlab("Treatment time (min)") +
    ylab("Microbial count (log CFU/g)") + 
    theme_bw() +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 16))

## Fit the model

# iso_Senf <- fit_isothermal_inactivation("Bigelow", data_Senf,
#                                           c(z = 5, D_R = 3),
#                                           list(temp_ref = 60))

iso_Ent <- fit_isothermal_inactivation("Bigelow", data_Ent,
                                         c(z = 5, D_R = .1),
                                         list(temp_ref = 57.5))

# summary(iso_Senf)
summary(iso_Ent)
sqrt(mean(residuals(iso_Ent$nls)^2))

# plot(iso_Senf, make_gg = TRUE)
# plot(iso_Ent, make_gg = TRUE)

## Mafart model

mafart_Senf <- fit_isothermal_inactivation("Mafart", data_Senf,
                                        c(z = 5, delta_ref = 3, p=1),
                                        list(temp_ref = 60))

mafart_Ent <- fit_isothermal_inactivation("Mafart", data_Ent,
                                       c(z = 5, delta_ref = .1, p=1),
                                       list(temp_ref = 57.5))

summary(mafart_Senf)
# plot(mafart_Senf, make_gg = TRUE)
sqrt(mean(residuals(mafart_Senf$nls)^2))

summary(mafart_Ent)
# plot(mafart_Ent, make_gg = TRUE)

## Peleg model

peleg_Senf <- fit_isothermal_inactivation("Peleg", 
                                        data_Senf,
                                        list(n=1, temp_crit = 59, k_b = 1),
                                        list()
                                        )

# iso_Ent <- fit_isothermal_inactivation("Mafart", data_Ent,
#                                        c(z = 5, delta_ref = .1, p=1),
#                                        list(temp_ref = 60))

summary(peleg_Senf)
sqrt(mean(residuals(peleg_Senf$nls)^2))
# plot(peleg_Senf, make_gg = TRUE)

# summary(iso_Ent)
# plot(iso_Ent, make_gg = TRUE)

## Figure 1

data_Senf %>%
    split(.$temp) %>%
    map(.,
        ~ tibble(time = seq(min(.$time), max(.$time), length = 1000),
                 temp = .$temp[1])
        ) %>%
    map(.,
        ~ mutate(., 
                 Peleg = predict(peleg_Senf$nls, newdata = .),
                 Mafart = predict(mafart_Senf$nls, newdata = .)
                 )
        ) %>%
    map(.,
        ~ pivot_longer(., c(-time, -temp), names_to = "model", values_to = "log_diff")
        ) %>%
    imap_dfr(., ~ mutate(.x, temp = paste0(.y, "ºC"))) %>%
    ggplot() +
    geom_line(aes(x = time, y = log_diff, colour = model, linetype = model),
              size = 1.5) +
    geom_point(aes(x = time, y = log_diff), 
               size = 3,
               data = mutate(data_Senf, temp = paste0(temp, "ºC")),
               inherit.aes = FALSE) +
    facet_wrap("temp", scales = "free") +
    theme_bw(base_size = 14) +
    xlab("Treatment time (min)") + ylab("Number of reductions (log CFU/g)") +
    scale_linetype_manual(values = 2:3)  +
    theme(legend.position = "none")

## Figure 3

data_Ent %>%
    split(.$temp) %>%
    map(.,
        ~ tibble(time = seq(min(.$time), max(.$time), length = 1000),
                 temp = .$temp[1])
    ) %>%
    map(.,
        ~ mutate(., 
                 Bigelow = predict(iso_Ent$nls, newdata = .),
                 Mafart = predict(mafart_Ent$nls, newdata = .)
        )
    ) %>%
    map(.,
        ~ pivot_longer(., c(-time, -temp), names_to = "model", values_to = "log_diff")
    ) %>%
    imap_dfr(., ~ mutate(.x, temp = paste0(.y, "ºC"))) %>%
    ggplot() +
    geom_line(aes(x = time, y = log_diff, colour = model, linetype = model),
              size = 1.5) +
    geom_point(aes(x = time, y = log_diff), 
               size = 3,
               data = mutate(data_Ent, temp = paste0(temp, "ºC")),
               inherit.aes = FALSE) +
    facet_wrap("temp", scales = "free") +
    theme_bw(base_size = 14) +
    xlab("Treatment time (min)") + ylab("Number of reductions (log CFU/g)") +
    scale_linetype_manual(values = 2:3)  +
    theme(legend.position = "none")





