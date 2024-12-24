
# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(ggsci)
library(gridExtra)
library(showtext)

font_add('Helvetica', 
         "./script/fonts/Helvetica.ttf",
         bold = "./script/fonts/Helvetica-Bold.ttf")

remove(list = ls())

load('./data/sars_2_cov.Rdata')

set.seed(202202)

# function ----------------------------------------------------------------

library(epitrix)

expose_date <- function(a, b){
     return(as.Date(seq.Date(a, b, by = 'day')))
}

expose_seq <- function(x, data){
     infector_expose_dates <- seq.Date(from = data[x, 'infector_exposedate1'],
                                       to = data[x, 'infector_exposedate2'],
                                       by = 'day')
     infectee_expose_dates <- seq.Date(from = data[x, 'infectee_exposedate1'],
                                       to = data[x, 'infectee_exposedate2'],
                                       by = 'day')
     datafile <- as.data.frame(expand.grid(infectee_expose_dates, infector_expose_dates)) %>% 
          mutate(median = as.numeric(Var1 - Var2),
                 median = ifelse(median<0, NA, median))
     gt_median <- median(datafile$median, na.rm = T)
     gt_max <- max(infectee_expose_dates) - min(infector_expose_dates)
     gt_min <- min(infectee_expose_dates) - max(infector_expose_dates)
     gt_min <- ifelse(gt_min<0, 0, gt_min)
     return(c(gt_min, gt_median, gt_max))
}

fit_best <- function(data, distribution.type = NULL){
     data[data == 0] <- 0.1
     fit.gamma <- try(fitdistr(data, "gamma"))
     if (class(fit.gamma) == "try-error") {
          fit.gamma$loglik <- NA
     }
     fit.weib <- try(fitdistr(data, "weibull"))
     if (class(fit.weib) == "try-error") {
          fit.weib$loglik <- NA
     }
     fit.lognorm <- try(fitdistr(data, "log-normal"))
     if (class(fit.lognorm) == "try-error") {
          fit.lognorm$loglik <- NA
     }
     fit.type <- c("gamma", "weibull", "lognormal")
     distribution.type <- ifelse(is.null(distribution.type),
                                 fit.type[which.max(c(fit.gamma$loglik, 
                                                      fit.weib$loglik, 
                                                      fit.lognorm$loglik))],
                                 distribution.type)
     x <- NULL
     rm(x)
     if (distribution.type == "gamma") {
          shape <- fit.gamma$estimate[1]
          rate <- fit.gamma$estimate[2]
          mean <- shape/rate
          sd <- sqrt(shape)/rate
          return(list(distr = 'gamma',
                      shape = shape,
                      rate = rate,
                      mean = mean,
                      sd = sd,
                      loglik = fit.gamma$loglik))
     }
     else if (distribution.type == "weibull") {
          shape <- fit.weib$estimate[1]
          scale <- fit.weib$estimate[2]
          mean <- scale * exp(lgamma(1 + 1/shape))
          sd <- sqrt(scale^2 * (exp(lgamma(1 + 2/shape)) - (exp(lgamma(1 + 1/shape)))^2))
          return(list(distr = 'weibull',
                      shape = shape,
                      scale = scale,
                      mean = mean,
                      sd = sd,
                      loglik = fit.weib$loglik))
     }
     else if (distribution.type == "lognormal") {
          meanlog <- fit.lognorm$estimate[1]
          sdlog <- fit.lognorm$estimate[2]
          mean <- exp(1/2 * sdlog^2 + meanlog)
          sd <- sqrt(exp(2 * meanlog + sdlog^2) * (exp(sdlog^2) - 1))
          return(list(distr = 'lognormal',
                      meanlog = meanlog,
                      sdlog = sdlog,
                      mean = mean,
                      sd = sd,
                      loglik = fit.lognorm$loglik))
     }
}

diff_dis <- function(shape1, rate1, shape2, rate2, n){
     a <- sort(rgamma(n, shape1, rate1))
     b <- sort(rgamma(n, shape2, rate2))
     return(a - b)
}

# estimate SI -------------------------------------------------------------

library(R0)

datafile_onset <- datafile_info_BA1 %>% 
     filter(type != 'Asymptomatic') %>% 
     dplyr::select(id, dateonset)

datafile_cont_1 <- datafile_chains_BA1 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_onset, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = dateonset) %>% 
     left_join(datafile_onset, by = c('infector' = 'id')) %>% 
     rename('infector_date' = dateonset) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

datafile_onset <- datafile_info_BA2 %>% 
     filter(type != 'Asymptomatic') %>% 
     dplyr::select(id, dateonset)

datafile_cont_2 <- datafile_chains_BA2 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_onset, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = dateonset) %>% 
     left_join(datafile_onset, by = c('infector' = 'id')) %>% 
     rename('infector_date' = dateonset) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

datafile_onset <- datafile_info_Delta %>% 
     filter(type != 'Asymptomatic') %>% 
     dplyr::select(id, dateonset)

datafile_cont_3 <- datafile_chains_Delta %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_onset, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = dateonset) %>% 
     left_join(datafile_onset, by = c('infector' = 'id')) %>% 
     rename('infector_date' = dateonset) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

si_ba1_gamma <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                         distribution.type = 'gamma')

si_ba2_gamma <- fit_best(data = as.numeric(datafile_cont_2$date_sep), 
                         distribution.type = 'gamma')

si_delta_gamma <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                           distribution.type = 'gamma')

si_ba1_lognormal <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                             distribution.type = 'lognormal')

si_ba2_lognormal <- fit_best(data = as.numeric(datafile_cont_2$date_sep), 
                             distribution.type = 'lognormal')

si_delta_lognormal <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                               distribution.type = 'lognormal')

si_ba1_weibull <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                           distribution.type = 'weibull')

si_ba2_weibull <- fit_best(data = as.numeric(datafile_cont_2$date_sep), 
                           distribution.type = 'weibull')

si_delta_weibull <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                             distribution.type = 'weibull')

# plot --------------------------------------------------------------------

fig_a <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_ba2_gamma$shape, rate = si_ba2_gamma$rate),
                   mapping = aes(colour = 'Gamma')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = si_ba2_lognormal$meanlog, sdlog = si_ba2_lognormal$sdlog),
                   mapping = aes(colour = 'Log-normal')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = si_ba2_weibull$shape, scale = si_ba2_weibull$scale),
                   mapping = aes(colour = 'Weibull')) +
     geom_histogram(data = datafile_cont_2,
                    mapping = aes(x = date_sep, y=..density..),
                    breaks = 0:15,
                    fill = 'gray',
                    alpha = 0.5)+
     annotate(geom = 'text',
              x = 7.5,
              y = 0.4,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Serial interval')+
     coord_cartesian(clip = "off")+
     scale_colour_npg()+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(limits = c(0, 0.4),
                        expand = c(0, 0),
                        labels = label_number(accuracy = 0.01))+
     labs(x = '',
          y = 'Relative frequency',
          colour = 'Fitted distribution',
          title = 'a')

table_a <- data.frame(
     Distribution = c('Gamma', 'Log-normal', 'Weibull'),
     para = c(paste0('Y ~ Gamma(', round(si_ba2_gamma$shape,2),
                     ', ', round(si_ba2_gamma$rate, 2), ')'),
              paste0('ln(Y) ~ Normal(', round(si_ba2_lognormal$mean, 2), 
                     ', ', round(si_ba2_lognormal$sd, 2), ')'),
              paste0('Y~Weibull(', round(si_ba2_weibull$shape,2), 
                     ', ', round(si_ba2_weibull$scale, 2), ')')),
     loglike = round(c(si_ba2_gamma$loglik,
                       si_ba2_lognormal$loglik,
                       si_ba2_weibull$loglik),
                     2))
names(table_a)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_a <- tableGrob(table_a,
                     rows = NULL,
                     theme = ttheme_default(
                          base_size = 6,
                          core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                          base_family = 'Helvetica'
                     ))

fig_a <- fig_a +
     inset_element(table_a, 0.95, 0.7, 0.4, 0.75)

fig_d <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_ba1_gamma$shape, rate = si_ba1_gamma$rate),
                   mapping = aes(colour = 'Gamma')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = si_ba1_lognormal$meanlog, sdlog = si_ba1_lognormal$sdlog),
                   mapping = aes(colour = 'Log-normal')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = si_ba1_weibull$shape, scale = si_ba1_weibull$scale),
                   mapping = aes(colour = 'Weibull')) +
     geom_histogram(data = datafile_cont_1,
                    mapping = aes(x = date_sep, y=..density..),
                    breaks = 0:15,
                    fill = 'gray',
                    alpha = 0.5)+
     scale_colour_npg()+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(limits = c(0, 0.4),
                        expand = c(0, 0),
                        labels = label_number(accuracy = 0.01))+
     labs(x = '',
          y = 'Relative frequency',
          colour = 'Fitted distribution',
          title = 'd')

table_d <- data.frame(
     Distribution = c('Gamma', 'Log-normal', 'Weibull'),
     para = c(paste0('Y ~ Gamma(', round(si_ba1_gamma$shape,2),
                     ', ', round(si_ba1_gamma$rate, 2), ')'),
              paste0('ln(Y) ~ Normal(', round(si_ba1_lognormal$mean, 2), 
                     ', ', round(si_ba1_lognormal$sd, 2), ')'),
              paste0('Y~Weibull(', round(si_ba1_weibull$shape,2), 
                     ', ', round(si_ba1_weibull$scale, 2), ')')),
     loglike = round(c(si_ba1_gamma$loglik,
                       si_ba1_lognormal$loglik,
                       si_ba1_weibull$loglik),
                     2))
names(table_d)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_d <- tableGrob(table_d,
                     rows = NULL,
                     theme = ttheme_default(
                          base_size = 6,
                          core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                          base_family = 'Helvetica'
                     ))

fig_d <- fig_d +
     inset_element(table_d, 0.95, 0.7, 0.4, 0.75)

fig_g <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_delta_gamma$shape, rate = si_delta_gamma$rate),
                   mapping = aes(colour = 'Gamma')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = si_delta_lognormal$meanlog, sdlog = si_delta_lognormal$sdlog),
                   mapping = aes(colour = 'Log-normal')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = si_delta_weibull$shape, scale = si_delta_weibull$scale),
                   mapping = aes(colour = 'Weibull')) +
     geom_histogram(data = datafile_cont_3,
                    mapping = aes(x = date_sep, y=..density..),
                    breaks = 0:15,
                    fill = 'gray',
                    alpha = 0.5)+
     scale_colour_npg()+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0),
                        limits = c(0, 15))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.3),
                        labels = label_number(accuracy = 0.01))+
     labs(x = 'Time (days)',
          y = 'Relative frequency',
          colour = 'Fitted distribution',
          title = 'g')

table_g <- data.frame(
     Distribution = c('Gamma', 'Log-normal', 'Weibull'),
     para = c(paste0('Y ~ Gamma(', round(si_delta_gamma$shape,2),
                     ', ', round(si_delta_gamma$rate, 2), ')'),
              paste0('ln(Y) ~ Normal(', round(si_delta_lognormal$mean, 2), 
                     ', ', round(si_delta_lognormal$sd, 2), ')'),
              paste0('Y~Weibull(', round(si_delta_weibull$shape,2), 
                     ', ', round(si_delta_weibull$scale, 2), ')')),
     loglike = round(c(si_delta_gamma$loglik,
                       si_delta_lognormal$loglik,
                       si_delta_weibull$loglik),
                     2))
names(table_g)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_g <- tableGrob(table_g,
                     rows = NULL,
                     theme = ttheme_default(
                          base_size = 6,
                          core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                          base_family = 'Helvetica'
                     ))

fig_g <- fig_g +
     inset_element(table_g, 0.95, 0.7, 0.4, 0.75)

fig_a / fig_d / fig_g

# estimate TG -------------------------------------------------------------

datafile_positive <- datafile_info_BA1 %>% 
     dplyr::select(id, datepositive)

datafile_cont_1 <- datafile_chains_BA1 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_positive, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = datepositive) %>% 
     left_join(datafile_positive, by = c('infector' = 'id')) %>% 
     rename('infector_date' = datepositive) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

datafile_positive <- datafile_info_BA2 %>% 
     dplyr::select(id, datepositive)

datafile_cont_2 <- datafile_chains_BA2 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_positive, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = datepositive) %>% 
     left_join(datafile_positive, by = c('infector' = 'id')) %>% 
     rename('infector_date' = datepositive) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

datafile_positive <- datafile_info_Delta %>% 
     dplyr::select(id, datepositive)

datafile_cont_3 <- datafile_chains_Delta %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_positive, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = datepositive) %>% 
     left_join(datafile_positive, by = c('infector' = 'id')) %>% 
     rename('infector_date' = datepositive) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

tg_ba1_gamma <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                         distribution.type = 'gamma')

tg_ba2_gamma <- fit_best(data = as.numeric(datafile_cont_2$date_sep), 
                         distribution.type = 'gamma')

tg_delta_gamma <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                           distribution.type = 'gamma')

tg_ba1_lognormal <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                             distribution.type = 'lognormal')

tg_ba2_lognormal <- fit_best(data = as.numeric(datafile_cont_2$date_sep), 
                             distribution.type = 'lognormal')

tg_delta_lognormal <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                               distribution.type = 'lognormal')

tg_ba1_weibull <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                           distribution.type = 'weibull')

tg_ba2_weibull <- fit_best(data = as.numeric(datafile_cont_2$date_sep), 
                           distribution.type = 'weibull')

tg_delta_weibull <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                             distribution.type = 'weibull')

# plot --------------------------------------------------------------------

# fig_b <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
#      stat_function(fun = dgamma, n = 100,
#                    args = list(shape = tg_ba2_gamma$shape, rate = tg_ba2_gamma$rate),
#                    mapping = aes(colour = 'Gamma')) +
#      stat_function(fun = dlnorm, n = 100,
#                    args = list(meanlog = tg_ba2_lognormal$meanlog, sdlog = tg_ba2_lognormal$sdlog),
#                    mapping = aes(colour = 'Log-normal')) +
#      stat_function(fun = dweibull, n = 100,
#                    args = list(shape = tg_ba2_weibull$shape, scale = tg_ba2_weibull$scale),
#                    mapping = aes(colour = 'Weibull')) +
#      geom_histogram(data = datafile_cont_2,
#                     mapping = aes(x = date_sep, y=..density..),
#                     breaks = 0:15,
#                     fill = 'gray',
#                     alpha = 0.5)+
#      scale_colour_npg()+
#      scale_x_continuous(breaks = seq(0, 15, 3),
#                         expand = c(0, 0))+
#      scale_y_continuous(expand = expansion(mult = c(0, .1)),
#                         labels = label_number(accuracy = 0.01))+
#      labs(x = 'Time (days)',
#           y = 'Relative frequency',
#           colour = 'Fitted distribution',
#           title = 'b   TG Omicron BA.2')

fig_b <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_ba2_gamma$shape, rate = tg_ba2_gamma$rate),
                   mapping = aes(colour = 'Gamma')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = tg_ba2_lognormal$meanlog, sdlog = tg_ba2_lognormal$sdlog),
                   mapping = aes(colour = 'Log-normal')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = tg_ba2_weibull$shape, scale = tg_ba2_weibull$scale),
                   mapping = aes(colour = 'Weibull')) +
     geom_histogram(data = datafile_cont_2,
                    mapping = aes(x = date_sep, y=..density..),
                    breaks = 0:15,
                    fill = 'gray',
                    alpha = 0.5)+
     annotate(geom = 'text',
              x = 7.5,
              y = 0.5,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Transmission generation (probability)')+
     coord_cartesian(clip = "off")+
     scale_colour_npg()+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(limits = c(0, 0.5),
                        expand = c(0, 0),
                        labels = label_number(accuracy = 0.01))+
     labs(x = '',
          y = '',
          colour = 'Fitted distribution',
          title = 'b')

table_b <- data.frame(
     Distribution = c('Gamma', 'Log-normal', 'Weibull'),
     para = c(paste0('Y ~ Gamma(', round(tg_ba2_gamma$shape,2),
                     ', ', round(tg_ba2_gamma$rate, 2), ')'),
              paste0('ln(Y) ~ Normal(', round(tg_ba2_lognormal$mean, 2), 
                     ', ', round(tg_ba2_lognormal$sd, 2), ')'),
              paste0('Y~Weibull(', round(tg_ba2_weibull$shape,2), 
                     ', ', round(tg_ba2_weibull$scale, 2), ')')),
     loglike = round(c(tg_ba2_gamma$loglik,
                       tg_ba2_lognormal$loglik,
                       tg_ba2_weibull$loglik),
                     2))
names(table_b)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_b <- tableGrob(table_b,
                     rows = NULL,
                     theme = ttheme_default(
                          base_size = 6,
                          core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                          base_family = 'Helvetica'
                     ))

fig_b <- fig_b +
     inset_element(table_b, 0.95, 0.7, 0.4, 0.75)

fig_e <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_ba1_gamma$shape, rate = tg_ba1_gamma$rate),
                   mapping = aes(colour = 'Gamma')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = tg_ba1_lognormal$meanlog, sdlog = tg_ba1_lognormal$sdlog),
                   mapping = aes(colour = 'Log-normal')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = tg_ba1_weibull$shape, scale = tg_ba1_weibull$scale),
                   mapping = aes(colour = 'Weibull')) +
     geom_histogram(data = datafile_cont_1,
                    mapping = aes(x = date_sep, y=..density..),
                    breaks = 0:15,
                    fill = 'gray',
                    alpha = 0.5)+
     scale_colour_npg()+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(limits = c(0, 0.4),
                        expand = c(0, 0),
                        labels = label_number(accuracy = 0.01))+
     labs(x = '',
          y = '',
          colour = 'Fitted distribution',
          title = 'e')

table_e <- data.frame(
     Distribution = c('Gamma', 'Log-normal', 'Weibull'),
     para = c(paste0('Y ~ Gamma(', round(tg_ba1_gamma$shape,2),
                     ', ', round(tg_ba1_gamma$rate, 2), ')'),
              paste0('ln(Y) ~ Normal(', round(tg_ba1_lognormal$mean, 2), 
                     ', ', round(tg_ba1_lognormal$sd, 2), ')'),
              paste0('Y~Weibull(', round(tg_ba1_weibull$shape,2), 
                     ', ', round(tg_ba1_weibull$scale, 2), ')')),
     loglike = round(c(tg_ba1_gamma$loglik,
                       tg_ba1_lognormal$loglik,
                       tg_ba1_weibull$loglik),
                     2))
names(table_e)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_e <- tableGrob(table_e,
                     rows = NULL,
                     theme = ttheme_default(
                          base_size = 6,
                          core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                          base_family = 'Helvetica'
                     ))

fig_e <- fig_e +
     inset_element(table_e, 0.95, 0.7, 0.4, 0.75)

fig_h <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_delta_gamma$shape, rate = tg_delta_gamma$rate),
                   mapping = aes(colour = 'Gamma')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = tg_delta_lognormal$meanlog, sdlog = tg_delta_lognormal$sdlog),
                   mapping = aes(colour = 'Log-normal')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = tg_delta_weibull$shape, scale = tg_delta_weibull$scale),
                   mapping = aes(colour = 'Weibull')) +
     geom_histogram(data = datafile_cont_3,
                    mapping = aes(x = date_sep, y=..density..),
                    breaks = 0:15,
                    fill = 'gray',
                    alpha = 0.5)+
     scale_colour_npg()+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.4),
                        labels = label_number(accuracy = 0.01))+
     labs(x = 'Time (days)',
          y = '',
          colour = 'Fitted distribution',
          title = 'h')

table_h <- data.frame(
     Distribution = c('Gamma', 'Log-normal', 'Weibull'),
     para = c(paste0('Y ~ Gamma(', round(tg_delta_gamma$shape,2),
                     ', ', round(tg_delta_gamma$rate, 2), ')'),
              paste0('ln(Y) ~ Normal(', round(tg_delta_lognormal$mean, 2), 
                     ', ', round(tg_delta_lognormal$sd, 2), ')'),
              paste0('Y~Weibull(', round(tg_delta_weibull$shape,2), 
                     ', ', round(tg_delta_weibull$scale, 2), ')')),
     loglike = round(c(tg_delta_gamma$loglik,
                       tg_delta_lognormal$loglik,
                       tg_delta_weibull$loglik),
                     2))
names(table_h)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_h <- tableGrob(table_h,
                     rows = NULL,
                     theme = ttheme_default(
                          base_size = 6,
                          core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                          base_family = 'Helvetica'
                     ))

fig_h <- fig_h +
     inset_element(table_h, 0.95, 0.7, 0.4, 0.75)

fig_b / fig_e / fig_h

# estimate GT -------------------------------------------------------------

data_expose <- datafile_info_BA1 %>% 
     dplyr::select(id, dateexpose1, dateexpose2)

datafile_cont_1 <- datafile_chains_BA1 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(data_expose, by = c('infectee' = 'id')) %>% 
     rename(c('infectee_exposedate1' = dateexpose1,
              'infectee_exposedate2' = dateexpose2))%>% 
     left_join(data_expose, by = c('infector' = 'id')) %>% 
     rename(c('infector_exposedate1' = dateexpose1,
              'infector_exposedate2' = dateexpose2))

datafile_expose_seq <- t(sapply(rownames(datafile_cont_1), expose_seq, 
                                data = datafile_cont_1))

datafile_cont_1 <- datafile_cont_1 %>% 
     mutate(gt_min = datafile_expose_seq[,1],
            gt_median = datafile_expose_seq[,2],
            gt_max = datafile_expose_seq[,3])

gt_ba1_gamma <- fit_best(as.numeric(datafile_cont_1$gt_median), "gamma")
gt_ba1_lognormal <- fit_best(as.numeric(datafile_cont_1$gt_median), "lognormal")
gt_ba1_weibull <- fit_best(as.numeric(datafile_cont_1$gt_median), "weibull")

data_expose <- datafile_info_BA2 %>% 
     dplyr::select(id, dateexpose1, dateexpose2)

datafile_cont_2 <- datafile_chains_BA2 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(data_expose, by = c('infectee' = 'id')) %>% 
     rename(c('infectee_exposedate1' = dateexpose1,
              'infectee_exposedate2' = dateexpose2))%>% 
     left_join(data_expose, by = c('infector' = 'id')) %>% 
     rename(c('infector_exposedate1' = dateexpose1,
              'infector_exposedate2' = dateexpose2))

datafile_expose_seq <- t(sapply(rownames(datafile_cont_2), expose_seq, 
                                data = datafile_cont_2))

datafile_cont_2 <- datafile_cont_2 %>% 
     mutate(gt_min = datafile_expose_seq[,1],
            gt_median = datafile_expose_seq[,2],
            gt_max = datafile_expose_seq[,3])

gt_ba2_gamma <- fit_best(as.numeric(datafile_cont_2$gt_median), "gamma")
gt_ba2_lognormal <- fit_best(as.numeric(datafile_cont_2$gt_median), "lognormal")
gt_ba2_weibull <- fit_best(as.numeric(datafile_cont_2$gt_median), "weibull")

data_expose <- datafile_info_Delta %>% 
     dplyr::select(id, dateexpose1, dateexpose2)

datafile_cont_3 <- datafile_chains_Delta %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(data_expose, by = c('infectee' = 'id')) %>% 
     rename(c('infectee_exposedate1' = dateexpose1,
              'infectee_exposedate2' = dateexpose2))%>% 
     left_join(data_expose, by = c('infector' = 'id')) %>% 
     rename(c('infector_exposedate1' = dateexpose1,
              'infector_exposedate2' = dateexpose2)) |> 
     filter(!is.na(infector_exposedate1) & !is.na(infectee_exposedate1)) |> 
     mutate(gt_median = as.numeric(infectee_exposedate1 - infector_exposedate1),
            gt_median = if_else(gt_median<0, 0, gt_median))

gt_delta_gamma <- fit_best(as.numeric(datafile_cont_3$gt_median), "gamma")
gt_delta_lognormal <- fit_best(as.numeric(datafile_cont_3$gt_median), "lognormal")
gt_delta_weibull <- fit_best(as.numeric(datafile_cont_3$gt_median), "weibull")

# plot --------------------------------------------------------------------

fig_c <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_ba2_gamma$shape, rate = gt_ba2_gamma$rate),
                   mapping = aes(colour = 'Gamma')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = gt_ba2_lognormal$meanlog, sdlog = gt_ba2_lognormal$sdlog),
                   mapping = aes(colour = 'Log-normal')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = gt_ba2_weibull$shape, scale = gt_ba2_weibull$scale),
                   mapping = aes(colour = 'Weibull')) +
     geom_histogram(data = datafile_cont_2,
                    mapping = aes(x = gt_median, y=..density..),
                    breaks = 0:15,
                    fill = 'gray',
                    alpha = 0.5)+
     annotate(geom = 'text',
              x = 7.5,
              y = 0.5,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Generation time')+
     annotate(geom = 'text',
              x = 15,
              y = 0.25,
              vjust = -2,
              size = 11*5/14,
              angle = -90,
              family = 'Helvetica',
              label = 'Omicron BA.2')+
     coord_cartesian(clip = "off")+
     scale_colour_npg()+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.5),
                        labels = label_number(accuracy = 0.01))+
     labs(x = '',
          y = '',
          colour = 'Fitted distribution',
          title = 'c')


table_c <- data.frame(
     Distribution = c('Gamma', 'Log-normal', 'Weibull'),
     para = c(paste0('Y ~ Gamma(', round(gt_ba2_gamma$shape,2),
                     ', ', round(gt_ba2_gamma$rate, 2), ')'),
              paste0('ln(Y) ~ Normal(', round(gt_ba2_lognormal$mean, 2), 
                     ', ', round(gt_ba2_lognormal$sd, 2), ')'),
              paste0('Y~Weibull(', round(gt_ba2_weibull$shape,2), 
                     ', ', round(gt_ba2_weibull$scale, 2), ')')),
     loglike = round(c(gt_ba2_gamma$loglik,
                       gt_ba2_lognormal$loglik,
                       gt_ba2_weibull$loglik),
                     2))
names(table_c)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_c <- tableGrob(table_c,
                     rows = NULL,
                     theme = ttheme_default(
                          base_size = 6,
                          core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                          base_family = 'Helvetica'
                     ))

fig_c <- fig_c +
     inset_element(table_c, 0.95, 0.7, 0.4, 0.75)

fig_f <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_ba1_gamma$shape, rate = gt_ba1_gamma$rate),
                   mapping = aes(colour = 'Gamma')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = gt_ba1_lognormal$meanlog, sdlog = gt_ba1_lognormal$sdlog),
                   mapping = aes(colour = 'Log-normal')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = gt_ba1_weibull$shape, scale = gt_ba1_weibull$scale),
                   mapping = aes(colour = 'Weibull')) +
     geom_histogram(data = datafile_cont_1,
                    mapping = aes(x = gt_median, y=..density..),
                    breaks = 0:15,
                    fill = 'gray',
                    alpha = 0.5)+
     annotate(geom = 'text',
              x = 15,
              y = 0.3,
              vjust = -2,
              size = 11*5/14,
              angle = -90,
              family = 'Helvetica',
              label = 'Omicron BA.1')+
     coord_cartesian(clip = "off")+
     scale_colour_npg()+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.6),
                        labels = label_number(accuracy = 0.01))+
     labs(x = '',
          y = '',
          colour = 'Fitted distribution',
          title = 'f')

table_f <- data.frame(
     Distribution = c('Gamma', 'Log-normal', 'Weibull'),
     para = c(paste0('Y ~ Gamma(', round(gt_ba1_gamma$shape,2),
                     ', ', round(gt_ba1_gamma$rate, 2), ')'),
              paste0('ln(Y) ~ Normal(', round(gt_ba1_lognormal$mean, 2), 
                     ', ', round(gt_ba1_lognormal$sd, 2), ')'),
              paste0('Y~Weibull(', round(gt_ba1_weibull$shape,2), 
                     ', ', round(gt_ba1_weibull$scale, 2), ')')),
     loglike = round(c(gt_ba1_gamma$loglik,
                       gt_ba1_lognormal$loglik,
                       gt_ba1_weibull$loglik),
                     2))
names(table_f)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_f <- tableGrob(table_f,
                     rows = NULL,
                     theme = ttheme_default(
                          base_size = 6,
                          core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                          base_family = 'Helvetica'
                     ))

fig_f <- fig_f +
     inset_element(table_f, 0.95, 0.7, 0.4, 0.75)

fig_i <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_delta_gamma$shape, rate = gt_delta_gamma$rate),
                   mapping = aes(colour = 'Gamma')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = gt_delta_lognormal$meanlog, sdlog = gt_delta_lognormal$sdlog),
                   mapping = aes(colour = 'Log-normal')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = gt_delta_weibull$shape, scale = gt_delta_weibull$scale),
                   mapping = aes(colour = 'Weibull')) +
     geom_histogram(data = datafile_cont_3,
                    mapping = aes(x = gt_median, y=..density..),
                    breaks = 0:15,
                    fill = 'gray',
                    alpha = 0.5)+
     annotate(geom = 'text',
              x = 15,
              y = 0.25,
              vjust = -2,
              size = 11*5/14,
              angle = -90,
              family = 'Helvetica',
              label = 'Delta')+
     coord_cartesian(clip = "off")+
     scale_colour_npg()+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.5),
                        labels = label_number(accuracy = 0.01))+
     labs(x = 'Time (days)',
          y = '',
          colour = 'Fitted distribution',
          title = 'i')

table_i <- data.frame(
     Distribution = c('Gamma', 'Log-normal', 'Weibull'),
     para = c(paste0('Y ~ Gamma(', round(gt_delta_gamma$shape,2),
                     ', ', round(gt_delta_gamma$rate, 2), ')'),
              paste0('ln(Y) ~ Normal(', round(gt_delta_lognormal$mean, 2), 
                     ', ', round(gt_delta_lognormal$sd, 2), ')'),
              paste0('Y~Weibull(', round(gt_delta_weibull$shape,2), 
                     ', ', round(gt_delta_weibull$scale, 2), ')')),
     loglike = round(c(gt_delta_gamma$loglik,
                       gt_delta_lognormal$loglik,
                       gt_delta_weibull$loglik),
                     2))
names(table_i)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_i <- tableGrob(table_i,
                     rows = NULL,
                     theme = ttheme_default(
                          base_size = 6,
                          core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                          base_family = 'Helvetica'
                     ))

fig_i <- fig_i +
     inset_element(table_i, 0.95, 0.7, 0.4, 0.75)

# combined plot -----------------------------------------------------------

fig <- fig_a + fig_b + fig_c + fig_d + fig_e + fig_f + fig_g + fig_h + fig_i+
     plot_layout(guides = 'collect', ncol = 3) &
     theme_classic(base_family = 'Helvetica')+
     theme(axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.text.y = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', color = 'black'),
           axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
           axis.title.y = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
           strip.text.x = element_text(size = 12, face = 'bold'),
           strip.text.y = element_text(size = 12, face = 'bold'),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           plot.margin = margin(1, 20, 1, 1),
           plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
           legend.spacing.y = unit(3, 'pt'),
           legend.text = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', 
                                      color = 'black'),
           legend.title =element_text(size = 12, hjust = 0, vjust = .5, face = 'bold'),
           legend.margin=margin(0,5,0,0),
           legend.background = element_rect(fill = "transparent", colour = 'transparent'),
           legend.position = 'bottom')

fig

ggsave(filename = './outcome/publish/extend/Figure S1.pdf', 
       height = 8, width = 12)

ggsave(filename = './outcome/publish/extend/Figure S1.tiff', 
       height = 8, width = 12)
