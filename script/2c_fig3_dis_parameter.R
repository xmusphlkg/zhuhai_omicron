
# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(gamlss)
library(ggsci)
library(gridExtra)

extrafont::font_import('./script/fonts', prompt = F)
extrafont::loadfonts(device = "win")

remove(list = ls())

load('data/data_case.Rdata')

set.seed(202202)

## R0 provide a function to estimate Generation Time distribution
## This method actually base on estimate serial intervals

library(R0)

# function ----------------------------------------------------------------

expose_seq <- function(x){
     infector_expose_dates <- seq.Date(from = datafile_cont_2[x, 'infector_exposedate1'],
                                       to = datafile_cont_2[x, 'infector_exposedate2'],
                                       by = 'day')
     infectee_expose_dates <- seq.Date(from = datafile_cont_2[x, 'infectee_exposedate1'],
                                       to = datafile_cont_2[x, 'infectee_exposedate2'],
                                       by = 'day')
     datafile <- as.data.frame(expand.grid(infectee_expose_dates, infector_expose_dates)) %>% 
          mutate(median = as.numeric(Var2 - Var1),
                 median = ifelse(median<0, 0, median))
     gt_median <- median(datafile$median)
     gt_max <- max(infectee_expose_dates) - min(infector_expose_dates)
     gt_min <- min(infectee_expose_dates) - max(infector_expose_dates)
     gt_min <- ifelse(gt_min<0, 0, gt_min)
     return(c(gt_min, gt_median, gt_max))
}

fit_best <- function(data, distribution.type = NULL){
     data[data == 0] <- 1/2
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

# estimate SI -------------------------------------------------------------
datafile_onset <- datafile_info_clean %>% 
     filter(type != 'Asymptomatic') %>% 
     dplyr::select(id, dateonset, age)

datafile_cont_1 <- datafile_cont %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_onset, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = dateonset) %>% 
     left_join(datafile_onset[,-3], by = c('infector' = 'id')) %>% 
     rename('infector_date' = dateonset) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

si_gamma <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                     distribution.type = 'gamma')
si_lognorm <- fit_best(data = as.numeric(datafile_cont_1$date_sep),
                       distribution.type = 'lognormal')
si_weibull <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                       distribution.type = 'weibull')

fig_a <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
        stat_function(fun = dgamma, n = 100,
                      args = list(shape = si_gamma$shape, rate = si_gamma$rate),
                      mapping = aes(colour = 'Gamma')) +
        stat_function(fun = dlnorm, n = 100,
                      args = list(meanlog = si_lognorm$meanlog, sdlog = si_lognorm$sdlog),
                      mapping = aes(colour = 'Log-normal')) +
        stat_function(fun = dweibull, n = 100,
                      args = list(shape = si_weibull$shape, scale = si_weibull$scale),
                      mapping = aes(colour = 'Weibull')) +
        geom_histogram(data = datafile_cont_1,
                       mapping = aes(x = date_sep, y=..density..),
                       breaks = 0:15,
                       fill = 'gray',
                       alpha = 0.5)+
        scale_colour_npg()+
        scale_x_continuous(breaks = seq(0, 15, 3),
                           expand = c(0, 0))+
        scale_y_continuous(expand = expansion(mult = c(0, .1)),
                           labels = label_number(accuracy = 0.01))+
        labs(x = 'Time (days)',
             y = 'Relative frequency',
             colour = 'Fitted distribution',
             title = 'a')

fig_a

# estimate TG -------------------------------------------------------------

datafile_positive <- datafile_info_clean %>% 
        dplyr::select(id, datepositive)

datafile_cont_3 <- datafile_cont %>% 
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

tg_gamma <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                     distribution.type = 'gamma')
tg_lognorm <- fit_best(data = as.numeric(datafile_cont_3$date_sep),
                       distribution.type = 'lognormal')
tg_weibull <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                       distribution.type = 'weibull')

fig_b <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
        stat_function(fun = dgamma, n = 100,
                      args = list(shape = tg_gamma$shape, rate = tg_gamma$rate),
                      mapping = aes(colour = 'Gamma')) +
        stat_function(fun = dlnorm, n = 100,
                      args = list(meanlog = tg_lognorm$meanlog, sdlog = tg_lognorm$sdlog),
                      mapping = aes(colour = 'Log-normal')) +
        stat_function(fun = dweibull, n = 100,
                      args = list(shape = tg_weibull$shape, scale = tg_weibull$scale),
                      mapping = aes(colour = 'Weibull')) +
        geom_histogram(data = datafile_cont_3,
                       mapping = aes(x = date_sep, y=..density..),
                       breaks = 0:15,
                       fill = 'gray',
                       alpha = 0.5)+
        scale_colour_npg()+
        scale_x_continuous(breaks = seq(0, 15, 3),
                           expand = c(0, 0))+
        scale_y_continuous(expand = expansion(mult = c(0, .1)),
                           labels = label_number(accuracy = 0.01))+
        labs(x = 'Time (days)',
             y = 'Relative frequency',
             colour = 'Fitted distribution',
             title = 'b')

fig_b

# plot --------------------------------------------------------------------

table_a <- data.frame(
        Distribution = c('Gamma', 'Log-normal', 'Weibull'),
        para = c(paste0('Y ~ Gamma(', round(si_gamma$shape,2),
                        ', ', round(si_gamma$rate, 2), ')'),
                 paste0('ln(Y) ~ Normal(', round(si_lognorm$mean, 2), 
                        ', ', round(si_lognorm$sd, 2), ')'),
                 paste0('Y~Weibull(', round(si_weibull$shape,2), 
                        ', ', round(si_weibull$scale, 2), ')')),
        loglike = round(c(si_gamma$loglik,
                          si_lognorm$loglik,
                          si_weibull$loglik),
                        2))
names(table_a)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_a <- tableGrob(table_a,
                     rows = NULL,
                     theme = ttheme_default(
                             base_size = 6,
                             core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                             base_family = 'Times New Roman'
                     ))

fig_a_1 <- fig_a +
        inset_element(table_a, 0.9, 0.7, 0.4, 0.95)


table_b <- data.frame(
        Distribution = c('Gamma', 'Log-normal', 'Weibull'),
        para = c(paste0('Y ~ Gamma(', round(tg_gamma$shape,2),
                        ', ', round(tg_gamma$rate, 2), ')'),
                 paste0('ln(Y) ~ Normal(', round(tg_lognorm$mean, 2), 
                        ', ', round(tg_lognorm$sd, 2), ')'),
                 paste0('Y~Weibull(', round(tg_weibull$shape,2), 
                        ', ', round(tg_weibull$scale, 2), ')')),
        loglike = round(c(tg_gamma$loglik,
                          tg_lognorm$loglik,
                          tg_weibull$loglik),
                        2))
names(table_b)[2:3] <- c('Parameter', 'Logarithmic\nlikelihood')

table_b <- tableGrob(table_b,
                     rows = NULL,
                     theme = ttheme_default(
                             base_size = 6,
                             core = list(fg_params = list(col = ggsci::pal_npg()(3))),
                             base_family = 'Times New Roman'
                     ))

fig_b_1 <- fig_b +
        inset_element(table_b, 0.9, 0.7, 0.4, 0.95)
fig_b_1

fig <- fig_a_1 + fig_b_1  &
     theme_bw(base_family = 'Helvetica')+
     theme(axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.text.y = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', color = 'black'),
           axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
           strip.text.x = element_text(size = 12, face = 'bold'),
           strip.text.y = element_text(size = 12, face = 'bold'),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           plot.margin = margin(1, 5, 1, 1),
           plot.title = element_text(size = 16, hjust = 0, vjust = .5, face = 'bold'),
           legend.spacing.y = unit(3, 'pt'),
           legend.text = element_text(size = 10, hjust = 0, vjust = .5, face = 'plain', 
                                      color = 'black'),
           legend.title =element_text(size = 12, hjust = 0, vjust = .5, face = 'bold'),
           legend.margin=margin(0,5,0,0),
           legend.background = element_rect(fill = "transparent", colour = 'transparent'),
           legend.position = 'none')

fig

ggsave(filename = './outcome/publish/extend/Figure S3.pdf', 
       device = cairo_pdf, height = 4, width = 8)
