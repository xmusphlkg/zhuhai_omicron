
# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(gamlss)

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
                                                      fit.weib$loglik, fit.lognorm$loglik))],
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
                      sd = sd))
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
                      sd = sd))
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
                      sd = sd))
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
     mutate(date_sep = infectee_date - infector_date,
            age_g = if_else(age <= 10, 'child', 'adult')) %>% 
     filter(date_sep >=0)

datafile_cont_c <- filter(datafile_cont_1, age_g == 'child')
datafile_cont_a <- filter(datafile_cont_1, age_g == 'adult')

si <- fit_best(data = as.numeric(datafile_cont_1$date_sep))
si_c <- fit_best(data = as.numeric(datafile_cont_c$date_sep), distribution.type = 'weibull')
si_a <- fit_best(data = as.numeric(datafile_cont_a$date_sep), distribution.type = 'weibull')

# estimate GT -------------------------------------------------------------

data_expose <- datafile_info_clean %>% 
     dplyr::select(id, dateexpose1, dateexpose2, age)

datafile_cont_2 <- datafile_cont %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(data_expose[,-4], by = c('infectee' = 'id')) %>% 
     rename(c('infectee_exposedate1' = dateexpose1,
              'infectee_exposedate2' = dateexpose2))%>% 
     left_join(data_expose, by = c('infector' = 'id')) %>% 
     rename(c('infector_exposedate1' = dateexpose1,
              'infector_exposedate2' = dateexpose2))

datafile_expose_seq <- t(sapply(rownames(datafile_cont_2), expose_seq))
datafile_cont_2 <- datafile_cont_2 %>% 
     mutate(gt_min = datafile_expose_seq[,1],
            age_g = if_else(age <= 10, 'child', 'adult'),
            gt_max = datafile_expose_seq[,3])

datafile_cont_c <- filter(datafile_cont_2, age_g == 'child')
datafile_cont_a <- filter(datafile_cont_2, age_g == 'adult')

fit_max <- fit_best(datafile_cont_c$gt_max)

gt_max_c <- fit_best(as.numeric(datafile_cont_c$gt_max), fit_max$distr)
gt_min_c <- fit_best(as.numeric(datafile_cont_c$gt_min + 0.002), fit_max$distr)

gt_max_a <- fit_best(as.numeric(datafile_cont_a$gt_max), fit_max$distr)
gt_min_a <- fit_best(as.numeric(datafile_cont_a$gt_min + 0.002), fit_max$distr)

gt_max <- fit_best(as.numeric(datafile_cont_2$gt_max), fit_max$distr)
gt_min <- fit_best(as.numeric(datafile_cont_2$gt_min + 0.002), fit_max$distr)

print(fit_max$distr)

# plot --------------------------------------------------------------------

library(ggsci)

fig_a <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = si$shape, scale = si$scale),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = si_a$shape, scale = si_a$scale),
                   mapping = aes(colour = 'B')) +
     stat_function(fun = dweibull, n = 100,
                   args = list(shape = si_c$shape, scale = si_c$scale),
                   mapping = aes(colour = 'C')) +
     scale_colour_npg(labels = c('Total', '> 10', '≤10'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = expansion(mult = c(0, .1)),
                        labels = label_number(accuracy = 0.01))+
     labs(x = 'Time (days)',
          y = 'Relative frequency',
          colour = 'Age of Infection',
          title = 'a')

fig_a

fig_b <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = gt_max$meanlog, sdlog = gt_max$sdlog),
                   mapping = aes(colour = 'A',
                                 linetype = 'GT (max)')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = gt_max_a$meanlog, sdlog = gt_max_a$sdlog),
                   mapping = aes(colour = 'B',
                                 linetype = 'GT (max)')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = gt_max_c$meanlog, sdlog = gt_max_c$sdlog),
                   mapping = aes(colour = 'C',
                                 linetype = 'GT (max)')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = gt_min$meanlog, sdlog = gt_min$sdlog),
                   mapping = aes(colour = 'A',
                                 linetype = 'GT (min)')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = gt_min_a$meanlog, sdlog = gt_min_a$sdlog),
                   mapping = aes(colour = 'B',
                                 linetype = 'GT (min)')) +
     stat_function(fun = dlnorm, n = 100,
                   args = list(meanlog = gt_min_c$meanlog, sdlog = gt_min_c$sdlog),
                   mapping = aes(colour = 'C',
                                 linetype = 'GT (min)')) +
     scale_colour_npg(labels = c('Total', '> 10', '≤10'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = expansion(mult = c(0, .1)),
                        labels = label_number(accuracy = 0.01))+
     labs(x = 'Time (days)',
          y = '',
          colour = 'Age of Infection',
          linetype = 'Generation Time',
          title = 'b')

fig_b

fig <- fig_a + fig_b  &
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
           legend.position = c(1,0.95),
           legend.justification = c(1, 1))

fig

ggsave(filename = './outcome/publish/extend/Figure S4.pdf', 
       device = cairo_pdf, height = 4, width = 8)
