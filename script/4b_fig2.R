
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

# estimate IP -------------------------------------------------------------

datafile_info <- datafile_info_BA1 %>% 
     dplyr::select(dateonset, dateexpose1, dateexpose2, type, age) %>% 
     mutate(date_seq_1 = dateonset - dateexpose1,
            date_seq_2 = dateonset - dateexpose2) %>% 
     filter(type != 'Asymptomatic')

datafile_info$expose_date <- mapply(expose_date, datafile_info$dateexpose1, datafile_info$dateexpose2)

ib_ba1 <- fit_gamma_incubation_dist(datafile_info, 
                                    dateonset, 
                                    dateexpose1,
                                    dateexpose2)
ib_ba1_c <- fit_gamma_incubation_dist(filter(datafile_info, age <= 18), 
                                      dateonset, 
                                      dateexpose1,
                                      dateexpose2)
ib_ba1_a <- fit_gamma_incubation_dist(filter(datafile_info, age > 18), 
                                      dateonset, 
                                      dateexpose1,
                                      dateexpose2)

datafile_BA2_ib <- datafile_info_BA2 %>% 
     dplyr::select(id, dateonset, dateexpose1, dateexpose2, type, age) %>% 
     mutate(date_seq_1 = dateonset - dateexpose1,
            date_seq_2 = dateonset - dateexpose2) %>% 
     filter(type != 'Asymptomatic' & !is.na(date_seq_1))

datafile_BA2_ib$expose_date <- mapply(expose_date, datafile_BA2_ib$dateexpose1, datafile_BA2_ib$dateexpose2)

ib_ba2 <- fit_gamma_incubation_dist(datafile_BA2_ib, 
                                    dateonset, 
                                    dateexpose1,
                                    dateexpose2) 
ib_ba2_c <- NA
ib_ba2_a <- fit_gamma_incubation_dist(filter(datafile_BA2_ib, age > 18), 
                                      dateonset, 
                                      dateexpose1,
                                      dateexpose2)


datafile_Delta_ib <- datafile_info_Delta %>% 
     dplyr::select(id, dateonset, dateexpose1, dateexpose2, type, age) %>% 
     mutate(date_seq_1 = dateonset - dateexpose1,
            date_seq_2 = dateonset - dateexpose2) %>% 
     filter(type != 'Asymptomatic' & !is.na(date_seq_1))

ib_delta <- fit_gamma_incubation_dist(datafile_Delta_ib,
                                      dateonset,
                                      dateexpose1,
                                      dateexpose2)
ib_delta_c <- fit_gamma_incubation_dist(filter(datafile_Delta_ib, age <= 18), 
                                        dateonset, 
                                        dateexpose1,
                                        dateexpose2)
ib_delta_a <- fit_gamma_incubation_dist(filter(datafile_Delta_ib, age > 18), 
                                        dateonset, 
                                        dateexpose1,
                                        dateexpose2)

# estimate SI -------------------------------------------------------------

library(R0)

datafile_onset <- datafile_info_BA1 %>% 
     filter(type != 'Asymptomatic') %>% 
     dplyr::select(id, dateonset, age)

datafile_cont_1 <- datafile_chains_BA1 %>% 
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

datafile_onset <- datafile_info_BA2 %>% 
     filter(type != 'Asymptomatic') %>% 
     dplyr::select(id, dateonset, age)

datafile_cont_2 <- datafile_chains_BA2 %>% 
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

datafile_onset <- datafile_info_Delta %>% 
     filter(type != 'Asymptomatic') %>% 
     dplyr::select(id, dateonset, age)

datafile_cont_3 <- datafile_chains_Delta %>% 
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

si_ba1_gamma_c <- fit_best(data = as.numeric(filter(datafile_cont_1, age <= 18)$date_sep), 
                           distribution.type = 'gamma')
si_ba1_gamma_a <- fit_best(data = as.numeric(filter(datafile_cont_1, age > 18)$date_sep), 
                           distribution.type = 'gamma')
si_ba1_gamma <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                         distribution.type = 'gamma')

si_ba2_gamma_c <- NA
si_ba2_gamma_a <- fit_best(data = as.numeric(filter(datafile_cont_2, age > 18)$date_sep), 
                           distribution.type = 'gamma')
si_ba2_gamma <- fit_best(data = as.numeric(datafile_cont_2$date_sep), 
                         distribution.type = 'gamma')

si_delta_gamma_c <- fit_best(data = as.numeric(filter(datafile_cont_3, age <= 18)$date_sep), 
                             distribution.type = 'gamma')
si_delta_gamma_a <- fit_best(data = as.numeric(filter(datafile_cont_3, age > 18)$date_sep), 
                             distribution.type = 'gamma')
si_delta_gamma <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                           distribution.type = 'gamma')

# estimate GT -------------------------------------------------------------

data_expose <- datafile_info_BA1 %>% 
     dplyr::select(id, dateexpose1, dateexpose2, age)

datafile_cont_1 <- datafile_chains_BA1 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(data_expose, by = c('infectee' = 'id')) %>% 
     rename(c('infectee_exposedate1' = dateexpose1,
              'infectee_exposedate2' = dateexpose2))%>% 
     left_join(data_expose[,-4], by = c('infector' = 'id')) %>% 
     rename(c('infector_exposedate1' = dateexpose1,
              'infector_exposedate2' = dateexpose2))

datafile_expose_seq <- t(sapply(rownames(datafile_cont_1), 
                                expose_seq, 
                                data = datafile_cont_1))

datafile_cont_1 <- datafile_cont_1 %>% 
     mutate(gt_min = datafile_expose_seq[,1],
            gt_median = datafile_expose_seq[,2],
            gt_max = datafile_expose_seq[,3])


data_expose <- datafile_info_BA2 %>% 
     dplyr::select(id, dateexpose1, dateexpose2, age)

datafile_cont_2 <- datafile_chains_BA2 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(data_expose, by = c('infectee' = 'id')) %>% 
     rename(c('infectee_exposedate1' = dateexpose1,
              'infectee_exposedate2' = dateexpose2))%>% 
     left_join(data_expose[,-4], by = c('infector' = 'id')) %>% 
     rename(c('infector_exposedate1' = dateexpose1,
              'infector_exposedate2' = dateexpose2))

datafile_expose_seq <- t(sapply(rownames(datafile_cont_2), 
                                expose_seq, 
                                data = datafile_cont_2))

datafile_cont_2 <- datafile_cont_2 %>% 
     mutate(gt_min = datafile_expose_seq[,1],
            gt_median = datafile_expose_seq[,2],
            gt_max = datafile_expose_seq[,3])

data_expose <- datafile_info_Delta %>% 
     dplyr::select(id, dateexpose1, dateexpose2, age)

datafile_cont_3 <- datafile_chains_Delta %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(data_expose, by = c('infectee' = 'id')) %>% 
     rename(c('infectee_exposedate1' = dateexpose1,
              'infectee_exposedate2' = dateexpose2))%>% 
     left_join(data_expose[,-4], by = c('infector' = 'id')) %>% 
     rename(c('infector_exposedate1' = dateexpose1,
              'infector_exposedate2' = dateexpose2)) |> 
     filter(!is.na(infector_exposedate1) & !is.na(infectee_exposedate1)) |> 
     mutate(gt_median = as.numeric(infectee_exposedate1 - infector_exposedate1),
            gt_median = if_else(gt_median<0, 0, gt_median))

datafile_expose_seq <- t(sapply(rownames(datafile_cont_3), 
                                expose_seq, 
                                data = datafile_cont_3))

datafile_cont_3 <- datafile_cont_3 %>% 
     mutate(gt_min = datafile_expose_seq[,1],
            gt_median = datafile_expose_seq[,2],
            gt_max = datafile_expose_seq[,3]) |> 
     filter(!is.na(gt_median))

gt_ba1_gamma_c <- fit_best(data = as.numeric(filter(datafile_cont_1, age <= 18)$gt_median), 
                           distribution.type = 'gamma')
gt_ba1_gamma_a <- fit_best(data = as.numeric(filter(datafile_cont_1, age > 18)$gt_median), 
                           distribution.type = 'gamma')
gt_ba1_gamma <- fit_best(data = as.numeric(datafile_cont_1$gt_median), 
                         distribution.type = 'gamma')

gt_ba2_gamma_c <- NA
gt_ba2_gamma_a <- fit_best(data = as.numeric(filter(datafile_cont_2, age > 18)$gt_median), 
                           distribution.type = 'gamma')
gt_ba2_gamma <- fit_best(data = as.numeric(datafile_cont_2$gt_median), 
                         distribution.type = 'gamma')

gt_delta_gamma_c <- fit_best(data = as.numeric(filter(datafile_cont_3, age <= 18)$gt_median), 
                             distribution.type = 'gamma')
gt_delta_gamma_a <- fit_best(data = as.numeric(filter(datafile_cont_3, age > 18)$gt_median), 
                             distribution.type = 'gamma')
gt_delta_gamma <- fit_best(data = as.numeric(datafile_cont_3$gt_median), 
                           distribution.type = 'gamma')

# estimate TG -------------------------------------------------------------

datafile_positive <- datafile_info_BA1 %>% 
     dplyr::select(id, datepositive, age)

datafile_cont_1 <- datafile_chains_BA1 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_positive, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = datepositive) %>% 
     left_join(datafile_positive[,-3], by = c('infector' = 'id')) %>% 
     rename('infector_date' = datepositive) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

datafile_positive <- datafile_info_BA2 %>% 
     dplyr::select(id, datepositive, age)

datafile_cont_2 <- datafile_chains_BA2 %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_positive, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = datepositive) %>% 
     left_join(datafile_positive[,-3], by = c('infector' = 'id')) %>% 
     rename('infector_date' = datepositive) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

datafile_positive <- datafile_info_Delta %>% 
     dplyr::select(id, datepositive, age)

datafile_cont_3 <- datafile_chains_Delta %>% 
     dplyr::select(to, from) %>% 
     rename(c(
          'infectee' = 'to',
          'infector' = 'from'
     )) %>% 
     left_join(datafile_positive, by = c('infectee' = 'id')) %>% 
     rename('infectee_date' = datepositive) %>% 
     left_join(datafile_positive[,-3], by = c('infector' = 'id')) %>% 
     rename('infector_date' = datepositive) %>% 
     mutate(date_sep = infectee_date - infector_date) %>% 
     filter(date_sep >=0)

tg_ba1_gamma_c <- fit_best(data = as.numeric(filter(datafile_cont_1, age <= 18)$date_sep), 
                           distribution.type = 'gamma')
tg_ba1_gamma_a <- fit_best(data = as.numeric(filter(datafile_cont_1, age > 18)$date_sep), 
                           distribution.type = 'gamma')
tg_ba1_gamma <- fit_best(data = as.numeric(datafile_cont_1$date_sep), distribution.type = 'gamma')

tg_ba2_gamma_c <- NA
tg_ba2_gamma_a <- fit_best(data = as.numeric(filter(datafile_cont_2, age > 18)$date_sep), 
                           distribution.type = 'gamma')
tg_ba2_gamma <- fit_best(data = as.numeric(datafile_cont_2$date_sep), distribution.type = 'gamma')

tg_delta_gamma_c <- fit_best(data = as.numeric(filter(datafile_cont_3, age <= 18)$date_sep), 
                             distribution.type = 'gamma')
tg_delta_gamma_a <- fit_best(data = as.numeric(filter(datafile_cont_3, age > 18)$date_sep), 
                             distribution.type = 'gamma')
tg_delta_gamma <- fit_best(data = as.numeric(datafile_cont_3$date_sep), distribution.type = 'gamma')


# plot IP -----------------------------------------------------------------

fig_a <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_ba2$distribution$parameters$shape, 
                               scale = ib_ba2$distribution$parameters$scale),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_ba2_a$distribution$parameters$shape, 
                               scale = ib_ba2_a$distribution$parameters$scale),
                   mapping = aes(colour = 'B')) +
     # stat_function(fun = dgamma, n = 100,
     #               args = list(shape = ib_ba2_c$distribution$parameters$shape, 
     #                           scale = ib_ba2_c$distribution$parameters$scale),
     #               mapping = aes(colour = 'C')) +
     annotate(geom = 'text',
              x = 7.5,
              y = 0.3,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Incubation period')+
     coord_cartesian(clip = "off")+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.3),
                        labels = label_number(accuracy = 0.01))+
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     labs(x = '',
          y = 'Relative frequency',
          colour = '',
          title = 'a')+
     guides(colour = 'none')

fig_e <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_ba1$distribution$parameters$shape, 
                               scale = ib_ba1$distribution$parameters$scale),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_ba1_a$distribution$parameters$shape, 
                               scale = ib_ba1_a$distribution$parameters$scale),
                   mapping = aes(colour = 'B')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_ba1_c$distribution$parameters$shape,
                               scale = ib_ba1_c$distribution$parameters$scale),
                   mapping = aes(colour = 'C')) +
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = expansion(mult = c(0, .1)),
                        labels = label_number(accuracy = 0.01))+
     labs(x = '',
          y = 'Relative frequency',
          colour = '',
          title = 'e')

fig_i <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_delta$distribution$parameters$shape, 
                               scale = ib_delta$distribution$parameters$scale),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_delta_a$distribution$parameters$shape, 
                               scale = ib_delta_a$distribution$parameters$scale),
                   mapping = aes(colour = 'B')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_delta_c$distribution$parameters$shape,
                               scale = ib_delta_c$distribution$parameters$scale),
                   mapping = aes(colour = 'C')) +
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = expansion(mult = c(0, .1)),
                        labels = label_number(accuracy = 0.01))+
     labs(x = 'Time (days)',
          y = 'Relative frequency',
          colour = '',
          title = 'i')

# plot SI -----------------------------------------------------------------

fig_b <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_ba2_gamma$shape, 
                               rate = si_ba2_gamma$rate),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_ba2_gamma_a$shape, 
                               rate = si_ba2_gamma_a$rate),
                   mapping = aes(colour = 'B')) +
     # stat_function(fun = dgamma, n = 100,
     #               args = list(shape = si_ba2_gamma_c$shape, 
     #                           rate = si_ba2_gamma_c$rate),
     #               mapping = aes(colour = 'C')) +
     annotate(geom = 'text',
              x = 7.5,
              y = 0.4,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Serial interval')+
     coord_cartesian(clip = "off")+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.4),
                        labels = label_number(accuracy = 0.01))+
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     labs(x = '',
          y = '',
          colour = '',
          title = 'b')+
     guides(colour = 'none')

fig_f <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_ba1_gamma$shape, 
                               rate = si_ba1_gamma$rate),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_ba1_gamma_a$shape, 
                               rate = si_ba1_gamma_a$rate),
                   mapping = aes(colour = 'B')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_ba1_gamma_c$shape,
                               rate = si_ba1_gamma_c$rate),
                   mapping = aes(colour = 'C')) +
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = expansion(mult = c(0, .1)),
                        labels = label_number(accuracy = 0.01))+
     labs(x = '',
          y = '',
          colour = '',
          title = 'f')

fig_j <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_delta_gamma$shape, 
                               rate = si_delta_gamma$rate),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_delta_gamma_a$shape, 
                               rate = si_delta_gamma_a$rate),
                   mapping = aes(colour = 'B')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_delta_gamma_c$shape,
                               rate = si_delta_gamma_c$rate),
                   mapping = aes(colour = 'C')) +
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = expansion(mult = c(0, .1)),
                        labels = label_number(accuracy = 0.01))+
     labs(x = 'Time (days)',
          y = '',
          colour = '',
          title = 'j')
# plot TG -----------------------------------------------------------------

fig_c <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_ba2_gamma$shape, 
                               rate = tg_ba2_gamma$rate),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_ba2_gamma_a$shape, 
                               rate = tg_ba2_gamma_a$rate),
                   mapping = aes(colour = 'B')) +
     # stat_function(fun = dgamma, n = 100,
     #               args = list(shape = tg_ba2_gamma_c$shape, 
     #                           rate = tg_ba2_gamma_c$rate),
     #               mapping = aes(colour = 'C')) +
     annotate(geom = 'text',
              x = 7.5,
              y = 0.5,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Transmission generation')+
     coord_cartesian(clip = "off")+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.5),
                        labels = label_number(accuracy = 0.01))+
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     labs(x = '',
          y = '',
          colour = '',
          title = 'c')+
     guides(colour = 'none')

fig_g <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_ba1_gamma$shape, 
                               rate = tg_ba1_gamma$rate),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_ba1_gamma_a$shape, 
                               rate = tg_ba1_gamma_a$rate),
                   mapping = aes(colour = 'B')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_ba1_gamma_c$shape,
                               rate = tg_ba1_gamma_c$rate),
                   mapping = aes(colour = 'C')) +
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = expansion(mult = c(0, .1)),
                        labels = label_number(accuracy = 0.01))+
     labs(x = '',
          y = '',
          colour = '',
          title = 'g')

fig_k <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_delta_gamma$shape, 
                               rate = tg_delta_gamma$rate),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_delta_gamma_a$shape, 
                               rate = tg_delta_gamma_a$rate),
                   mapping = aes(colour = 'B')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_delta_gamma_c$shape,
                               rate = tg_delta_gamma_c$rate),
                   mapping = aes(colour = 'C')) +
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = expansion(mult = c(0, .1)),
                        labels = label_number(accuracy = 0.01))+
     labs(x = 'Time (days)',
          y = '',
          colour = '',
          title = 'k')

# plot GT -----------------------------------------------------------------

fig_d <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_ba2_gamma$shape, 
                               rate = gt_ba2_gamma$rate),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_ba2_gamma_a$shape, 
                               rate = gt_ba2_gamma_a$rate),
                   mapping = aes(colour = 'B')) +
     # stat_function(fun = dgamma, n = 100,
     #               args = list(shape = gt_ba2_gamma_c$shape, 
     #                           rate = gt_ba2_gamma_c$rate),
     #               mapping = aes(colour = 'C')) +
     annotate(geom = 'text',
              x = 7.5,
              y = 0.4,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Generation time(probability)')+
     annotate(geom = 'text',
              x = 15,
              y = 0.2,
              vjust = -0.7,
              size = 11*5/14,
              angle = -90,
              family = 'Helvetica',
              label = 'Omicron BA.2')+
     coord_cartesian(clip = "off")+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.4),
                        labels = label_number(accuracy = 0.01))+
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     labs(x = '',
          y = '',
          colour = '',
          title = 'd')+
     guides(colour = 'none')

fig_h <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_ba1_gamma$shape, 
                               rate = gt_ba1_gamma$rate),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_ba1_gamma_a$shape, 
                               rate = gt_ba1_gamma_a$rate),
                   mapping = aes(colour = 'B')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_ba1_gamma_c$shape,
                               rate = gt_ba1_gamma_c$rate),
                   mapping = aes(colour = 'C')) +
     annotate(geom = 'text',
              x = 15,
              y = 0.35,
              vjust = -0.7,
              size = 11*5/14,
              angle = -90,
              family = 'Helvetica',
              label = 'Omicron BA.1')+
     coord_cartesian(clip = "off")+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.7),
                        labels = label_number(accuracy = 0.01))+
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     labs(x = '',
          y = '',
          colour = '',
          title = 'h')

fig_l <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_delta_gamma$shape, 
                               rate = gt_delta_gamma$rate),
                   mapping = aes(colour = 'A')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_delta_gamma_a$shape, 
                               rate = gt_delta_gamma_a$rate),
                   mapping = aes(colour = 'B')) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_delta_gamma_c$shape,
                               rate = gt_delta_gamma_c$rate),
                   mapping = aes(colour = 'C')) +
     annotate(geom = 'text',
              x = 15,
              y = 0.3,
              vjust = -0.7,
              size = 11*5/14,
              angle = -90,
              family = 'Helvetica',
              label = 'Delta')+
     coord_cartesian(clip = "off")+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.6),
                        labels = label_number(accuracy = 0.01))+
     scale_colour_nejm(labels = c('Total', '19-', '0-18'))+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     labs(x = 'Time (days)',
          y = '',
          colour = '',
          title = 'l')

# combined plot -----------------------------------------------------------

fig_a + fig_b +fig_c + fig_d + fig_e + fig_f + fig_g + fig_h + fig_i + fig_j + fig_k + fig_l +
     plot_layout(guides = 'collect', ncol = 4)&
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

ggsave(filename = './outcome/publish/extend/Figure S2.pdf', 
       height = 8, width = 12)

ggsave(filename = './outcome/publish/extend/Figure S2.tiff', 
       height = 8, width = 12)

