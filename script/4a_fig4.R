
# packages ----------------------------------------------------------------
# devtools::install_github("reconhub/epitrix")
# devtools::install_github('tobadia/R0')
library(tidyverse)
library(patchwork)
library(scales)
library(openxlsx)
library(ggsci)

library(showtext)
font_add('Helvetica', 
         "./script/fonts/Helvetica.ttf",
         bold = "./script/fonts/Helvetica-Bold.ttf")

remove(list = ls())

load('data/sars_2_cov.Rdata')

set.seed(202202)

# function ----------------------------------------------------------------

library(epitrix)
library(MASS)

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

# diff_dis <- function(shape1, rate1, shape2, rate2, n){
#      a <- rgamma(n, shape1, rate1)
#      b <- rgamma(n, shape2, rate2)
#      return(a - b)
# }

# fit incubation time -----------------------------------------------------

datafile_BA1_ib <- datafile_info_BA1 %>% 
     dplyr::select(id, dateonset, dateexpose1, dateexpose2, type) %>% 
     mutate(date_seq_1 = dateonset - dateexpose1,
            date_seq_2 = dateonset - dateexpose2) %>% 
     filter(type != 'Asymptomatic')

datafile_BA1_ib$expose_date <- mapply(expose_date, datafile_BA1_ib$dateexpose1, datafile_BA1_ib$dateexpose2)

ib_ba1 <- fit_gamma_incubation_dist(datafile_BA1_ib, 
                                    dateonset, 
                                    dateexpose1,
                                    dateexpose2)

datafile_BA2_ib <- datafile_info_BA2 %>% 
     dplyr::select(id, dateonset, dateexpose1, dateexpose2, type) %>% 
     mutate(date_seq_1 = dateonset - dateexpose1,
            date_seq_2 = dateonset - dateexpose2) %>% 
     filter(type != 'Asymptomatic' & !is.na(date_seq_1))

datafile_BA2_ib$expose_date <- mapply(expose_date, datafile_BA2_ib$dateexpose1, datafile_BA2_ib$dateexpose2)

ib_ba2 <- fit_gamma_incubation_dist(datafile_BA2_ib, 
                                    dateonset, 
                                    dateexpose1,
                                    dateexpose2) 


datafile_Delta_ib <- datafile_info_Delta %>% 
     dplyr::select(id, dateonset, dateexpose1, dateexpose2, type) %>% 
     mutate(date_seq_1 = dateonset - dateexpose1,
            date_seq_2 = dateonset - dateexpose2) %>% 
     filter(type != 'Asymptomatic' & !is.na(date_seq_1))

ib_delta <- fit_gamma_incubation_dist(datafile_Delta_ib,
                                      dateonset,
                                      dateexpose1,
                                      dateexpose2)

# estimate SI -------------------------------------------------------------
## R0 provide a function to estimate Generation Time distribution
## This method actually base on estimate serial intervals

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

si_ba1 <- fit_best(data = as.numeric(datafile_cont_1$date_sep), 
                   distribution.type = 'gamma')

si_ba2 <- fit_best(data = as.numeric(datafile_cont_2$date_sep), 
                   distribution.type = 'gamma')

si_delta <- fit_best(data = as.numeric(datafile_cont_3$date_sep), 
                     distribution.type = 'gamma')

# estimate differ ---------------------------------------------------------

datafile_ib <- datafile_BA1_ib |> dplyr::select(id, date_seq_1, date_seq_2)
datafile_ib$ib <- mapply(seq, from = datafile_ib$date_seq_2, to = datafile_ib$date_seq_1)
datafile_ib <- datafile_ib |> 
        mutate(freq = 1/n()) |> 
        unnest(cols = ib) |> 
        group_by(id) |> 
        mutate(freq = freq/n()) |> 
        dplyr::select(id, ib, freq)

datafile_si <- datafile_cont_1 |> 
        dplyr::select(infectee, date_sep) |> 
        group_by(infectee) |> 
        summarise(si = median(date_sep),
                  .groups = 'drop') |> 
        rename(id = 'infectee')

datafile_differ_BA1 <- datafile_ib |> 
        left_join(datafile_si) |> 
        mutate(diff = round(si - ib)) |> 
        ungroup() |> 
        filter(!is.na(diff)) |> 
        mutate(freq = as.numeric(freq/sum(freq))) |> 
        group_by(diff) |> 
        summarise(freq = sum(freq),
                  .groups = 'drop')

datafile_ib <- datafile_BA2_ib |> dplyr::select(id, date_seq_1, date_seq_2)
datafile_ib$ib <- mapply(seq, from = datafile_ib$date_seq_2, to = datafile_ib$date_seq_1)
datafile_ib <- datafile_ib |> 
        mutate(freq = 1/n()) |> 
        unnest(cols = ib) |> 
        group_by(id) |> 
        mutate(freq = freq/n()) |> 
        dplyr::select(id, ib, freq)

datafile_si <- datafile_cont_2 |> 
        dplyr::select(infectee, date_sep) |> 
        group_by(infectee) |> 
        summarise(si = median(date_sep),
                  .groups = 'drop') |> 
        rename(id = 'infectee')

datafile_differ_BA2 <- datafile_ib |> 
        left_join(datafile_si) |> 
        mutate(diff = round(si - ib)) |> 
        ungroup() |> 
        filter(!is.na(diff)) |> 
        mutate(freq = as.numeric(freq/sum(freq))) |> 
        group_by(diff) |> 
        summarise(freq = sum(freq),
                  .groups = 'drop')

datafile_ib <- datafile_Delta_ib |> dplyr::select(id, date_seq_1, date_seq_2)
datafile_ib$ib <- mapply(seq, from = datafile_ib$date_seq_2, to = datafile_ib$date_seq_1)
datafile_ib <- datafile_ib |> 
        mutate(freq = 1/n()) |> 
        unnest(cols = ib) |> 
        group_by(id) |> 
        mutate(freq = freq/n()) |> 
        dplyr::select(id, ib, freq)

datafile_si <- datafile_cont_3 |> 
        dplyr::select(infectee, date_sep) |> 
        group_by(infectee) |> 
        summarise(si = median(date_sep),
                  .groups = 'drop') |> 
        rename(id = 'infectee')

datafile_differ_Delta <- datafile_ib |> 
        left_join(datafile_si) |> 
        mutate(diff = round(si - ib)) |> 
        ungroup() |> 
        filter(!is.na(diff)) |> 
        mutate(freq = as.numeric(freq/sum(freq))) |> 
        group_by(diff) |> 
        summarise(freq = sum(freq),
                  .groups = 'drop')

# estimate R0 -------------------------------------------------------------

data_onset_1 <- datafile_info_BA1 %>%
     group_by(dateonset) %>%
     count() %>%
     rename(c(
          'cases' = 'n',
          'date' = 'dateonset'
     )) %>%
     as.data.frame() %>%
     complete(
          date = seq.Date(from = min(date), to = max(date), by = 'day'),
          fill = list(cases = 0)
     )

data_onset_2 <- datafile_info_BA2 %>%
     group_by(dateonset) %>%
     count() %>%
     rename(c(
          'cases' = 'n',
          'date' = 'dateonset'
     )) %>%
     as.data.frame() %>%
     complete(
          date = seq.Date(from = min(date), to = max(date), by = 'day'),
          fill = list(cases = 0)
     )

data_onset_3 <- datafile_info_Delta %>%
     group_by(dateonset) %>%
     count() %>%
     rename(c(
          'cases' = 'n',
          'date' = 'dateonset'
     )) %>%
     as.data.frame() %>%
     complete(
          date = seq.Date(from = min(date), to = max(date), by = 'day'),
          fill = list(cases = 0)
     )

r0_ba1 <- est.R0.ML(data_onset_1$cases, 
                    generation.time('gamma', c(si_ba1$mean, si_ba1$sd)), 
                    begin = 1, end = 8)
# r0_out <- capture.output(r0)
r0_ba2 <- est.R0.ML(data_onset_2$cases, 
                    generation.time('gamma', c(si_ba2$mean, si_ba2$sd)), 
                    begin = 1, end = 5)

r0_delta <- est.R0.ML(data_onset_3$cases, 
                      generation.time('gamma', c(si_ba2$mean, si_ba2$sd)), 
                      begin = 1, end = 11)

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

gt_ba1 <- fit_best(as.numeric(datafile_cont_1$gt_median), "gamma")

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

gt_ba2 <- fit_best(as.numeric(datafile_cont_2$gt_median), "gamma")

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

gt_delta <- fit_best(as.numeric(datafile_cont_3$gt_median), "gamma")

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

tg_ba1 <- fit_best(as.numeric(datafile_cont_1$date_sep), "gamma")
tg_ba2 <- fit_best(as.numeric(datafile_cont_2$date_sep), "gamma")
tg_delta <- fit_best(as.numeric(datafile_cont_3$date_sep), "gamma")

# write data --------------------------------------------------------------

save(si_ba1, si_ba2, si_delta, 
     ib_ba1, ib_ba2, ib_delta, 
     r0_ba1, r0_ba2, r0_delta, 
     file = './outcome/base/para.RData')

# outcome <- data.frame(
#         name = c('Incubation period', 'Serial interval',
#                  'Transmission generation',
#                  'Generation time (probability)',
#                  'Effective reproductive number'),
#         shape = c(time_ib_shape, time_si_shape, tg$estimate[1],
#                   gt_max$estimate[1], gt_min$estimate[1]),
#         rate = c(time_ib_scale, time_si_rate, tg$estimate[2],
#                  gt_max$estimate[2], gt_min$estimate[2])
# )

# write.csv(outcome, paste0('./outcome/base/parameter distribution.csv'),
#           quote = F, row.names = F)

# plot --------------------------------------------------------------------

library(cowplot)

fill_color <- rev(ggsci::pal_nejm()(4))[-3]

## plot ib -----------------------------------------------------------------

datafile <- data.frame(
     parameter = 'Incubation Period',
     variant = c('Delta', 'BA.1', 'BA.2'),
     mean = c(ib_delta$mu, ib_ba1$mu, ib_ba2$mu),
     sd = c(ib_delta$sd, ib_ba1$sd, ib_ba2$sd),
     shape = c(ib_delta$distribution$parameters$shape,
               ib_ba1$distribution$parameters$shape,
               ib_ba2$distribution$parameters$shape),
     rate = 1/c(ib_delta$distribution$parameters$scale,
               ib_ba1$distribution$parameters$scale,
               ib_ba2$distribution$parameters$scale)
)

datafile_table <- datafile

fig_ib_inside <- ggplot(data = datafile,
                        mapping = aes(x = variant,
                                      y = mean,
                                      fill = variant))+
     geom_bar(stat="identity", color="black", 
              position=position_dodge())+
     geom_errorbar(mapping = aes(ymin=mean-sd, ymax=mean+sd), 
                   width=.2,
                   position=position_dodge())+
     scale_fill_manual(values = fill_color,
                       breaks = datafile$variant)+
     scale_x_discrete(limits = datafile$variant)+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 12))+
     theme_bw(base_family = 'Helvetica')+
     theme(axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, face = 'plain', color = 'black'))+
     labs(x = '',
          y = '')+
     guides(fill = 'none')

fig_ib <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_delta$distribution$parameters$shape, 
                               scale = ib_delta$distribution$parameters$scale),
                   mapping = aes(colour = 'Delta'))+
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_ba1$distribution$parameters$shape, 
                               scale = ib_ba1$distribution$parameters$scale),
                   mapping = aes(colour = 'BA.1'))+
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = ib_ba2$distribution$parameters$shape, 
                               scale = ib_ba2$distribution$parameters$scale),
                   mapping = aes(colour = 'BA.2'))+
     annotate(geom = 'text',
              x = 7.5,
              y = 0.4,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Incubation period')+
     coord_cartesian(clip = "off")+
     scale_color_manual(values = fill_color,
                        breaks = datafile$variant)+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.4),
                        labels = label_number(accuracy = 0.1))+
     theme_classic(base_family = 'Helvetica')+
     labs(x = '',
          y = 'Relative frequency',
          colour = '',
          title = 'b')+
     guides(color = 'none')

fig_ib <- fig_ib + inset_element(fig_ib_inside, 0.4, 0.2, 1, 1)

## plot si -----------------------------------------------------------------

datafile <- data.frame(
     parameter = 'Serial interval',
     variant = c('Delta', 'BA.1', 'BA.2'),
     mean = c(si_delta$mean, si_ba1$mean, si_ba2$mean),
     sd = c(si_delta$sd, si_ba1$sd, si_ba2$sd),
     shape = c(si_delta$shape,
               si_ba1$shape,
               si_ba2$shape),
     rate = c(si_delta$rate,
               si_ba1$rate,
               si_ba2$rate)
)

datafile_table <- rbind(datafile, datafile_table)

fig_si_inside <- ggplot(data = datafile,
                        mapping = aes(x = variant,
                                      y = mean,
                                      fill = variant))+
     geom_bar(stat="identity", color="black", 
              position=position_dodge())+
     geom_errorbar(mapping = aes(ymin=mean-sd, ymax=mean+sd), 
                   width=.2,
                   position=position_dodge())+
     scale_fill_manual(values = fill_color,
                       breaks = datafile$variant)+
     scale_x_discrete(limits = datafile$variant)+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 12))+
     theme_bw(base_family = 'Helvetica')+
     theme(axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, face = 'plain', color = 'black'))+
     labs(x = '',
          y = '')+
     guides(fill = 'none')

fig_si <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_delta$shape, 
                               rate = si_delta$rate),
                   mapping = aes(colour = 'Delta'))+
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_ba1$shape, 
                               rate = si_ba1$rate),
                   mapping = aes(colour = 'BA.1'))+
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = si_ba2$shape, 
                               rate = si_ba2$rate),
                   mapping = aes(colour = 'BA.2'))+
     annotate(geom = 'text',
              x = 7.5,
              y = 0.4,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Serial interval')+
     coord_cartesian(clip = "off")+
     scale_color_manual(values = fill_color,
                        breaks = datafile$variant)+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.4),
                        labels = label_number(accuracy = 0.1))+
     labs(x = 'Time (days)',
          y = 'Relative frequency',
          colour = 'Variants',
          title = 'a')+
     theme_classic(base_family = 'Helvetica')+
     theme(legend.position = 'bottom')+ 
     guides(color = 'none')

fig_si <- fig_si + 
     inset_element(fig_si_inside, 0.4, 0.2, 1, 1)

## plot tg -----------------------------------------------------------------

datafile <- data.frame(
     parameter = 'Transmission generation',
     variant = c('Delta', 'BA.1', 'BA.2'),
     shape = c(tg_delta$shape, tg_ba1$shape, tg_ba2$shape),
     rate = c(tg_delta$rate, tg_ba1$rate, tg_ba2$rate)
) |> 
     mutate(mean = shape/rate,
            sd = sqrt(shape)/rate)

datafile_table <- rbind(datafile, datafile_table)

fig_tg_inside <- ggplot(data = datafile,
                        mapping = aes(x = variant,
                                      y = mean,
                                      fill = variant))+
     geom_bar(stat="identity", color="black", 
              position=position_dodge())+
     geom_errorbar(mapping = aes(ymin=ifelse(mean-sd < 0, 0.01, mean-sd), ymax=mean+sd), 
                   width=.2,
                   position=position_dodge())+
     scale_fill_manual(values = fill_color,
                       breaks = datafile$variant)+
     scale_x_discrete(limits = datafile$variant)+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 12))+
     theme_bw(base_family = 'Helvetica')+
     theme(axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, face = 'plain', color = 'black'))+
     labs(x = '',
          y = '')+
     guides(fill = 'none')

fig_tg <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_delta$shape, 
                               rate = tg_delta$rate),
                   mapping = aes(colour = 'Delta'))+
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_ba1$shape, 
                               rate = tg_ba1$rate),
                   mapping = aes(colour = 'BA.1'))+
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = tg_ba2$shape, 
                               rate = tg_ba2$rate),
                   mapping = aes(colour = 'BA.2'))+
     annotate(geom = 'text',
              x = 7.5,
              y = 0.5,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Transmission generation')+
     coord_cartesian(clip = "off")+
     scale_color_manual(values = fill_color,
                        breaks = datafile$variant)+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.5),
                        labels = label_number(accuracy = 0.1))+
     theme_classic(base_family = 'Helvetica')+
     labs(x = 'Time (days)',
          y = '',
          colour = '',
          title = 'e')+
     guides(color = 'none')

fig_tg <- fig_tg + 
     inset_element(fig_tg_inside, 0.4, 0.2, 1, 1)

## plot gt -----------------------------------------------------------------

datafile <- data.frame(
     parameter = 'Generation time',
     variant = c('Delta', 'BA.1', 'BA.2'),
     shape = c(gt_delta$shape, gt_ba1$shape, gt_ba2$shape),
     rate = c(gt_delta$rate, gt_ba1$rate, gt_ba2$rate)
) |> 
     mutate(mean = shape/rate,
            sd = sqrt(shape)/rate)

datafile_table <- rbind(datafile, datafile_table)

fig_gt_inside <- ggplot(data = datafile,
                        mapping = aes(x = variant,
                                      y = mean,
                                      fill = variant))+
     geom_bar(stat="identity", color="black", 
              position=position_dodge())+
     geom_errorbar(mapping = aes(ymin=ifelse(mean-sd < 0, 0.01, mean-sd), ymax=mean+sd), 
                   width=.2,
                   position=position_dodge())+
     scale_fill_manual(values = fill_color,
                       breaks = datafile$variant)+
     scale_x_discrete(limits = datafile$variant)+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 12))+
     theme_bw(base_family = 'Helvetica')+
     theme(axis.text.x = element_text(size = 10, hjust = .5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.text.y = element_text(size = 10, hjust = 1, vjust = 0.5, face = 'plain', color = 'black'))+
     labs(x = '',
          y = '')+
     guides(fill = 'none')

fig_gt <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_delta$shape, 
                               rate = gt_delta$rate),
                   mapping = aes(colour = 'Delta'))+
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_ba1$shape, 
                               rate = gt_ba1$rate),
                   mapping = aes(colour = 'BA.1'))+
     stat_function(fun = dgamma, n = 100,
                   args = list(shape = gt_ba2$shape, 
                               rate = gt_ba2$rate),
                   mapping = aes(colour = 'BA.2'))+
     annotate(geom = 'text',
              x = 7.5,
              y = 0.5,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Generation time (probability)')+
     coord_cartesian(clip = "off")+
     scale_color_manual(values = fill_color, 
                        breaks = datafile$variant)+
     scale_x_continuous(breaks = seq(0, 15, 3),
                        expand = c(0, 0))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.5),
                        labels = label_number(accuracy = 0.1))+
     theme_classic(base_family = 'Helvetica')+
     labs(x = '',
          y = '',
          colour = '',
          title = 'd')+
     guides(color = 'none')

fig_gt <- fig_gt + 
     inset_element(fig_gt_inside, 0.4, 0.2, 1, 1)


## plot differ -----------------------------------------------------------

# diff_delta <- diff_dis(si_delta$shape,si_delta$rate,
#                        ib_delta$distribution$parameters$shape,
#                        1 / ib_delta$distribution$parameters$scale,
#                        100000)
# diff_ba1 <- diff_dis(si_ba1$shape,si_ba1$rate,
#                      ib_ba1$distribution$parameters$shape,
#                      1 / ib_ba1$distribution$parameters$scale,
#                      100000)
# diff_ba2 <- diff_dis(si_ba2$shape,
#                      si_ba2$rate,
#                      ib_ba2$distribution$parameters$shape,
#                      1 / ib_ba2$distribution$parameters$scale,
#                      100000)

datafile_differ_BA1$lineage <- 'BA.1'
datafile_differ_BA2$lineage <- 'BA.2'
datafile_differ_Delta$lineage <- 'Delta'

datafile <- rbind(datafile_differ_BA1,
                  datafile_differ_BA2,
                  datafile_differ_Delta) |> 
        group_by(lineage) |> 
        mutate(lineage = factor(lineage,
                                levels = c('Delta', 'BA.1', 'BA.2')),
               diff = as.numeric(diff),
               cfreq = cumsum(freq))



fig_diff <- ggplot(datafile)+
     geom_col(mapping = aes(x = diff,
                            y = freq,
                            fill = lineage),
              color = 'white',
              position = position_dodge2(width = 1, preserve = "single"),
              show.legend = T)+
     scale_fill_manual(values = fill_color)+
     scale_x_continuous(expand = c(0, 0),
                        limits = c(-12, 6))+
     scale_y_continuous(expand = c(0, 0),
                        limits = c(0, 0.4))+
     annotate(geom = 'text',
              x = -3,
              y = Inf,
              vjust = -0.7,
              size = 11*5/14,
              family = 'Helvetica',
              label = 'Differ between SI and IP')+
     coord_cartesian(clip = "off")+
     theme_classic(base_family = 'Helvetica')+
     labs(x = 'Time (days)',
          y = 'Relative frequency',
          colour = '',
          title = 'c')+
     guides(fill = 'none')

# plot Reff ---------------------------------------------------------------

datafile <- data.frame(
     variant = c('Delta', 'BA.1', 'BA.2'),
     mean = c(r0_delta$R, r0_ba1$R, r0_ba2$R),
     low_ci = c(r0_delta$conf.int[1], r0_ba1$conf.int[1], r0_ba2$conf.int[1]),
     top_ci = c(r0_delta$conf.int[2], r0_ba1$conf.int[2], r0_ba2$conf.int[2])
) |> 
     mutate_if(is.numeric, round, digits = 2) |> 
     mutate(labels = paste0(mean, 
                            " (", low_ci, "-", top_ci, ")"))

fig_reff <- ggplot(data = datafile,
                   mapping = aes(y = variant,
                                 x = mean,
                                 fill = variant))+
     geom_bar(stat="identity", color="black", 
              position=position_dodge())+
     geom_errorbar(mapping = aes(xmin=low_ci, xmax=top_ci), 
                   width=.2,
                   position=position_dodge())+
     geom_text(aes(label = labels),
               color = fill_color,
               position = position_dodge(width = 0.7),
               size = 10*5/14,
               vjust = -1,
               hjust = -0.1
     )+
     # annotate(geom = 'text',
     #          x = 5,
     #          y = Inf,
     #          vjust = -0.7,
     #          size = 11*5/14,
     #          family = 'Helvetica',
     #          label = 'Effective reproductive number')+
     coord_cartesian(clip = "off")+
     scale_fill_manual(values = fill_color,
                       breaks = datafile$variant)+
     scale_y_discrete(limits = datafile$variant)+
     scale_x_continuous(expand = c(0, 0),
                        limits = c(0, 10),
                        breaks = seq(0, 10, 2))+
     theme_classic(base_family = 'Helvetica')+
     labs(x = 'Effective reproductive number',
          y = '',
          colour = '',
          title = 'f')+
     guides(fill = 'none')

# multiplot ---------------------------------------------------------------

fig <- fig_si + fig_ib + fig_diff + fig_gt + fig_tg + fig_reff+
     plot_layout(ncol = 2, byrow = F) &
     theme(axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 0.5, face = 'plain', color = 'black'),
           axis.text.y = element_text(size = 10, hjust = 1, vjust = .5, face = 'plain', color = 'black'),
           axis.title.x = element_text(size = 12, hjust = .5, vjust = .5, face = 'bold'),
           axis.title.y = element_text(size = 12, hjust = .5, vjust = 0, face = 'bold'),
           strip.text.x = element_text(size = 12, face = 'bold'),
           strip.text.y = element_text(size = 12, face = 'bold'),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin(1, 5, 1, 1),
           plot.background = element_blank(),
           plot.title = element_text(size = 16, hjust = 0, vjust = 0, face = 'bold'))

fig

ggsave(filename = './outcome/publish/Figure 4.pdf', 
       height = 7, width = 7)

# write table -------------------------------------------------------------

datafile_table <- datafile_table |> 
     mutate_if(is.numeric, round, digits = 2)

write.csv(datafile_table,
          file = './outcome/base/Table 3.csv',
          row.names = F)
