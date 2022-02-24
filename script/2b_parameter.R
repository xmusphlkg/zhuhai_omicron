
# packages ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
library(scales)
library(ghibli)

extrafont::font_import('./script/fonts', prompt = F)
extrafont::loadfonts(device = "win")

remove(list = ls())

load('data/data_case.Rdata')

set.seed(202202)

# fit incubation time -----------------------------------------------------

datafile_info <- datafile_info_clean %>% 
  select(dateonset, dateexpose1, dateexpose2, type) %>% 
  mutate(date_seq_1 = dateonset - dateexpose1,
         date_seq_2 = dateonset - dateexpose2) %>% 
  filter(type != 'Asymptomatic')

library(epitrix)

expose_date <- function(a, b){
  return(as.Date(seq.Date(a, b, by = 'day')))
}

datafile_info$expose_date <- mapply(expose_date, datafile_info$dateexpose1, datafile_info$dateexpose2)

ib <- fit_gamma_incubation_dist(datafile_info, 
                                 dateonset, 
                                 dateexpose1,
                                 dateexpose2)

ib_out <- capture.output(ib)

# estimate SI -------------------------------------------------------------
## R0 provide a function to estimate Generation Time distribution
## This method actually base on estimate serial intervals

library(R0)

datafile_onset <- datafile_info_clean %>% 
  filter(type != 'Asymptomatic') %>% 
  dplyr::select(id, dateonset)

datafile_cont_1 <- datafile_cont %>% 
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

si <- est.GT(serial.interval = as.numeric(datafile_cont_1$date_sep))
si_out <- capture.output(si)

# estimate R0 -------------------------------------------------------------

data_onset <- datafile_info_clean %>% 
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

r0 <- est.R0.ML(data_onset$cases, si, begin = 1, end = 8)

r0_out <- capture.output(r0)

# estimate GT -------------------------------------------------------------

data_expose <- datafile_info_clean %>% 
  dplyr::select(id, dateexpose1, dateexpose2)

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

datafile_cont_2 <- datafile_cont %>% 
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

datafile_expose_seq <- t(sapply(rownames(datafile_cont_2), expose_seq))
datafile_cont_2 <- datafile_cont_2 %>% 
  mutate(gt_min = datafile_expose_seq[,1],
         # gt_median = datafile_expose_seq[,2],
         gt_max = datafile_expose_seq[,3])

gt_max <- fitdistr(as.numeric(datafile_cont_2$gt_max), "gamma")
gt_min <- fitdistr(as.numeric(datafile_cont_2$gt_min + 0.002), "gamma")

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

tg <- fitdistr(as.numeric(datafile_cont_3$date_sep), "gamma")

# write data --------------------------------------------------------------

## Serial Interval
time_si <- si
time_si_shape <- (time_si$mean/time_si$sd)^2
time_si_rate <- time_si_shape / time_si$mean

## Incubation Period
time_ib <- ib
time_ib_shape <- time_ib$distribution$parameters$shape
time_ib_scale <- time_ib$distribution$parameters$scale

outcome <- data.frame(
  name = c('Incubation Time', 'Serial Tnterval', 'Transmission Generation',
           'Generation Time (max)', 'Generation Time (min)'),
  shape = c(time_ib_shape, time_si_shape, tg$estimate[1],
            gt_max$estimate[1], gt_min$estimate[1]),
  rate = c(time_ib_scale, time_si_rate, tg$estimate[2],
           gt_max$estimate[2], gt_min$estimate[2])
) %>% 
  mutate(mean = shape/rate,
         sd = sqrt(shape)/rate)

write.csv(outcome, paste0('./outcome/share/Figure distribution.csv'),
          quote = F, row.names = F)

save(si, ib, r0, file = './outcome/para.RData')


# plot --------------------------------------------------------------------

library(cowplot)

fill_color <- ghibli::ghibli_palette('PonyoMedium', 5)

fig_a <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
  stat_function(fun = dgamma, n = 100,
                args = list(shape = time_ib_shape, scale = time_ib_scale),
                mapping = aes(colour = 'IP')) +
  stat_function(fun = dgamma, n = 100,
                args = list(shape = time_si_shape, rate = time_si_rate),
                mapping = aes(colour = 'SI')) +
  stat_function(fun = dgamma, n = 100,
                args = list(shape = tg$estimate[1], rate = tg$estimate[2]),
                mapping = aes(colour = 'TG')) +
  scale_color_manual(values = fill_color[1:3])+
  scale_x_continuous(breaks = seq(0, 15, 3),
                     expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     labels = label_number(accuracy = 0.01))+
  labs(x = '',
       y = 'Relative frequency',
       colour = '',
       title = 'A')

fig_b <- ggplot(data = data.frame(x = c(0, 15)), aes(x)) +
  stat_function(fun = dgamma, n = 100,
                args = list(shape = gt_max$estimate[1], rate = gt_max$estimate[2]),
                mapping = aes(colour = 'GT (max)')) +
  stat_function(fun = dgamma, n = 100,
                args = list(shape = gt_min$estimate[1], rate = gt_min$estimate[2]),
                mapping = aes(colour = 'GT (min)')) +
  scale_color_manual(values = fill_color[4:5])+
  scale_x_continuous(breaks = seq(0, 15, 3),
                     expand = c(0, 0))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)),
                     labels = label_number(accuracy = 0.01))+
  labs(x = 'Time (days)',
       y = 'Relative frequency',
       colour = '',
       title = 'B')

fig <- fig_a / fig_b  &
  theme_bw(base_family = 'Helvetica', base_line_size = 0.2)+
  theme(axis.title.x = element_text(face = 'bold', size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 5, color = 'black'),
        axis.text.y = element_text(size = 5, color = 'black'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.title = element_text(face = 'bold', size = 6),
        plot.title = element_text(size = 9, hjust = 0, vjust = .5, face = 'bold'),
        plot.margin = margin(1, 5, 1, 1),
        legend.spacing.y = unit(0, 'pt'),
        legend.text = element_text(margin = margin(r = 10, unit = "pt"),
                                   size = 5),
        legend.background = element_rect(fill = "transparent", colour = 'transparent'),
        legend.position = c(0.85,0.7))

fig

ggsave(filename = './outcome/science/Figure 2.pdf', 
       device = cairo_pdf, height = 10, width = 5.5, units = 'cm')

ggsave(filename = './outcome/science/Figure 2.tiff', 
       width = 5.5, height = 10, dpi = 300, units = 'cm')

ggsave(filename = './outcome/science/Figure 2.png', 
       width = 5.5, height = 10, dpi = 300, units = 'cm')
